//! Module for handling I/O, data loading, and sample alignment.
// ADDED: Import IndexOrder
use polars::prelude::{IndexOrder, *};
use std::path::Path;
use ndarray::{Array1, Array2, s}; // <-- Added 's' macro import
use thiserror::Error;

#[derive(Error, Debug)]
pub enum IoError {
    #[error("Polars error: {0}")]
    Polars(#[from] PolarsError),
    #[error("File not found: {0}")]
    NotFound(String),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    // CORRECTED: Added the missing Alignment variant
    #[error("Data alignment error: {0}")]
    Alignment(String),
}

/// Holds the aligned phenotype (y) and covariate (X) data.
pub struct AlignedData {
    /// Phenotype vector (y)
    pub y: Array1<f64>,
    /// Covariate matrix (X), including intercept
    pub x: Array2<f64>,
    /// Aligned sample IDs, in order
    pub sample_ids: Vec<String>,
}

/// Loads phenotype and covariate files, aligns them to a given
/// list of master samples (e.g., from a PLINK .fam file).
pub fn load_aligned_data(
    pheno_file: &Path,
    covar_file: &Path,
    trait_name: &str,
    sample_id_pheno: &str,
    sample_id_covar: &str,
    master_sample_ids: &[String],
) -> Result<AlignedData, IoError> {
    
    // Create a DataFrame for master sample IDs to join against
    let master_df = DataFrame::new(vec![
        Series::new("MASTER_SAMPLES", master_sample_ids)
    ])?;

    // 1. Load Phenotype
    // In Polars 0.41.x, CsvReadOptions is used to configure the reader
    let mut pheno_df = CsvReadOptions::default()
        .with_has_header(true)
        .with_parse_options(
            CsvParseOptions::default()
                .with_separator(b'\t')
        )
        .try_into_reader_with_file_path(Some(pheno_file.into()))?
        .finish()?;
    
    // Cast sample ID column to string to handle both numeric and string IDs
    let pheno_id_col = pheno_df.column(sample_id_pheno)?.clone();
    let pheno_id_as_str = pheno_id_col.cast(&DataType::String)?;
    pheno_df.replace(sample_id_pheno, pheno_id_as_str)?;

    // 2. Load Covariates
    // In Polars 0.41.x, CsvReadOptions is used to configure the reader
    let mut covar_df = CsvReadOptions::default()
        .with_has_header(true)
        .with_parse_options(
            CsvParseOptions::default()
                .with_separator(b'\t')
        )
        .try_into_reader_with_file_path(Some(covar_file.into()))?
        .finish()?;
    
    // Cast sample ID column to string to handle both numeric and string IDs
    let covar_id_col = covar_df.column(sample_id_covar)?.clone();
    let covar_id_as_str = covar_id_col.cast(&DataType::String)?;
    covar_df.replace(sample_id_covar, covar_id_as_str)?;

    // 3. Get final sample list (those present in all files)
    let final_samples_df = master_df
        .lazy()
        .join(
            pheno_df.clone().lazy(),
            [col("MASTER_SAMPLES")],
            [col(sample_id_pheno)],
            JoinType::Inner.into(),
        )
        .join(
            covar_df.clone().lazy(),
            [col("MASTER_SAMPLES")],
            [col(sample_id_covar)],
            JoinType::Inner.into(),
        )
        .select([col("MASTER_SAMPLES")])
        .collect()?;

    let sample_ids: Vec<String> = final_samples_df
        .column("MASTER_SAMPLES")?
        .str()?
        .into_iter()
        .map(|opt_s| opt_s.map(String::from))
        .collect::<Option<Vec<String>>>()
        .ok_or(IoError::Alignment("Failed to unwrap sample IDs".into()))?;

    // 4. Now filter pheno and covar to only include the final sample intersection
    let final_samples_series = Series::new("FINAL_SAMPLES", &sample_ids);
    let final_samples_for_filter = DataFrame::new(vec![final_samples_series])?;
    
    let filtered_pheno = final_samples_for_filter.clone()
        .lazy()
        .join(
            pheno_df.clone().lazy(),
            [col("FINAL_SAMPLES")],
            [col(sample_id_pheno)],
            JoinType::Inner.into(),
        )
        .drop(["FINAL_SAMPLES"])
        .collect()?;
    
    let filtered_covar = final_samples_for_filter
        .lazy()
        .join(
            covar_df.clone().lazy(),
            [col("FINAL_SAMPLES")],
            [col(sample_id_covar)],
            JoinType::Inner.into(),
        )
        .drop(["FINAL_SAMPLES"])
        .collect()?;

    // 5. Convert to ndarray
    let y_vec: Vec<f64> = filtered_pheno
        .column(trait_name)?
        .f64()?
        .into_iter()
        .collect::<Option<Vec<f64>>>()
        .ok_or(IoError::Alignment("Phenotype column contains nulls".into()))?;
    let y: Array1<f64> = Array1::from_vec(y_vec);

    let n_samples = sample_ids.len();
    let n_covars = filtered_covar.width();
    
    // Add an intercept term
    let intercept = Array1::ones(n_samples);
    let mut x = Array2::zeros((n_samples, n_covars + 1));
    x.column_mut(0).assign(&intercept);
    
    let covar_matrix_no_intercept = filtered_covar.to_ndarray::<Float64Type>(IndexOrder::Fortran)?;
    x.slice_mut(s![.., 1..]).assign(&covar_matrix_no_intercept);

    if y.len() != n_samples || x.nrows() != n_samples {
        // CORRECTED: Use `IoError::Alignment`
        return Err(IoError::Alignment(format!(
            "Dimension mismatch after alignment: y ({}), X ({}), samples ({})",
            y.len(), x.nrows(), n_samples
        )));
    }

    Ok(AlignedData { y, x, sample_ids })
}

/// Reads sample IDs from a .fam file
pub fn get_fam_samples(plink_file_no_ext: &Path) -> Result<Vec<String>, IoError> {
    let fam_path = plink_file_no_ext.with_extension("fam");
    if !fam_path.exists() {
        return Err(IoError::NotFound(fam_path.to_string_lossy().into()));
    }
    
    // In Polars 0.41.x, CsvReadOptions is used to configure the reader
    // .fam files have 6 columns: FID IID PID MID SEX PHENO
    let fam_df = CsvReadOptions::default()
        .with_has_header(false)
        .with_parse_options(
            CsvParseOptions::default()
                .with_separator(b'\t')
        )
        .try_into_reader_with_file_path(Some(fam_path.clone().into()))?
        .finish()?;

    // Get the second column (IID) by index
    // Cast to string to handle both numeric and string IDs
    let iid_column = fam_df.get_columns()[1].clone(); // Get the second column (0-indexed)
    let iid_as_str = iid_column.cast(&DataType::String)?;
    
    let sample_ids: Vec<String> = iid_as_str
        .str()?
        .into_iter()
        .map(|opt_s| opt_s.map(String::from))
        .collect::<Option<Vec<String>>>()
        .ok_or(IoError::Alignment("Failed to parse .fam sample IDs".into()))?;

    Ok(sample_ids)
}