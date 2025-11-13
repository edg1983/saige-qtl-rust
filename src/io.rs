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
    let pheno_df = CsvReadOptions::default()
        .with_has_header(true)
        .with_parse_options(
            CsvParseOptions::default()
                .with_separator(b'\t')
        )
        .try_into_reader_with_file_path(Some(pheno_file.into()))?
        .finish()?;

    // 2. Load Covariates
    // In Polars 0.41.x, CsvReadOptions is used to configure the reader
    let covar_df = CsvReadOptions::default()
        .with_has_header(true)
        .with_parse_options(
            CsvParseOptions::default()
                .with_separator(b'\t')
        )
        .try_into_reader_with_file_path(Some(covar_file.into()))?
        .finish()?;

    // 3. Align pheno and covar to master list
    // This join ensures we are in the correct order and subset
    let aligned_pheno = master_df.clone()
        .lazy()
        .join(
            pheno_df.clone().lazy(),
            [col("MASTER_SAMPLES")],
            [col(sample_id_pheno)],
            JoinType::Inner.into(),
        )
        // CORRECTED: `drop_columns` is now `drop`
        .drop(["MASTER_SAMPLES"])
        .collect()?;

    // CORRECTED: Re-added the missing definition for `aligned_covar`
    let aligned_covar = master_df.clone()
        .lazy()
        .join(
            covar_df.clone().lazy(),
            [col("MASTER_SAMPLES")],
            [col(sample_id_covar)],
            JoinType::Inner.into(),
        )
        .drop(["MASTER_SAMPLES"])
        .collect()?;

    // 4. Get final sample list (those present in all files)
    let final_samples_df = master_df
        .lazy()
        .join(
            pheno_df.lazy(),
            [col("MASTER_SAMPLES")],
            [col(sample_id_pheno)],
            JoinType::Inner.into(),
        )
        .join(
            covar_df.lazy(),
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
        // CORRECTED: Use `ok_or` for `Option`
        .ok_or(IoError::Alignment("Failed to unwrap sample IDs".into()))?;

    // 5. Convert to ndarray
    // CORRECTED: `to_ndarray` is gone. We must collect into a Vec first.
    let y_vec: Vec<f64> = aligned_pheno
        .column(trait_name)?
        .f64()?
        .into_iter()
        .collect::<Option<Vec<f64>>>()
        // CORRECTED: `map_err` is not for `Option`. Use `ok_or`.
        .ok_or(IoError::Alignment("Phenotype column contains nulls".into()))?;
    let y: Array1<f64> = Array1::from_vec(y_vec);

    let n_samples = sample_ids.len();
    // CORRECTED: `aligned_covar` is now in scope
    let n_covars = aligned_covar.width();
    
    // Add an intercept term
    let intercept = Array1::ones(n_samples);
    let mut x = Array2::zeros((n_samples, n_covars + 1));
    x.column_mut(0).assign(&intercept);
    
    // CORRECTED: `aligned_covar` is now in scope
    // FIXED: `IndexOrder::F` is now `IndexOrder::Fortran`.
    let covar_matrix_no_intercept = aligned_covar.to_ndarray::<Float64Type>(IndexOrder::Fortran)?;
    x.slice_mut(s![.., 1..]).assign(&covar_matrix_no_intercept); // This line now works

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