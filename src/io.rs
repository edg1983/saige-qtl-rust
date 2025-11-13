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

/// Loads a single file containing both phenotypes and covariates, then subsets
/// to samples present in the master sample list (e.g., from a PLINK .fam file).
pub fn load_aligned_data(
    pheno_covar_file: &Path,
    trait_name: &str,
    sample_id_col: &str,
    covariate_cols: &[String],
    master_sample_ids: &[String],
) -> Result<AlignedData, IoError> {
    
    log::info!("Loading phenotype/covariate file: {:?}", pheno_covar_file);
    
    // 1. Load the combined phenotype/covariate file
    let mut data_df = CsvReadOptions::default()
        .with_has_header(true)
        .with_parse_options(
            CsvParseOptions::default()
                .with_separator(b'\t')
        )
        .try_into_reader_with_file_path(Some(pheno_covar_file.into()))?
        .finish()?;
    
    log::info!("Loaded data file with {} rows and {} columns", data_df.height(), data_df.width());
    
    // Cast sample ID column to string to handle both numeric and string IDs
    let sample_id_series = data_df.column(sample_id_col)?.clone();
    let sample_id_as_str = sample_id_series.cast(&DataType::String)?;
    data_df.replace(sample_id_col, sample_id_as_str)?;
    
    // 2. Filter to only rows where sample ID is in the FAM file
    log::info!("Filtering to samples present in {} FAM samples", master_sample_ids.len());
    
    let master_df = DataFrame::new(vec![
        Series::new("MASTER_SAMPLES", master_sample_ids)
    ])?;
    
    // Semi-join: keep only rows from data_df where sample_id is in master_df
    let filtered_data = data_df
        .lazy()
        .join(
            master_df.lazy(),
            [col(sample_id_col)],
            [col("MASTER_SAMPLES")],
            JoinType::Semi.into(),
        )
        .collect()?;
    
    log::info!("After filtering: {} rows remaining", filtered_data.height());

    // 3. Extract phenotype column
    log::info!("Extracting phenotype column '{}'", trait_name);
    
    // The joined dataframe will have both FAM_SAMPLES and the original sample_id_col
    // We'll use the phenotype from the joined data
    let y_vec: Vec<f64> = filtered_data
        .column(trait_name)?
        .f64()?
        .into_iter()
        .collect::<Option<Vec<f64>>>()
        .ok_or(IoError::Alignment(format!("Phenotype column '{}' contains nulls or missing values", trait_name)))?;
    let y: Array1<f64> = Array1::from_vec(y_vec);
    
    log::info!("Phenotype vector has {} values", y.len());

    // 4. Extract sample IDs from filtered data
    let sample_ids: Vec<String> = filtered_data
        .column(sample_id_col)?
        .str()?
        .into_iter()
        .map(|opt_s| opt_s.map(String::from))
        .collect::<Option<Vec<String>>>()
        .ok_or(IoError::Alignment("Failed to extract sample IDs from filtered data".into()))?;

    // 5. Extract covariate columns and create design matrix with intercept
    let n_samples = filtered_data.height();
    let n_covars = covariate_cols.len();
    
    log::info!("Building design matrix: {} rows x {} covariates (+ intercept)", n_samples, n_covars);
    
    let mut x = Array2::zeros((n_samples, n_covars + 1));
    
    // Add intercept as first column
    x.column_mut(0).fill(1.0);
    
    // Add each covariate column
    for (i, covar_name) in covariate_cols.iter().enumerate() {
        log::debug!("Extracting covariate column '{}' (column {})", covar_name, i + 1);
        
        let covar_vec: Vec<f64> = filtered_data
            .column(covar_name)?
            .f64()?
            .into_iter()
            .collect::<Option<Vec<f64>>>()
            .ok_or(IoError::Alignment(format!("Covariate column '{}' contains nulls", covar_name)))?;
        
        if covar_vec.len() != n_samples {
            return Err(IoError::Alignment(format!(
                "Covariate '{}' has {} values but expected {} samples",
                covar_name, covar_vec.len(), n_samples
            )));
        }
        
        for (j, &val) in covar_vec.iter().enumerate() {
            x[[j, i + 1]] = val;
        }
    }
    
    log::info!("Design matrix created: shape [{}, {}]", x.nrows(), x.ncols());

    if y.len() != n_samples || x.nrows() != n_samples {
        return Err(IoError::Alignment(format!(
            "Dimension mismatch after alignment: y ({}), X ({}), samples ({})",
            y.len(), x.nrows(), n_samples
        )));
    }
    
    log::info!("Data alignment completed successfully: {} samples with {} covariates", n_samples, n_covars);

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