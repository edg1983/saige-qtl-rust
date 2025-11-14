//! Module for handling I/O, data loading, and sample alignment.
use polars::prelude::*;
use std::path::Path;
use ndarray::{Array1, Array2};
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
    /// Aligned sample IDs, in order (CELL IDs for single-cell data)
    pub sample_ids: Vec<String>,
    /// Aligned donor IDs, in order (for VCF matching in single-cell data)
    /// For bulk data, this is identical to sample_ids
    pub donor_ids: Vec<String>,
}

/// Loads a single file containing both phenotypes and covariates, then subsets
/// to samples present in the master sample list (e.g., from a PLINK .fam file).
pub fn load_aligned_data(
    pheno_covar_file: &Path,
    trait_name: &str,
    sample_id_col: &str,      // Donor/sample ID column (matches GRM)
    cell_id_col: Option<&str>, // Optional cell ID column (for single-cell data)
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
    
    // Inner join: keep only rows from data_df where sample_id is in master_df
    // This will keep all rows from data_df that have a matching sample ID in master_df
    let joined_data = data_df
        .lazy()
        .join(
            master_df.lazy(),
            [col(sample_id_col)],
            [col("MASTER_SAMPLES")],
            JoinType::Inner.into(),
        )
        .collect()?;
    
    // Drop the MASTER_SAMPLES column if it exists
    let filtered_data = if joined_data.get_column_names().contains(&"MASTER_SAMPLES") {
        joined_data.drop("MASTER_SAMPLES")?
    } else {
        joined_data
    };
    
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

    // 4. Extract donor IDs (from sample_id_col - these match the GRM)
    let donor_ids: Vec<String> = filtered_data
        .column(sample_id_col)?
        .str()?
        .into_iter()
        .map(|opt_s| opt_s.map(String::from))
        .collect::<Option<Vec<String>>>()
        .ok_or(IoError::Alignment("Failed to extract donor IDs from filtered data".into()))?;

    // 5. Extract cell IDs if specified (for single-cell data), otherwise use donor IDs
    let sample_ids = if let Some(cell_col) = cell_id_col {
        log::info!("Extracting cell IDs from column '{}' (single-cell mode)", cell_col);
        
        // Cast cell ID column to string
        let cell_id_series = filtered_data.column(cell_col)?.clone();
        let cell_id_as_str = cell_id_series.cast(&DataType::String)?;
        
        cell_id_as_str
            .str()?
            .into_iter()
            .map(|opt_s| opt_s.map(String::from))
            .collect::<Option<Vec<String>>>()
            .ok_or(IoError::Alignment(format!("Cell ID column '{}' contains nulls", cell_col)))?
    } else {
        // For bulk data, sample IDs = donor IDs
        log::info!("No cell ID column specified - using donor IDs as sample IDs (bulk data mode)");
        donor_ids.clone()
    };

    // 6. Extract covariate columns and create design matrix with intercept
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
    
    // Validate that donor_ids and sample_ids have the same length
    if donor_ids.len() != n_samples || sample_ids.len() != n_samples {
        return Err(IoError::Alignment(format!(
            "Dimension mismatch: donor_ids ({}), sample_ids ({}), n_samples ({})",
            donor_ids.len(), sample_ids.len(), n_samples
        )));
    }
    
    log::info!("Data alignment completed successfully: {} samples ({} unique donors) with {} covariates", 
              n_samples, donor_ids.iter().collect::<std::collections::HashSet<_>>().len(), n_covars);

    Ok(AlignedData { y, x, sample_ids, donor_ids })
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

/// Gets unique sample IDs from a phenotype/covariate file
pub fn get_unique_sample_ids(
    pheno_covar_file: &Path,
    sample_id_col: &str,
) -> Result<Vec<String>, IoError> {
    log::info!("Extracting unique sample IDs from {:?}", pheno_covar_file);
    
    let mut data_df = CsvReadOptions::default()
        .with_has_header(true)
        .with_parse_options(
            CsvParseOptions::default()
                .with_separator(b'\t')
        )
        .try_into_reader_with_file_path(Some(pheno_covar_file.into()))?
        .finish()?;
    
    // Cast sample ID column to string
    let sample_id_series = data_df.column(sample_id_col)?.clone();
    let sample_id_as_str = sample_id_series.cast(&DataType::String)?;
    data_df.replace(sample_id_col, sample_id_as_str)?;
    
    // Get unique sample IDs
    let unique_ids = data_df
        .column(sample_id_col)?
        .str()?
        .unique()?
        .into_iter()
        .map(|opt_s| opt_s.map(String::from))
        .collect::<Option<Vec<String>>>()
        .ok_or(IoError::Alignment("Failed to extract unique sample IDs".into()))?;
    
    log::info!("Found {} unique samples in data file", unique_ids.len());
    
    Ok(unique_ids)
}

/// Loads sample-level covariates from a separate file (one row per sample/donor)
/// Returns a matrix with one row per sample, aligned to the provided sample_ids
pub fn load_sample_covariates(
    sample_covar_file: &Path,
    sample_id_col: &str,
    sample_covariate_cols: &[String],
    sample_ids: &[String], // Master list of samples to align to
) -> Result<Array2<f64>, IoError> {
    log::info!("Loading sample-level covariates from {:?}", sample_covar_file);
    log::info!("Expected {} samples", sample_ids.len());
    
    // Load the sample covariate file
    let mut data_df = CsvReadOptions::default()
        .with_has_header(true)
        .with_parse_options(
            CsvParseOptions::default()
                .with_separator(b'\t')
        )
        .try_into_reader_with_file_path(Some(sample_covar_file.into()))?
        .finish()?;
    
    log::info!("Loaded sample covariate file with {} rows", data_df.height());
    
    // Cast sample ID column to string
    let sample_id_series = data_df.column(sample_id_col)?.clone();
    let sample_id_as_str = sample_id_series.cast(&DataType::String)?;
    data_df.replace(sample_id_col, sample_id_as_str)?;
    
    // Check that all required samples are present
    let available_samples: std::collections::HashSet<String> = data_df
        .column(sample_id_col)?
        .str()?
        .into_iter()
        .filter_map(|opt_s| opt_s.map(String::from))
        .collect();
    
    let missing_samples: Vec<&String> = sample_ids
        .iter()
        .filter(|id| !available_samples.contains(*id))
        .collect();
    
    if !missing_samples.is_empty() {
        return Err(IoError::Alignment(format!(
            "Sample covariate file is missing {} samples. First few: {:?}",
            missing_samples.len(),
            missing_samples.iter().take(5).collect::<Vec<_>>()
        )));
    }
    
    log::info!("All {} required samples found in covariate file", sample_ids.len());
    
    let n_samples = sample_ids.len();
    let n_covars = sample_covariate_cols.len();
    
    // Build design matrix with intercept
    let mut x_sample = Array2::zeros((n_samples, n_covars + 1));
    x_sample.column_mut(0).fill(1.0); // Intercept
    
    // For each sample in our master list, extract their covariates
    for (sample_idx, sample_id) in sample_ids.iter().enumerate() {
        // Find this sample in the data_df
        let sample_row = data_df
            .column(sample_id_col)?
            .str()?
            .into_iter()
            .position(|id_opt| id_opt == Some(sample_id.as_str()))
            .ok_or_else(|| IoError::Alignment(format!("Sample {} not found in data", sample_id)))?;
        
        // Extract covariates for this sample
        for (cov_idx, covar_name) in sample_covariate_cols.iter().enumerate() {
            let val = data_df
                .column(covar_name)?
                .f64()?
                .get(sample_row)
                .ok_or_else(|| IoError::Alignment(format!("Missing covariate {} for sample {}", covar_name, sample_id)))?;
            
            x_sample[[sample_idx, cov_idx + 1]] = val;
        }
    }
    
    log::info!("Sample covariate matrix: {} samples x {} covariates (+ intercept)", n_samples, n_covars + 1);
    
    Ok(x_sample)
}

/// Loads donor-level covariates (covariates that vary by donor, not by cell)
/// Returns a matrix with one row per unique donor
pub fn load_donor_covariates(
    pheno_covar_file: &Path,
    sample_id_col: &str,
    donor_covariate_cols: &[String],
    donor_sample_ids: &[String],
) -> Result<Array2<f64>, IoError> {
    log::info!("Loading donor-level covariates from {:?}", pheno_covar_file);
    
    // Load the data
    let mut data_df = CsvReadOptions::default()
        .with_has_header(true)
        .with_parse_options(
            CsvParseOptions::default()
                .with_separator(b'\t')
        )
        .try_into_reader_with_file_path(Some(pheno_covar_file.into()))?
        .finish()?;
    
    // Cast sample ID column to string
    let sample_id_series = data_df.column(sample_id_col)?.clone();
    let sample_id_as_str = sample_id_series.cast(&DataType::String)?;
    data_df.replace(sample_id_col, sample_id_as_str)?;
    
    // Group by sample ID and take the first row for each donor
    // (donor-level covariates should be the same for all cells from a donor)
    let donor_df = data_df
        .lazy()
        .group_by([col(sample_id_col)])
        .agg([
            col(sample_id_col).first().alias(sample_id_col),
        ]
        .into_iter()
        .chain(donor_covariate_cols.iter().map(|c| col(c).first().alias(c)))
        .collect::<Vec<_>>())
        .collect()?;
    
    let n_donors = donor_sample_ids.len();
    let n_covars = donor_covariate_cols.len();
    
    // Build design matrix with intercept
    let mut x_donor = Array2::zeros((n_donors, n_covars + 1));
    x_donor.column_mut(0).fill(1.0); // Intercept
    
    // For each donor in our master list, extract their covariates
    for (donor_idx, donor_id) in donor_sample_ids.iter().enumerate() {
        // Find this donor in the donor_df
        let donor_row = donor_df
            .column(sample_id_col)?
            .str()?
            .into_iter()
            .position(|id_opt| id_opt == Some(donor_id.as_str()))
            .ok_or_else(|| IoError::Alignment(format!("Donor {} not found in data", donor_id)))?;
        
        // Extract covariates for this donor
        for (cov_idx, covar_name) in donor_covariate_cols.iter().enumerate() {
            let val = donor_df
                .column(covar_name)?
                .f64()?
                .get(donor_row)
                .ok_or_else(|| IoError::Alignment(format!("Missing covariate {} for donor {}", covar_name, donor_id)))?;
            
            x_donor[[donor_idx, cov_idx + 1]] = val;
        }
    }
    
    log::info!("Donor covariate matrix: {} donors x {} covariates (+ intercept)", n_donors, n_covars + 1);
    
    Ok(x_donor)
}