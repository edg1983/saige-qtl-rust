//! Module for calculating the Genetic Relationship Matrix (GRM)
use bed_reader::Bed;
use ndarray::{Array2, Axis};
use rayon::prelude::*;
use std::path::Path;
use std::fs::File;
use std::io::BufReader;
use serde::{Serialize, Deserialize};

/// Structure to hold GRM along with sample IDs
#[derive(Serialize, Deserialize)]
pub struct GrmData {
    /// The genetic relationship matrix
    pub grm: Array2<f64>,
    /// Sample IDs in the order corresponding to GRM rows/columns
    pub sample_ids: Vec<String>,
}

/// Loads a pre-computed GRM from a binary file
pub fn load_grm_from_file(grm_file: &Path) -> Result<GrmData, Box<dyn std::error::Error>> {
    log::info!("Loading pre-computed GRM from {:?}", grm_file);
    
    let file = File::open(grm_file)?;
    let reader = BufReader::new(file);
    let grm_data: GrmData = bincode::deserialize_from(reader)?;
    
    log::info!("GRM loaded: {} x {} matrix with {} sample IDs", 
        grm_data.grm.nrows(), grm_data.grm.ncols(), grm_data.sample_ids.len());
    
    Ok(grm_data)
}

/// Saves GRM and sample IDs to a binary file
pub fn save_grm_to_file(
    grm_file: &Path, 
    grm: &Array2<f64>, 
    sample_ids: &[String]
) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("Saving GRM to {:?}", grm_file);
    
    let grm_data = GrmData {
        grm: grm.clone(),
        sample_ids: sample_ids.to_vec(),
    };
    
    let file = File::create(grm_file)?;
    let writer = std::io::BufWriter::new(file);
    bincode::serialize_into(writer, &grm_data)?;
    
    log::info!("GRM saved successfully");
    
    Ok(())
}

/// Subsets and reorders a GRM to match a specific list of sample IDs
pub fn subset_grm(
    grm_data: &GrmData,
    target_sample_ids: &[String],
) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
    log::info!("Subsetting GRM from {} to {} samples", 
        grm_data.sample_ids.len(), target_sample_ids.len());
    
    // Create a mapping from sample ID to index in the original GRM
    let sample_to_idx: std::collections::HashMap<&str, usize> = grm_data
        .sample_ids
        .iter()
        .enumerate()
        .map(|(idx, id)| (id.as_str(), idx))
        .collect();
    
    // Find indices for target samples
    let mut indices = Vec::with_capacity(target_sample_ids.len());
    for target_id in target_sample_ids {
        match sample_to_idx.get(target_id.as_str()) {
            Some(&idx) => indices.push(idx),
            None => {
                return Err(format!(
                    "Sample '{}' not found in GRM. GRM contains {} samples.",
                    target_id, grm_data.sample_ids.len()
                ).into());
            }
        }
    }
    
    // Subset the GRM matrix
    let n = indices.len();
    let mut subset_grm = Array2::zeros((n, n));
    
    for (i, &idx_i) in indices.iter().enumerate() {
        for (j, &idx_j) in indices.iter().enumerate() {
            subset_grm[[i, j]] = grm_data.grm[[idx_i, idx_j]];
        }
    }
    
    log::info!("GRM subset complete: {} x {} matrix", subset_grm.nrows(), subset_grm.ncols());
    
    Ok(subset_grm)
}

/// Expands a donor-level GRM to cell-level using a mapping matrix
/// For single-cell data where multiple cells belong to the same donor
pub fn expand_grm_to_cell_level(
    grm: &Array2<f64>,
    unique_donor_ids: &[String],
    cell_donor_ids: &[String],
) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
    log::info!("Expanding GRM from {} donors to {} cells", 
        unique_donor_ids.len(), cell_donor_ids.len());
    
    // Create mapping from donor ID to index in the GRM
    let donor_to_idx: std::collections::HashMap<&str, usize> = unique_donor_ids
        .iter()
        .enumerate()
        .map(|(idx, id)| (id.as_str(), idx))
        .collect();
    
    // Create mapping from each cell to its donor's GRM index
    let mut cell_to_donor_idx = Vec::with_capacity(cell_donor_ids.len());
    for cell_donor in cell_donor_ids {
        match donor_to_idx.get(cell_donor.as_str()) {
            Some(&idx) => cell_to_donor_idx.push(idx),
            None => {
                return Err(format!(
                    "Cell belongs to donor '{}' which is not in the GRM",
                    cell_donor
                ).into());
            }
        }
    }
    
    // Expand the GRM: cell_grm[i,j] = donor_grm[donor_of_cell_i, donor_of_cell_j]
    let n_cells = cell_donor_ids.len();
    let mut expanded_grm = Array2::zeros((n_cells, n_cells));
    
    for i in 0..n_cells {
        for j in 0..n_cells {
            let donor_i = cell_to_donor_idx[i];
            let donor_j = cell_to_donor_idx[j];
            expanded_grm[[i, j]] = grm[[donor_i, donor_j]];
        }
    }
    
    log::info!("GRM expansion complete: {} x {} matrix", 
        expanded_grm.nrows(), expanded_grm.ncols());
    
    Ok(expanded_grm)
}

/// Calculates the GRM (A = XX^T / M) from a PLINK .bed file.
/// Assumes the .bed file samples are in the *exact* order as `master_sample_ids`.
/// This alignment must be guaranteed by the caller.
/// Note: The global rayon thread pool should be initialized by the caller.
pub fn build_grm_from_plink(
    plink_file_no_ext: &Path,
    n_threads: usize,
) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
    
    // Note: Thread pool is initialized globally by the caller (main)
    // We just use the n_threads parameter for logging purposes
    log::info!("Building GRM using {} threads", n_threads);

    let bed_path = plink_file_no_ext.with_extension("bed");
    let mut bed = Bed::new(bed_path)?;
    
    let n_samples = bed.iid_count()?;
    let n_variants = bed.sid_count()?;

    // Read all genotypes into memory
    // This is (n_variants, n_samples)
    let genotypes_flat = bed.read()?;
    let mut genotype_matrix = genotypes_flat.into_shape((n_variants, n_samples))?;
    
    log::info!(
        "Read genotype matrix: {} variants x {} samples",
        n_variants, n_samples
    );

    // 1. Center and standardize genotypes
    // We parallelize over variants (rows)
    let genotype_means_vars: Vec<(f64, f64)> = genotype_matrix
        .axis_iter(Axis(0)) // Iterate over rows (variants)
        .into_par_iter()
        .map(|row| {
            let mut sum = 0.0;
            let mut sum_sq = 0.0;
            let mut count = 0;
            
            for &g_code in row.iter() {
                if g_code != i8::MIN { // i8::MIN (-128) is missing in bed-reader
                    let g = g_code as f64;
                    sum += g;
                    sum_sq += g * g;
                    count += 1;
                }
            }

            if count == 0 {
                return (0.0, 1.0); // No data, mean 0, variance 1 (will be zeroed)
            }
            
            let mean = sum / count as f64;
            let var = (sum_sq / count as f64) - (mean * mean);
            // Handle invariant sites
            if var < 1e-6 {
                (mean, 1.0) // Set std_dev to 1.0 to avoid division by zero
            } else {
                (mean, var)
            }
        })
        .collect();

    genotype_matrix.axis_iter_mut(Axis(0)) // Iterate over rows
        .into_par_iter()
        .enumerate()
        .for_each(|(i, mut row)| {
            let (mean, var) = genotype_means_vars[i];
            let std_dev = var.sqrt();
            
            for g_code in row.iter_mut() {
                if *g_code == i8::MIN { // Missing
                    *g_code = 0; // Impute to mean (which is 0 after centering)
                } else {
                    let g_centered = (*g_code as f64 - mean) / std_dev;
                    // Store as i8 for memory efficiency before final cast
                    // This is a slight precision loss, but common
                    *g_code = (g_centered * 32.0).round() as i8; 
                }
            }
        });

    // 2. Compute GRM: A = X * X^T / M
    // X is now the standardized (n_variants, n_samples) matrix
    // We want (n_samples, n_samples)
    // Let's transpose X to (n_samples, n_variants)
    let x_t = genotype_matrix.t().mapv(|v| v as f64 / 32.0); // Back to f64
    
    log::info!("Computing GRM ({} x {})", n_samples, n_samples);
    
    // A = (X_t * X_t^T) / M
    let grm = x_t.dot(&x_t.t()) / (n_variants as f64);

    Ok(grm)
}