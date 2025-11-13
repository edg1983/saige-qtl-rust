//! Module for calculating the Genetic Relationship Matrix (GRM)
use bed_reader::Bed;
use ndarray::{Array2, Axis}; // Removed unused 's'
use rayon::prelude::*;
use std::path::Path;

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