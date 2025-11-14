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
/// Uses parallel processing for efficient expansion of large matrices
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
    // Use parallel processing for rows
    let n_cells = cell_donor_ids.len();
    
    // Collect rows in parallel
    let rows: Vec<Vec<f64>> = (0..n_cells)
        .into_par_iter()
        .map(|i| {
            let donor_i = cell_to_donor_idx[i];
            (0..n_cells)
                .map(|j| {
                    let donor_j = cell_to_donor_idx[j];
                    grm[[donor_i, donor_j]]
                })
                .collect()
        })
        .collect();
    
    // Flatten and create array
    let flat_data: Vec<f64> = rows.into_iter().flatten().collect();
    let expanded_grm = Array2::from_shape_vec((n_cells, n_cells), flat_data)
        .map_err(|e| format!("Failed to create expanded GRM: {}", e))?;
    
    log::info!("GRM expansion complete: {} x {} matrix", 
        expanded_grm.nrows(), expanded_grm.ncols());
    
    Ok(expanded_grm)
}

/// Calculates the GRM (A = XX^T / M) from a PLINK .bed file with filtering options.
/// Assumes the .bed file samples are in the *exact* order as `master_sample_ids`.
/// This alignment must be guaranteed by the caller.
/// Note: The global rayon thread pool should be initialized by the caller.
/// 
/// # Arguments
/// * `plink_file_no_ext` - Path to PLINK file prefix (without .bed/.bim/.fam extension)
/// * `n_threads` - Number of threads to use (for logging, actual parallelism uses global rayon pool)
/// * `min_maf` - Minimum MAF threshold (default: 0.001)
/// * `max_missing_rate` - Maximum missing rate threshold (default: 0.05)
/// * `num_random_markers` - If Some(n), randomly select n markers; if None, use all passing filters
/// * `relatedness_cutoff` - If Some(cutoff), set GRM values < cutoff to 0 (sparsification)
pub fn build_grm_from_plink(
    plink_file_no_ext: &Path,
    n_threads: usize,
) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
    build_grm_from_plink_filtered(plink_file_no_ext, n_threads, 0.0, 1.0, None, None)
}

/// Calculates the GRM with full filtering options
pub fn build_grm_from_plink_filtered(
    plink_file_no_ext: &Path,
    n_threads: usize,
    min_maf: f64,
    max_missing_rate: f64,
    num_random_markers: Option<usize>,
    relatedness_cutoff: Option<f64>,
) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
    
    // Note: Thread pool is initialized globally by the caller (main)
    // We just use the n_threads parameter for logging purposes
    log::info!("Building GRM using {} threads", n_threads);
    log::info!("  Filtering: minMAF={}, maxMissingRate={}", min_maf, max_missing_rate);
    if let Some(n) = num_random_markers {
        log::info!("  Random marker selection: {} markers", n);
    }
    if let Some(cutoff) = relatedness_cutoff {
        log::info!("  Relatedness cutoff (sparsification): {}", cutoff);
    }

    let bed_path = plink_file_no_ext.with_extension("bed");
    let mut bed = Bed::new(bed_path)?;
    
    let n_samples = bed.iid_count()?;
    let n_variants = bed.sid_count()?;

    // Read all genotypes into memory
    // This is (n_variants, n_samples)
    let genotypes_flat = bed.read::<i8>()?;
    let genotype_matrix = genotypes_flat.into_shape((n_variants, n_samples))?;
    
    log::info!(
        "Read genotype matrix: {} variants x {} samples",
        n_variants, n_samples
    );

    // 1. Calculate MAF and missing rate for each variant, and filter
    log::info!("Filtering variants by MAF and missing rate...");
    let variant_stats: Vec<(f64, f64, f64, f64)> = genotype_matrix
        .axis_iter(Axis(0)) // Iterate over rows (variants)
        .into_par_iter()
        .map(|row| {
            let mut sum = 0.0;
            let mut sum_sq = 0.0;
            let mut count = 0;
            let mut missing_count = 0;
            
            for &g_code in row.iter() {
                if g_code != i8::MIN { // i8::MIN (-128) is missing in bed-reader
                    let g = g_code as f64;
                    sum += g;
                    sum_sq += g * g;
                    count += 1;
                } else {
                    missing_count += 1;
                }
            }

            let total = count + missing_count;
            let missing_rate = if total > 0 { missing_count as f64 / total as f64 } else { 1.0 };
            
            if count == 0 {
                return (0.0, 1.0, 1.0, 0.0); // No data: mean=0, var=1, missing_rate=1, maf=0
            }
            
            let mean = sum / count as f64;
            let var = (sum_sq / count as f64) - (mean * mean);
            
            // Calculate MAF: mean / 2 gives allele frequency, MAF is min(freq, 1-freq)
            let freq = mean / 2.0;
            let maf = freq.min(1.0 - freq);
            
            (mean, var.max(1e-6), missing_rate, maf)
        })
        .collect();

    // Filter variants based on MAF and missing rate
    let mut valid_indices: Vec<usize> = variant_stats
        .iter()
        .enumerate()
        .filter(|(_, &(_, _, missing_rate, maf))| {
            maf >= min_maf && missing_rate <= max_missing_rate
        })
        .map(|(idx, _)| idx)
        .collect();
    
    log::info!("Variants passing filters: {} / {}", valid_indices.len(), n_variants);
    
    if valid_indices.is_empty() {
        return Err("No variants passed MAF and missing rate filters".into());
    }
    
    // 2. Random marker selection if requested
    if let Some(n_random) = num_random_markers {
        if n_random < valid_indices.len() {
            use rand::seq::SliceRandom;
            use rand::SeedableRng;
            let mut rng = rand::rngs::StdRng::seed_from_u64(42); // Fixed seed for reproducibility
            valid_indices.shuffle(&mut rng);
            valid_indices.truncate(n_random);
            log::info!("Selected {} random markers from {} passing filters", n_random, valid_indices.len());
        } else {
            log::info!("Requested {} random markers but only {} available, using all", 
                      n_random, valid_indices.len());
        }
    }
    
    let n_markers_used = valid_indices.len();
    log::info!("Using {} markers for GRM calculation", n_markers_used);

    // 3. Center and standardize genotypes (only for selected markers)
    // We parallelize over variants (rows)
    let genotype_means_vars: Vec<(f64, f64)> = valid_indices
        .par_iter()
        .map(|&idx| variant_stats[idx])
        .map(|(mean, var, _, _)| {
            if var < 1e-6 {
                (mean, 1.0) // Set std_dev to 1.0 to avoid division by zero
            } else {
                (mean, var)
            }
        })
        .collect();

    // Now center and standardize each selected variant
    // We only process the selected variants, not all
    let mut selected_genotypes = Array2::<i8>::zeros((n_markers_used, n_samples));
    
    selected_genotypes
        .axis_iter_mut(Axis(0))
        .into_par_iter()
        .zip(valid_indices.par_iter())
        .zip(genotype_means_vars.par_iter())
        .for_each(|((mut out_row, &var_idx), &(mean, var))| {
            let std_dev = var.sqrt();
            let in_row = genotype_matrix.row(var_idx);
            for (out_g, &in_g) in out_row.iter_mut().zip(in_row.iter()) {
                if in_g != i8::MIN {
                    let g_f = in_g as f64;
                    let standardized = (g_f - mean) / std_dev;
                    // Scale by 32 for i8 storage efficiency
                    *out_g = (standardized * 32.0).round() as i8;
                } else {
                    *out_g = 0; // Impute missing to 0 after standardization (mean = 0)
                }
            }
        });

    // 4. Convert to f64 matrix for GRM calculation
    // Shape: (n_markers_used, n_samples)
    let genotypes_f64: Array2<f64> = selected_genotypes.mapv(|x| x as f64 / 32.0);

    // 5. Compute GRM = (1/M) * X^T * X
    // X is (n_markers_used, n_samples), so X^T is (n_samples, n_markers_used)
    // GRM is (n_samples, n_samples)
    log::info!("Computing GRM ({} x {})...", n_samples, n_samples);
    let xt = genotypes_f64.t(); // (n_samples, n_markers_used)
    let mut grm = xt.dot(&genotypes_f64) / (n_markers_used as f64);

    // 6. Apply relatedness cutoff (sparsification) if requested
    if let Some(cutoff) = relatedness_cutoff {
        log::info!("Applying relatedness cutoff: setting values < {} to 0", cutoff);
        let mut n_zeroed = 0;
        for i in 0..n_samples {
            for j in 0..n_samples {
                if i != j && grm[[i, j]].abs() < cutoff {
                    grm[[i, j]] = 0.0;
                    n_zeroed += 1;
                }
            }
        }
        log::info!("Set {} off-diagonal elements to 0 ({}% of off-diagonal)", 
                  n_zeroed, 100.0 * n_zeroed as f64 / (n_samples * (n_samples - 1)) as f64);
    }

    log::info!("GRM computed with shape {:?}", grm.dim());
    Ok(grm)
}