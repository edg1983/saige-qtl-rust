//! Module for running association tests (Score test + SPA)
use crate::{NullModelFit, TraitType};
use rust_htslib::{bcf, bcf::Read};
use ndarray::{Array1, Array2};
use std::path::Path;
// CORRECTED: Import `Continuous` trait for `.pdf()`
use statrs::distribution::{Continuous, ContinuousCDF, Normal};
// REMOVED bad dependency
use rayon::prelude::*;
use csv::Writer;
use std::io::Write;
use flate2::write::GzEncoder;
// use flate2::Compression; // <-- Removed unused import

#[derive(serde::Serialize)]
struct AssocResult {
    chr: String,
    pos: i64,
    rsid: String,
    gene: String,
    beta: f64,
    se: f64,
    pval: f64,
    mac: f64,
    n_samples: usize,
}

/// Helper to parse region strings like "chr1:1-100000"
fn parse_region(
    region_str: &str,
    header: &bcf::header::HeaderView,
) -> Result<(u32, u64, Option<u64>), Box<dyn std::error::Error>> {
    let parts: Vec<&str> = region_str.split(':').collect();
    if parts.len() != 2 {
        return Err("Invalid region format. Use 'chr:start-end' or 'chr:start'".into());
    }
    let rid = header.name2rid(parts[0].as_bytes())?;
    
    let range_parts: Vec<&str> = parts[1].split('-').collect();
    let start = range_parts[0].parse::<u64>()?.saturating_sub(1); // 1-based to 0-based
    let end = if range_parts.len() > 1 {
        Some(range_parts[1].parse::<u64>()?)
    } else {
        None
    };
    Ok((rid, start, end))
}


/// Main entry point for Step 2.
/// Iterates VCF and runs tests in parallel.
pub fn run_parallel_tests(
    vcf_file: &Path,
    vcf_field: &str,
    region: Option<&str>,
    null_model: &NullModelFit,
    output_writer: &mut Writer<GzEncoder<impl Write>>,
    min_mac: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    
    log::info!("Starting association testing");
    log::info!("VCF file: {:?}", vcf_file);
    log::info!("VCF field: {}", vcf_field);
    log::info!("Minimum MAC: {}", min_mac);
    log::info!("Region: {:?}", region);
    
    let mut vcf = bcf::IndexedReader::from_path(vcf_file)?;
    let header = vcf.header().clone(); // Clone the header to avoid borrowing issues
    
    // Align VCF samples (DONOR IDs) to null model donor IDs
    let vcf_samples = header.samples();
    let model_donor_ids = &null_model.donor_ids;  // For VCF matching
    let model_sample_ids = &null_model.sample_ids;  // Cell IDs
    
    log::info!("VCF has {} samples (donors)", vcf_samples.len());
    log::info!("Null model has {} donors", model_donor_ids.len());
    log::info!("Null model has {} cells", model_sample_ids.len());
    
    // Align VCF donors to model donors
    let (vcf_indices, donor_indices) = align_samples(&vcf_samples, model_donor_ids);
    
    log::info!("Sample alignment: {} overlapping donors", donor_indices.len());
    
    if donor_indices.len() < 10 {
        return Err(format!("Only {} donors overlap between VCF and null model. Stopping.", donor_indices.len()).into());
    }
    
    log::debug!("First 5 overlapping donors (VCF indices): {:?}", &vcf_indices.iter().take(5).collect::<Vec<_>>());
    log::debug!("First 5 overlapping donors (Model indices): {:?}", &donor_indices.iter().take(5).collect::<Vec<_>>());
    
    // Build donor-to-cell mapping
    // For each cell, record which donor it belongs to (using donor indices after VCF alignment)
    let mut cell_to_donor_idx: Vec<Option<usize>> = vec![None; model_sample_ids.len()];
    for (cell_idx, donor_id) in model_donor_ids.iter().enumerate() {
        // Find which position this donor_id is in the aligned donor list
        if let Some(pos) = donor_indices.iter().position(|&idx| &model_donor_ids[idx] == donor_id) {
            cell_to_donor_idx[cell_idx] = Some(pos);
        }
    }
    
    // Count how many cells we can analyze (those with aligned donors)
    let n_cells_with_donors = cell_to_donor_idx.iter().filter(|x| x.is_some()).count();
    log::info!("{} cells have donors in the VCF", n_cells_with_donors);
    
    if n_cells_with_donors < 10 {
        return Err(format!("Only {} cells have donors in VCF. Stopping.", n_cells_with_donors).into());
    }
    
    // No subsetting - use full cell-level matrices
    let p_x = &null_model.p_x_matrix;  // n_cells x n_cells
    let res = &null_model.residuals;  // n_cells
    let mu = &null_model.mu;  // n_cells
    let y = &null_model.y;  // n_cells

    log::debug!("Using full cell-level matrices: {} cells", model_sample_ids.len());
    log::debug!("P_X matrix shape: [{}, {}]", p_x.nrows(), p_x.ncols());

    // Pre-calculate P_X * res
    let p_x_res = p_x.dot(res);  // No & needed - both are references
    log::debug!("Pre-calculated P_X * residuals");

    if let Some(r) = region {
        // CORRECTED: `fetch_str` is gone. Parse region and use `fetch`.
        let (rid, start, end) = parse_region(r, &header)?;
        log::info!("Fetching region: rid={}, start={}, end={:?}", rid, start, end);
        vcf.fetch(rid, start, end)?;
    }
    
    // Buffer records for parallel processing and extract metadata that requires header access
    // This avoids capturing the non-Sync header in the parallel closure
    let mut record_data: Vec<_> = Vec::new();
    let mut total_variants: usize = 0;
    
    log::info!("Reading variants from VCF...");
    for record_result in vcf.records() {
        let record = record_result?;
        total_variants += 1;
        
        // Extract information that requires header access
        let chr_bytes = header.rid2name(record.rid().unwrap()).unwrap();
        let chr = String::from_utf8_lossy(chr_bytes).to_string();
        let pos = record.pos();
        let id_bytes = record.id();
        let rsid = String::from_utf8_lossy(&id_bytes).to_string();
        
        record_data.push((record, chr, pos, rsid));
        
        if total_variants % 10000 == 0 {
            log::debug!("Read {} variants so far...", total_variants);
        }
    }
    
    log::info!("Read {} total variants from VCF", total_variants);

    if total_variants == 0 {
        log::warn!("No variants found in VCF file!");
        return Ok(());
    }

    log::info!("Processing variants in parallel...");

    let results: Vec<AssocResult> = record_data
        .into_par_iter()
        .filter_map(|(mut record, chr, pos, rsid)| {

            // 1. Extract Genotypes at DONOR level and calculate MAC
            let (g_donor, mac) = match extract_genotypes(&mut record, vcf_field, &vcf_indices) {
                Ok(result) => result,
                Err(e) => {
                    log::trace!("Failed to extract genotypes for {}:{}:{} - {}", chr, pos, rsid, e);
                    return None; // Skip variant
                }
            };
            
            // 2. Expand donor-level genotypes to cell-level
            // For each cell, copy the genotype from its corresponding donor
            let mut g_cells = Array1::zeros(model_sample_ids.len());
            let mut n_cells_analyzed = 0;
            
            for (cell_idx, donor_idx_opt) in cell_to_donor_idx.iter().enumerate() {
                if let Some(donor_idx) = donor_idx_opt {
                    // This cell has a donor in the VCF
                    g_cells[cell_idx] = g_donor[*donor_idx];
                    n_cells_analyzed += 1;
                } else {
                    // This cell's donor is not in the VCF - set to 0 (will be centered anyway)
                    g_cells[cell_idx] = 0.0;
                }
            }
            
            if n_cells_analyzed == 0 {
                log::trace!("Variant {}:{}:{} has no cells with genotypes", chr, pos, rsid);
                return None;
            }
            
            // 3. Filter by MAC (MAC is already calculated at donor level)
            if mac < min_mac {
                log::trace!("Variant {}:{}:{} filtered by MAC: {:.2} < {}", chr, pos, rsid, mac, min_mac);
                return None; // Skip
            }

            // 4. Center cell-level genotypes for association testing
            let mean_g = g_cells.sum() / n_cells_analyzed as f64;
            g_cells.mapv_inplace(|v| v - mean_g);
            
            // 5. Run Score Test with cell-level genotypes
            let (score, var2) = score_test(&g_cells, &p_x, &p_x_res);
            
            // Apply variance ratio correction (matches R SAIGE: var1 = var2 * varRatio)
            let var = var2 * null_model.var_ratio;
            
            if var <= 1e-6 { 
                log::trace!("Variant {}:{}:{} has insufficient variance: {:.2e}", chr, pos, rsid, var);
                return None; 
            }

            // 6. Calculate effect size (BETA) and standard error (SE)
            // SAIGE formula: beta = Score / var1, se = |beta| / sqrt(|stat|), where stat = Score^2 / var1
            let beta = score / var;
            let se = 1.0 / var.sqrt();

            // 7. Calculate P-value
            let pval = if null_model.trait_type == TraitType::Quantitative {
                // Standard score test (Chi-sq 1-dof)
                let chi_sq = score * score / var;
                let normal = Normal::new(0.0, 1.0).unwrap();
                // Two-sided p-value
                2.0 * normal.cdf(-chi_sq.sqrt())
            } else {
                // Use Saddlepoint Approximation (SPA)
                run_spa(score, &g_cells, &mu, &y).unwrap_or(1.0)
            };
            
            log::trace!("Variant {}:{}:{} passed: MAC={:.2}, beta={:.4}, se={:.4}, pval={:.2e}, var2={:.4}, var_ratio={:.4}, n_cells={}", 
                       chr, pos, rsid, mac, beta, se, pval, var2, null_model.var_ratio, n_cells_analyzed);
            
            Some(AssocResult {
                chr, // This is now a String
                pos,
                rsid,
                gene: null_model.gene_name.clone(),
                beta,
                se,
                pval,
                mac,
                n_samples: n_cells_analyzed,  // Now reports actual cell count
            })
        })
        .collect();

    log::info!("Association testing complete:");
    log::info!("  Total variants processed: {}", total_variants);
    log::info!("  Variants passing filters: {}", results.len());
    log::info!("  Filtered by genotype extraction: estimated ~{}", total_variants.saturating_sub(results.len()));
    
    if results.is_empty() {
        log::warn!("WARNING: No variants passed filters!");
        log::warn!("  This could be due to:");
        log::warn!("  1. All variants have MAC < {} (min_mac)", min_mac);
        log::warn!("  2. Genotype extraction failed (wrong --vcf-field?)", );
        log::warn!("  3. No genetic variance in variants");
        log::warn!("  Enable RUST_LOG=trace to see per-variant filtering reasons");
    } else {
        log::info!("Writing {} results to output file", results.len());
    }

    // Write results
    for res in results {
        output_writer.serialize(res)?;
    }
    output_writer.flush()?;
    
    log::info!("Results written successfully");
    Ok(())
}

/// Aligns VCF samples and Model samples
fn align_samples(vcf_samples: &[&[u8]], model_samples: &[String]) -> (Vec<usize>, Vec<usize>) {
    use std::collections::HashMap;
    let model_sample_map: HashMap<&str, usize> = model_samples
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let mut vcf_indices = Vec::new();
    let mut model_indices = Vec::new();

    log::debug!("Aligning VCF samples to null model samples...");
    
    for (i, vcf_sample_bytes) in vcf_samples.iter().enumerate() {
        let vcf_sample_str = String::from_utf8_lossy(vcf_sample_bytes);
        if let Some(&model_idx) = model_sample_map.get(vcf_sample_str.as_ref()) {
            vcf_indices.push(i);
            model_indices.push(model_idx);
            
            if vcf_indices.len() <= 5 {
                log::debug!("  Sample '{}': VCF index {} -> Model index {}", vcf_sample_str, i, model_idx);
            }
        } else {
            if i < 5 {
                log::trace!("  Sample '{}' (VCF index {}) not found in null model", vcf_sample_str, i);
            }
        }
    }
    
    if vcf_indices.is_empty() {
        log::error!("NO SAMPLES OVERLAP between VCF and null model!");
        log::error!("VCF samples (first 5): {:?}", 
                   vcf_samples.iter().take(5).map(|s| String::from_utf8_lossy(s)).collect::<Vec<_>>());
        log::error!("Model samples (first 5): {:?}", 
                   model_samples.iter().take(5).collect::<Vec<_>>());
    }
    
    (vcf_indices, model_indices)
}

/// Extracts a genotype vector 'g' aligned to the model
/// Returns (centered_genotypes, mac)
fn extract_genotypes(
    record: &mut bcf::Record,
    vcf_field: &str,
    vcf_indices: &[usize], // Indices of VCF samples that overlap with model
) -> Result<(Array1<f64>, f64), Box<dyn std::error::Error>> {
    
    let n_samples_out = vcf_indices.len();
    let mut g = Array1::zeros(n_samples_out);

    if vcf_field == "DS" {
        let dosages = record.format(b"DS").float()?;
        
        // Check if we have dosage data
        if dosages.is_empty() {
            return Err("No DS (dosage) field found in VCF record".into());
        }
        
        for (i, &vcf_idx) in vcf_indices.iter().enumerate() {
            if vcf_idx >= dosages.len() {
                return Err(format!("VCF index {} out of bounds (dosages len: {})", vcf_idx, dosages.len()).into());
            }
            if dosages[vcf_idx].is_empty() {
                return Err(format!("Empty dosage for sample index {}", vcf_idx).into());
            }
            g[i] = dosages[vcf_idx][0] as f64; // DS is usually Float, take first value
        }
    } else if vcf_field == "GT" {
        // CORRECTED: `genotypes()` is a method on `record`, not `format`.
        let genotypes = record.genotypes()?;
        for (i, &vcf_idx) in vcf_indices.iter().enumerate() {
            let gt = genotypes.get(vcf_idx);
            let allele1 = gt[0].index().unwrap_or(0);
            let allele2 = gt[1].index().unwrap_or(0);
            let g_val = (allele1 + allele2) as f64;
            g[i] = g_val;
        }
    } else {
        return Err(format!("VCF field '{}' not supported. Use 'DS' or 'GT'.", vcf_field).into());
    }
    
    // Calculate MAC BEFORE centering
    // MAC = Minor Allele Count
    // For DS: sum of dosages (already represents allele count)
    // For GT: sum of allele counts
    let allele_sum = g.sum();
    let allele_count = n_samples_out as f64 * 2.0; // Diploid
    
    // MAC is the minimum of (alt allele count, ref allele count)
    let mac = allele_sum.min(allele_count - allele_sum);
    
    // Calculate stats for debugging
    let mean_before = g.mean().unwrap_or(0.0);
    
    log::trace!("Genotype extraction (DONOR level): mean={:.4}, allele_sum={:.4}, MAC={:.4}, n_donors={}", 
               mean_before, allele_sum, mac, n_samples_out);
    
    // DO NOT center here - we'll center at cell level after expansion
    
    Ok((g, mac))
}

/// Performs the score test
#[inline]
fn score_test(g: &Array1<f64>, p_x: &Array2<f64>, p_x_res: &Array1<f64>) -> (f64, f64) {
    // Numerator: t(G) %*% P_X %*% res = t(G) %*% (P_X_res)
    let score = g.dot(p_x_res);
    
    // Denominator (Variance): t(G) %*% P_X %*% G
    let p_x_g = p_x.dot(g);
    let var = g.dot(&p_x_g);
    
    (score, var)
}

/// Implements the Saddlepoint Approximation (SPA)
fn run_spa(
    score: f64,
    g: &Array1<f64>, // Genotype vector (centered)
    mu: &Array1<f64>, // Fitted means (n x 1)
    y: &Array1<f64>,  // Phenotypes (n x 1)
) -> Result<f64, Box<dyn std::error::Error>> {
    
    // CGF = Cumulant Generating Function
    // K(t) = sum(log( (1-mu_i)*exp(-t*g_i*y_i) + mu_i*exp(t*g_i*(1-y_i)) ))
    // This is for BINARY traits. Count traits (Poisson/NB) have a different CGF.
    // This highlights the complexity. We'll use the binary CGF.
    
    // K'(t) = score. We need to find t_hat that solves this.
    // We define the functions `f` and `f_prime` for the solver.
    let k_prime = |t: f64, g: &Array1<f64>, mu: &Array1<f64>, y: &Array1<f64>| -> f64 {
        (g * y * (1.0 - mu) * (-t * g * y).mapv(f64::exp)
            + g * (1.0 - y) * mu * (t * g * (1.0 - y)).mapv(f64::exp))
            .sum()
    };
    
    // CORRECTED: `f664` typo
    let k_double_prime = |t: f64, g: &Array1<f64>, mu: &Array1<f64>, y: &Array1<f64>| -> f64 {
        (g.mapv(|v| v*v) * y.mapv(|v| v*v) * (1.0 - mu) * (-t * g * y).mapv(f64::exp)
            + g.mapv(|v| v*v) * (1.0 - y).mapv(|v| v*v) * mu * (t * g * (1.0 - y)).mapv(f64::exp))
            .sum()
    };

    // `f` is the function we want to find the root of (K'(t) - score = 0)
    let f = |t: f64| -> f64 {
        k_prime(t, g, mu, y) - score
    };
    
    // `f_prime` is its derivative (K''(t))
    let f_prime = |t: f64| -> f64 {
        k_double_prime(t, g, mu, y)
    };
    
    // Find t_hat
    let t_hat = match manual_newton_raphson(0.0, f, f_prime, 1e-6, 100) {
        Ok(root) => root,
        Err(_) => return Ok(1.0), // Failed to converge
    };
    
    // Now calculate K(t_hat) and K''(t_hat)
    let k_t_hat = ((1.0 - mu) * (-t_hat * g * y).mapv(f64::exp)
                  + mu * (t_hat * g * (1.0 - y)).mapv(f64::exp))
                  .mapv(f64::ln)
                  .sum();
    let k_double_prime_val = k_double_prime(t_hat, g, mu, y);

    if k_double_prime_val < 1e-6 { return Ok(1.0); } // Numerical issue

    // SPA p-value formula
    let w = (2.0 * (t_hat * score - k_t_hat)).sqrt() * t_hat.signum();
    let v = t_hat * k_double_prime_val.sqrt();
    
    let normal = Normal::new(0.0, 1.0).unwrap();
    // CORRECTED: `pdf` is available now that `Continuous` trait is in scope.
    let pval = normal.cdf(w) + normal.pdf(w) * (1.0/w - 1.0/v);
    
    Ok(pval)
}

/// A simple, self-contained Newton-Raphson solver.
/// Replaces the external `rust-roots` dependency.
fn manual_newton_raphson<F, FPrime>(
    mut x0: f64,
    f: F,
    f_prime: FPrime,
    tolerance: f64,
    max_iter: u32,
) -> Result<f64, &'static str>
where
    F: Fn(f64) -> f64,
    FPrime: Fn(f64) -> f64,
{
    for _ in 0..max_iter {
        let y = f(x0);
        let y_prime = f_prime(x0);

        if y_prime.abs() < 1e-10 {
            // Avoid division by zero
            return Err("Derivative is zero");
        }

        let x1 = x0 - y / y_prime;

        if (x1 - x0).abs() < tolerance {
            return Ok(x1);
        }

        x0 = x1;
    }

    Err("Newton-RDaphson failed to converge")
}