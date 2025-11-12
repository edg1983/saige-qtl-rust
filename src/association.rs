//! Module for running association tests (Score test + SPA)
use crate::{NullModelFit, TraitType};
use rust_htslib::{bcf, bcf::Read};
use ndarray::{Array1, Array2, Axis}; // <-- Added Axis
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
    // TODO: Add other fields like BETA, SE
    pval: f64,
    mac: f64,
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
    
    let mut vcf = bcf::IndexedReader::from_path(vcf_file)?;
    let header = vcf.header();
    
    // Align VCF samples to null model samples
    let vcf_samples = header.samples();
    let model_samples = &null_model.sample_ids;
    
    // CORRECTED: `vcf_samples` is Vec<&[u8]>, function expects &[&[u8]]. Add borrow.
    let (vcf_indices, model_indices) = align_samples(&vcf_samples, model_samples);
    
    if model_indices.len() < 10 {
        return Err(format!("Only {} samples overlap between VCF and null model. Stopping.", model_indices.len()).into());
    }
    
    // Subset the null model components to match the *VCF sample order*
    let p_x = null_model.p_x_matrix.select(Axis(0), &model_indices)
                                 .select(Axis(1), &model_indices);
    let res = null_model.residuals.select(Axis(0), &model_indices);
    let mu = null_model.mu.select(Axis(0), &model_indices);
    let y = null_model.y.select(Axis(0), &model_indices);

    // Pre-calculate P_X * res
    let p_x_res = p_x.dot(&res);

    if let Some(r) = region {
        // CORRECTED: `fetch_str` is gone. Parse region and use `fetch`.
        let (rid, start, end) = parse_region(r, header)?;
        vcf.fetch(rid, start, end)?;
    }
    
    // Buffer records for parallel processing
    let records: Vec<_> = vcf.records().collect::<Result<_, _>>()?;

    let results: Vec<AssocResult> = records
        .into_par_iter()
        .filter_map(|mut record| {
            // CORRECTED: `chr` is `&[u8]`, must be converted.
            let chr_bytes = header.rid2name(record.rid().unwrap()).unwrap();
            let chr = String::from_utf8_lossy(chr_bytes).to_string();
            let pos = record.pos();
            // CORRECTED: Use `as_deref` to fix type ambiguity (Errors 8, 9)
            let rsid = record.id()
                .as_deref() // Converts Option<Vec<u8>> -> Option<&[u8]>
                .map(|v_slice| String::from_utf8_lossy(v_slice).to_string())
                .unwrap_or_else(|| ".".to_string());

            // 1. Extract Genotypes
            let g = match extract_genotypes(&mut record, header, vcf_field, &vcf_indices) {
                Ok(g) => g,
                Err(_) => return None, // Skip variant
            };
            
            // 2. Compute MAC
            let mac = g.sum();
            if mac < min_mac {
                return None; // Skip
            }

            // 3. Run Score Test
            let (score, var) = score_test(&g, &p_x, &p_x_res);
            
            if var <= 1e-6 { return None; } // No variance

            // 4. Calculate P-value
            let pval = if null_model.trait_type == TraitType::Quantitative {
                // Standard score test (Chi-sq 1-dof)
                let chi_sq = score * score / var;
                let normal = Normal::new(0.0, 1.0).unwrap();
                // Two-sided p-value
                2.0 * normal.cdf(-chi_sq.sqrt())
            } else {
                // Use Saddlepoint Approximation (SPA)
                run_spa(score, &g, &mu, &y).unwrap_or(1.0)
            };
            
            Some(AssocResult {
                chr, // This is now a String
                pos,
                rsid,
                gene: null_model.gene_name.clone(),
                pval,
                mac,
            })
        })
        .collect();

    // Write results
    for res in results {
        output_writer.serialize(res)?;
    }
    output_writer.flush()?;
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

    for (i, vcf_sample_bytes) in vcf_samples.iter().enumerate() {
        let vcf_sample_str = String::from_utf8_lossy(vcf_sample_bytes);
        if let Some(&model_idx) = model_sample_map.get(vcf_sample_str.as_ref()) {
            vcf_indices.push(i);
            model_indices.push(model_idx);
        }
    }
    (vcf_indices, model_indices)
}

/// Extracts a genotype vector 'g' aligned to the model
fn extract_genotypes(
    record: &mut bcf::Record,
    // CORRECTED: Mark `header` as unused
    _header: &rust_htslib::bcf::header::HeaderView, // <-- Fixed to full, public path
    vcf_field: &str,
    vcf_indices: &[usize], // Indices of VCF samples that overlap with model
) -> Result<Array1<f64>, Box<dyn std::error::Error>> {
    
    let n_samples_out = vcf_indices.len();
    let mut g = Array1::zeros(n_samples_out);

    if vcf_field == "DS" {
        let dosages = record.format(b"DS").float()?;
        for (i, &vcf_idx) in vcf_indices.iter().enumerate() {
            g[i] = dosages[vcf_idx][0] as f64; // DS is usually Float, take first value
        }
    } else if vcf_field == "GT" {
        // CORRECTED: `genotypes()` is a method on `record`, not `format`.
        let genotypes = record.genotypes()?;
        for (i, &vcf_idx) in vcf_indices.iter().enumerate() {
            let gt = genotypes.get(vcf_idx);
            let g_val = (gt[0].index().unwrap_or(0) + gt[1].index().unwrap_or(0)) as f64;
            g[i] = g_val;
        }
    } else {
        return Err(format!("VCF field '{}' not supported.", vcf_field).into());
    }
    
    // Center genotypes
    let mean = g.mean().unwrap_or(0.0);
    g.mapv_inplace(|v| v - mean);
    
    Ok(g)
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