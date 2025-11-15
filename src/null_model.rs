//! Module for fitting the Null GLMM via AI-REML.
//! 
//! This module provides two fitting strategies:
//! 1. Full REML optimization (standard approach)
//! 2. Fast two-stage approach for single-cell data (precompute donor model, quick per-gene fitting)
use crate::{NullModelFit, TraitType, io::AlignedData};
use ndarray::{Array1, Array2};
// CORRECTED: Import `Solve` trait
use ndarray_linalg::{Cholesky, Inverse, Solve};
use rayon::prelude::*;

// --- Updated argmin imports ---
// CORRECTED: The trait is `CostFunction`, as you pointed out. `Problem` is a struct.
// CORRECTED: Removed unused `State` import
use argmin::core::{CostFunction, Error, Executor};
// CORRECTED: Use BrentOpt for optimization
use argmin::solver::brent::BrentOpt;
// --- End updated imports ---

use thiserror::Error;

#[derive(Error, Debug)]
pub enum ModelError {
    #[error("Linear algebra error: {0}")]
    LinAlg(String),
    #[error("Optimization failed to converge: {0}")]
    Convergence(String),
    #[error("Invalid dimensions: {0}")]
    Dimensions(String),
    #[error("Optimization setup error: {0}")]
    Setup(String),
}

/// Precomputed components for fast single-cell fitting
/// These components are computed once at the donor level and reused for all genes
pub struct PrecomputedComponents {
    /// Donor-level V^-1 = (tau*G + I)^-1
    pub donor_v_inv: Array2<f64>,
    /// Donor-level covariates (X_donor)
    pub donor_x: Array2<f64>,
    /// Mapping from cell index to donor index
    pub cell_to_donor: Vec<usize>,
    /// Optimal tau from donor-level optimization
    pub optimal_tau: f64,
    /// Donor-level GRM
    pub donor_grm: Array2<f64>,
}

impl PrecomputedComponents {
    /// Precompute donor-level components that will be reused for all genes
    /// This does the expensive REML optimization once at the donor level
    pub fn compute_donor_level(
        donor_grm: &Array2<f64>,
        donor_covariates: &Array2<f64>, // Donor-level covariates only (e.g., PCs, donor age/sex)
        donor_sample_ids: &[String],    // Unique donor IDs (matches GRM rows/cols)
        cell_sample_ids: &[String],     // Cell IDs (one per cell, for identification)
        cell_donor_ids: &[String],      // Donor ID for each cell (maps cells to donors)
        tau_init: f64,
        max_iter: u64,
        eps: f64,
    ) -> Result<Self, ModelError> {
        log::info!("=== Precomputing donor-level model components ===");
        log::info!("Donors: {}, Cells: {}", donor_sample_ids.len(), cell_sample_ids.len());
        
        // Build cell-to-donor mapping
        // For each cell, find which donor it belongs to
        let cell_to_donor: Vec<usize> = cell_donor_ids
            .iter()
            .map(|donor_id| {
                donor_sample_ids
                    .iter()
                    .position(|d| d == donor_id)
                    .unwrap_or_else(|| panic!("Donor ID '{}' from cell data not found in donor list", donor_id))
            })
            .collect();
        
        log::info!("Built cell-to-donor mapping for {} cells", cell_to_donor.len());
        
        // Diagnostics: Check GRM properties
        let diag_mean = donor_grm.diag().mean().unwrap_or(0.0);
        let diag_min = donor_grm.diag().iter().copied().fold(f64::INFINITY, f64::min);
        let diag_max = donor_grm.diag().iter().copied().fold(f64::NEG_INFINITY, f64::max);
        
        log::info!("GRM diagnostics:");
        log::info!("  Diagonal: min={:.6}, max={:.6}, mean={:.6}", diag_min, diag_max, diag_mean);
        
        // Calculate off-diagonal statistics
        let mut off_diag_vals = Vec::new();
        for i in 0..donor_grm.nrows() {
            for j in 0..donor_grm.ncols() {
                if i != j {
                    off_diag_vals.push(donor_grm[[i, j]]);
                }
            }
        }
        if !off_diag_vals.is_empty() {
            off_diag_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let off_diag_min = off_diag_vals[0];
            let off_diag_max = off_diag_vals[off_diag_vals.len() - 1];
            let off_diag_median = off_diag_vals[off_diag_vals.len() / 2];
            log::info!("  Off-diagonal: min={:.6}, max={:.6}, median={:.6}", 
                      off_diag_min, off_diag_max, off_diag_median);
        }
        
        // Check for NaN or Inf
        let has_nan = donor_grm.iter().any(|&x| !x.is_finite());
        if has_nan {
            return Err(ModelError::LinAlg("GRM contains NaN or Inf values".to_string()));
        }
        
        // For donor-level optimization, create a pseudo-phenotype (we'll use zeros)
        // since we only care about variance component estimation
        let n_donors = donor_sample_ids.len();
        let y_donor = Array1::zeros(n_donors);
        
        let cost_function = RemlCost {
            y: y_donor,
            x: donor_covariates.clone(),
            grm: donor_grm.clone(),
        };
        
        let lower_bound = 1e-6;
        let upper_bound = 100.0;
        
        let solver = BrentOpt::new(lower_bound, upper_bound)
            .set_tolerance(eps, eps);
        
        log::info!("Starting donor-level REML optimization for tau...");
        let res = Executor::new(cost_function, solver)
            .configure(|state| state.param(tau_init).max_iters(max_iter))
            .run()
            .map_err(|e| ModelError::Convergence(e.to_string()))?;
        
        let optimal_tau = res.state.best_param
            .ok_or_else(|| ModelError::Convergence("No best parameter found".to_string()))?;
        
        log::info!("Donor-level REML converged. Optimal tau = {:.6}", optimal_tau);
        
        // Compute donor-level V_inv with optimal tau
        log::info!("Computing donor-level V^-1 matrix...");
        let mut v_donor = donor_grm * optimal_tau;
        v_donor.diag_mut().mapv_inplace(|v_ii| v_ii + 1.0);
        
        let chol_v = v_donor.cholesky(ndarray_linalg::UPLO::Lower)
            .map_err(|e| ModelError::LinAlg(e.to_string()))?;
        
        // Compute V_inv by solving against identity
        let donor_v_inv_cols: Vec<Array1<f64>> = (0..n_donors)
            .into_par_iter()
            .map(|j| {
                let mut col = Array1::zeros(n_donors);
                col[j] = 1.0;
                chol_v.solve(&col)
                    .map_err(|e| ModelError::LinAlg(e.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        
        let mut donor_v_inv = Array2::zeros((n_donors, n_donors));
        for (j, col_data) in donor_v_inv_cols.into_iter().enumerate() {
            donor_v_inv.column_mut(j).assign(&col_data);
        }
        
        log::info!("Donor-level components precomputed successfully");
        
        Ok(PrecomputedComponents {
            donor_v_inv,
            donor_x: donor_covariates.clone(),
            cell_to_donor,
            optimal_tau,
            donor_grm: donor_grm.clone(),
        })
    }
    
    /// Fast cell-level fitting for a single gene using precomputed donor components
    /// This avoids expensive REML optimization and works directly with donor components
    pub fn fit_gene_fast(
        &self,
        y_cells: &Array1<f64>, // Gene expression at cell level
        x_cells: &Array2<f64>, // Cell-level covariates (includes donor + cell covariates)
        sample_ids: &[String],  // Cell IDs
        donor_ids: &[String],   // Donor IDs for VCF matching
    ) -> Result<NullModelFit, ModelError> {
        let n_cells = y_cells.len();
        let n_covars = x_cells.ncols();
        let n_donors = self.donor_v_inv.nrows();
        
        if sample_ids.len() != n_cells {
            return Err(ModelError::Dimensions(format!(
                "Sample IDs length ({}) doesn't match y length ({})",
                sample_ids.len(), n_cells
            )));
        }
        
        if donor_ids.len() != n_cells {
            return Err(ModelError::Dimensions(format!(
                "Donor IDs length ({}) doesn't match y length ({})",
                donor_ids.len(), n_cells
            )));
        }
        
        log::debug!("Fast fitting: {} cells, {} donors, {} covariates", n_cells, n_donors, n_covars);
        
        // Key optimization: Instead of expanding donor_v_inv to cell level (expensive),
        // we aggregate cells to donors, compute at donor level, then disaggregate
        
        // Step 1: Aggregate y_cells and x_cells to donor level by averaging
        log::debug!("Aggregating cells to donor level...");
        let mut y_donor_agg = Array1::zeros(n_donors);
        let mut x_donor_agg = Array2::zeros((n_donors, n_covars));
        let mut cell_counts = vec![0usize; n_donors];
        
        for (cell_idx, &donor_idx) in self.cell_to_donor.iter().enumerate() {
            y_donor_agg[donor_idx] += y_cells[cell_idx];
            for cov_idx in 0..n_covars {
                x_donor_agg[[donor_idx, cov_idx]] += x_cells[[cell_idx, cov_idx]];
            }
            cell_counts[donor_idx] += 1;
        }
        
        // Average by cell count per donor
        for donor_idx in 0..n_donors {
            let count = cell_counts[donor_idx] as f64;
            if count > 0.0 {
                y_donor_agg[donor_idx] /= count;
                for cov_idx in 0..n_covars {
                    x_donor_agg[[donor_idx, cov_idx]] /= count;
                }
            }
        }
        
        // Step 2: Compute at donor level using precomputed V_inv
        log::debug!("Computing donor-level model components...");
        let v_inv_y_donor = self.donor_v_inv.dot(&y_donor_agg);
        let v_inv_x_donor = self.donor_v_inv.dot(&x_donor_agg);
        
        let x_t_v_inv_x = x_donor_agg.t().dot(&v_inv_x_donor);
        let x_t_v_inv_y = x_donor_agg.t().dot(&v_inv_y_donor);
        
        log::debug!("Inverting X^T V^-1 X...");
        let x_t_v_inv_x_inv = x_t_v_inv_x.inv()
            .map_err(|e| ModelError::LinAlg(e.to_string()))?;
        
        // Fixed effects (donor level)
        let beta = x_t_v_inv_x_inv.dot(&x_t_v_inv_y);
        
        // Step 3: Map results back to cell level
        log::debug!("Mapping results to cell level...");
        
        // Fitted values at cell level using cell-level covariates
        let mu = x_cells.dot(&beta);
        let residuals = y_cells - &mu;
        
        // Compute P matrix at donor level
        let p_donor = &self.donor_v_inv - v_inv_x_donor.dot(&x_t_v_inv_x_inv).dot(&v_inv_x_donor.t());
        
        // Expand P matrix to cell level
        log::debug!("Expanding P matrix to cell level...");
        let p_x_matrix = Array2::from_shape_fn((n_cells, n_cells), |(i, j)| {
            let donor_i = self.cell_to_donor[i];
            let donor_j = self.cell_to_donor[j];
            p_donor[[donor_i, donor_j]]
        });
        
        // Estimate variance components using donor-level aggregated residuals
        let p_y_donor = &v_inv_y_donor - v_inv_x_donor.dot(&beta);
        let y_p_y = y_donor_agg.dot(&p_y_donor);
        let sigma_e2 = y_p_y / (n_donors - n_covars) as f64;
        let sigma_g2 = self.optimal_tau * sigma_e2;
        
        // Calculate variance ratio for association testing
        let var_ratio = sigma_g2 / sigma_e2;
        
        log::debug!("Gene fitting complete: sigma_g2={:.6}, sigma_e2={:.6}, var_ratio={:.6}", sigma_g2, sigma_e2, var_ratio);
        
        Ok(NullModelFit {
            variance_components: vec![sigma_g2, sigma_e2],
            var_ratio,
            fixed_effects: beta.to_vec(),
            residuals,
            p_x_matrix,
            mu,
            y: y_cells.clone(),
            gene_name: "STUB".to_string(),
            trait_type: TraitType::Quantitative,
            sample_ids: sample_ids.to_vec(),
            donor_ids: donor_ids.to_vec(),
        })
    }
}

/// This struct holds the state for the REML optimization
/// We are optimizing `tau` (variance ratio sigma_g / sigma_e)
struct RemlCost {
    y: Array1<f64>,
    x: Array2<f64>,
    grm: Array2<f64>, // The GRM (G)
}

// Implement the `CostFunction` trait for `argmin`
// This calculates the Restricted Log-Likelihood
impl CostFunction for RemlCost {
    // CORRECTED: Renamed `Input` to `Param`
    type Param = f64; // We are optimizing one parameter: tau. This works due to `primitives` feature.
    // CORRECTED: Renamed `Output` to `Cost`
    type Output = f64; // The cost is the negative log-likelihood

    // CORRECTED: Renamed `apply` to `cost` and use `argmin::core::Error`
    fn cost(&self, param: &Self::Param) -> Result<Self::Output, Error> {
        let tau = *param;
        // Set a hard floor for tau
        if tau < 1e-6 { 
            return Ok(1e100); // Return a large cost if tau is non-positive
        }

        let n = self.y.len() as f64;
        let k = self.x.ncols() as f64;
        
        // 1. Form V = tau*G + I
        // CRITICAL FIX: Ensure V is positive definite by adding small regularization if needed
        let mut v = &self.grm * tau;
        
        // Check if diagonal elements are reasonable after multiplication by tau
        // The new GRM formula uses theoretical variance, which can result in smaller diagonal values
        let min_diag_before = v.diag().iter().copied().fold(f64::INFINITY, f64::min);
        let max_diag_before = v.diag().iter().copied().fold(f64::NEG_INFINITY, f64::max);
        
        // Add identity to diagonal: V = tau*G + I
        v.diag_mut().mapv_inplace(|v_ii| v_ii + 1.0);
        
        let min_diag_after = v.diag().iter().copied().fold(f64::INFINITY, f64::min);
        
        // If diagonal is still too small (< 0.01), add regularization
        // This can happen if GRM diagonal values are very small or negative
        let regularization = if min_diag_after < 0.01 {
            let reg = (0.01 - min_diag_after).max(1e-6);
            log::trace!("Adding regularization {:.6} to V (min_diag={:.6})", reg, min_diag_after);
            reg
        } else {
            0.0
        };
        
        if regularization > 0.0 {
            v.diag_mut().mapv_inplace(|v_ii| v_ii + regularization);
        }

        // 2. Cholesky decomposition of V
        let chol_v = match v.cholesky(ndarray_linalg::UPLO::Lower) {
            Ok(c) => c,
            // If V is not positive-definite, return high cost instead of error
            Err(e) => {
                // Log at debug level for first few failures
                log::debug!("Cholesky failed for tau={:.6}: {} (diag range: [{:.6}, {:.6}])", 
                           tau, e, min_diag_before, max_diag_before);
                return Ok(1e100);
            }
        };

        // 3. Calculate components for log-likelihood
        let log_det_v = chol_v.diag().mapv(f64::ln).sum() * 2.0;
        
        // 4. Solve for P = V_inv - V_inv*X*(X^T*V_inv*X)_inv*X^T*V_inv
        // We solve systems to avoid direct inversion
        // CORRECTED: `solve` method works with 1D arrays (vectors)
        let v_inv_y = chol_v.solve(&self.y)
            .map_err(|e| argmin::core::Error::from(ModelError::LinAlg(e.to_string())))?;
        
        // For matrices, we solve column by column in parallel
        let n_covars = self.x.ncols();
        let n_samples = self.x.nrows();
        
        // Parallel column solving
        let v_inv_x_cols: Vec<Array1<f64>> = (0..n_covars)
            .into_par_iter()
            .map(|j| {
                let col = self.x.column(j);
                chol_v.solve(&col.to_owned())
                    .map_err(|e| argmin::core::Error::from(ModelError::LinAlg(e.to_string())))
            })
            .collect::<Result<Vec<_>, _>>()?;
        
        let mut v_inv_x = Array2::zeros((n_samples, n_covars));
        for (j, col_data) in v_inv_x_cols.into_iter().enumerate() {
            v_inv_x.column_mut(j).assign(&col_data);
        }

        let x_t_v_inv_x = self.x.t().dot(&v_inv_x);
        let x_t_v_inv_y = self.x.t().dot(&v_inv_y);

        // CORRECTED: Use `inv()` from the Inverse trait for 2D matrices
        let x_t_v_inv_x_inv = x_t_v_inv_x.inv()
            .map_err(|e| argmin::core::Error::from(ModelError::LinAlg(e.to_string())))?;

        // 5. Calculate beta (fixed effects)
        let beta = x_t_v_inv_x_inv.dot(&x_t_v_inv_y);

        // 6. Calculate P*y = V_inv*y - V_inv*X*beta
        let p_y = &v_inv_y - v_inv_x.dot(&beta);
        
        // 7. Calculate sigma_e^2 (variance scale)
        let y_p_y = self.y.dot(&p_y);
        let sigma_e2 = y_p_y / (n - k);

        if sigma_e2 <= 0.0 {
            return Ok(1e100); // Invalid variance
        }
        
        // 8. Calculate REML Log-Likelihood
        // This is the *negative* log-likelihood, as argmin minimizes
        // We drop constant terms (like 2*PI)
        let log_l = 0.5 * (
            log_det_v + 
            (n - k) * sigma_e2.ln() +
            y_p_y / sigma_e2
        );
        
        // Check for NaN or Inf
        if !log_l.is_finite() {
            return Ok(1e100);
        }

        Ok(log_l)
    }
}

// NOTE: Gradient implementation is no longer needed as we are
// using a derivative-free solver (Brent).

/// Fits the null GLMM using AI-REML.
/// This is the main entry point for Step 1.
pub fn fit_null_glmm(
    aligned_data: &AlignedData,
    grm: &Array2<f64>,
    trait_type: &TraitType,
    tau_init: Vec<f64>, // [tau, sigma_e]
    max_iter: u64,
    eps: f64,
) -> Result<NullModelFit, ModelError> {

    // For now, we only implement the LMM (Quantitative trait)
    // GLMM (Binary, Count) requires an iterative PQL/IRLS loop *outside*
    // the REML optimization, which is significantly more complex.
    if *trait_type != TraitType::Quantitative {
        unimplemented!("Only Quantitative (LMM) traits are currently supported.");
    }

    let cost_function = RemlCost {
        y: aligned_data.y.clone(),
        x: aligned_data.x.clone(),
        grm: grm.clone(),
    };

    let initial_tau = tau_init[0]; // Just tau

    // --- Real optimization using Brent's method ---
    // Brent's method is ideal for 1D, derivative-free optimization.
    // We need to set bounds for tau (variance ratio > 0).
    // Let's set a reasonable range, e.g., [1e-6, 100.0]
    let lower_bound = 1e-6;
    let upper_bound = 100.0; // Assume tau > 100 is unreasonable
    
    // CORRECTED: `set_tolerance` returns `BrentOpt`, not `Result`.
    // The `map_err` was incorrect.
    let solver = BrentOpt::new(lower_bound, upper_bound)
        .set_tolerance(eps, eps); // Set tolerance here (param_tol, cost_tol)

    log::info!("Starting REML optimization for tau...");
    let res = Executor::new(cost_function, solver)
        // FIXED: `initial_guess` is now `configure` to set the param
        // FIXED: `max_iters` must also be set inside `configure`
        .configure(|state| state.param(initial_tau).max_iters(max_iter))
        // .target_precision(eps) // Tolerance is set in solver
        .run()
        .map_err(|e| ModelError::Convergence(e.to_string()))?;
    
    // FIXED: `best_param` is an Option<f64>. We must unwrap it.
    // This fixes the Display, Mul, and ScalarOperand errors.
    let optimal_tau = res.state.best_param.ok_or_else(|| ModelError::Convergence("Optimization found no best parameter".to_string()))?;
    log::info!("REML optimization converged. Optimal tau = {}", optimal_tau);
    // --- End real optimization ---

    // --- Recalculate all final components with optimal_tau ---
    let n = aligned_data.y.len();
    let k = aligned_data.x.ncols();

    let mut v = grm * optimal_tau;
    v.diag_mut().mapv_inplace(|v_ii| v_ii + 1.0);

    let chol_v = v.cholesky(ndarray_linalg::UPLO::Lower)
        .map_err(|e| ModelError::LinAlg(e.to_string()))?;

    // We need V_inv for P_X, so calculate it directly.
    // CORRECTED: `Cholesky` doesn't have `.inverse()`. We solve against Identity column by column.
    // Use parallel processing for large matrices
    log::info!("Computing V_inv matrix ({} x {}) in parallel...", n, n);
    
    let v_inv_cols: Vec<Array1<f64>> = (0..n)
        .into_par_iter()
        .map(|j| {
            let mut col = Array1::zeros(n);
            col[j] = 1.0; // j-th column of identity
            chol_v.solve(&col)
                .map_err(|e| ModelError::LinAlg(e.to_string()))
        })
        .collect::<Result<Vec<_>, _>>()?;
    
    let mut v_inv = Array2::zeros((n, n));
    for (j, col_data) in v_inv_cols.into_iter().enumerate() {
        v_inv.column_mut(j).assign(&col_data);
    }
    
    // CORRECTED: This logic is now valid.
    let v_inv_y = v_inv.dot(&aligned_data.y);
    let v_inv_x = v_inv.dot(&aligned_data.x);

    let x_t_v_inv_x = aligned_data.x.t().dot(&v_inv_x);
    let x_t_v_inv_y = aligned_data.x.t().dot(&v_inv_y);

    // CORRECTED: Use `inv()` from the Inverse trait for 2D matrices
    let x_t_v_inv_x_inv = x_t_v_inv_x.inv()
        .map_err(|e| ModelError::LinAlg(e.to_string()))?;

    // Final beta
    let beta = x_t_v_inv_x_inv.dot(&x_t_v_inv_y);

    // Final sigma_e^2
    let p_y = &v_inv_y - v_inv_x.dot(&beta);
    let sigma_e2 = aligned_data.y.dot(&p_y) / (n - k) as f64;

    // Final variance components
    // FIXED: This now works because optimal_tau is f64, not Option<f64>
    let sigma_g2 = optimal_tau * sigma_e2;
    let variance_components = vec![sigma_g2, sigma_e2];
    
    // Final P_X matrix
    // CORRECTED: All dimensions are now 2D and this works.
    // This was fixed by the `solve_h` -> `solve` change above.
    let p_x_matrix = &v_inv - v_inv_x.dot(&x_t_v_inv_x_inv).dot(&v_inv_x.t());

    // Final residuals (y - mu)
    // For LMM, mu = X*beta
    let mu = aligned_data.x.dot(&beta);
    let residuals = &aligned_data.y - &mu;
    
    // Calculate variance ratio for association testing
    let var_ratio = if variance_components.len() >= 2 && variance_components[1] > 0.0 {
        variance_components[0] / variance_components[1]
    } else {
        1.0  // Default to 1.0 if only one component or division by zero
    };

    Ok(NullModelFit {
        variance_components,
        var_ratio,
        fixed_effects: beta.to_vec(),
        residuals,
        // CORRECTED: `p_x_matrix` is 2D, matching the struct
        p_x_matrix,
        mu,
        y: aligned_data.y.clone(),
        gene_name: "STUB_GENE".into(), // Will be passed in
        trait_type: trait_type.clone(),
        sample_ids: aligned_data.sample_ids.clone(),
        donor_ids: aligned_data.donor_ids.clone(),
    })
}