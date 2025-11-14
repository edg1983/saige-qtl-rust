//! Module for fitting the Null GLMM via AI-REML.
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
        let mut v = &self.grm * tau;
        v.diag_mut().mapv_inplace(|v_ii| v_ii + 1.0); // V = tau*G + I*sigma_e^2 (we factor out sigma_e^2)

        // 2. Cholesky decomposition of V
        let chol_v = match v.cholesky(ndarray_linalg::UPLO::Lower) {
            Ok(c) => c,
            // If V is not positive-definite, optimization is bad
            Err(e) => return Err(Error::from(ModelError::LinAlg(e.to_string()))),
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

    Ok(NullModelFit {
        variance_components,
        fixed_effects: beta.to_vec(),
        residuals,
        // CORRECTED: `p_x_matrix` is 2D, matching the struct
        p_x_matrix,
        mu,
        y: aligned_data.y.clone(),
        gene_name: "STUB_GENE".into(), // Will be passed in
        trait_type: trait_type.clone(),
        sample_ids: aligned_data.sample_ids.clone(),
    })
}