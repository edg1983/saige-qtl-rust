//! # SAIGE-QTL-RUST Crate
//!
//! This library contains the core statistical logic for the `saige-qtl-rust` package.
//! The binaries `step1-fit-null` and `step2-run-tests` call functions from this library.

// Re-export key modules
pub mod grm;
pub mod io;
pub mod null_model;
pub mod association;

use ndarray::{Array1, Array2};

/// Represents the fitted null model from Step 1.
/// This struct will be serialized to disk (e.g., using bincode)
/// and read by Step 2.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct NullModelFit {
    /// Estimated variance components (e.g., [sigma_g, sigma_e])
    pub variance_components: Vec<f64>,
    /// Variance ratio: tau_genetic / tau_residual, used for score test variance adjustment
    pub var_ratio: f64,
    /// Fixed-effect coefficients (beta)
    pub fixed_effects: Vec<f64>,
    /// Residuals at convergence (y - mu)
    pub residuals: Array1<f64>,
    /// Projection matrix P_X = V_inv - V_inv*X*(X^T*V_inv*X)_inv*X^T*V_inv
    /// This is a dense n x n matrix, where n = number of samples.
    pub p_x_matrix: Array2<f64>,
    /// Fitted mean values (mu)
    pub mu: Array1<f64>,
    /// Original phenotype vector (y)
    pub y: Array1<f64>,
    /// The gene/trait this model was fit for
    pub gene_name: String,
    /// Trait type used for the model
    pub trait_type: TraitType,
    /// Sample IDs used in the fit (CELL IDs for single-cell data)
    pub sample_ids: Vec<String>,
    /// Donor IDs for single-cell data (for VCF matching)
    /// For bulk data, this will be identical to sample_ids
    /// For single-cell data, this maps cells to donors (length = sample_ids.len())
    pub donor_ids: Vec<String>,
    // Note: The projection matrix P_X and residuals are the core
    // components needed by Step 2. Other components like `mu` and `y`
    // are included to support SPA calculations, matching the R
    // script's functionality.
}

#[derive(serde::Serialize, serde::Deserialize, Debug, Clone, PartialEq)]
pub enum TraitType {
    Quantitative,
    Binary,
    Count,
}

impl std::str::FromStr for TraitType {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "quantitative" | "q" => Ok(TraitType::Quantitative),
            "binary" | "b" => Ok(TraitType::Binary),
            "count" | "c" => Ok(TraitType::Count),
            _ => Err("Unknown trait type. Use 'quantitative', 'binary', or 'count'."),
        }
    }
}
