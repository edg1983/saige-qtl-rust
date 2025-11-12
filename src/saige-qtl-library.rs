//! # SAIGE-QTL-RUST Crate
//!
//! This library contains the core statistical logic for the `saige-qtl-rust` package.
//! The binaries `step1-fit-null` and `step2-run-tests` call functions from this library.

// Re-export key modules
pub mod grm;
pub mod io;
pub mod null_model;
pub mod association;

/// Represents the fitted null model from Step 1.
/// This struct will be serialized to disk (e.g., using bincode)
/// and read by Step 2.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct NullModelFit {
    /// Estimated variance components (e.g., [sigma_g, sigma_e])
    pub variance_components: Vec<f64>,
    /// Fixed-effect coefficients (beta)
    pub fixed_effects: Vec<f64>,
    /// Residuals (y - X*beta)
    pub residuals: Vec<f64>,
    /// Trait type used for the model
    pub trait_type: TraitType,
    /// Sample IDs used in the fit
    pub sample_ids: Vec<String>,
    // TODO: Add other necessary components like the inverse of the
    // variance-covariance matrix, or its Cholesky decomposition,
    // which are needed for the Step 2 score test.
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
