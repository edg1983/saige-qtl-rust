//! Step 1: Fit Null GLMM
//!
//! This binary is the Rust equivalent of `extdata/step1_fitNULLGLMM_qtl.R`.
//! It reads a single file containing both phenotypes and covariates along with 
//! genotype data to fit a null Generalized Linear Mixed Model (GLMM) and 
//! estimates the variance components.
//!
//! The output is a serialized `NullModelFit` struct.

use clap::Parser;
use saige_qtl_rust::{
    io::{get_fam_samples, get_unique_sample_ids, load_aligned_data, load_sample_covariates},
    grm::{build_grm_from_plink, load_grm_from_file, subset_grm},
    null_model::{fit_null_glmm, PrecomputedComponents},
    TraitType
};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "step1-fit-null",
    version,
    about = "Fits a null GLMM for SAIGEQTL analysis (Rust port)"
)]
struct Cli {
    /// Path to the PLINK file prefix (.bed/.bim/.fam) for sample list
    /// Required only if --grm-file is not provided
    #[arg(long)]
    plink_file: Option<PathBuf>,

    /// Path to pre-computed GRM file (use compute-grm to generate this)
    /// If provided, --plink-file is not required
    #[arg(long)]
    grm_file: Option<PathBuf>,

    /// Path to the file containing phenotypes (and optionally cell-level covariates)
    #[arg(long, required = true)]
    pheno_covar_file: PathBuf,
    
    /// Path to sample-level covariate file (one row per sample/donor)
    /// This file should contain sample IDs and covariates that vary by sample (e.g., PCs, donor age/sex)
    /// Required for single-cell mode with --sample-covariate-cols
    #[arg(long)]
    sample_covar_file: Option<PathBuf>,

    /// Column name for the trait (gene) to analyze
    #[arg(long, required = true)]
    trait_name: String,

    /// Column name for sample IDs in phenotype file
    #[arg(long, required = true)]
    sample_id_col: String,
    
    /// Column name for sample IDs in sample covariate file
    #[arg(long, default_value = "IID")]
    sample_id_col_covar: String,

    /// Comma-separated list of covariate column names from pheno_covar_file
    /// These are cell-level covariates (e.g., cell_type, batch)
    #[arg(long, value_delimiter = ',')]
    covariate_cols: Option<Vec<String>>,
    
    /// Comma-separated list of sample-level covariate names from sample_covar_file
    /// These are covariates that vary by sample/donor (e.g., PCs, donor_age, donor_sex)
    /// Used for the donor-level model in single-cell optimization
    #[arg(long, value_delimiter = ',')]
    sample_covariate_cols: Option<Vec<String>>,

    /// Type of trait ('quantitative', 'binary', or 'count')
    #[arg(long, required = true)]
    trait_type: TraitType,

    /// Prefix for output files
    #[arg(long, default_value = "saige_step1")]
    output_prefix: String,

    /// Number of threads to use
    #[arg(long, default_value_t = 1)]
    n_threads: usize,

    /// Initial values for tau (variance ratio) and sigma_e, comma-separated (e.g., "0.1,0.5")
    #[arg(long, default_value = "0.1,0.5", value_delimiter = ',')]
    tau_init: Vec<f64>,

    /// Maximum iterations for REML optimization
    #[arg(long, default_value_t = 100)]
    max_iter: u64,

    /// Convergence epsilon for REML optimization
    #[arg(long, default_value_t = 1e-4)]
    eps: f64,
    
    /// Cutoff for variance ratio (tau). Not used in this implementation.
    #[arg(long, default_value = "0.01")]
    variance_ratio_cutoff: Option<String>, 
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init(); // Initialize logger
    let cli = Cli::parse();

    log::info!("Starting Step 1: Fitting Null GLMM for {}", cli.trait_name);
    log::info!("Trait type: {:?}", cli.trait_type);
    log::info!("Using {} threads", cli.n_threads);

    // Validate arguments
    if cli.grm_file.is_none() && cli.plink_file.is_none() {
        return Err("Either --grm-file or --plink-file must be provided".into());
    }

    // Set the global thread pool for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.n_threads)
        .build_global()?;

    // ===================================================================
    // 1. Determine master sample list
    // ===================================================================
    let master_sample_ids = if let Some(grm_file) = &cli.grm_file {
        // Using pre-computed GRM: get unique samples from phenotype file
        log::info!("Using pre-computed GRM mode");
        log::info!("Getting unique sample IDs from phenotype/covariate file...");
        get_unique_sample_ids(&cli.pheno_covar_file, &cli.sample_id_col)?
    } else if let Some(plink_file) = &cli.plink_file {
        // Computing GRM from PLINK: use FAM file as master list
        log::info!("Using PLINK mode (computing GRM on-the-fly)");
        log::info!("Reading sample list from .fam file...");
        let fam_samples = get_fam_samples(plink_file)?;
        log::info!("Found {} samples in .fam file.", fam_samples.len());
        fam_samples
    } else {
        unreachable!()
    };

    log::info!("Master sample list: {} samples", master_sample_ids.len());

    // ===================================================================
    // 2. Load Data (Pheno + optional cell-level covariates)
    // ===================================================================
    log::info!("Loading and aligning phenotype data...");
    let cell_covariate_cols = cli.covariate_cols.as_ref().map(|v| v.as_slice()).unwrap_or(&[]);
    let aligned_data = load_aligned_data(
        &cli.pheno_covar_file,
        &cli.trait_name,
        &cli.sample_id_col,
        cell_covariate_cols,
        &master_sample_ids,
    )?;
    log::info!("{} rows in aligned data.", aligned_data.sample_ids.len());

    // ===================================================================
    // 3. Load or Build GRM, then subset/align to unique donors
    // ===================================================================
    let donor_grm = if let Some(grm_file) = &cli.grm_file {
        // Load pre-computed GRM
        log::info!("Loading pre-computed GRM from {:?}", grm_file);
        let grm_data = load_grm_from_file(grm_file)?;
        
        // Subset GRM to match unique samples from phenotype file
        log::info!("Subsetting GRM to {} unique samples", master_sample_ids.len());
        subset_grm(&grm_data, &master_sample_ids)?
    } else if let Some(plink_file) = &cli.plink_file {
        // Compute GRM from PLINK files
        log::info!("Computing GRM from PLINK file: {:?}", plink_file);
        log::warn!("Computing GRM on-the-fly. For better performance, pre-compute with 'compute-grm' and use --grm-file");
        build_grm_from_plink(plink_file, cli.n_threads)?
    } else {
        unreachable!()
    };
    
    log::info!("Donor-level GRM ready: {} x {} matrix", donor_grm.nrows(), donor_grm.ncols());
    
    // ===================================================================
    // 4. Detect single-cell data and use optimized two-stage approach
    // ===================================================================
    let n_cells = aligned_data.sample_ids.len();
    let n_donors = master_sample_ids.len();
    
    let null_model_fit = if n_cells > n_donors {
        log::info!("=== SINGLE-CELL MODE DETECTED ===");
        log::info!("Cells: {}, Donors: {}", n_cells, n_donors);
        log::info!("Using optimized two-stage fitting approach");
        
        // Validate that sample covariate file and columns are provided
        let sample_covar_file = cli.sample_covar_file.as_ref()
            .ok_or("Single-cell mode requires --sample-covar-file")?;
        let sample_covars = cli.sample_covariate_cols.as_ref()
            .ok_or("Single-cell mode requires --sample-covariate-cols")?;
        
        log::info!("Sample-level covariates: {:?}", sample_covars);
        
        // Load sample-level covariates (one row per donor)
        log::info!("Loading sample-level covariates from {:?}", sample_covar_file);
        let x_sample = load_sample_covariates(
            sample_covar_file,
            &cli.sample_id_col_covar,
            sample_covars,
            &master_sample_ids,
        )?;
        
        log::info!("Loaded sample covariate matrix: {} samples x {} covariates (+ intercept)", 
                   x_sample.nrows(), x_sample.ncols());
        
        // Precompute donor-level components using sample covariates (this is the expensive part, done once)
        log::info!("=== STAGE 1: Precomputing donor-level model with sample covariates ===");
        let precomputed = PrecomputedComponents::compute_donor_level(
            &donor_grm,
            &x_sample,
            &master_sample_ids,
            &aligned_data.sample_ids,
            cli.tau_init[0],
            cli.max_iter,
            cli.eps,
        )?;
        
        log::info!("Donor-level components precomputed successfully");
        log::info!("Optimal tau from donor model: {:.6}", precomputed.optimal_tau);
        
        // Fast per-gene fitting (this should be very quick)
        log::info!("=== STAGE 2: Fitting gene {} at cell level ===", cli.trait_name);
        let mut fit = precomputed.fit_gene_fast(
            &aligned_data.y,
            &aligned_data.x,
            &aligned_data.sample_ids,
        )?;
        
        fit.gene_name = cli.trait_name.clone();
        fit.trait_type = cli.trait_type.clone();
        
        log::info!("Gene fitting completed");
        fit
    } else {
        log::info!("=== STANDARD MODE (one observation per sample) ===");
        log::info!("Using full REML optimization");
        
        let mut fit = fit_null_glmm(
            &aligned_data,
            &donor_grm,
            &cli.trait_type,
            cli.tau_init,
            cli.max_iter,
            cli.eps,
        )?;
        
        fit.gene_name = cli.trait_name.clone();
        fit
    };

    // ===================================================================
    // 5. Save Output
    // ===================================================================
    let output_file = PathBuf::from(format!("{}.{}.model.bin", cli.output_prefix, cli.trait_name));
    log::info!("Saving fitted null model to {:?}", &output_file);

    let file = std::fs::File::create(&output_file)?;
    let writer = std::io::BufWriter::new(file);
    bincode::serialize_into(writer, &null_model_fit)?;

    log::info!("Step 1 for {} completed successfully.", cli.trait_name);
    Ok(())
}