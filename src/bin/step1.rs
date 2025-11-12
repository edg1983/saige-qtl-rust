//! Step 1: Fit Null GLMM
//!
//! This binary is the Rust equivalent of `extdata/step1_fitNULLGLMM_qtl.R`.
//! It reads phenotype, covariate, and genotype data to fit a
//! null Generalized Linear Mixed Model (GLMM) and estimates
//! the variance components.
//!
//! The output is a serialized `NullModelFit` struct.

use clap::Parser;
use saige_qtl_rust::{
    io::{get_fam_samples, load_aligned_data},
    grm::build_grm_from_plink,
    null_model::fit_null_glmm,
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
    /// Path to the PLINK file prefix (.bed/.bim/.fam) for GRM calculation
    #[arg(long, required = true)]
    plink_file: PathBuf,

    /// Path to the phenotype file
    #[arg(long, required = true)]
    pheno_file: PathBuf,

    /// Column name in the phenotype file for the trait (gene)
    #[arg(long, required = true)]
    trait_name: String,

    /// Path to the covariate file
    #[arg(long, required = true)]
    covar_file: PathBuf,

    /// Column name in the phenotype file for sample IDs
    #[arg(long, required = true)]
    sample_id_in_pheno: String,

    /// Column name in the covariate file for sample IDs
    #[arg(long, required = true)]
    sample_id_in_covar: String,

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

    // Set the global thread pool for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.n_threads)
        .build_global()?;

    // ===================================================================
    // 1. Get sample list from PLINK .fam file
    // ===================================================================
    log::info!("Reading sample list from .fam file...");
    let master_sample_ids = get_fam_samples(&cli.plink_file)?;
    log::info!("Found {} samples in .fam file.", master_sample_ids.len());

    // ===================================================================
    // 2. Load Data (Pheno, Covar) and align to .fam samples
    // ===================================================================
    log::info!("Loading and aligning phenotype and covariate data...");
    let aligned_data = load_aligned_data(
        &cli.pheno_file,
        &cli.covar_file,
        &cli.trait_name,
        &cli.sample_id_in_pheno,
        &cli.sample_id_in_covar,
        &master_sample_ids,
    )?;
    log::info!("{} samples aligned between .fam, pheno, and covar.", aligned_data.sample_ids.len());

    // ===================================================================
    // 3. Build GRM
    // ===================================================================
    // NOTE: This assumes the .bed file is ordered identically to the .fam file,
    // which bed-rs guarantees.
    log::info!("Building GRM from PLINK file: {:?}", cli.plink_file);
    let grm = build_grm_from_plink(&cli.plink_file, cli.n_threads)?;
    
    // We must subset the GRM to match the *aligned* samples
    // This is complex. For now, we assume `load_aligned_data`
    // *only* returns samples present in the .fam, and `build_grm_from_plink`
    // builds the GRM for *all* samples.
    // TODO: Add subsetting of GRM based on `aligned_data.sample_ids`
    log::info!("(Stub) Assuming GRM samples match aligned data.");

    // ===================================================================
    // 4. Fit Null GLMM
    // ===================================================================
    log::info!("Fitting GLMM using AI-REML...");
    let mut null_model_fit = fit_null_glmm(
        &aligned_data,
        &grm,
        &cli.trait_type,
        cli.tau_init,
        cli.max_iter,
        cli.eps,
    )?;

    // Update gene name from stub
    null_model_fit.gene_name = cli.trait_name.clone();

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
