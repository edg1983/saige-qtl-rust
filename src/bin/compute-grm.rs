//! Compute GRM: Pre-compute and save Genetic Relationship Matrix
//!
//! This binary computes the GRM from PLINK files once and saves it
//! to disk for reuse in Step 1 analyses. This avoids re-computing
//! the GRM for each gene/trait.

use clap::Parser;
use saige_qtl_rust::{
    grm::{build_grm_from_plink_filtered, save_grm_to_file},
    io::get_fam_samples,
};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "compute-grm",
    version,
    about = "Pre-computes and saves the Genetic Relationship Matrix (GRM) from PLINK files"
)]
struct Cli {
    /// Path to the PLINK file prefix (.bed/.bim/.fam)
    #[arg(long, required = true)]
    plink_file: PathBuf,

    /// Output file path for the GRM (will be saved in binary format)
    #[arg(long, required = true)]
    output_file: PathBuf,

    /// Number of threads to use
    #[arg(long, default_value_t = 1)]
    n_threads: usize,

    /// Minimum minor allele frequency (MAF) for variants to include in GRM
    #[arg(long, default_value_t = 0.001)]
    min_maf_for_grm: f64,

    /// Maximum missing rate for variants to include in GRM
    #[arg(long, default_value_t = 0.05)]
    max_missing_rate_for_grm: f64,

    /// Number of random markers to select for sparse GRM (optional)
    /// If specified, randomly selects this many markers from those passing MAF/missing filters
    #[arg(long)]
    num_random_marker_for_sparse_kin: Option<usize>,

    /// Relatedness cutoff for sparsification (optional)
    /// If specified, sets GRM[i,j] = 0 for all i≠j where |GRM[i,j]| < cutoff
    #[arg(long)]
    relatedness_cutoff: Option<f64>,

    /// Whether to set diagonal elements of GRM to 1.0
    /// Matches R SAIGE-QTL --isDiagofKinSetAsOne parameter
    /// If false (default), diagonal = sum(standardized_genotype^2) / M (computed from data)
    /// If true, diagonal is forced to 1.0 (assumes perfect self-relatedness)
    #[arg(long, default_value_t = false)]
    is_diag_of_kin_set_as_one: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let cli = Cli::parse();

    log::info!("Starting GRM computation");
    log::info!("PLINK file: {:?}", cli.plink_file);
    log::info!("Output file: {:?}", cli.output_file);
    log::info!("Using {} threads", cli.n_threads);
    log::info!("Filtering parameters:");
    log::info!("  minMAFforGRM = {}", cli.min_maf_for_grm);
    log::info!("  maxMissingRateforGRM = {}", cli.max_missing_rate_for_grm);
    if let Some(n) = cli.num_random_marker_for_sparse_kin {
        log::info!("  numRandomMarkerforSparseKin = {}", n);
    }
    if let Some(cutoff) = cli.relatedness_cutoff {
        log::info!("  relatednessCutoff = {}", cutoff);
    }
    log::info!("  isDiagofKinSetAsOne = {}", cli.is_diag_of_kin_set_as_one);

    // Set the global thread pool for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.n_threads)
        .build_global()?;

    // Get sample IDs from FAM file
    log::info!("Reading sample IDs from .fam file...");
    let sample_ids = get_fam_samples(&cli.plink_file)?;
    log::info!("Found {} samples in .fam file", sample_ids.len());

    // Compute GRM with filtering
    log::info!("Computing GRM from PLINK file...");
    let mut grm = build_grm_from_plink_filtered(
        &cli.plink_file,
        cli.n_threads,
        cli.min_maf_for_grm,
        cli.max_missing_rate_for_grm,
        cli.num_random_marker_for_sparse_kin,
        cli.relatedness_cutoff,
    )?;
    
    log::info!("GRM computed: {} x {} matrix", grm.nrows(), grm.ncols());

    // Apply isDiagofKinSetAsOne if requested
    // MATCHES R SAIGE-QTL: get_DiagofKin() returns 1.0 if isDiagofKinSetAsOne=TRUE
    if cli.is_diag_of_kin_set_as_one {
        log::info!("Setting GRM diagonal to 1.0 (isDiagofKinSetAsOne=TRUE)");
        for i in 0..grm.nrows() {
            grm[[i, i]] = 1.0;
        }
        log::info!("Diagonal set to 1.0 for all {} samples", grm.nrows());
    }

    // Save GRM with sample IDs to file
    log::info!("Saving GRM to {:?}", cli.output_file);
    save_grm_to_file(&cli.output_file, &grm, &sample_ids)?;

    log::info!("GRM saved successfully");
    log::info!("File size: {} MB", 
        std::fs::metadata(&cli.output_file)?.len() as f64 / 1_048_576.0);

    Ok(())
}
