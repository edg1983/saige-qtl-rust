//! Compute GRM: Pre-compute and save Genetic Relationship Matrix
//!
//! This binary computes the GRM from PLINK files once and saves it
//! to disk for reuse in Step 1 analyses. This avoids re-computing
//! the GRM for each gene/trait.

use clap::Parser;
use saige_qtl_rust::{
    grm::{build_grm_from_plink, save_grm_to_file},
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
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let cli = Cli::parse();

    log::info!("Starting GRM computation");
    log::info!("PLINK file: {:?}", cli.plink_file);
    log::info!("Output file: {:?}", cli.output_file);
    log::info!("Using {} threads", cli.n_threads);

    // Set the global thread pool for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.n_threads)
        .build_global()?;

    // Get sample IDs from FAM file
    log::info!("Reading sample IDs from .fam file...");
    let sample_ids = get_fam_samples(&cli.plink_file)?;
    log::info!("Found {} samples in .fam file", sample_ids.len());

    // Compute GRM
    log::info!("Computing GRM from PLINK file...");
    let grm = build_grm_from_plink(&cli.plink_file, cli.n_threads)?;
    
    log::info!("GRM computed: {} x {} matrix", grm.nrows(), grm.ncols());

    // Save GRM with sample IDs to file
    log::info!("Saving GRM to {:?}", cli.output_file);
    save_grm_to_file(&cli.output_file, &grm, &sample_ids)?;

    log::info!("GRM saved successfully");
    log::info!("File size: {} MB", 
        std::fs::metadata(&cli.output_file)?.len() as f64 / 1_048_576.0);

    Ok(())
}
