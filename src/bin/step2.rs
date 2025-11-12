//! Step 2: Run Association Tests
//!
//! This binary is the Rust equivalent of `extdata/step2_tests_qtl.R`.
//! It loads the fitted null model from Step 1, iterates through
//! a VCF file, and performs a score test (with SPA) for each variant.

use clap::Parser;
use saige_qtl_rust::{
    association::run_parallel_tests,
    NullModelFit
};
use std::path::PathBuf;
use std::fs::File;
use std::io::BufReader;
use flate2::write::GzEncoder;
use flate2::Compression;
use csv::WriterBuilder;

#[derive(Parser, Debug)]
#[command(
    name = "step2-run-tests",
    version,
    about = "Runs association tests using a fitted SAIGEQTL null model (Rust port)"
)]
struct Cli {
    /// Path to the VCF/BCF file with genotypes
    #[arg(long, required = true)]
    vcf_file: PathBuf,
    
    /// Path to the VCF/BCF index file (e.g., .tbi or .csi)
    #[arg(long)]
    vcf_file_index: Option<PathBuf>, // Handled by rust-htslib automatically if named conventionally

    /// Field in VCF to use for genotypes (e.g., 'DS' for dosage, 'GT' for genotype)
    #[arg(long, default_value = "DS")]
    vcf_field: String,

    /// Path to the phenotype file (for sample alignment, not strictly needed if model is good)
    #[arg(long, required = true)]
    pheno_file: PathBuf,
    
    /// Column name in the phenotype file for sample IDs
    #[arg(long, required = true)]
    sample_id_in_pheno: String,

    /// Column name in the phenotype file for the trait (gene)
    #[arg(long, required = true)]
    trait_name: String,

    /// Path to the covariate file (for sample alignment)
    #[arg(long, required = true)]
    covar_file: PathBuf,

    /// Column name in the covariate file for sample IDs
    #[arg(long, required = true)]
    sample_id_in_covar: String,

    /// Path to the fitted null model file from Step 1
    #[arg(long, required = true)]
    step1_output_file: PathBuf,

    /// Genomic region to test (e.g., 'chr1:1-1000000')
    #[arg(long)]
    region: Option<String>,

    /// Path for the output association results file (gzipped)
    #[arg(long, default_value = "saige_step2.results.txt.gz")]
    output_file: PathBuf,

    /// Minimum Minor Allele Count (MAC) to test
    #[arg(long, default_value_t = 0.5)]
    min_mac: f64,

    /// Number of lines to buffer from VCF (not used in this parallel model)
    #[arg(long)]
    num_lines: Option<usize>,

    /// Number of threads to use
    #[arg(long, default_value_t = 1)]
    n_threads: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let cli = Cli::parse();

    log::info!("Starting Step 2: Running Association Tests for {}", cli.trait_name);
    log::info!("Using {} threads", cli.n_threads);

    // Set thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.n_threads)
        .build_global()?;

    // ===================================================================
    // 1. Load Null Model
    // ===================================================================
    log::info!("Loading null model from {:?}", &cli.step1_output_file);
    let file = File::open(&cli.step1_output_file)?;
    let reader = BufReader::new(file);
    let null_model: NullModelFit = bincode::deserialize_from(reader)?;
    
    if null_model.gene_name != cli.trait_name {
        log::warn!(
            "Model gene name '{}' does not match CLI trait name '{}'",
            null_model.gene_name, cli.trait_name
        );
    }
    log::info!("Loaded model for trait type: {:?}", null_model.trait_type);

    // ===================================================================
    // 2. Prepare Output File
    // ===================================================================
    let out_file = File::create(&cli.output_file)?;
    let out_writer_gz = GzEncoder::new(out_file, Compression::default());
    let mut csv_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(out_writer_gz);
        
    // Write header
    csv_writer.write_record(&["chr", "pos", "rsid", "gene", "pval", "mac"])?;

    // ===================================================================
    // 3. Run Association Tests
    // ===================================================================
    // The R scripts re-load pheno/covar to align with VCF samples.
    // Our `run_parallel_tests` handles this internally by aligning
    // VCF samples to the null model's stored samples.
    log::info!("Opening VCF file: {:?}", cli.vcf_file);
    
    run_parallel_tests(
        &cli.vcf_file,
        &cli.vcf_field,
        cli.region.as_deref(),
        &null_model,
        &mut csv_writer,
        cli.min_mac,
    )?;

    log::info!("Step 2 for {} completed successfully.", cli.trait_name);
    Ok(())
}
