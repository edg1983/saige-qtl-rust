# saige-qtl-rust

A high-performance Rust implementation of SAIGE-QTL for single-cell eQTL analysis. This reimplementation provides improved performance and memory safety while maintaining compatibility with the original SAIGE-QTL workflow.

## Features

- **Fast Mixed Model Fitting**: AI-REML algorithm for variance component estimation
- **Parallel Association Testing**: Multi-threaded score tests with Saddlepoint Approximation (SPA)
- **Memory Efficient**: Optimized data structures and lazy evaluation
- **Type Safe**: Leverages Rust's type system to prevent common errors
- **Simplified Data Format**: Single file for phenotypes and covariates for faster loading

## Installation

### System Dependencies

Load required system libraries before compilation:

```bash
module load llvm-clang/20.0.0git openblas/0.3.24
```

### Compilation

Build the release binaries using Cargo:

```bash
cargo build --release
```

Executables will be available in `target/release/`:
- `compute-grm` - Pre-computes the GRM (one-time setup)
- `step1-fit-null` - Fits the null model
- `step2-run-tests` - Runs association tests

### Docker Installation

For easier deployment across different systems, use the provided Docker container:

#### Build the Docker Image

```bash
docker build -t saige-qtl-rust:latest .
```

#### Pull from Registry (if available)

```bash
docker pull ghcr.io/edg1983/saige-qtl-rust:latest
```

#### Quick Test

```bash
# Test that the container works
docker run --rm saige-qtl-rust:latest compute-grm --help
docker run --rm saige-qtl-rust:latest step1-fit-null --help
docker run --rm saige-qtl-rust:latest step2-run-tests --help
```

## Usage

### Compute GRM (One-time Setup)

Pre-compute the Genetic Relationship Matrix (GRM) once for your dataset. This significantly speeds up Step 1 when analyzing multiple traits/genes.

```bash
target/release/compute-grm \
  --plink-file data/genotypes \
  --output-file data/genotypes.grm.bin \
  --n-threads 16
```

#### Arguments

| Argument        | Description                                |
| --------------- | ------------------------------------------ |
| `--plink-file`  | Path to PLINK file prefix (.bed/.bim/.fam) |
| `--output-file` | Output path for the GRM (binary format)    |
| `--n-threads`   | Number of threads to use (default: 1)      |

**Note**: The GRM file can be reused for all traits analyzed on the same set of samples.

### Step 1: Fit Null Model

Fits a null Generalized Linear Mixed Model (GLMM) to estimate variance components and fixed effects.

**With pre-computed GRM (recommended):**

```bash
target/release/step1-fit-null \
  --grm-file data/genotypes.grm.bin \
  --pheno-covar-file data/phenotypes_and_covariates.tsv \
  --trait-name GENE1 \
  --sample-id-col IID \
  --covariate-cols PC1,PC2,PC3,age,sex \
  --trait-type quantitative \
  --output-prefix output/step1 \
  --n-threads 8
```

**Without pre-computed GRM (slower):**

```bash
target/release/step1-fit-null \
  --plink-file data/genotypes \
  --pheno-covar-file data/phenotypes_and_covariates.tsv \
  --trait-name GENE1 \
  --sample-id-col IID \
  --covariate-cols PC1,PC2,PC3,age,sex \
  --trait-type quantitative \
  --output-prefix output/step1 \
  --n-threads 8
```
  --covariate-cols PC1,PC2,PC3,age,sex \
  --trait-type quantitative \
  --output-prefix output/step1 \
  --n-threads 8 \
  --tau-init 0.1,0.5 \
  --max-iter 100 \
  --eps 1e-4
```

#### Required Arguments

| Argument             | Description                                                  |
| -------------------- | ------------------------------------------------------------ |
| `--pheno-covar-file` | Path to file containing both phenotypes and covariates (TSV) |
| `--trait-name`       | Column name for the trait/gene to analyze                    |
| `--sample-id-col`    | Column name for sample IDs                                   |
| `--covariate-cols`   | Comma-separated list of covariate column names               |
| `--trait-type`       | Type of trait: `quantitative`, `binary`, or `count`          |

**Note**: Either `--grm-file` OR `--plink-file` must be provided:
- `--grm-file`: Path to pre-computed GRM (faster, recommended for multiple traits)
- `--plink-file`: Path to PLINK files (will compute GRM on-the-fly, slower)

#### Optional Arguments

| Argument          | Default       | Description                                         |
| ----------------- | ------------- | --------------------------------------------------- |
| `--output-prefix` | `saige_step1` | Prefix for output files                             |
| `--n-threads`     | `1`           | Number of threads for parallel processing           |
| `--tau-init`      | `0.1,0.5`     | Initial values for tau (variance ratio) and sigma_e |
| `--max-iter`      | `100`         | Maximum iterations for REML optimization            |
| `--eps`           | `1e-4`        | Convergence tolerance for REML                      |

#### Output

Generates a binary file: `{output_prefix}.{trait_name}.model.bin` containing the fitted null model.

### Step 2: Run Association Tests

Performs association testing for variants in a VCF/BCF file using the fitted null model.

```bash
target/release/step2-run-tests \
  --vcf-file data/genotypes.vcf.gz \
  --vcf-field DS \
  --pheno-file data/phenotypes.tsv \
  --sample-id-in-pheno IID \
  --trait-name GENE1 \
  --covar-file data/covariates.tsv \
  --sample-id-in-covar IID \
  --step1-output-file output/step1.GENE1.model.bin \
  --region chr1:1000000-2000000 \
  --output-file output/step2.GENE1.results.txt.gz \
  --min-mac 5 \
  --n-threads 8
```

#### Required Arguments

| Argument               | Description                                             |
| ---------------------- | ------------------------------------------------------- |
| `--vcf-file`           | Path to VCF/BCF file with genotypes (indexed)           |
| `--pheno-file`         | Path to phenotype/covariate file (for sample alignment) |
| `--sample-id-in-pheno` | Column name for sample IDs in phenotype file            |
| `--trait-name`         | Trait/gene name (must match Step 1)                     |
| `--covar-file`         | Path to covariate file (for sample alignment)           |
| `--sample-id-in-covar` | Column name for sample IDs in covariate file            |
| `--step1-output-file`  | Path to fitted null model from Step 1                   |

**Note**: For Step 2, `--pheno-file` and `--covar-file` should point to the same combined file used in Step 1, with the same `--sample-id-in-pheno` and `--sample-id-in-covar` both set to your sample ID column name.

#### Optional Arguments

| Argument        | Default                      | Description                                           |
| --------------- | ---------------------------- | ----------------------------------------------------- |
| `--vcf-field`   | `DS`                         | VCF field to use (`DS` for dosage, `GT` for genotype) |
| `--region`      | All variants                 | Genomic region to test (e.g., `chr1:1-1000000`)       |
| `--output-file` | `saige_step2.results.txt.gz` | Path for output file (gzipped)                        |
| `--min-mac`     | `0.5`                        | Minimum Minor Allele Count threshold                  |
| `--n-threads`   | `1`                          | Number of threads for parallel processing             |

#### Output

Generates a tab-delimited gzipped file with columns:
- `chr`: Chromosome
- `pos`: Position
- `rsid`: Variant ID
- `gene`: Gene/trait name
- `beta`: Effect size estimate
- `se`: Standard error of the effect size
- `pval`: Association p-value
- `mac`: Minor Allele Count

## Example Workflow

```bash
# Step 0: Pre-compute GRM (one-time setup)
target/release/compute-grm \
  --plink-file data/genotypes \
  --output-file data/genotypes.grm.bin \
  --n-threads 16

# Step 1: Fit null model for a gene (using pre-computed GRM)
target/release/step1-fit-null \
  --grm-file data/genotypes.grm.bin \
  --pheno-covar-file data/expression_with_covariates.tsv \
  --trait-name ENSG00000000003 \
  --sample-id-col individual_id \
  --covariate-cols PC1,PC2,PC3,PC4,PC5,age,sex,batch \
  --trait-type quantitative \
  --output-prefix results/null_model \
  --n-threads 16

# Step 2: Test variants in cis-region
target/release/step2-run-tests \
  --vcf-file data/chr1.vcf.gz \
  --vcf-field DS \
  --pheno-file data/expression_with_covariates.tsv \
  --sample-id-in-pheno individual_id \
  --trait-name ENSG00000000003 \
  --covar-file data/expression_with_covariates.tsv \
  --sample-id-in-covar individual_id \
  --step1-output-file results/null_model.ENSG00000000003.model.bin \
  --region chr1:50000000-51000000 \
  --output-file results/ENSG00000000003.associations.txt.gz \
  --min-mac 5 \
  --n-threads 16
```

## Docker Usage

### Running with Docker

Mount your data directory and run the analysis inside the container:

```bash
# Step 1: Fit null model
docker run --rm \
  -v /path/to/data:/data:ro \
  -v /path/to/output:/output \
  -u $(id -u):$(id -g) \
  saige-qtl-rust:latest \
  step1-fit-null \
    --plink-file /data/genotypes \
    --pheno-covar-file /data/expression_with_covariates.tsv \
    --trait-name ENSG00000000003 \
    --sample-id-col individual_id \
    --covariate-cols PC1,PC2,PC3,PC4,PC5,age,sex,batch \
    --trait-type quantitative \
    --output-prefix /output/null_model \
    --n-threads 16

# Step 2: Run association tests
docker run --rm \
  -v /path/to/data:/data:ro \
  -v /path/to/output:/output \
  -u $(id -u):$(id -g) \
  saige-qtl-rust:latest \
  step2-run-tests \
    --vcf-file /data/chr1.vcf.gz \
    --vcf-field DS \
    --pheno-file /data/expression_with_covariates.tsv \
    --sample-id-in-pheno individual_id \
    --trait-name ENSG00000000003 \
    --covar-file /data/expression_with_covariates.tsv \
    --sample-id-in-covar individual_id \
    --step1-output-file /output/null_model.ENSG00000000003.model.bin \
    --region chr1:50000000-51000000 \
    --output-file /output/ENSG00000000003.associations.txt.gz \
    --min-mac 5 \
    --n-threads 16
```

### Singularity

For HPC environments that use Singularity:

```bash
# Convert Docker image to Singularity
singularity pull saige-qtl-rust.sif docker://saige-qtl-rust:latest

# Or build directly from Dockerfile
singularity build saige-qtl-rust.sif docker-daemon://saige-qtl-rust:latest

# Run Step 1
singularity exec \
  --bind /path/to/data:/data \
  --bind /path/to/output:/output \
  saige-qtl-rust.sif \
  step1-fit-null \
    --plink-file /data/genotypes \
    --pheno-covar-file /data/expression_with_covariates.tsv \
    --trait-name ENSG00000000003 \
    --sample-id-col individual_id \
    --covariate-cols PC1,PC2,PC3,PC4,PC5,age,sex,batch \
    --trait-type quantitative \
    --output-prefix /output/null_model \
    --n-threads 16

# Run Step 2
singularity exec \
  --bind /path/to/data:/data \
  --bind /path/to/output:/output \
  saige-qtl-rust.sif \
  step2-run-tests \
    --vcf-file /data/chr1.vcf.gz \
    --pheno-file /data/expression_with_covariates.tsv \
    --sample-id-in-pheno individual_id \
    --trait-name ENSG00000000003 \
    --covar-file /data/expression_with_covariates.tsv \
    --sample-id-in-covar individual_id \
    --step1-output-file /output/null_model.ENSG00000000003.model.bin \
    --region chr1:50000000-51000000 \
    --output-file /output/ENSG00000000003.associations.txt.gz \
    [other options...]
```

## Input File Formats

### Phenotype and Covariate File (Combined)
Tab-delimited file with header containing both phenotypes (gene expression) and covariates:
```
IID         GENE1    GENE2    GENE3    PC1    PC2    AGE    SEX
sample1     0.5      1.2      -0.3     0.1    -0.5   45     1
sample2     1.1      0.8       0.9    -0.3     0.2   38     0
```

**Note**: You specify which column is the phenotype using `--trait-name` and which columns are covariates using `--covariate-cols`. The intercept is automatically added.

### PLINK Files
Standard PLINK binary format (.bed/.bim/.fam) for genetic relationship matrix calculation.

### VCF/BCF File
Indexed VCF or BCF file with dosage (DS) or genotype (GT) information.

## Performance Tips

### Multi-threading Configuration

The software uses both Rust parallelism (rayon) and BLAS multi-threading for optimal performance:

```bash
# Set BLAS threads (should match --n-threads)
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16  # If using Intel MKL
export BLIS_NUM_THREADS=16  # If using BLIS

# Run with matching thread count
target/release/step1-fit-null \
  --n-threads 16 \
  [other options...]
```

**Important**: Always set BLAS threading environment variables to match `--n-threads` for maximum efficiency during:
- REML optimization (tau estimation)
- Matrix inversions and Cholesky decompositions
- GRM expansion for single-cell data

### Additional Optimization Tips

- **Parallel Region Testing**: For large datasets, test regions in parallel using `--region` across multiple jobs
- **Pre-filter VCF**: Filter VCF files to relevant variants before testing to reduce I/O
- **BCF Format**: Use BCF format for faster I/O compared to compressed VCF
- **Pre-compute GRM**: Use `compute-grm` once and reuse for all traits to avoid repeated computation
- **Single-cell Data**: For very large cell counts (>50k cells), monitor memory during GRM expansion as it creates NxN matrices

## Logging

Set the `RUST_LOG` environment variable to control logging verbosity:

```bash
# Info level (default)
export RUST_LOG=info

# Debug level (verbose)
export RUST_LOG=debug

# Run with logging
RUST_LOG=info target/release/step1-fit-null [options...]
```

## Current Limitations

- Only quantitative traits (LMM) are fully implemented
- Binary and count traits require PQL/IRLS implementation
- GRM subsetting assumes all samples are present

## License

See LICENSE file for details.

## Citation

If you use this software, please cite the original SAIGE-QTL paper:
- Zhou et al. (2020) - SAIGE-QTL: scalable and accurate inference of expression quantitative trait loci

