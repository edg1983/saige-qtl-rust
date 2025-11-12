# saige-qtl-rust

A high-performance Rust implementation of SAIGE-QTL for single-cell eQTL analysis. This reimplementation provides improved performance and memory safety while maintaining compatibility with the original SAIGE-QTL workflow.

## Features

- **Fast Mixed Model Fitting**: AI-REML algorithm for variance component estimation
- **Parallel Association Testing**: Multi-threaded score tests with Saddlepoint Approximation (SPA)
- **Memory Efficient**: Optimized data structures and lazy evaluation
- **Type Safe**: Leverages Rust's type system to prevent common errors

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
- `step1-fit-null` - Fits the null model
- `step2-run-tests` - Runs association tests

## Usage

### Step 1: Fit Null Model

Fits a null Generalized Linear Mixed Model (GLMM) to estimate variance components and fixed effects.

```bash
target/release/step1-fit-null \
  --plink-file data/genotypes \
  --pheno-file data/phenotypes.tsv \
  --trait-name GENE1 \
  --covar-file data/covariates.tsv \
  --sample-id-in-pheno IID \
  --sample-id-in-covar IID \
  --trait-type quantitative \
  --output-prefix output/step1 \
  --n-threads 8 \
  --tau-init 0.1,0.5 \
  --max-iter 100 \
  --eps 1e-4
```

#### Required Arguments

| Argument               | Description                                                    |
| ---------------------- | -------------------------------------------------------------- |
| `--plink-file`         | Path to PLINK file prefix (.bed/.bim/.fam) for GRM calculation |
| `--pheno-file`         | Path to phenotype file (TSV format)                            |
| `--trait-name`         | Column name in phenotype file for the trait/gene               |
| `--covar-file`         | Path to covariate file (TSV format)                            |
| `--sample-id-in-pheno` | Column name for sample IDs in phenotype file                   |
| `--sample-id-in-covar` | Column name for sample IDs in covariate file                   |
| `--trait-type`         | Type of trait: `quantitative`, `binary`, or `count`            |

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

| Argument               | Description                                   |
| ---------------------- | --------------------------------------------- |
| `--vcf-file`           | Path to VCF/BCF file with genotypes (indexed) |
| `--pheno-file`         | Path to phenotype file (for sample alignment) |
| `--sample-id-in-pheno` | Column name for sample IDs in phenotype file  |
| `--trait-name`         | Trait/gene name (must match Step 1)           |
| `--covar-file`         | Path to covariate file (for sample alignment) |
| `--sample-id-in-covar` | Column name for sample IDs in covariate file  |
| `--step1-output-file`  | Path to fitted null model from Step 1         |

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
# Step 1: Fit null model for a gene
target/release/step1-fit-null \
  --plink-file data/genotypes \
  --pheno-file data/expression.tsv \
  --trait-name ENSG00000000003 \
  --covar-file data/covariates.tsv \
  --sample-id-in-pheno individual_id \
  --sample-id-in-covar individual_id \
  --trait-type quantitative \
  --output-prefix results/null_model \
  --n-threads 16

# Step 2: Test variants in cis-region
target/release/step2-run-tests \
  --vcf-file data/chr1.vcf.gz \
  --vcf-field DS \
  --pheno-file data/expression.tsv \
  --sample-id-in-pheno individual_id \
  --trait-name ENSG00000000003 \
  --covar-file data/covariates.tsv \
  --sample-id-in-covar individual_id \
  --step1-output-file results/null_model.ENSG00000000003.model.bin \
  --region chr1:50000000-51000000 \
  --output-file results/ENSG00000000003.associations.txt.gz \
  --min-mac 5 \
  --n-threads 16
```

## Input File Formats

### Phenotype File
Tab-delimited file with header:
```
IID    GENE1    GENE2    GENE3
sample1    0.5      1.2      -0.3
sample2    1.1      0.8       0.9
```

### Covariate File
Tab-delimited file with header:
```
IID    PC1    PC2    AGE    SEX
sample1    0.1    -0.5    45    1
sample2    -0.3    0.2    38    0
```

### PLINK Files
Standard PLINK binary format (.bed/.bim/.fam) for genetic relationship matrix calculation.

### VCF/BCF File
Indexed VCF or BCF file with dosage (DS) or genotype (GT) information.

## Performance Tips

- Use `--n-threads` to leverage multiple CPU cores
- For large datasets, test regions in parallel using `--region`
- Pre-filter VCF files to relevant variants before testing
- Use BCF format for faster I/O compared to VCF

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

