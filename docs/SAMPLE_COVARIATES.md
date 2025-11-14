# Sample-Level Covariate File

## Overview

For single-cell eQTL analysis, the optimized two-stage fitting requires a **sample-level covariate file** separate from the cell-level phenotype file. This separation allows for efficient donor-level model fitting.

## File Format

### Sample Covariate File

A tab-delimited file with **one row per sample/donor**:

```
IID      PC1       PC2       PC3       donor_age  donor_sex
donor1   0.123     -0.456    0.789     45         1
donor2   -0.234    0.567     -0.123    38         0
donor3   0.345     0.678     0.234     52         1
```

**Requirements:**
- Tab-delimited (TSV)
- Header row with column names
- One row per unique sample/donor
- Sample ID column (default: `IID`, can be changed with `--sample-id-col-covar`)
- All samples from GRM must be present

### Phenotype File

A tab-delimited file with **one row per cell** (for single-cell data):

```
donor_id    GENE1     GENE2     cell_type    batch
donor1      2.3       1.5       CD4          1
donor1      2.1       1.8       CD4          1
donor1      3.2       0.9       CD8          1
donor2      1.9       2.1       CD4          2
donor2      2.5       1.3       CD8          2
```

**Requirements:**
- Tab-delimited (TSV)
- Header row with column names
- One row per cell (multiple rows per donor)
- Sample ID column identifying the donor
- Gene expression columns
- Optional: cell-level covariates (e.g., cell_type, batch)

## Usage

### Single-Cell Mode (New Approach)

```bash
./target/release/step1-fit-null \
  --grm-file data/genotypes.grm.bin \
  --pheno-covar-file data/expression.tsv \
  --sample-covar-file data/donor_covariates.tsv \
  --trait-name GENE1 \
  --sample-id-col donor_id \
  --sample-id-col-covar IID \
  --covariate-cols cell_type,batch \
  --sample-covariate-cols PC1,PC2,PC3,PC4,PC5,donor_age,donor_sex \
  --trait-type quantitative \
  --output-prefix results/null_model \
  --n-threads 16
```

**Key Arguments:**

| Argument                  | Description                               | Example                           |
| ------------------------- | ----------------------------------------- | --------------------------------- |
| `--pheno-covar-file`      | Cell-level phenotype file                 | `expression.tsv`                  |
| `--sample-covar-file`     | Sample-level covariate file               | `donor_covariates.tsv`            |
| `--sample-id-col`         | Sample ID column in phenotype file        | `donor_id`                        |
| `--sample-id-col-covar`   | Sample ID column in sample covariate file | `IID` (default)                   |
| `--covariate-cols`        | Optional cell-level covariates            | `cell_type,batch`                 |
| `--sample-covariate-cols` | Sample-level covariates for donor model   | `PC1,PC2,PC3,donor_age,donor_sex` |

### Standard Mode (One Observation per Sample)

For bulk RNA-seq or when each sample has one observation:

```bash
./target/release/step1-fit-null \
  --grm-file data/genotypes.grm.bin \
  --pheno-covar-file data/expression_with_covariates.tsv \
  --trait-name GENE1 \
  --sample-id-col individual_id \
  --covariate-cols PC1,PC2,PC3,age,sex \
  --trait-type quantitative \
  --output-prefix results/null_model \
  --n-threads 16
```

In standard mode:
- No `--sample-covar-file` needed (all covariates in phenotype file)
- No `--sample-covariate-cols` needed
- Uses full REML optimization (no two-stage approach)

## Two-Stage Fitting Logic

### Stage 1: Donor-Level Model

Uses the **sample covariate file** with GRM:

1. Loads sample-level covariates (one row per donor)
2. Verifies all GRM samples are present
3. Runs REML optimization at donor level
4. Computes and caches V⁻¹ matrix

**This stage uses sample covariates ONLY** (not cell-level covariates).

### Stage 2: Cell-Level Gene Fitting

Uses both sample and cell covariates:

1. Loads gene expression at cell level
2. Aggregates cells to donor level
3. Uses precomputed donor V⁻¹ from Stage 1
4. Incorporates cell-level covariates (if provided)
5. Maps results back to cell level

## Covariate Selection Guide

### Sample-Level Covariates (--sample-covariate-cols)

Include covariates that vary **by donor**:
- ✅ Genotype PCs (PC1, PC2, PC3, ...)
- ✅ Donor demographics (age, sex)
- ✅ Donor-level technical factors (sequencing batch for donor)
- ✅ Donor disease status, treatment

These are used in Stage 1 for genetic variance component estimation.

### Cell-Level Covariates (--covariate-cols)

Include covariates that vary **by cell**:
- ✅ Cell type
- ✅ Cell QC metrics (nUMI, percent_mito)
- ✅ Cell-level batch effects
- ✅ Pseudotime, cell cycle phase

These are incorporated in Stage 2 during gene fitting.

## Data Validation

The program performs these checks:

1. **GRM samples**: All samples in GRM must be in sample covariate file
2. **Phenotype samples**: Sample IDs in phenotype file must be subset of GRM samples
3. **Single-cell detection**: Automatically detects if n_cells > n_donors
4. **Required arguments**: Single-cell mode requires both `--sample-covar-file` and `--sample-covariate-cols`

## Example Workflow

### 1. Prepare Sample Covariate File

Extract unique donors with their covariates:

```bash
# From bulk genotype data
awk 'NR>1 {print $2}' genotypes.fam | sort -u > donor_list.txt

# Merge with covariate data
join donor_list.txt covariates.txt > donor_covariates.tsv
```

### 2. Prepare Phenotype File

Your single-cell expression matrix with cell-level covariates:

```bash
# Format: donor_id, gene1, gene2, ..., cell_type, batch
# Multiple rows per donor (one per cell)
```

### 3. Pre-compute GRM (once)

```bash
./target/release/compute-grm \
  --plink-file genotypes \
  --output-file genotypes.grm.bin \
  --n-threads 16
```

### 4. Fit Null Models (per gene)

```bash
# Set BLAS threads
export OPENBLAS_NUM_THREADS=16

# Fit gene models
for gene in GENE1 GENE2 GENE3; do
  ./target/release/step1-fit-null \
    --grm-file genotypes.grm.bin \
    --pheno-covar-file expression.tsv \
    --sample-covar-file donor_covariates.tsv \
    --trait-name $gene \
    --sample-id-col donor_id \
    --sample-id-col-covar IID \
    --covariate-cols cell_type,batch \
    --sample-covariate-cols PC1,PC2,PC3,PC4,PC5,donor_age,donor_sex \
    --trait-type quantitative \
    --output-prefix results/null_model \
    --n-threads 16
done
```

## Performance Benefits

With sample covariate file separation:

1. **Clear separation**: Sample-level vs cell-level effects
2. **Efficient Stage 1**: Donor model uses only relevant covariates
3. **Reusable**: Sample covariate file used for all genes
4. **Validated**: Ensures all GRM samples have covariate data

## Troubleshooting

### Error: "Single-cell mode requires --sample-covar-file"

**Cause**: Detected multiple cells per donor but no sample covariate file provided.

**Solution**: Provide `--sample-covar-file` and `--sample-covariate-cols`.

### Error: "Sample covariate file is missing X samples"

**Cause**: Some samples in GRM are not in the sample covariate file.

**Solution**: Ensure sample covariate file has one row for each sample in the GRM.

### Warning: Not using two-stage optimization

**Cause**: Number of cells equals number of donors (not single-cell data).

**Solution**: This is expected for bulk data. Two-stage optimization only activates when n_cells > n_donors.

## Migration from Previous Version

**Old approach** (no separate files):
```bash
--pheno-covar-file combined.tsv \
--covariate-cols PC1,PC2,PC3,donor_age,cell_type \
--donor-covariate-cols PC1,PC2,PC3,donor_age
```

**New approach** (separate files):
```bash
--pheno-covar-file expression.tsv \
--sample-covar-file donor_covariates.tsv \
--covariate-cols cell_type \
--sample-covariate-cols PC1,PC2,PC3,donor_age
```

**Benefits of new approach:**
- Clearer data organization
- Explicit validation that all GRM samples have covariates
- More efficient loading (sample covariates loaded once)
- Easier to reuse sample covariates across analyses
