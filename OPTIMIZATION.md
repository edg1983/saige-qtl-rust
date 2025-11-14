# Performance Optimization Guide

## Two-Stage Fitting for Single-Cell eQTL

### Problem

The original implementation was slow for single-cell data because:

1. **Large GRM Expansion**: For 26,593 cells from 1,611 donors, we expanded the donor-level GRM (1,611×1,611) to cell-level (26,593×26,593), creating a matrix with ~700M elements
2. **Repeated Optimization**: REML optimization with Cholesky decomposition and matrix inversions was performed for each gene
3. **Redundant Computation**: The genetic variance structure (GRM) doesn't change across genes, yet we recomputed expensive operations

### Solution: Two-Stage Approach

The key insight is that for single-cell eQTL analysis:
- **Genetic effects** operate at the donor level (shared across all cells from a donor)
- **Gene expression** varies per gene but the covariate structure remains the same

#### Stage 1: Donor-Level Precomputation (Once per Dataset)

Performed once when fitting the first gene:

```rust
PrecomputedComponents::compute_donor_level(
    donor_grm,              // 1,611 × 1,611 donor GRM
    donor_covariates,       // Donor-level covariates (PCs, donor age/sex)
    donor_sample_ids,       // List of unique donors
    cell_sample_ids,        // All cell IDs (with duplicates for multiple cells/donor)
    tau_init,
    max_iter,
    eps,
)
```

This stage:
1. Runs REML optimization **once** at donor level to find optimal τ (variance ratio)
2. Computes and caches **V⁻¹** at donor level using the optimal τ
3. Builds a mapping from cell indices to donor indices

**Time complexity**: O(n_donors³) for Cholesky decomposition - done once
**Memory**: Stores donor-level V⁻¹ (1,611×1,611 = 2.6M elements)

#### Stage 2: Fast Per-Gene Fitting (Seconds per Gene)

For each gene, we use the precomputed components:

```rust
precomputed.fit_gene_fast(
    y_cells,           // Gene expression at cell level
    x_cells,           // Cell-level covariates
    sample_ids,        // Cell sample IDs
)
```

This stage:
1. **Aggregates** cell-level data to donor level (fast averaging)
2. **Reuses** precomputed donor V⁻¹ for all calculations (no Cholesky!)
3. **Computes** fixed effects at donor level
4. **Maps** results back to cell level for fitted values and residuals
5. **Expands** only the final P matrix (needed for Step 2) to cell level

**Time complexity**: 
- Aggregation: O(n_cells)
- Donor-level computation: O(n_donors² × n_covars)
- P matrix expansion: O(n_cells² × n_donors²) but using simple indexing

### Performance Gains

| Operation         | Old Approach       | New Approach       | Speedup   |
| ----------------- | ------------------ | ------------------ | --------- |
| REML optimization | Per gene           | Once (donor level) | 100-1000× |
| GRM expansion     | 26,593² matrix     | Not needed         | ∞         |
| Cholesky decomp   | Per gene (26,593³) | Once (1,611³)      | ~15,000×  |
| Per-gene time     | Minutes            | Seconds            | ~100×     |

### Usage

#### Standard Mode (One Sample per Observation)

```bash
./target/release/step1-fit-null \
  --grm-file data/genotypes.grm.bin \
  --pheno-covar-file data/expression.tsv \
  --trait-name GENE1 \
  --sample-id-col individual_id \
  --covariate-cols PC1,PC2,PC3,age,sex \
  --trait-type quantitative
```

Uses traditional full REML optimization (no optimization needed).

#### Single-Cell Mode (Multiple Cells per Donor)

```bash
./target/release/step1-fit-null \
  --grm-file data/genotypes.grm.bin \
  --pheno-covar-file data/expression_with_covariates.tsv \
  --trait-name GENE1 \
  --sample-id-col donor_id \
  --covariate-cols PC1,PC2,PC3,donor_age,donor_sex,cell_type,batch \
  --donor-covariate-cols PC1,PC2,PC3,donor_age,donor_sex \
  --trait-type quantitative
```

Key differences:
- `sample-id-col` identifies the donor (multiple rows with same ID)
- `donor-covariate-cols` specifies which covariates vary by donor (used for Stage 1)
- Automatically detects single-cell data (n_cells > n_donors) and uses two-stage approach

**Note**: If `--donor-covariate-cols` is not specified, all covariates are assumed to be donor-level.

### BLAS Threading

For maximum performance, set BLAS threads to match your CPU cores:

```bash
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16  # If using Intel MKL

./target/release/step1-fit-null --n-threads 16 [options...]
```

Or use the provided helper script:

```bash
./run_optimized.sh 16 target/release/step1-fit-null [options...]
```

### Implementation Details

#### Cell-to-Donor Aggregation

For gene expression `y` at cell level:
- Cells from the same donor are averaged: `y_donor[d] = mean(y_cell[c] for c where donor[c] = d)`
- Same aggregation for covariates

#### Variance Component Estimation

Uses donor-level aggregated data:
```
σ²_e = y_donor^T P y_donor / (n_donors - n_covars)
σ²_g = τ * σ²_e
```

Where τ is the optimal variance ratio from Stage 1.

#### P Matrix Expansion

Only the P matrix needs cell-level expansion (for Step 2 association testing):
```
P_cell[i, j] = P_donor[donor_of_cell_i, donor_of_cell_j]
```

This is fast because it's just indexing, not matrix operations.

### Mathematical Justification

For single-cell data, the variance structure is:
```
V = τ * G ⊗ I_cells + I_total
```

Where:
- `G` is the donor-level GRM
- `⊗` is the Kronecker product
- `I_cells` handles multiple cells per donor

Key property: `V⁻¹` can be efficiently computed at donor level and then mapped to cells through aggregation, avoiding explicit expansion of the full cell-level matrix.

### Limitations

1. Assumes genetic effects are constant across cells from the same donor (reasonable for eQTL)
2. Cell-level covariates (e.g., cell type) are included but don't interact with genetic effects
3. P matrix expansion still creates large matrices - memory scales as O(n_cells²)

For datasets with >100k cells, consider:
- Processing genes in batches
- Using sparse P matrix representation (future optimization)
- Running on high-memory nodes

## Multi-Threading Configuration

### Rayon (Rust Parallelism)

Set via `--n-threads` argument:
- Controls parallel processing of GRM operations
- Parallel column-wise matrix solves
- Parallel cell aggregation

### BLAS (Linear Algebra)

Set via environment variables:
- `OPENBLAS_NUM_THREADS`: For OpenBLAS
- `MKL_NUM_THREADS`: For Intel MKL  
- `BLIS_NUM_THREADS`: For BLIS

**Important**: BLAS threads and Rayon threads should match:
```bash
export OPENBLAS_NUM_THREADS=16
./target/release/step1-fit-null --n-threads 16 [options...]
```

### Optimal Thread Count

- Use physical cores (not hyperthreads) for compute-bound tasks
- For AMD EPYC: cores within same NUMA node for best memory bandwidth
- Monitor with `htop` or `top` to verify all threads are active

## Profiling

To profile and verify optimizations:

```bash
# Install perf
sudo apt install linux-tools-generic

# Record performance
perf record -g ./target/release/step1-fit-null [options...]

# View report
perf report

# Or use flamegraph
cargo install flamegraph
cargo flamegraph --bin step1-fit-null -- [options...]
```

Look for:
- REML optimization should only appear once in single-cell mode
- Most time in Stage 2 should be in matrix multiplications (BLAS)
- Cell aggregation should be fast (< 1% of time)
