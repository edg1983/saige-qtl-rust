# Performance Improvements Summary

## What Was Changed

I've completely redesigned the null model fitting to dramatically speed up single-cell eQTL analysis:

### Key Changes

1. **Two-Stage Optimization Strategy** (`src/null_model.rs`)
   - **Stage 1** (once): Fit donor-level model, optimize τ, compute V⁻¹ at donor level
   - **Stage 2** (per gene): Fast fitting using precomputed donor components
   - Avoids repeated REML optimization and massive GRM expansions

2. **Donor-Level Precomputation** (New: `PrecomputedComponents`)
   - Caches donor-level V⁻¹ matrix (1,611×1,611 instead of 26,593×26,593)
   - Stores optimal τ from donor-level optimization
   - Builds cell-to-donor mapping once

3. **Fast Gene Fitting** (New: `fit_gene_fast()`)
   - Aggregates cell data to donor level (fast averaging)
   - Uses precomputed V⁻¹ for all calculations (no Cholesky decomposition!)
   - Maps results back to cell level
   - Only expands P matrix (needed for Step 2) at the end

4. **Donor Covariate Loading** (`src/io.rs`)
   - New `load_donor_covariates()` function
   - Extracts donor-level covariates for Stage 1 optimization
   - Handles multiple cells per donor gracefully

5. **Updated CLI** (`src/bin/step1.rs`)
   - New `--donor-covariate-cols` argument to specify which covariates are donor-level
   - Automatically detects single-cell mode (n_cells > n_donors)
   - Routes to appropriate fitting strategy

## Expected Performance Gains

| Metric                 | Before            | After               | Improvement            |
| ---------------------- | ----------------- | ------------------- | ---------------------- |
| REML optimization      | Every gene        | Once (donor-level)  | **100-1000× faster**   |
| GRM size               | 26,593² = 700M    | 1,611² = 2.6M       | **270× smaller**       |
| Cholesky decomposition | O(n³) per gene    | O(n³) once          | **Per gene: ∞ faster** |
| Per-gene fitting time  | **Minutes**       | **Seconds**         | **~100× speedup**      |
| Memory usage           | ~5-10 GB per gene | ~50-100 MB per gene | **50-100× less**       |

## How to Test

### 1. Rebuild the Code

```bash
cd /ssu/gassu/GAU_tools/saige-qtl-rust
cargo build --release
```

### 2. Run with Your Single-Cell Data

For single-cell eQTL with multiple cells per donor:

```bash
# Set BLAS threads for maximum performance
export OPENBLAS_NUM_THREADS=16

# Run step1 (automatically uses two-stage optimization)
./target/release/step1-fit-null \
  --grm-file /path/to/genotypes.grm.bin \
  --pheno-covar-file /path/to/expression_with_covariates.tsv \
  --trait-name YOUR_GENE_NAME \
  --sample-id-col donor_id \
  --covariate-cols PC1,PC2,PC3,PC4,PC5,donor_age,donor_sex,cell_type,batch \
  --donor-covariate-cols PC1,PC2,PC3,PC4,PC5,donor_age,donor_sex \
  --trait-type quantitative \
  --output-prefix results/null_model \
  --n-threads 16
```

**Key arguments**:
- `--donor-covariate-cols`: Specify which covariates vary by donor (not by cell)
  - Include: PCs, donor age/sex, donor batch
  - Exclude: cell type, cell-specific QC metrics
- If omitted, all covariates are assumed donor-level

### 3. Monitor Performance

Watch CPU usage while running:
```bash
htop
```

You should now see:
- **Stage 1**: Brief period of high CPU usage (all cores) during donor-level REML optimization
- **Stage 2**: Fast per-gene fitting (should complete in seconds)

### 4. Check Logs

Look for these log messages:

```
=== SINGLE-CELL MODE DETECTED ===
Cells: 26593, Donors: 1611
Using optimized two-stage fitting approach

=== STAGE 1: Precomputing donor-level model ===
Starting donor-level REML optimization for tau...
Donor-level REML converged. Optimal tau = 0.XXXXXX
Donor-level components precomputed successfully

=== STAGE 2: Fitting gene ENSG00000000003 at cell level ===
Aggregating cells to donor level...
Computing donor-level model components...
Mapping results to cell level...
Gene fitting completed
```

### 5. Benchmark Timing

Time a single gene:

```bash
time ./target/release/step1-fit-null [your options] --trait-name TEST_GENE
```

Expected result:
- **Before**: Several minutes per gene
- **After**: A few seconds per gene (after Stage 1 precomputation)

## Architecture Overview

```
Single-Cell Data Flow:

Input:
├── 26,593 cells
├── 1,611 unique donors  
└── Donor-level GRM (1,611×1,611)

Stage 1 (Once):
├── Fit donor-level model
├── REML optimization → optimal τ
├── Compute V⁻¹ at donor level
└── Cache for reuse

Stage 2 (Per Gene):
├── Aggregate cells → donors
├── Use cached V⁻¹
├── Compute β (fixed effects)
├── Map back to cells
└── Expand P matrix

Output:
└── NullModelFit with cell-level results
```

## Troubleshooting

### Build Fails with hts-sys Error

This is a system dependency issue. Load required modules:

```bash
module load llvm-clang/20.0.0git openblas/0.3.24
cargo clean
cargo build --release
```

### Still Slow Performance

Check:
1. **BLAS threads set?**
   ```bash
   echo $OPENBLAS_NUM_THREADS  # Should equal --n-threads
   ```

2. **Using pre-computed GRM?**
   - Use `--grm-file` (fast) not `--plink-file` (slow)

3. **Donor covariates specified?**
   - Must use `--donor-covariate-cols` for single-cell optimization

4. **Log shows single-cell mode?**
   - Check logs for "SINGLE-CELL MODE DETECTED"
   - If not detected, check that your data has n_cells > n_donors

### Memory Issues

For very large datasets (>100k cells):
- The P matrix expansion is O(n_cells²)
- Consider processing on high-memory nodes
- Future: sparse P matrix implementation

## What's Not Changed

- **Step 2** (association testing): Still uses the same algorithm
- **Standard mode** (one obs per sample): Uses original full REML
- **Binary/Count traits**: Still unimplemented (require PQL/IRLS)
- **Output format**: Identical to before (backward compatible)

## Next Steps

1. **Benchmark** your specific dataset
2. **Compare results** with original implementation (if available)
3. **Profile** if still slow to identify remaining bottlenecks
4. **Consider** sparse P matrix for very large datasets

## Technical Details

See `OPTIMIZATION.md` for:
- Mathematical justification
- Implementation details
- Profiling instructions
- Memory/performance tradeoffs

## Questions?

If per-gene fitting still takes more than ~10 seconds:
1. Check that single-cell mode is detected (look at logs)
2. Verify donor covariates are specified correctly
3. Ensure BLAS threading is configured
4. Profile to see where time is spent

The two-stage approach should make per-gene fitting nearly instantaneous compared to full REML optimization!
