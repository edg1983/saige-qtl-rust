# Debugging Step 2 Empty Results

## Quick Diagnostic Commands

### 1. Check with INFO level logging (recommended first)
```bash
RUST_LOG=info ./target/release/step2-run-tests \
  --vcf-file data.vcf.gz \
  --vcf-field DS \
  [other args...]
```

Look for these key messages:
- `VCF has X samples` - Number of samples in VCF
- `Null model has Y samples` - Number of samples in null model
- `Sample alignment: Z overlapping samples` - **Critical**: Must be > 10
- `Read X total variants from VCF` - Total variants loaded
- `Variants passing filters: Y` - Final count

### 2. Check with DEBUG level for more details
```bash
RUST_LOG=debug ./target/release/step2-run-tests \
  [args...]
```

Additional information:
- First 5 sample alignments (VCF → Model mapping)
- Null model matrix dimensions
- Per-10k variant progress

### 3. Check with TRACE level for per-variant details
```bash
RUST_LOG=trace ./target/release/step2-run-tests \
  [args...] 2>&1 | tee step2_trace.log
```

**Warning**: Very verbose! Shows:
- Why each variant was filtered
- Genotype extraction details
- MAC values for each variant

## Common Issues and Solutions

### Issue 1: No Sample Overlap

**Symptom**:
```
ERROR: NO SAMPLES OVERLAP between VCF and null model!
Error: Only 0 samples overlap
```

**Cause**: Sample IDs in VCF don't match sample IDs in null model

**Solution**:
1. Check VCF sample names:
   ```bash
   bcftools query -l data.vcf.gz | head
   ```

2. Check null model sample IDs (from Step 1):
   ```bash
   # The sample IDs used in Step 1 (--sample-id-col)
   ```

3. Ensure they match **exactly** (case-sensitive)

### Issue 2: All Variants Filtered by MAC

**Symptom**:
```
WARNING: No variants passed filters!
  This could be due to:
  1. All variants have MAC < 1 (min_mac)
```

**Cause**: Minor Allele Count too low

**Debug**:
```bash
RUST_LOG=trace ./target/release/step2-run-tests [args...] 2>&1 | grep "filtered by MAC"
```

**Solutions**:
- Lower `--min-mac` (try `--min-mac 0.5`)
- Check if VCF has variants with sufficient allele frequency
- Verify genotype field is correct (`--vcf-field DS` or `GT`)

### Issue 3: Genotype Extraction Failed

**Symptom**:
```
WARNING: No variants passed filters!
  This could be due to:
  2. Genotype extraction failed (wrong --vcf-field?)
```

**Cause**: Wrong `--vcf-field` specified or format issue

**Debug**:
```bash
RUST_LOG=trace ./target/release/step2-run-tests [args...] 2>&1 | grep "Failed to extract"
```

Common errors:
- `No DS (dosage) field found` - Use `--vcf-field GT` instead
- `VCF field 'XX' not supported` - Only DS and GT are supported

**Solutions**:
- Check VCF format:
  ```bash
  bcftools view data.vcf.gz | grep -m 1 "^#CHROM" -A 1
  ```
- If you see `GT:DP:GQ`, use `--vcf-field GT`
- If you see `GT:DS:GP`, use `--vcf-field DS`

### Issue 4: No Genetic Variance

**Symptom**:
```
WARNING: No variants passed filters!
  This could be due to:
  3. No genetic variance in variants
```

**Cause**: All variants are monomorphic or have very low variance after centering

**Debug**:
```bash
RUST_LOG=trace ./target/release/step2-run-tests [args...] 2>&1 | grep "insufficient variance"
```

**Solutions**:
- Check VCF has polymorphic variants in the region
- Verify samples in VCF have genotype calls (not all missing)

### Issue 5: No Variants in Region

**Symptom**:
```
Read 0 total variants from VCF
WARNING: No variants found in VCF file!
```

**Cause**: 
- Wrong region specified
- VCF not indexed
- No variants in the specified region

**Solutions**:
1. Check VCF is indexed:
   ```bash
   ls -lh data.vcf.gz.tbi  # Should exist
   ```

2. Check region exists:
   ```bash
   bcftools view data.vcf.gz --regions chr1:1000000-2000000 | wc -l
   ```

3. Verify chromosome naming (chr1 vs 1):
   ```bash
   bcftools view -H data.vcf.gz | head -1 | cut -f1
   ```

4. Try without `--region` to test all variants

## Verification Checklist

Before running Step 2, verify:

- [ ] Step 1 completed successfully
- [ ] VCF file is indexed (.tbi or .csi file exists)
- [ ] Sample IDs match between:
  - VCF samples
  - Phenotype file (used in Step 1)
  - Null model (output of Step 1)
- [ ] VCF has the genotype field you're requesting (DS or GT)
- [ ] Region syntax is correct (if using --region)
- [ ] Minimum MAC is reasonable (--min-mac, try 1 for testing)

## Example Debug Session

```bash
# 1. Quick check with INFO
RUST_LOG=info ./target/release/step2-run-tests \
  --vcf-file chr1.vcf.gz \
  --vcf-field DS \
  --pheno-file expression.tsv \
  --sample-id-in-pheno donor_id \
  --trait-name GENE1 \
  --covar-file expression.tsv \
  --sample-id-in-covar donor_id \
  --step1-output-file null_model.GENE1.model.bin \
  --region chr1:1000000-2000000 \
  --output-file results.txt.gz \
  --min-mac 1

# Look for:
# - "Sample alignment: X overlapping samples" (should be > 10)
# - "Read X total variants from VCF" (should be > 0)
# - "Variants passing filters: Y" (should be > 0)

# 2. If sample alignment is 0, check sample IDs
bcftools query -l chr1.vcf.gz | head
# Compare with samples in your phenotype file

# 3. If no variants, check region
bcftools view chr1.vcf.gz --regions chr1:1000000-2000000 | grep -v "^#" | wc -l

# 4. If variants filtered, enable TRACE
RUST_LOG=trace ./target/release/step2-run-tests [same args] 2>&1 | head -500
# Look for specific filtering reasons
```

## Expected Output

**Successful run**:
```
[INFO] Starting association testing
[INFO] VCF has 1611 samples
[INFO] Null model has 1611 samples
[INFO] Sample alignment: 1611 overlapping samples
[INFO] Reading variants from VCF...
[INFO] Read 5432 total variants from VCF
[INFO] Processing variants in parallel...
[INFO] Association testing complete:
[INFO]   Total variants processed: 5432
[INFO]   Variants passing filters: 4891
[INFO] Writing 4891 results to output file
[INFO] Results written successfully
```

**Failed run (no overlap)**:
```
[INFO] VCF has 1611 samples
[INFO] Null model has 1611 samples
[ERROR] NO SAMPLES OVERLAP between VCF and null model!
[ERROR] VCF samples (first 5): ["sample_1", "sample_2", ...]
[ERROR] Model samples (first 5): ["donor1", "donor2", ...]
Error: Only 0 samples overlap between VCF and null model
```

## Performance Note

With TRACE logging, output can be very large (MB/GB for many variants). Use:
- `INFO` for normal runs
- `DEBUG` for troubleshooting
- `TRACE` only for specific problem variants (small regions)

You can also grep the output:
```bash
RUST_LOG=trace ./target/release/step2-run-tests [args] 2>&1 | grep -E "(filtered|WARNING|ERROR)"
```
