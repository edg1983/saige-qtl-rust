# MAC (Minor Allele Count) Calculation

## Definition

**MAC (Minor Allele Count)** is the count of the less frequent allele in the dataset.

For a biallelic variant:
- Reference allele count = sum of ref alleles
- Alternate allele count = sum of alt alleles  
- MAC = min(ref_count, alt_count)

## Calculation in saige-qtl-rust

### For DS (Dosage) Field

```rust
// DS represents the dosage of the ALTERNATE allele (0.0 to 2.0)
let allele_sum = g.sum();              // Sum of alternate allele dosages
let allele_count = n_samples * 2.0;    // Total alleles (diploid)

// MAC is minimum of alt and ref counts
let mac = allele_sum.min(allele_count - allele_sum);
```

**Example with 3 samples:**
```
DS values: [0.1, 0.2, 0.9]
allele_sum = 0.1 + 0.2 + 0.9 = 1.2 (alternate alleles)
allele_count = 3 * 2 = 6 (total alleles)
ref_count = 6 - 1.2 = 4.8 (reference alleles)
MAC = min(1.2, 4.8) = 1.2
```

### For GT (Genotype) Field

```rust
// GT is converted to 0/1/2 (count of alternate alleles)
let allele_sum = g.sum();              // Sum of alternate alleles
let allele_count = n_samples * 2.0;    // Total alleles (diploid)

// MAC is minimum of alt and ref counts
let mac = allele_sum.min(allele_count - allele_sum);
```

**Example with 3 samples:**
```
Genotypes: 0/0, 0/1, 1/1
Converted: [0, 1, 2]
allele_sum = 0 + 1 + 2 = 3 (alternate alleles)
allele_count = 3 * 2 = 6 (total alleles)
ref_count = 6 - 3 = 3 (reference alleles)
MAC = min(3, 3) = 3
```

## Important: Calculation Before Centering

**Critical**: MAC must be calculated **before** centering genotypes for association testing.

### Why?

The association test requires centered genotypes (mean = 0), but MAC needs the original allele counts:

```rust
// CORRECT ORDER:
1. Extract genotypes: [0.1, 0.2, 0.9]
2. Calculate MAC: 1.2 (from original values)
3. Center genotypes: [-0.3, -0.2, 0.5] (mean = 0.4)
4. Run association test with centered genotypes

// WRONG ORDER (previous bug):
1. Extract genotypes: [0.1, 0.2, 0.9]
2. Center genotypes: [-0.3, -0.2, 0.5]
3. Calculate MAC from centered: wrong!
```

## Common MAC Thresholds

| Threshold | Use Case                  | Example                   |
| --------- | ------------------------- | ------------------------- |
| MAC ≥ 1   | Very liberal, for testing | Allows singletons         |
| MAC ≥ 5   | Common for eQTL           | ~0.3% MAF in 1000 samples |
| MAC ≥ 10  | More stringent            | ~0.6% MAF in 1000 samples |
| MAC ≥ 20  | Conservative              | ~1.2% MAF in 1000 samples |

## Relationship to MAF

**MAF (Minor Allele Frequency)** = MAC / (2 × n_samples)

For 1000 samples:
- MAC = 1 → MAF = 0.05% (very rare)
- MAC = 5 → MAF = 0.25%
- MAC = 10 → MAF = 0.5%
- MAC = 20 → MAF = 1.0%

## Implementation Details

### Current Implementation (Fixed)

```rust
fn extract_genotypes(
    record: &mut bcf::Record,
    vcf_field: &str,
    vcf_indices: &[usize],
) -> Result<(Array1<f64>, f64), Box<dyn std::error::Error>> {
    
    // 1. Extract raw genotypes (DS or GT)
    let mut g = Array1::zeros(n_samples_out);
    // ... fill g with values ...
    
    // 2. Calculate MAC from ORIGINAL values
    let allele_sum = g.sum();
    let allele_count = n_samples_out as f64 * 2.0;
    let mac = allele_sum.min(allele_count - allele_sum);
    
    // 3. Center genotypes for association testing
    let mean = g.mean().unwrap_or(0.0);
    g.mapv_inplace(|v| v - mean);
    
    // 4. Return both centered genotypes and MAC
    Ok((g, mac))
}
```

### Previous Implementation (Bug)

```rust
// OLD (INCORRECT):
// 1. Extract and center genotypes
let g = extract_genotypes(...)?;  // Already centered!

// 2. Calculate MAC from centered genotypes (WRONG!)
let mac = g.iter().map(|&x| x.abs()).sum();  // This is meaningless!
```

**Why the old approach was wrong:**
- Centered genotypes have mean = 0
- Absolute values of centered genotypes don't represent allele counts
- For rare variants, could give artificially high MAC
- For common variants, could give artificially low MAC

## Logging

With `RUST_LOG=trace`, you'll see for each variant:

```
Genotype extraction: mean=0.4000, allele_sum=1.2000, MAC=1.2000, n_samples=3
```

- `mean`: Average genotype (before centering)
- `allele_sum`: Total alternate allele count
- `MAC`: Minor Allele Count (used for filtering)
- `n_samples`: Number of samples

## Testing MAC Filtering

To see which variants pass MAC filtering:

```bash
# Very permissive (allow rare variants)
--min-mac 1

# Standard for eQTL
--min-mac 5

# Conservative
--min-mac 10
```

With `RUST_LOG=trace`, filtered variants will show:
```
Variant chr1:12345:rs123 filtered by MAC: 0.80 < 5
```

This tells you the variant had MAC=0.8, which is below your threshold of 5.

## Edge Cases

### Monomorphic Variants
```
All samples: 0/0
allele_sum = 0
MAC = min(0, 6) = 0
→ Filtered (MAC < min_mac)
```

### Fixed Alternate
```
All samples: 1/1
allele_sum = 6
MAC = min(6, 0) = 0
→ Filtered (MAC < min_mac)
```

### Perfect 50/50
```
Samples: 0/0, 0/1, 1/1
allele_sum = 3
MAC = min(3, 3) = 3
→ Passes if min_mac ≤ 3
```

## Summary

**Key Points:**
1. MAC = Minor Allele Count (not frequency)
2. Calculated from **original** genotypes before centering
3. MAC = min(alt_count, ref_count)
4. For DS: sum of dosages represents alt allele count
5. For GT: sum of 0/1/2 codes represents alt allele count
6. Filtering happens before association testing
7. Use `--min-mac` to control quality threshold
