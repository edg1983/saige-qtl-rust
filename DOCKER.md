# Docker Quick Reference

## Building the Image

```bash
# Build locally
docker build -t saige-qtl-rust:latest .

# Or use the build script
./build-docker.sh
```

## Running Examples

### Step 1: Fit Null Model

```bash
docker run --rm \
  -v $(pwd)/data:/data:ro \
  -v $(pwd)/output:/output \
  -u $(id -u):$(id -g) \
  saige-qtl-rust:latest \
  step1-fit-null \
    --plink-file /data/genotypes \
    --pheno-file /data/phenotypes.tsv \
    --trait-name GENE1 \
    --covar-file /data/covariates.tsv \
    --sample-id-in-pheno IID \
    --sample-id-in-covar IID \
    --trait-type quantitative \
    --output-prefix /output/step1 \
    --n-threads 8
```

### Step 2: Run Association Tests

```bash
docker run --rm \
  -v $(pwd)/data:/data:ro \
  -v $(pwd)/output:/output \
  -u $(id -u):$(id -g) \
  saige-qtl-rust:latest \
  step2-run-tests \
    --vcf-file /data/genotypes.vcf.gz \
    --vcf-field DS \
    --pheno-file /data/phenotypes.tsv \
    --sample-id-in-pheno IID \
    --trait-name GENE1 \
    --covar-file /data/covariates.tsv \
    --sample-id-in-covar IID \
    --step1-output-file /output/step1.GENE1.model.bin \
    --region chr1:1000000-2000000 \
    --output-file /output/step2.GENE1.results.txt.gz \
    --min-mac 5 \
    --n-threads 8
```

## Docker Compose Usage

```bash
# Run with docker-compose
UID=$(id -u) GID=$(id -g) docker-compose run --rm saige-qtl step1-fit-null --help

# Run predefined step1 service
UID=$(id -u) GID=$(id -g) docker-compose up step1

# Run predefined step2 service
UID=$(id -u) GID=$(id -g) docker-compose up step2
```

## Singularity/Apptainer for HPC

```bash
# Build from Docker image
singularity build saige-qtl-rust.sif docker-daemon://saige-qtl-rust:latest

# Or pull from registry (when published)
singularity pull saige-qtl-rust.sif docker://ghcr.io/edg1983/saige-qtl-rust:latest

# Run Step 1
singularity exec \
  --bind ./data:/data,./output:/output \
  saige-qtl-rust.sif \
  step1-fit-null \
    --plink-file /data/genotypes \
    --pheno-file /data/phenotypes.tsv \
    --trait-name GENE1 \
    --covar-file /data/covariates.tsv \
    --sample-id-in-pheno IID \
    --sample-id-in-covar IID \
    --trait-type quantitative \
    --output-prefix /output/step1 \
    --n-threads 8

# Run Step 2
singularity exec \
  --bind ./data:/data,./output:/output \
  saige-qtl-rust.sif \
  step2-run-tests \
    --vcf-file /data/genotypes.vcf.gz \
    --vcf-field DS \
    --pheno-file /data/phenotypes.tsv \
    --sample-id-in-pheno IID \
    --trait-name GENE1 \
    --covar-file /data/covariates.tsv \
    --sample-id-in-covar IID \
    --step1-output-file /output/step1.GENE1.model.bin \
    --region chr1:1000000-2000000 \
    --output-file /output/step2.GENE1.results.txt.gz \
    --min-mac 5 \
    --n-threads 8
```

## Environment Variables

- `RUST_LOG`: Set logging level (`debug`, `info`, `warn`, `error`)
- `OPENBLAS_NUM_THREADS`: Control OpenBLAS threading (default: 1)

```bash
docker run --rm \
  -e RUST_LOG=debug \
  -e OPENBLAS_NUM_THREADS=4 \
  -v $(pwd)/data:/data:ro \
  -v $(pwd)/output:/output \
  saige-qtl-rust:latest \
  step1-fit-null [options...]
```

## Volume Mounts

- `/data`: Input data directory (mount read-only with `:ro`)
- `/output`: Output directory (mount read-write)

## Tips

1. **Use `-u $(id -u):$(id -g)`** to run as your user and avoid permission issues
2. **Mount data as read-only** (`:ro`) for safety
3. **Set `OPENBLAS_NUM_THREADS=1`** when using `--n-threads` for better parallelization
4. **Use absolute paths** or `$(pwd)` when mounting volumes

## Troubleshooting

### GLIBC Version Mismatch in Singularity

**Error**: `version 'GLIBC_2.38' not found`

**Cause**: The builder and runtime images use different GLIBC versions.

**Solution**: The Dockerfile now uses `rust:1-bookworm` (Debian 12) for both build and runtime stages, ensuring GLIBC compatibility. Rebuild the image:

```bash
docker build -t saige-qtl-rust:latest .
singularity build saige-qtl-rust.sif docker-daemon://saige-qtl-rust:latest
```

**Verification**: Test the container with the provided script:
```bash
./test-container-glibc.sh
```

This will verify:
- Docker image builds correctly
- All binaries run in Docker
- Singularity conversion works
- All binaries run in Singularity

### Permission Issues

If you get permission errors, ensure you're running with your user ID:
```bash
docker run --rm -u $(id -u):$(id -g) ...
```

### Memory Issues

Increase Docker memory limit in Docker Desktop settings or add `--memory` flag:
```bash
docker run --rm --memory=16g ...
```

### Debugging

Run with debug logging:
```bash
docker run --rm -e RUST_LOG=debug ...
```

Or get a shell inside the container:
```bash
docker run --rm -it \
  -v $(pwd)/data:/data \
  -v $(pwd)/output:/output \
  --entrypoint /bin/bash \
  saige-qtl-rust:latest
```
