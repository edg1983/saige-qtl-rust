# GLIBC Compatibility Fix

## Problem

When running the Singularity container, the binaries fail with:
```
/usr/local/bin/step2-run-tests: /lib/x86_64-linux-gnu/libc.so.6: version `GLIBC_2.38' not found
```

## Root Cause

The multi-stage Docker build was using different base images:
- **Builder stage**: `rust:slim` (based on Debian Sid/unstable with GLIBC 2.38+)
- **Runtime stage**: `debian:bookworm-slim` (Debian 12 with GLIBC 2.36)

When Rust compiles the binaries, they link against the GLIBC version in the builder image (2.38). However, the runtime image only has GLIBC 2.36, causing the mismatch.

## Solution

Changed the builder image from `rust:slim` to `rust:1-bookworm`:

```dockerfile
# Before (mismatched)
FROM rust:slim AS builder              # GLIBC 2.38+ (Debian Sid)
...
FROM debian:bookworm-slim              # GLIBC 2.36 (Debian 12)

# After (matched)
FROM rust:1-bookworm AS builder        # GLIBC 2.36 (Debian 12)
...
FROM debian:bookworm-slim              # GLIBC 2.36 (Debian 12)
```

Both stages now use Debian 12 (Bookworm) with the same GLIBC version (2.36).

## Why rust:1-bookworm?

- `rust:1-bookworm`: Official Rust image based on Debian 12 Bookworm
- `rust:slim`: Based on Debian Sid (unstable), has newer packages but can cause compatibility issues
- Bookworm is stable and widely compatible with HPC systems using Debian-based distributions

## Verification

Rebuild the Docker image and test:

```bash
# Rebuild
docker build -t saige-qtl-rust:latest .

# Test with provided script
chmod +x test-container-glibc.sh
./test-container-glibc.sh
```

The test script verifies:
1. Docker image builds successfully
2. All three binaries run in Docker
3. Singularity conversion works
4. All binaries run in Singularity container

## Manual Testing

```bash
# Build Docker image
docker build -t saige-qtl-rust:latest .

# Check GLIBC version in container
docker run --rm saige-qtl-rust:latest bash -c "ldd --version | head -1"
# Should show: ldd (Debian GLIBC 2.36-9+deb12u...) 2.36

# Check binary requirements
docker run --rm saige-qtl-rust:latest bash -c "ldd /usr/local/bin/step2-run-tests | grep libc"
# Should show: libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x...)

# Test binaries work
docker run --rm saige-qtl-rust:latest step2-run-tests --help
docker run --rm saige-qtl-rust:latest step1-fit-null --help
docker run --rm saige-qtl-rust:latest compute-grm --help

# Convert to Singularity
singularity build saige-qtl-rust.sif docker-daemon://saige-qtl-rust:latest

# Test in Singularity
singularity exec saige-qtl-rust.sif step2-run-tests --help
singularity exec saige-qtl-rust.sif step1-fit-null --help
singularity exec saige-qtl-rust.sif compute-grm --help
```

## For HPC Users

After rebuilding the container:

```bash
# Pull new image (when published to registry)
singularity pull saige-qtl-rust.sif docker://ghcr.io/edg1983/saige-qtl-rust:latest

# Or build from local Docker
singularity build saige-qtl-rust.sif docker-daemon://saige-qtl-rust:latest

# Test
singularity exec saige-qtl-rust.sif step2-run-tests --help
```

The error should be resolved.

## Alternative: Static Linking

If GLIBC issues persist on very old systems, consider using `musl` for static linking:

```dockerfile
FROM rust:1-alpine AS builder
# Alpine uses musl libc instead of glibc
# Produces fully static binaries
```

However, this requires:
- All dependencies support musl compilation
- May have different performance characteristics
- More complex build process for scientific libraries (OpenBLAS, etc.)

The current solution (matching Debian versions) is simpler and sufficient for most HPC environments.

## GLIBC Version Compatibility

Debian/Ubuntu GLIBC versions:
- **Debian 10 (Buster)**: GLIBC 2.28
- **Debian 11 (Bullseye)**: GLIBC 2.31
- **Debian 12 (Bookworm)**: GLIBC 2.36 ← Current choice
- **Debian Sid (Unstable)**: GLIBC 2.38+
- **Ubuntu 20.04**: GLIBC 2.31
- **Ubuntu 22.04**: GLIBC 2.35
- **Ubuntu 24.04**: GLIBC 2.39

Our choice (Debian 12 Bookworm with GLIBC 2.36) is compatible with:
✅ Ubuntu 22.04 and newer
✅ Debian 12 and newer
✅ Most modern HPC systems
✅ Singularity/Apptainer environments

For older systems (Ubuntu 20.04, Debian 11), consider using `rust:1-bullseye` instead.
