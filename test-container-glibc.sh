#!/bin/bash
# Test script to verify GLIBC compatibility in Docker image

set -e

echo "Building Docker image..."
docker build -t saige-qtl-rust:glibc-test .

echo ""
echo "Checking GLIBC version requirements..."
docker run --rm saige-qtl-rust:glibc-test bash -c "ldd --version | head -1"

echo ""
echo "Checking binary GLIBC requirements..."
docker run --rm saige-qtl-rust:glibc-test bash -c "ldd /usr/local/bin/step2-run-tests | grep libc"

echo ""
echo "Testing step2-run-tests --help..."
docker run --rm saige-qtl-rust:glibc-test step2-run-tests --help | head -5

echo ""
echo "Testing compute-grm --help..."
docker run --rm saige-qtl-rust:glibc-test compute-grm --help | head -5

echo ""
echo "Testing step1-fit-null --help..."
docker run --rm saige-qtl-rust:glibc-test step1-fit-null --help | head -5

echo ""
echo "✓ All binaries work correctly in container!"
echo ""
echo "Converting to Singularity (if singularity is available)..."
if command -v singularity &> /dev/null; then
    singularity build saige-qtl-rust-test.sif docker-daemon://saige-qtl-rust:glibc-test
    echo ""
    echo "Testing in Singularity container..."
    singularity exec saige-qtl-rust-test.sif step2-run-tests --help | head -5
    echo ""
    echo "✓ Singularity container works!"
    rm -f saige-qtl-rust-test.sif
else
    echo "Singularity not available, skipping Singularity test"
fi

echo ""
echo "All tests passed!"
