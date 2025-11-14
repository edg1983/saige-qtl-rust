#!/bin/bash
# Helper script to run SAIGE-QTL with optimized threading settings
# This ensures both Rayon (Rust parallelism) and OpenBLAS (linear algebra) 
# use all available CPU cores

# Detect number of CPU cores
if command -v nproc &> /dev/null; then
    NUM_CORES=$(nproc)
elif [ -f /proc/cpuinfo ]; then
    NUM_CORES=$(grep -c ^processor /proc/cpuinfo)
else
    NUM_CORES=4  # Default fallback
fi

echo "Detected $NUM_CORES CPU cores"
echo "Setting threading environment variables for optimal performance..."

# Set OpenBLAS threading (used by ndarray-linalg)
export OPENBLAS_NUM_THREADS=$NUM_CORES
export OPENBLAS_MAIN_FREE=1  # Allow freeing memory from any thread

# Set generic BLAS threading (fallback for other BLAS implementations)
export OMP_NUM_THREADS=$NUM_CORES

# Set Rayon thread pool (Rust parallelism)
export RAYON_NUM_THREADS=$NUM_CORES

echo "Environment configured:"
echo "  OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
echo "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "  RAYON_NUM_THREADS=$RAYON_NUM_THREADS"
echo ""

# Pass all arguments to the command
exec "$@"
