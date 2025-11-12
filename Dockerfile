# Multi-stage build for saige-qtl-rust
# Stage 1: Build the Rust binaries
# Use Rust 1.84+ which supports edition 2024 dependencies
FROM rust:slim AS builder

# Install system dependencies required for compilation
RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    libssl-dev \
    libopenblas-dev \
    libgfortran5 \
    clang \
    llvm \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /build

# Copy only dependency files first for better caching
COPY Cargo.toml Cargo.lock ./

# Create dummy source to build dependencies
RUN mkdir src && \
    echo "fn main() {}" > src/main.rs && \
    mkdir -p src/bin && \
    echo "fn main() {}" > src/bin/step1.rs && \
    echo "fn main() {}" > src/bin/step2.rs

# Build dependencies (this layer will be cached)
RUN cargo build --release && \
    rm -rf src

# Copy actual source code
COPY src ./src

# Build the actual binaries
# Touch files to ensure rebuild
RUN touch src/lib.rs src/bin/step1.rs src/bin/step2.rs && \
    cargo build --release

# Verify binaries were built
RUN ls -lh target/release/step1-fit-null target/release/step2-run-tests

# Stage 2: Create minimal runtime image
FROM debian:bookworm-slim

# Install runtime dependencies only
RUN apt-get update && apt-get install -y \
    libopenblas0 \
    libgfortran5 \
    libgomp1 \
    libssl3 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for running the application
RUN useradd -m -u 1000 -s /bin/bash saige && \
    mkdir -p /data /output && \
    chown -R saige:saige /data /output

# Copy binaries from builder
COPY --from=builder /build/target/release/step1-fit-null /usr/local/bin/
COPY --from=builder /build/target/release/step2-run-tests /usr/local/bin/

# Make binaries executable
RUN chmod +x /usr/local/bin/step1-fit-null /usr/local/bin/step2-run-tests

# Set up environment
ENV RUST_LOG=info
ENV OPENBLAS_NUM_THREADS=1

# Switch to non-root user
USER saige
WORKDIR /data

# Default command shows help
CMD ["step1-fit-null", "--help"]

# Labels for metadata
LABEL maintainer="Edoardo Giacopuzzi <edoardo.giacopuzzi@fht.org>"
LABEL description="SAIGE-QTL Rust implementation for single-cell eQTL analysis"
LABEL version="0.1.0"
LABEL org.opencontainers.image.source="https://github.com/edg1983/saige-qtl-rust"
