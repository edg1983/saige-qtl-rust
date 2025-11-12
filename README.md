# saige-qtl-rust

Reimplements key components of SAIGE-QTL in Rust for improved performance and safety.

## Usage

Main executables are in target/release/ after compilation.

## Compilation Environment

First load needed system dependencies:

```bash
module load llvm-clang/20.0.0git openblas/0.3.24
```

Then compile using cargo:

```bash
cargo build --release
```
