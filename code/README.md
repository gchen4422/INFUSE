# code/

This directory contains helper scripts and utility functions used by the INFUSE analysis pipeline.

## Contents

- `simulation/` — Simulation scripts for the benchmarking analysis
  - `run_simulation.R` — Main simulation runner
  - `simulate_sumstats.R` — Generates simulated GWAS summary statistics from 1KG LD
  - `evaluate_methods.R` — Computes calibration, coverage, and power metrics

- `preprocessing/` — Data preprocessing utilities
  - `gwas_qc.R` — GWAS quality control and harmonization functions
  - `ld_construction.R` — LD matrix construction from plink2 genotype files
  - `ld_mismatch.R` — Kriging RSS-based LD mismatch detection

*(Scripts to be released with the manuscript.)*
