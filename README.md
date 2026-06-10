# INFUSE

**INtegrating Functional annotations into mUlti-ancestry fine-mapping via Sum of single Effect model**

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr

---

## Overview

INFUSE is an R package that performs multi-ancestry genetic fine-mapping by incorporating functional genomic annotations into a Sum of Single Effects (SuSiE) framework. It accepts GWAS summary statistics and ancestry-specific LD matrices from multiple populations, and optionally integrates SNP-level functional annotations (e.g., open chromatin, transcription factor binding, conservation scores) to improve causal variant prioritization.

The core function `infuse_core()` jointly models shared and ancestry-specific causal signals across diverse populations, providing:

- **Posterior Inclusion Probabilities (PIPs)** per SNP per ancestry
- **95% Credible Sets** for causal variant discovery
- **Annotation-informed priors** via sparse PCA (`mlk`) or GLM-based methods

---

## Installation

Install the development version from GitHub:

```r
install.packages("devtools")
devtools::install_github("gchen4422/INFUSE")
library(INFUSE)
```

---

## Quick Start

```r
library(INFUSE)

# Provide per-ancestry LD matrices and summary statistics
result <- infuse_core(
  R_mat_list        = list(EUR = eur_ld, AFR = afr_ld),
  summary_stat_list = list(EUR = eur_gwas, AFR = afr_gwas),
  L                 = 10,
  annot             = snp_annot_matrix,
  annot_method      = "mlk"
)

# Inspect credible sets
result$cs

# Inspect PIPs
head(result$pip)
```

---

## Website

The companion website documents the INFUSE methodology, installation, usage, and worked examples:

- [Home](docs/index.html) — INFUSE overview and model description
- [Installation](docs/installation.html) — Installation and motivating example
- [Run INFUSE](docs/Run_INFUSE.html) — Input format and running `infuse_core()`
- [Annotation Integration](docs/annotation_integration.html) — Incorporating functional annotations
- [Simulation Benchmark](docs/simulation_benchmark.html) — Performance evaluation under simulated data
- [Real Data Example](docs/real_data_example.html) — Applied analysis with multi-ancestry GWAS data
- [About](docs/about.html) — Authors and citation
- [License](docs/license.html)

---

## Citation

If you use INFUSE in your research, please cite:

> Chen G, et al. INFUSE improves multi-ancestry fine-mapping through high-dimensional functional annotation integration and robust LD discrepancy detection. *In preparation*.

---

## License

This project is licensed under the MIT License. See `LICENSE` for details.
