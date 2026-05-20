# INFUSE

**INFUSE**: **IN**tegrating **F**unctional annotations into m**U**lti-ancestry fine-mapping via **S**um of single **E**ffect model.

![INFUSE logo](INFUSE.png)

## Overview

INFUSE is an R package for statistical fine-mapping of causal genetic variants
using GWAS summary statistics from multiple ancestry groups. It makes two key
methodological contributions:

**1. High-dimensional functional annotation integration.**
INFUSE parameterizes per-SNP causal priors via a multinomial logit model
defined on the *sparse principal components* of a functional annotation matrix
(e.g., BaselineLF v2.2). This allows it to incorporate hundreds of correlated
annotations while learning a separate, signal-specific functional prior for
each independent association signal at a locus — improving fine-mapping
resolution, particularly at loci with multiple causal variants.

**2. Joint multi-ancestry LD-discrepancy test (`kriging_rss_joint`).**
INFUSE introduces a novel joint diagnostic that simultaneously tests for LD
mismatches and allele flips across *both* ancestries. It stacks z-scores as
$(z_1, z_2)$ and constructs a block-diagonal LD model
$R = \mathrm{blockdiag}(R_1, R_2)$, estimates a shared regularization
parameter $s$ by maximum likelihood, and computes per-variant log-likelihood
ratios (logLR) using a mixture-of-normals model on the conditional residuals.
A joint logLR that sums evidence across ancestries provides greater sensitivity
than single-ancestry diagnostics, and is used to flag and remove problematic
variants before fine-mapping. In simulations, INFUSE with joint LD-discrepancy
correction achieved power of 0.81 under allele-flipping conditions, close to
the oracle (0.84), while uncorrected methods ranged from 0.05 to 0.27.

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install from GitHub
devtools::install_github("gchen4422/INFUSE")

# Load the package
library(INFUSE)
```

## Example data

A ready-to-use example dataset is bundled with the package:

```r
data(infuse_example)

R_mat_list           <- infuse_example$R_mat_list           # LD matrices (EUR, AFR)
summary_stat_sd_list <- infuse_example$summary_stat_sd_list # GWAS summary stats
annot_file_all       <- infuse_example$annot_file_all       # 1432 x 186 annotation matrix
ancestry_weight      <- infuse_example$ancestry_weight      # prior weights c(3,3,1)/7
```

## Core functions

### `infuse_core()`

```r
infuse_core <- function(
  R_mat_list,              # named list of LD correlation matrices by ancestry (SNP x SNP)
  summary_stat_list,       # named list of data.frames with SNP, A1, A2, BETA, SE (harmonized, same SNP order)
  L,                       # integer: max number of causal signals
  residual_variance = NULL,
  prior_weights = NULL,    # optional numeric vector length m (per-SNP prior multipliers)
  ancestry_weight = NULL,  # optional named numeric vector (e.g., c(EUR = 0.6, AFR = 0.4))
  optim_method = "optim",
  estimate_residual_variance = FALSE,
  max_iter = 100,
  cor_method = "min.abs.corr",
  cor_threshold = 0.5,
  annot = NULL,            # SNP x feature matrix/data.frame aligned to SNP order
  annot_method = NULL,     # e.g., "mlk" for sparsePCA, "glm", or NULL
  est_annot_prior = "fixed"
)
```

#### Arguments

- **`R_mat_list`**  
  Named list of ancestry-specific **LD correlation matrices** (symmetric; identical SNP IDs and identical SNP order across ancestries).

- **`summary_stat_list`**  
  Named list of GWAS summary-stat data.frames aligned to LD. Required columns: `SNP`, `A1`, `A2`, `BETA`, `SE` (optionally `CHR`, `BP`). Alleles should be uppercase and harmonized across ancestries.

- **`L`**  
  Integer: maximum number of causal signals to model for the region.

- **`residual_variance`**  
  Numeric or `NULL`. If `NULL`, use default or estimate when `estimate_residual_variance = TRUE`.

- **`prior_weights`**  
  Optional numeric vector of length *m* (per-SNP prior multipliers), aligned to SNP order.

- **`ancestry_weight`**  
  Optional **named** numeric vector weighting ancestries in the joint objective (e.g., `c(EUR = 0.6, AFR = 0.4)`).

- **`optim_method`**  
  Optimizer choice (default `"optim"`).

- **`estimate_residual_variance`**  
  Logical. If `TRUE`, estimate residual variance.

- **`max_iter`**  
  Integer: maximum optimization iterations.

- **`cor_method`, `cor_threshold`**  
  Internal correlation/LD-consistency handling (defaults `"min.abs.corr"` and `0.5`).

- **`annot`**  
  Annotation matrix/data.frame (SNP × features), rows aligned to SNP order.

- **`annot_method`**  
  How to use annotations (e.g., `"mlk"` for sparse PCA, `"glm"`, or `NULL`).

- **`est_annot_prior`**  
  Annotation-prior mode (e.g., `"fixed"`).

---

### `kriging_rss_joint()` — Joint multi-ancestry LD-discrepancy test

A novel diagnostic that simultaneously detects LD mismatches and allele flips
across two ancestries. This is a key novelty of INFUSE: unlike single-ancestry
diagnostics (e.g., `susieR::kriging_rss`), `kriging_rss_joint` stacks both
ancestry z-vectors and constructs a joint conditional model, yielding a combined
log-likelihood ratio with higher sensitivity for catching problematic variants.

```r
kriging_rss_joint(
  z1, z2,          # z-score vectors for ancestry 1 and 2 (same SNP order)
  R1, R2,          # LD correlation matrices for ancestry 1 and 2
  n1 = NULL,       # sample size for ancestry 1 (enables finite-sample shrinkage)
  n2 = NULL,       # sample size for ancestry 2
  r_tol = 1e-8,    # eigenvalue floor for PSD enforcement
  s = NULL         # regularization in [0,1]; estimated by MLE if NULL
)
```

**Returns** a list with:
- `$plot`: ggplot2 scatter of observed vs. expected z-scores colored by ancestry;
  allele-flip candidates (logLR > 2 and |z_std_diff| > 4) highlighted in red.
- `$conditional_dist`: data frame with columns `ancestry`, `z`, `condmean`,
  `condvar`, `z_std_diff`, `logLR` (per-variant), `logLR_joint` (summed across
  both ancestries for the same SNP).

**Typical usage** (run before `infuse_core()` to clean the input data):

```r
diag <- kriging_rss_joint(z_eur, z_afr, R_eur, R_afr, n1 = 193593, n2 = 60760)
diag$plot   # inspect; red points = candidate mismatches

# Remove flagged SNPs
bad_snps <- diag$conditional_dist$logLR_joint > 2 &
            abs(diag$conditional_dist$z_std_diff) > 4
clean_idx <- which(!bad_snps[seq_along(z_eur)])  # ancestry-1 positions
```

---

## Runtime and reproducibility notes

INFUSE may rely on numerical linear algebra routines that use multi-threaded BLAS, LAPACK, or OpenMP backends. On shared computing environments, such as HPC clusters, uncontrolled multi-threading can lead to CPU oversubscription, unstable runtime, and small numerical differences across runs.

To improve runtime stability and numerical reproducibility, we recommend limiting the number of threads used by common linear algebra backends before running INFUSE.

For shell-based workflows or Slurm jobs, add the following lines before launching R:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

For example, in a Slurm job script:

```bash
#!/bin/bash
#SBATCH --job-name=infuse
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

Rscript run_infuse.R
```

For RStudio, these variables can be set before loading INFUSE:

```r
Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

library(INFUSE)
```

For repeated analyses, we recommend adding these settings to an `.Renviron` file so they are applied automatically when R starts:

```bash
OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

These settings are especially useful when running many loci in parallel, where each job should typically use a single computational thread unless the user explicitly allocates additional CPUs per task.

