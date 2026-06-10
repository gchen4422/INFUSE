# ============================================================
# Build script for the INFUSE workflowr website
# Run interactively or via: Rscript build_site.R
# ============================================================

library(workflowr)

# ------------------------------------------------------------
# 1. Git identity (only needed once per machine)
# ------------------------------------------------------------
wflow_git_config(
  user.name  = "Guangxin Chen",
  user.email = "chen4422@purdue.edu",
  overwrite  = TRUE
)

# ------------------------------------------------------------
# 2. Initialise workflowr (safe to re-run; preserves existing files)
# ------------------------------------------------------------
setwd("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/website")
wflow_start("INFUSE_website", existing = TRUE, overwrite = FALSE, change_wd = TRUE)

# ------------------------------------------------------------
# 3. Build pages (in navbar order; eval=FALSE chunks are skipped)
# ------------------------------------------------------------
pages <- c(
  "analysis/index.Rmd",
  "analysis/installation.Rmd",
  "analysis/Data_prepare.Rmd",
  "analysis/Run_INFUSE.Rmd",
  "analysis/annotation_integration.Rmd",
  "analysis/real_data_example.Rmd",
  "analysis/about.Rmd",
  "analysis/license.Rmd"
)
wflow_build(pages)


wflow_view()
# ------------------------------------------------------------
# 4. Publish: commit rendered HTML + source Rmd together
# ------------------------------------------------------------
wflow_publish(pages, message = "Rebuild INFUSE website")

# ------------------------------------------------------------
# 5. Status check
# ------------------------------------------------------------
wflow_status()
