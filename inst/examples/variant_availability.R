library(INFUSE)
data(infuse_availability_example)
x <- infuse_availability_example

# Observed input sets differ: EUR has 1,556 variants and AFR has 2,733.
vapply(x$summary_stat_list, nrow, integer(1))
dim(x$annot) # 2,857 union variants and 186 annotations

prep <- prepare_infuse_union(
  summary_stat_list = x$summary_stat_list,
  R_mat_list = x$R_mat_list,
  annot = x$annot,
  annot_ids = x$annot_ids,
  union_ids = x$union_ids,
  allele_ref = x$allele_ref,
  ld_allele_coding = "ancestry_effect",
  require_complete_ld_dimnames = TRUE
)

table(prep$variant_map$availability_pattern)
# EUR+AFR: 1,432; EUR: 124; AFR: 1,301.

fit <- infuse_core(
  R_mat_list = prep$R_mat_list,
  summary_stat_list = prep$summary_stat_list,
  L = 10,
  variant_available = prep$variant_available,
  annot = prep$annot,
  annot_method = "mlk"
)

fit$pip
fit$cs
fit$availability$config_names
