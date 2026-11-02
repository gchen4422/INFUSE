% Please edit documentation in R/data.R
\docType{data}
\name{infuse_availability_example}
\alias{infuse_availability_example}
\title{Two-ancestry variant-availability example for INFUSE}
\format{
A list with the following elements:
\describe{
\item{summary_stat_list}{Named EUR and AFR summary-statistic data frames with
1,556 and 2,733 observed variants. Columns include SNP, Beta, Se, Z, N, A1 and A2.}
\item{R_mat_list}{Named LD matrices with 1,556 EUR and 2,733 AFR rows and columns,
aligned to the observed summary statistics.}
\item{annot}{A numeric matrix of 2,857 union variants and 186 annotations.}
\item{annot_ids, union_ids}{Variant identifiers in union order.}
\item{allele_ref}{Reference alleles used for harmonization.}
\item{variant_available}{A logical 2,857 by 2 matrix: 1,432 shared, 124 EUR-only,
and 1,301 AFR-only variants.}
\item{ancestry_names}{EUR and AFR.}
\item{source}{Trait, chromosome, locus and annotation-matching information.}
}}
\usage{infuse_availability_example}
\description{
Observed EUR and AFR inputs for a chromosome-1 BMI region. Its common variants
have the same SNP identifiers as the original \code{\link{infuse_example}} dataset.
}
\details{
Unavailable variants are absent from each ancestry's input. Availability describes
the supplied data and does not establish biological absence. Unknown annotation
records were excluded, and annotations were matched by rsID and unordered allele
pair rather than genomic position. Use \code{\link{prepare_infuse_union}} before
fitting the union model.
}
\source{Existing two-ancestry BMI locus inputs used in the INFUSE
variant-availability analysis, with baselineLF2.2.UKB annotations.}
\examples{
data(infuse_availability_example)
x <- infuse_availability_example
vapply(x$summary_stat_list, nrow, integer(1))
table(rowSums(x$variant_available))
dim(x$annot)
\dontrun{
prep <- prepare_infuse_union(
  summary_stat_list = x$summary_stat_list,
  R_mat_list = x$R_mat_list,
  annot = x$annot, annot_ids = x$annot_ids,
  union_ids = x$union_ids, allele_ref = x$allele_ref,
  ld_allele_coding = "ancestry_effect",
  require_complete_ld_dimnames = TRUE
)
fit <- infuse_core(
  R_mat_list = prep$R_mat_list,
  summary_stat_list = prep$summary_stat_list,
  L = 10, variant_available = prep$variant_available,
  annot = prep$annot, annot_method = "mlk"
)
}}
\keyword{datasets}
