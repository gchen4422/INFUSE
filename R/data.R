#' Example dataset for INFUSE
#'
#' A list of ready-to-use inputs for \code{\link{infuse_core}}, drawn from a
#' 1,432-SNP locus on chromosome 1 (GRCh38: ~7.9–9.5 Mb) with EUR and AFR
#' GWAS summary statistics from the \emph{All of Us} Research Program.
#'
#' @format A list with four elements:
#' \describe{
#'   \item{R_mat_list}{Named list of two 1432 × 1432 LD correlation matrices
#'     (\code{EUR} and \code{AFR}), derived from 1000 Genomes Project (GRCh38).
#'     Row and column names are rsIDs.}
#'   \item{summary_stat_sd_list}{Named list of two data frames (\code{EUR} and
#'     \code{AFR}), each with 1,432 rows. Columns include \code{SNP}, \code{CHR},
#'     \code{POS}, \code{Z}, \code{Beta}, \code{Se}, \code{MAF}, \code{N},
#'     \code{INFO}, \code{ALT}, \code{REF}, \code{PVAL}, \code{CHR_POS}.}
#'   \item{annot_file_all}{1432 × 186 numeric matrix of BaselineLF v2.2
#'     functional annotations from PolyFun (Weissbrod et al., 2022), covering
#'     coding, regulatory, chromatin, conservation, and CpG features.}
#'   \item{ancestry_weight}{Numeric vector of length 3 giving the prior weight
#'     on EUR-specific, AFR-specific, and shared causal signals respectively.
#'     Set to \code{c(3, 3, 1) / 7} to encourage detection of ancestry-specific
#'     signals.}
#' }
#' @source
#' GWAS summary statistics: All of Us Research Program
#' (\url{https://allofus.nih.gov}).
#' LD reference: 1000 Genomes Project, GRCh38.
#' Functional annotations: BaselineLF v2.2 from PolyFun
#' (Weissbrod et al., \emph{Nature Genetics}, 2022,
#' \doi{10.1038/s41588-020-00735-5}).
#' @examples
#' data(infuse_example)
#'
#' # Inspect dimensions
#' lapply(infuse_example$R_mat_list, dim)          # 1432 x 1432 for EUR and AFR
#' lapply(infuse_example$summary_stat_sd_list, dim) # 1432 x 13
#' dim(infuse_example$annot_file_all)               # 1432 x 186
#' print(infuse_example$ancestry_weight)            # c(0.429, 0.429, 0.143)
"infuse_example"

#' Two-ancestry variant-availability example for INFUSE
#'
#' Observed EUR and AFR inputs for a chromosome-1 BMI region, including variants
#' available in only one ancestry. Its 1,432 common variants have the same SNP
#' identifiers as the original \code{\link{infuse_example}} dataset.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{summary_stat_list}{Named EUR and AFR summary-statistic data frames with
#'     1,556 and 2,733 observed variants, respectively. Columns include SNP, Beta,
#'     Se, Z, N, A1 and A2. Unavailable variants are not present in these frames.}
#'   \item{R_mat_list}{Named LD matrices of dimensions 1,556 by 1,556 for EUR and
#'     2,733 by 2,733 for AFR, aligned to the observed summary statistics.}
#'   \item{annot}{A numeric 2,857 by 186 functional-annotation matrix, with unknown
#'     annotation records excluded before constructing the analysis union.}
#'   \item{annot_ids, union_ids}{Variant identifiers in union order.}
#'   \item{allele_ref}{Reference alleles used for harmonization.}
#'   \item{variant_available}{A logical 2,857 by 2 matrix describing availability
#'     in EUR and AFR: 1,432 shared, 124 EUR-only, and 1,301 AFR-only variants.}
#'   \item{ancestry_names}{The ancestry labels EUR and AFR.}
#'   \item{source}{Trait, chromosome, locus and annotation-matching information.}
#' }
#' @details Availability describes the supplied data, rather than establishing
#'   biological absence from another ancestry. Annotations were matched by rsID
#'   and unordered allele pair; genomic position was not used for that matching.
#'   Use \code{\link{prepare_infuse_union}} before fitting the union model.
#' @source The existing two-ancestry BMI locus inputs used for the INFUSE
#'   variant-availability analysis, with baselineLF2.2.UKB functional annotations.
#' @examples
#' data(infuse_availability_example)
#' x <- infuse_availability_example
#' vapply(x$summary_stat_list, nrow, integer(1))
#' table(rowSums(x$variant_available))
#' dim(x$annot)
#' \dontrun{
#' prep <- prepare_infuse_union(
#'   summary_stat_list = x$summary_stat_list,
#'   R_mat_list = x$R_mat_list,
#'   annot = x$annot, annot_ids = x$annot_ids,
#'   union_ids = x$union_ids, allele_ref = x$allele_ref,
#'   ld_allele_coding = "ancestry_effect",
#'   require_complete_ld_dimnames = TRUE
#' )
#' fit <- infuse_core(
#'   R_mat_list = prep$R_mat_list,
#'   summary_stat_list = prep$summary_stat_list,
#'   L = 10, variant_available = prep$variant_available,
#'   annot = prep$annot, annot_method = "mlk"
#' )
#' }
"infuse_availability_example"
