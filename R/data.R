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
