############################################################################
#
#		This file defines GWAS QC, reference SNP alignment, and union
#		preprocessing for multi-ancestry fine-mapping
#
############################################################################
##############################################
#	Step 1. QC of GWAS for each ancestry
###############################################

#' GWAS Quality Control Function
#'
#' This function performs several quality control steps on GWAS data:
#' 1. Removes strand ambiguous variants.
#' 2. Removes multi-allelic variants.
#' 3. Excludes the MHC complex region.
#' 4. Filters SNPs based on Minor Allele Frequency (MAF).
#'
#' @param dat A data frame containing the GWAS data. Columns "CHR", "POS", "CHR_POS",
#' "REF", "ALT", and (optionally) "MAF" are required.
#' @param MAF_threshold The threshold for Minor Allele Frequency filtering. SNPs
#' with MAF below this threshold or above (1 - threshold) will be removed.
#'
#' @return A filtered data frame after performing the quality control steps.
#' @import dplyr
#' @export

GWAS_QC <- function(dat, MAF_threshold, rm_strand = TRUE, rm_multi = TRUE, rm_indel = TRUE, rm_mhc = TRUE,rm_maf = TRUE) {
  # Column name to uppers
  colnames(dat) <- toupper(colnames(dat))
  # Check if required columns are present and to upper
  if(!all(c("CHR", "POS", "REF", "ALT") %in% colnames(dat))) {
    stop("Required columns 'CHR', 'POS', 'REF', and 'ALT' are missing.")
  }

  dat <- dat %>%
    mutate(CHR_POS = paste0(CHR, "_", POS))
  # CHR, POS, REF, ALT, and MAF are required for QC
  # QC Step 1: Remove Strand Ambiguous Variants
  # These are SNPs that can't be differentiated based on the strand and may introduce errors.
  if(rm_strand){
    dat <- dat %>%
      filter(!(
        (REF == "G" & ALT == "C") |
          (REF == "C" & ALT == "G") |
          (REF == "T" & ALT == "A") |
          (REF == "A" & ALT == "T")
      ))
  }
  # QC Step 2: Remove Multi-allelic Variants
  if(rm_multi){
    multi_allelic <- dat$CHR_POS[duplicated(dat$CHR_POS)]
    if(length(multi_allelic)>0){
      dat<-dat %>%filter(!(CHR_POS%in%multi_allelic))
    }
  }

  # QC Step 3: Remove Indel if necessary
  if(rm_indel){
    dat <- dat %>%
      filter(!(nchar(REF) + nchar(ALT) > 2))
  }

  # QC Step 4: Exclude MHC Complex Region
  # This region can be particularly complex and may need special handling.
  if(rm_mhc){
    dat <- dat %>%
      filter(!(CHR == 6 & POS > 25e6 & POS < 34e6))
  }


  # QC Step 5: Minor Allele Frequency (MAF) Filter
  # Remove SNPs with MAF below the threshold or above (1 - threshold).
  if(rm_maf){
    if("MAF"%in%colnames(dat)){
      dat <- dat %>%
        filter(MAF > MAF_threshold, MAF < 1 - MAF_threshold)
    }
  }
  return(dat)
}
##############################################
#	Step 2. Find common set of SNPs, and align
#	the reference allele, switch the sign of beta
#	or z score based on the allele
###############################################

#' Identify Common SNPs Between Two Datasets
#'
#' This function identifies common SNPs (based on "CHR_POS") between a summary statistics
#' dataset (`sumstat`) and a reference dataset (`ref`). It then ensures that the
#' alleles match between the two datasets, either in the original or flipped orientation.
#'
#' @param sumstat A data frame containing the summary statistics data with at least
#' the columns "CHR_POS", "REF", and "ALT".
#' @param ref A data frame serving as the reference, with at least the columns
#' "CHR_POS", "REF", and "ALT".
#'
#' @return A vector of "CHR_POS" values corresponding to the SNPs common between
#' `sumstat` and `ref` with matching allele information.
#' @import dplyr
#' @export

find_common_snps <- function(sumstat, ref) {

  # Identifying common SNPs between sumstat and ref datasets
  common_SNP <- intersect(sumstat$CHR_POS, ref$CHR_POS)

  # Subsetting data frames based on common SNPs and arranging by order of common_SNP
  sumstat <- sumstat %>%
    filter(CHR_POS %in% common_SNP) %>%
    arrange(match(CHR_POS, common_SNP))

  ref <- ref %>%
    filter(CHR_POS %in% common_SNP) %>%
    arrange(match(CHR_POS, common_SNP))

  matched_pos<-which((sumstat$REF==ref$REF&sumstat$ALT==ref$ALT)|(sumstat$REF==ref$ALT&sumstat$ALT==ref$REF))

  common_SNP<-sumstat[matched_pos,]%>%pull(CHR_POS)

  return(common_SNP)
}
# Match the reference allele with the sumstat used,
# flip the beta and zscore if the reference allele is flipped.

#' Allele Flip Function for SNP Matching
#'
#' This function modifies the reference dataset (`ref`) based on the summary
#' statistics dataset (`sumstat`) to ensure that the alleles match. If the alleles
#' do not match, the function will flip the alleles and corresponding statistics
#' in the reference dataset.
#'
#' @param sumstat A data frame containing the summary statistics data with at least
#' the columns "REF".
#' @param ref A data frame serving as the reference with at least the columns
#' "REF", "ALT", "BETA", and "SE".
#'
#' @return A modified version of `ref` where alleles have been flipped to match
#' `sumstat` and corresponding statistics (BETA and Z) have been updated accordingly.
#' @import dplyr
#' @export

allele_flip<-function(sumstat,ref){
  ref<-ref %>% mutate(flip_index = ifelse(ref$REF==sumstat$REF,1,-1))%>%
    mutate(BETA = BETA * flip_index)%>%
    mutate(Z = BETA/SE, REF = sumstat$REF, ALT = sumstat$ALT)%>%
    dplyr::select(-flip_index)
  return(ref)
}
#################################################
#	Step 3. Match up the reference panel with
#	the GWAS reference and alternative
#	alleles
###############################################
# Here we use plink to do the referernce allele flip
# Write out the reference SNP for the genotype data
# plink --bfile ref_geno --ref-allele ref.txt --make-bed --out ref_geno_flip

####################################################################################
#
#   Union-of-variants preprocessing for availability-aware INFUSE
#
#   This layer is deliberately SEPARATE from the statistical layer.  Its only job
#   is bookkeeping:
#
#       * allele-align each ancestry's summary statistics to a reference
#       * build the union of variants and remember where each one came from
#       * record the m x K availability matrix H
#       * expand each ancestry's summary statistics and LD matrix to the union
#       * align annotations to the union once
#
#   The statistical layer (infuse_core) then uses `variant_available` to
#   decide what is computable.  It never inspects the padded values.
#
#   Padding convention
#   ------------------
#   For a variant j that is unavailable in ancestry k we write
#
#       Beta = 0,  Z = 0,  Se = 1,  N = median(N in ancestry k)
#       LD row/col = e_j    (1 on the diagonal, 0 elsewhere)
#
#   This is *neutral*, not informative:
#
#     - XtX.diag[k][j] is then a finite O(1) number, so shat2 never becomes 0/Inf
#       and no division by zero can occur even if some code path touches it;
#     - Xty[k][j] = XtX.diag[k][j] * 0 = 0;
#     - the identity LD row means XtX[k][j, ] has no off-diagonal entries, so a
#       padded variant cannot leak into any other variant's residual;
#     - because M_{js}=0 for every configuration s containing ancestry k,
#       alpha_{lj}^{(s)}=0 for every effect l, hence EB1[j,k]=EB2[j,k]=0 and the
#       padded cell contributes exactly zero to Xr, the ELBO, and the
#       residual-variance update.
#
#   Setting the padded LD *diagonal* to 1 rather than 0 matters: XtX.diag comes
#   from Se, so a zero diagonal in R would make XtX[j,j] and XtX.diag[j]
#   inconsistent between compute_Xr() and compute_SER_posterior_loglik().
#
####################################################################################


.validate_ids <- function(ids, label, allow_empty = FALSE) {
  ids <- trimws(as.character(ids))
  if (!allow_empty && length(ids) == 0L) {
    stop(label, " is empty.", call. = FALSE)
  }
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop(label, " must contain non-missing, non-empty identifiers.", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop(label, " contains duplicated variant identifiers.", call. = FALSE)
  }
  ids
}


.normalize_alleles <- function(x, label) {
  x <- toupper(trimws(as.character(x)))
  if (anyNA(x) || any(!nzchar(x))) {
    stop(label, " contains missing or empty allele values.", call. = FALSE)
  }
  x
}


.validate_allele_ref <- function(ref) {
  if (!is.data.frame(ref)) {
    stop("allele_ref must be a data frame with columns SNP, A1, and A2.",
         call. = FALSE)
  }
  miss <- setdiff(c("SNP", "A1", "A2"), names(ref))
  if (length(miss)) {
    stop("allele_ref is missing column(s): ", paste(miss, collapse = ", "), ".",
         call. = FALSE)
  }
  ref$SNP <- .validate_ids(ref$SNP, "allele_ref$SNP")
  ref$A1 <- .normalize_alleles(ref$A1, "allele_ref$A1")
  ref$A2 <- .normalize_alleles(ref$A2, "allele_ref$A2")
  if (any(ref$A1 == ref$A2)) {
    stop("allele_ref contains variant(s) with identical A1 and A2.", call. = FALSE)
  }
  ref
}


#' Align annotations to a variant union
#'
#' @param annot Data frame or matrix of annotations.
#' @param annot_ids Character vector of variant identifiers for the rows of
#'   `annot`.  If `NULL`, `rownames(annot)` is used.
#' @param union_ids Character vector of union variant identifiers.
#' @param fill Finite numeric value used for union variants absent from `annot`
#'   (default 0). Use zero only when it is a scientifically valid annotation value;
#'   if absence means unknown annotation, filter those variants before calling this
#'   function rather than encoding unknown values as zero.
#' @return A numeric matrix with `length(union_ids)` rows in union order.
#' @export
align_annotations_to_union <- function(annot, annot_ids = NULL, union_ids, fill = 0) {
  if (is.null(annot)) return(NULL)
  A <- as.matrix(annot)
  if (is.null(annot_ids)) annot_ids <- rownames(A)
  if (is.null(annot_ids)) {
    stop("Cannot align annotations: supply annot_ids or set rownames(annot).", call. = FALSE)
  }
  annot_ids <- .validate_ids(annot_ids, "annot_ids", allow_empty = TRUE)
  union_ids <- .validate_ids(union_ids, "union_ids", allow_empty = TRUE)
  if (length(annot_ids) != nrow(A)) {
    stop(sprintf("length(annot_ids) is %d but annot has %d rows.",
                 length(annot_ids), nrow(A)), call. = FALSE)
  }
  if (ncol(A) == 0L) {
    stop("annot must contain at least one annotation column.", call. = FALSE)
  }
  if (is.null(colnames(A)) || anyNA(colnames(A)) ||
      any(!nzchar(trimws(colnames(A)))) || anyDuplicated(colnames(A))) {
    stop("annot must have unique, non-missing, non-empty column names.",
         call. = FALSE)
  }
  if (!(is.numeric(A) || is.logical(A))) {
    stop("annot must contain only numeric or logical values.", call. = FALSE)
  }
  if (anyNA(A) || any(!is.finite(A))) {
    stop("annot contains missing or non-finite values.", call. = FALSE)
  }
  if (length(fill) != 1L || !is.numeric(fill) || !is.finite(fill)) {
    stop("fill must be one finite numeric value.", call. = FALSE)
  }

  out <- matrix(fill, nrow = length(union_ids), ncol = ncol(A),
                dimnames = list(union_ids, colnames(A)))
  hit <- match(union_ids, annot_ids)
  ok <- !is.na(hit)
  out[ok, ] <- A[hit[ok], , drop = FALSE]
  storage.mode(out) <- "double"
  attr(out, "n_missing") <- sum(!ok)
  out
}


#' Allele-align one ancestry's summary statistics to a reference
#'
#' Variants whose alleles are swapped relative to the reference have `Beta` and
#' `Z` sign-flipped, and their returned allele labels are set to the reference
#' order. Re-aligning the returned summary statistics to that reference therefore
#' does not flip them again. Variants whose alleles are genuinely incompatible (not the
#' same unordered pair under exact allele labels) are reported as `mismatch` -- they are
#' **never** silently recoded as unavailable.
#'
#' @param ss Data frame with at least `SNP`, `Beta`, `Se`, `Z`, `N`; allele
#'   columns named by `a1_col`/`a2_col` if allele alignment is wanted.
#' @param ref Data frame with `SNP`, `A1`, `A2` giving the reference alleles.
#' @param a1_col,a2_col Names of the effect / other allele columns in `ss`.
#' @return List with `ss` (aligned), `flipped` and `mismatch` (logical vectors
#'   over `nrow(ss)`).
#' @keywords internal
allele_align_one <- function(ss, ref, a1_col = "A1", a2_col = "A2") {
  if (!is.data.frame(ss)) {
    stop("ss must be a data frame.", call. = FALSE)
  }
  n <- nrow(ss)
  flipped <- rep(FALSE, n)
  mismatch <- rep(FALSE, n)

  present <- c(a1_col, a2_col) %in% names(ss)
  if (xor(present[1L], present[2L])) {
    stop("summary statistics must contain both ", a1_col, " and ", a2_col,
         ", or neither.", call. = FALSE)
  }
  if (is.null(ref) || !all(present)) {
    return(list(ss = ss, flipped = flipped, mismatch = mismatch))
  }

  req <- c("SNP", "Beta", "Z", a1_col, a2_col)
  miss <- setdiff(req, names(ss))
  if (length(miss)) {
    stop("summary statistics are missing column(s): ", paste(miss, collapse = ", "),
         ".", call. = FALSE)
  }
  ss$SNP <- .validate_ids(ss$SNP, "summary-statistic SNP identifiers",
                          allow_empty = TRUE)
  ref <- .validate_allele_ref(ref)

  a1 <- .normalize_alleles(ss[[a1_col]], paste0("summary statistics$", a1_col))
  a2 <- .normalize_alleles(ss[[a2_col]], paste0("summary statistics$", a2_col))
  if (any(a1 == a2)) {
    stop("summary statistics contain variant(s) with identical effect and other alleles.",
         call. = FALSE)
  }
  if (!is.numeric(ss$Beta) || !is.numeric(ss$Z) ||
      any(!is.finite(ss$Beta)) || any(!is.finite(ss$Z))) {
    stop("summary-statistic Beta and Z must be numeric and finite.", call. = FALSE)
  }

  hit <- match(as.character(ss$SNP), as.character(ref$SNP))
  have <- !is.na(hit)

  r1 <- ref$A1[hit]
  r2 <- ref$A2[hit]

  same  <- have & a1 == r1 & a2 == r2
  swap  <- have & a1 == r2 & a2 == r1
  mismatch <- have & !same & !swap

  ss$Beta[swap] <- -ss$Beta[swap]
  ss$Z[swap]    <- -ss$Z[swap]
  # Keep the returned effect-allele labels consistent with the flipped effects.
  # Unmatched and incompatible variants retain their original normalized alleles.
  a1[swap] <- r1[swap]
  a2[swap] <- r2[swap]
  ss[[a1_col]] <- a1
  ss[[a2_col]] <- a2
  flipped <- swap

  list(ss = ss, flipped = flipped, mismatch = mismatch)
}


#' Build availability-aware INFUSE inputs on the union of ancestry variant sets
#'
#' @param summary_stat_list Named list of length `K`. Element `k` holds ancestry
#'   `k`'s summary statistics -- only the variants actually available in that
#'   ancestry -- with columns `SNP`, `Beta`, `Se`, `Z`, `N` and optionally `A1`,
#'   `A2`.
#' @param R_mat_list Named list of length `K` of LD correlation matrices. Element
#'   `k` must be square with `dimnames` (or at least `colnames`) matching
#'   `summary_stat_list[[k]]$SNP`.
#' @param annot Optional annotation matrix / data frame.
#' @param annot_ids Optional variant identifiers for `annot` rows.
#' @param union_ids Optional explicit union ordering.  Defaults to first-seen
#'   order across ancestries, which keeps ancestry 1's variants in their original
#'   order followed by the variants unique to later ancestries.
#' @param allele_ref Optional data frame with `SNP`, `A1`, `A2` used as the allele
#'   reference. For two ancestries, this defaults to ancestry 1's alleles when
#'   available. When `K > 2` and allele columns are supplied, an explicit reference
#'   is required and must cover every variant observed in more than one ancestry;
#'   a full-union reference is recommended.
#' @param on_allele_mismatch `"error"` (default) or `"drop"`.  Mismatching
#'   variants are never treated as merely unavailable.
#' @param ld_allele_coding **The LD allele contract — you must state which applies.**
#'   Allele alignment sign-flips `Beta` and `Z` for swapped variants, and the LD
#'   matrix must end up in the *same* allele coding as the flipped statistics.
#'   \describe{
#'     \item{`"reference"` (default)}{`R_mat_list[[k]]` is already oriented to the
#'       same reference alleles used for alignment (the usual case when LD comes
#'       from a reference panel keyed on a fixed allele order). No LD change is
#'       made.}
#'     \item{`"ancestry_effect"`}{`R_mat_list[[k]]` follows ancestry `k`'s own
#'       effect-allele coding. Row *and* column `j` of `R_k` are then sign-flipped
#'       for every variant `j` whose statistics were flipped (the diagonal is
#'       unaffected, since \eqn{(-1)(-1)=1}).}
#'   }
#'   Getting this wrong silently corrupts every cross-variant correlation involving
#'   a flipped variant, so there is no "auto" option.
#' @param drop_unavailable_everywhere If `TRUE` (default) variants absent from all
#'   ancestries are dropped with a message; if `FALSE` an error is raised.
#' @param check_ld_psd If `TRUE`, perform an exact eigenvalue check that each LD
#'   matrix is positive semidefinite. This is `FALSE` by default because it is
#'   cubic in the number of variants; it is recommended for input audits and small
#'   or moderate loci.
#' @param ld_tolerance Numeric tolerance for LD symmetry, unit diagonal, and the
#'   correlation bounds.
#' @param ld_psd_tolerance Allowed magnitude of a negative LD eigenvalue when
#'   `check_ld_psd = TRUE`.
#' @param require_complete_ld_dimnames If `TRUE`, require both row and column SNP
#'   names on every LD matrix. The default is `FALSE` for compatibility with legacy
#'   drivers, but strict mode is recommended for new analyses because the order of
#'   an unnamed LD axis cannot be verified.
#' @return List with `summary_stat_list`, `R_mat_list`, `variant_available`,
#'   `annot`, `variant_map` and `dropped`.  The first three are ready to pass
#'   straight to [infuse_core()].
#' @export
prepare_infuse_union <- function(summary_stat_list, R_mat_list, annot = NULL,
                                annot_ids = NULL, union_ids = NULL,
                                allele_ref = NULL,
                                on_allele_mismatch = c("error", "drop"),
                                ld_allele_coding = c("reference", "ancestry_effect"),
                                drop_unavailable_everywhere = TRUE,
                                verbose = TRUE,
                                check_ld_psd = FALSE,
                                ld_tolerance = 1e-8,
                                ld_psd_tolerance = 1e-8,
                                require_complete_ld_dimnames = FALSE) {

  on_allele_mismatch <- match.arg(on_allele_mismatch)
  ld_allele_coding <- match.arg(ld_allele_coding)
  for (nm in c("drop_unavailable_everywhere", "verbose", "check_ld_psd",
               "require_complete_ld_dimnames")) {
    val <- get(nm, inherits = FALSE)
    if (length(val) != 1L || is.na(val) || !is.logical(val)) {
      stop(nm, " must be one non-missing logical value.", call. = FALSE)
    }
  }
  if (length(ld_tolerance) != 1L || !is.numeric(ld_tolerance) ||
      !is.finite(ld_tolerance) || ld_tolerance < 0) {
    stop("ld_tolerance must be one finite non-negative number.", call. = FALSE)
  }
  if (length(ld_psd_tolerance) != 1L || !is.numeric(ld_psd_tolerance) ||
      !is.finite(ld_psd_tolerance) || ld_psd_tolerance < 0) {
    stop("ld_psd_tolerance must be one finite non-negative number.", call. = FALSE)
  }
  ancestry_names <- resolve_ancestry_names(summary_stat_list, R_mat_list)
  K <- length(ancestry_names)

  req <- c("SNP", "Beta", "Se", "Z", "N")
  allele_presence <- vapply(ancestry_names, function(k) {
    ss <- summary_stat_list[[k]]
    if (!is.data.frame(ss)) {
      stop(sprintf("summary_stat_list[['%s']] must be a data frame.", k),
           call. = FALSE)
    }
    present <- c("A1", "A2") %in% names(ss)
    if (xor(present[1L], present[2L])) {
      stop(sprintf("summary_stat_list[['%s']] must contain both A1 and A2, or neither.",
                   k), call. = FALSE)
    }
    all(present)
  }, logical(1))
  if (any(allele_presence) && !all(allele_presence)) {
    stop("A1/A2 allele metadata are present for only some ancestries. Supply both columns ",
         "for every ancestry, or remove them from every ancestry to disable allele alignment.",
         call. = FALSE)
  }
  if (!is.null(allele_ref) && !all(allele_presence)) {
    stop("allele_ref was supplied, but A1/A2 are not present for every ancestry.",
         call. = FALSE)
  }

  for (k in ancestry_names) {
    miss <- setdiff(req, names(summary_stat_list[[k]]))
    if (length(miss)) {
      stop(sprintf("summary_stat_list[['%s']] is missing column(s): %s.",
                   k, paste(miss, collapse = ", ")), call. = FALSE)
    }
    if (nrow(summary_stat_list[[k]]) == 0L) {
      stop(sprintf("summary_stat_list[['%s']] has no variants.", k), call. = FALSE)
    }
    summary_stat_list[[k]]$SNP <- .validate_ids(
      summary_stat_list[[k]]$SNP,
      sprintf("summary_stat_list[['%s']]$SNP", k)
    )
    for (nm in c("Beta", "Se", "Z", "N")) {
      x <- summary_stat_list[[k]][[nm]]
      if (!is.numeric(x) || any(!is.finite(x))) {
        stop(sprintf("summary_stat_list[['%s']]$%s must be numeric and finite.", k, nm),
             call. = FALSE)
      }
    }
    if (any(summary_stat_list[[k]]$Se <= 0)) {
      stop(sprintf("summary_stat_list[['%s']]$Se must be strictly positive.", k),
           call. = FALSE)
    }
    if (any(summary_stat_list[[k]]$N <= 2)) {
      stop(sprintf("summary_stat_list[['%s']]$N must be greater than 2.", k),
           call. = FALSE)
    }
    # Validate and normalize BOTH LD axes.  The union expansion indexes
    # rows and columns with the same positions, so checking colnames alone
    # would silently corrupt LD when named rows are in a different order.
    Rk <- R_mat_list[[k]]
    if (!is.matrix(Rk) || !is.numeric(Rk)) {
      stop(sprintf("R_mat_list[['%s']] must be a numeric matrix.", k), call. = FALSE)
    }
    if (nrow(Rk) != ncol(Rk)) {
      stop(sprintf("R_mat_list[['%s']] must be square; got %d x %d.",
                   k, nrow(Rk), ncol(Rk)), call. = FALSE)
    }
    ss_ids <- summary_stat_list[[k]]$SNP
    rn <- rownames(Rk); cn <- colnames(Rk)

    if (require_complete_ld_dimnames && (is.null(rn) || is.null(cn))) {
      stop(sprintf("R_mat_list[['%s']] must have both row and column SNP names when require_complete_ld_dimnames = TRUE.",
                   k), call. = FALSE)
    }

    if (is.null(rn) && is.null(cn)) {
      if (nrow(Rk) != length(ss_ids)) {
        stop(sprintf("R_mat_list[['%s']] has no dimnames and its dimension (%d) does ",
                     k, nrow(Rk)),
             sprintf("not match summary_stat_list[['%s']] (%d rows).",
                     k, length(ss_ids)), call. = FALSE)
      }
      if (verbose) {
        warning(sprintf("R_mat_list[['%s']] has no dimnames; assuming both axes follow summary-statistic row order.",
                        k), call. = FALSE)
      }
      dimnames(Rk) <- list(ss_ids, ss_ids)
    } else {
      # Preserve the one-named-axis convention used by legacy drivers, but make
      # the uncheckable positional assumption explicit.
      if (is.null(rn) || is.null(cn)) {
        if (verbose) {
          warning(sprintf("R_mat_list[['%s']] has SNP identifiers on only one axis; assuming the unnamed axis has the same positional order.",
                          k), call. = FALSE)
        }
        if (is.null(rn)) rn <- cn
        if (is.null(cn)) cn <- rn
      }
      rn <- .validate_ids(rn, sprintf("rownames(R_mat_list[['%s']])", k))
      cn <- .validate_ids(cn, sprintf("colnames(R_mat_list[['%s']])", k))
      if (!setequal(rn, ss_ids) || !setequal(cn, ss_ids)) {
        stop(sprintf("row/col names of R_mat_list[['%s']] do not match summary_stat_list[['%s']]$SNP.",
                     k, k), call. = FALSE)
      }
      dimnames(Rk) <- list(rn, cn)
      Rk <- Rk[ss_ids, ss_ids, drop = FALSE]
    }
    if (any(!is.finite(Rk))) {
      stop(sprintf("R_mat_list[['%s']] contains non-finite values.", k), call. = FALSE)
    }
    if (!isTRUE(all.equal(Rk, t(Rk), tolerance = ld_tolerance,
                          check.attributes = FALSE))) {
      stop(sprintf("R_mat_list[['%s']] is not symmetric after SNP alignment.", k),
           call. = FALSE)
    }
    if (any(abs(diag(Rk) - 1) > ld_tolerance)) {
      stop(sprintf("R_mat_list[['%s']] must have unit diagonal within tolerance %g.",
                   k, ld_tolerance), call. = FALSE)
    }
    if (any(abs(Rk) > 1 + ld_tolerance)) {
      stop(sprintf("R_mat_list[['%s']] contains values outside the correlation range [-1, 1].",
                   k), call. = FALSE)
    }
    if (check_ld_psd) {
      min_eigen <- min(eigen((Rk + t(Rk)) / 2, symmetric = TRUE,
                             only.values = TRUE)$values)
      if (min_eigen < -ld_psd_tolerance) {
        stop(sprintf("R_mat_list[['%s']] is not positive semidefinite: minimum eigenvalue %.6g.",
                     k, min_eigen), call. = FALSE)
      }
    }
    R_mat_list[[k]] <- Rk
  }

  ## ------------------------------------------------------- allele alignment ----
  if (all(allele_presence) && is.null(allele_ref) && K > 2L) {
    stop("For K > 2 with A1/A2 columns, supply an explicit allele_ref. Ancestry 1 ",
         "need not contain variants shared only by later ancestries.", call. = FALSE)
  }
  if (all(allele_presence) && is.null(allele_ref)) {
    ss1 <- summary_stat_list[[1]]
    allele_ref <- data.frame(SNP = ss1$SNP, A1 = ss1$A1, A2 = ss1$A2,
                             stringsAsFactors = FALSE)
  }
  if (!is.null(allele_ref)) {
    allele_ref <- .validate_allele_ref(allele_ref)
    observed_ids <- unlist(lapply(summary_stat_list, `[[`, "SNP"), use.names = FALSE)
    id_count <- table(observed_ids)
    multi_ancestry_ids <- names(id_count)[id_count > 1L]
    missing_ref <- setdiff(multi_ancestry_ids, allele_ref$SNP)
    if (length(missing_ref)) {
      stop(sprintf("allele_ref omits %d variant(s) observed in multiple ancestries (first few: %s).",
                   length(missing_ref),
                   paste(utils::head(missing_ref, 5), collapse = ", ")),
           call. = FALSE)
    }
  }

  flip_flag  <- vector("list", K); names(flip_flag)  <- ancestry_names
  mismatch_flag <- vector("list", K); names(mismatch_flag) <- ancestry_names
  # Allele mismatch is ancestry-specific: a variant whose alleles cannot be
  # reconciled in ancestry k is unusable *in ancestry k only*.  Dropping it from
  # every ancestry would silently discard usable data.
  drop_mismatch <- vector("list", K); names(drop_mismatch) <- ancestry_names

  for (k in ancestry_names) {
    al <- allele_align_one(summary_stat_list[[k]], allele_ref)
    summary_stat_list[[k]] <- al$ss
    flip_flag[[k]] <- al$flipped
    mismatch_flag[[k]] <- al$mismatch
    drop_mismatch[[k]] <- character(0)
    if (any(al$mismatch)) {
      bad <- summary_stat_list[[k]]$SNP[al$mismatch]
      if (on_allele_mismatch == "error") {
        stop(sprintf("%d variant(s) in ancestry '%s' have alleles incompatible with the ",
                     length(bad), k),
             "reference (first few: ", paste(utils::head(bad, 5), collapse = ", "),
             "). Resolve the allele coding explicitly, or pass ",
             "on_allele_mismatch = \"drop\". Allele mismatch is NOT the same as ",
             "variant unavailability.", call. = FALSE)
      }
      drop_mismatch[[k]] <- bad
    }
  }
  n_mismatch <- sum(vapply(drop_mismatch, length, integer(1)))
  if (n_mismatch && verbose) {
    message(sprintf("prepare_infuse_union(): marking %d allele-mismatched variant/ancestry pair(s) unavailable (%s).",
                    n_mismatch,
                    paste(sprintf("%s=%d", ancestry_names,
                                  vapply(drop_mismatch, length, integer(1))),
                          collapse = ", ")))
  }

  ## ------------------------------------------------------------------ union ----
  per_anc_ids <- lapply(ancestry_names, function(k)
    setdiff(summary_stat_list[[k]]$SNP, drop_mismatch[[k]]))
  names(per_anc_ids) <- ancestry_names

  if (is.null(union_ids)) {
    union_ids <- unique(unlist(per_anc_ids, use.names = FALSE))
    union_ids <- .validate_ids(union_ids, "automatically constructed union_ids")
  } else {
    union_ids <- .validate_ids(union_ids, "union_ids")
    extra <- setdiff(unique(unlist(per_anc_ids, use.names = FALSE)), union_ids)
    if (length(extra)) {
      stop(sprintf("union_ids omits %d variant(s) present in the input.", length(extra)),
           call. = FALSE)
    }
  }

  variant_available <- vapply(ancestry_names, function(k) union_ids %in% per_anc_ids[[k]],
                              logical(length(union_ids)))
  variant_available <- matrix(variant_available, nrow = length(union_ids),
                              dimnames = list(union_ids, ancestry_names))

  dead <- which(!apply(variant_available, 1, any))
  dropped <- character(0)
  if (length(dead)) {
    if (!drop_unavailable_everywhere) {
      stop(sprintf("%d variant(s) are unavailable in every ancestry: %s.",
                   length(dead), paste(utils::head(union_ids[dead], 5), collapse = ", ")),
           call. = FALSE)
    }
    dropped <- union_ids[dead]
    if (verbose) {
      message(sprintf("prepare_infuse_union(): dropping %d variant(s) unavailable in every ancestry.",
                      length(dropped)))
    }
    union_ids <- union_ids[-dead]
    variant_available <- variant_available[-dead, , drop = FALSE]
  }
  p <- length(union_ids)
  if (p == 0L) stop("The variant union is empty.", call. = FALSE)

  ## ------------------------------------------------------ expand to the union ---
  ss_out <- vector("list", K); names(ss_out) <- ancestry_names
  R_out  <- vector("list", K); names(R_out)  <- ancestry_names
  orig_index <- matrix(NA_integer_, nrow = p, ncol = K,
                       dimnames = list(union_ids, ancestry_names))
  flip_out <- matrix(FALSE, nrow = p, ncol = K,
                     dimnames = list(union_ids, ancestry_names))

  for (k in ancestry_names) {
    ss <- summary_stat_list[[k]]
    hit <- match(union_ids, ss$SNP)
    avail_k <- variant_available[, k]
    hit[!avail_k] <- NA_integer_          # mismatched-and-dropped stay unavailable
    orig_index[, k] <- hit
    flip_out[avail_k, k] <- flip_flag[[k]][hit[avail_k]]

    Nmed <- stats::median(ss$N, na.rm = TRUE)

    new_ss <- data.frame(
      SNP  = union_ids,
      Beta = 0,
      Se   = 1,
      Z    = 0,
      N    = Nmed,
      stringsAsFactors = FALSE
    )
    ok <- !is.na(hit)
    new_ss$Beta[ok] <- ss$Beta[hit[ok]]
    new_ss$Se[ok]   <- ss$Se[hit[ok]]
    new_ss$Z[ok]    <- ss$Z[hit[ok]]
    new_ss$N[ok]    <- ss$N[hit[ok]]

    # guard against degenerate observed values that would break XtX.diag
    bad_se <- ok & (!is.finite(new_ss$Se) | new_ss$Se <= 0)
    if (any(bad_se)) {
      stop(sprintf("ancestry '%s': %d observed variant(s) have non-positive or ",
                   k, sum(bad_se)),
           "non-finite Se.", call. = FALSE)
    }
    ss_out[[k]] <- new_ss

    # ---- LD allele orientation ------------------------------------------
    # If the LD follows this ancestry's own effect-allele coding, then flipping
    # Beta/Z for a swapped variant must be mirrored in the LD: negate row j and
    # column j of R_k.  The diagonal is unchanged because (-1)*(-1) = 1.
    if (ld_allele_coding == "ancestry_effect" && any(flip_flag[[k]])) {
      Rk0 <- R_mat_list[[k]]
      fl <- as.character(ss$SNP)[flip_flag[[k]]]
      pos_fl <- match(fl, as.character(colnames(Rk0)))
      pos_fl <- pos_fl[!is.na(pos_fl)]
      if (length(pos_fl)) {
        sgn <- rep(1, ncol(Rk0)); sgn[pos_fl] <- -1
        Rk0 <- Rk0 * outer(sgn, sgn)          # negates row j and col j, keeps diag
        R_mat_list[[k]] <- Rk0
        if (verbose) {
          message(sprintf("prepare_infuse_union(): sign-flipped %d LD row/column pair(s) in '%s' (ld_allele_coding = \"ancestry_effect\").",
                          length(pos_fl), k))
        }
      }
    }

    Rk <- R_mat_list[[k]][ss$SNP, ss$SNP, drop = FALSE]
    Rnew <- diag(1, nrow = p, ncol = p)          # identity padding
    dimnames(Rnew) <- list(union_ids, union_ids)
    src <- match(union_ids[ok], as.character(colnames(Rk)))
    Rnew[ok, ok] <- Rk[src, src, drop = FALSE]
    # Input diagonals were validated above; this normalizes only tolerated
    # numerical drift and establishes the padded identity cells exactly.
    diag(Rnew) <- 1
    R_out[[k]] <- Rnew
  }

  annot_out <- if (is.null(annot)) NULL else
    align_annotations_to_union(annot, annot_ids, union_ids)

  variant_map <- data.frame(
    union_index = seq_len(p),
    snp_id = union_ids,
    stringsAsFactors = FALSE
  )
  for (k in ancestry_names) {
    variant_map[[paste0("index_", k)]]     <- orig_index[, k]
    variant_map[[paste0("available_", k)]] <- variant_available[, k]
    variant_map[[paste0("allele_flipped_", k)]] <- flip_out[, k]
  }
  variant_map$n_ancestry_available <- rowSums(variant_available)
  variant_map$availability_pattern <- apply(variant_available, 1, function(r)
    paste(ancestry_names[r], collapse = "+"))

  if (verbose) {
    tab <- table(variant_map$availability_pattern)
    message("prepare_infuse_union(): union size ", p, "; availability patterns: ",
            paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", "))
  }

  list(
    summary_stat_list = ss_out,
    R_mat_list        = R_out,
    variant_available = variant_available,
    annot             = annot_out,
    variant_map       = variant_map,
    union_ids         = union_ids,
    ancestry_names    = ancestry_names,
    dropped           = list(unavailable_everywhere = dropped,
                             allele_mismatch = drop_mismatch,
                             allele_mismatch_any = unique(unlist(drop_mismatch,
                                                                 use.names = FALSE)))
  )
}
