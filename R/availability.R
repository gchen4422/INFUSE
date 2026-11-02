####################################################################################
#
#   Variant-availability layer for INFUSE
#
#   This file contains everything that turns the m x K availability matrix H
#   into (i) the configuration feasibility mask M, (ii) the availability-aware
#   configuration probabilities phi^H_{s,j}, and (iii) the restricted SNP index
#   sets that the likelihood / posterior-moment code must respect.  In the code,
#   the matrix whose (j,s) entry is phi^H_{s,j} is named `Qstar`.
#
#   Statistical target
#   ------------------
#   Let k = 1,...,K index ancestries, j = 1,...,m index variants in the
#   union, H_{jk} indicate whether variant j is available in ancestry k, and
#   s be a nonzero K-vector indicating an ancestry configuration.  Then
#
#       M_{js} = 1  iff H_{jk} = 1 for every ancestry k with s_k = 1,
#
#       phi^H_{s,j} = phi_s M_{js} / sum_{s'} phi_{s'} M_{js'}.
#
#   Thus, sum_s phi^H_{s,j} = 1 for every retained variant j.  Missingness is
#   a *configuration feasibility constraint*, not evidence against causality:
#   a variant keeps its full SNP-level prior mass regardless of how many
#   ancestries it happens to be observed in.
#
#   Combined with the original component-specific functional-annotation prior
#   pi_{lj} = Pr(gamma_{lj}=1 | A,w_l), where sum_j pi_{lj}=1, the joint
#   variant-configuration prior is
#
#       Pr(gamma_{lj}=1, s_l=s | A,w_l,H) = pi_{lj} phi^H_{s,j},
#
#   with sum_{j,s} pi_{lj} phi^H_{s,j}=1.  This avoids using pi for both the
#   SNP-level prior and the joint prior; pi_{lj} retains exactly the meaning used
#   in the original INFUSE Supplementary Note.
#
#   NOTE ON NORMALIZATION.  The software object `Qstar`, which stores
#   phi^H_{s,j}, is row-stochastic by construction, and the joint prior is
#   *never* globally renormalized.  A bug in the masking logic must surface as a
#   failed assertion, not be silently absorbed by a rescaling step.
#
#   NOTE ON ORDERING.  Nothing here assumes that ancestry 1 is EUR or AFR, or
#   that the configuration columns are in any particular order.  Configurations
#   are always identified by ancestry *membership* derived from ancestry names.
#
####################################################################################


#' Build the INFUSE ancestry-configuration list
#'
#' Reproduces the `column_config` construction used by the upstream INFUSE
#' `meSuSieObject`: all non-empty subsets of `seq_len(n_ancestry)`, ordered by
#' increasing subset size and, within a size, by `utils::combn` column order.
#'
#' @param n_ancestry Total number of ancestries, `K`; `k` indexes ancestry.
#' @return A list of length `2^K - 1`; each element is a one-column matrix of
#'   ancestry indices.
#' @export
build_column_config <- function(n_ancestry) {
  stopifnot(length(n_ancestry) == 1L, n_ancestry >= 1L)
  Reduce(append, lapply(seq_len(n_ancestry), function(x) {
    poss_config <- utils::combn(seq_len(n_ancestry), x)
    lapply(seq_len(ncol(poss_config)), function(y) matrix(poss_config[, y]))
  }))
}


#' Ancestry-name membership of each configuration
#'
#' @param column_config Output of [build_column_config()] (or the `column_config`
#'   field of an INFUSE model object).
#' @param ancestry_names Character vector of length `K`, in the same order as
#'   `summary_stat_list` / `R_mat_list`.
#' @return A list of character vectors giving the ancestries required by each
#'   configuration.
#' @export
config_ancestry_sets <- function(column_config, ancestry_names) {
  lapply(column_config, function(idx) ancestry_names[as.integer(idx)])
}


#' Stable human-readable configuration labels
#'
#' Single-ancestry configurations are labelled by the ancestry name; multi-ancestry
#' configurations are labelled by the `+`-joined ancestry names, e.g. `"AFR+EUR"`.
#'
#' @inheritParams config_ancestry_sets
#' @return Character vector of configuration labels.
#' @export
build_config_labels <- function(column_config, ancestry_names) {
  vapply(config_ancestry_sets(column_config, ancestry_names),
         function(x) paste(x, collapse = "+"),
         character(1))
}


#' Resolve the ancestry names of an INFUSE input
#'
#' @param summary_stat_list Named list of per-ancestry summary statistics.
#' @param R_mat_list Optional named list of per-ancestry LD matrices; if supplied
#'   its names must agree with `summary_stat_list`.
#' @return Character vector of ancestry names.
#' @export
resolve_ancestry_names <- function(summary_stat_list, R_mat_list = NULL) {
  if (!is.list(summary_stat_list)) {
    stop("summary_stat_list must be a list.", call. = FALSE)
  }
  nm <- names(summary_stat_list)
  nm_trim <- if (is.null(nm)) NULL else trimws(as.character(nm))
  if (is.null(nm_trim) || anyNA(nm_trim) || any(!nzchar(nm_trim)) ||
      anyDuplicated(nm_trim)) {
    stop("summary_stat_list must be a named list with unique, non-empty ancestry names. ",
         "Ancestry identity must be carried by names, never by list position.",
         call. = FALSE)
  }
  if (any(nm != nm_trim)) {
    stop("summary_stat_list ancestry names must not contain surrounding whitespace.",
         call. = FALSE)
  }
  if (!is.null(R_mat_list)) {
    if (!is.list(R_mat_list)) {
      stop("R_mat_list must be a list.", call. = FALSE)
    }
    nm_r <- names(R_mat_list)
    nm_r_trim <- if (is.null(nm_r)) NULL else trimws(as.character(nm_r))
    if (is.null(nm_r_trim) || anyNA(nm_r_trim) || any(!nzchar(nm_r_trim)) ||
        anyDuplicated(nm_r_trim) || !identical(nm_r, nm)) {
      stop("names(R_mat_list) must be identical to names(summary_stat_list); got (",
           paste(nm_r, collapse = ", "), ") vs (", paste(nm, collapse = ", "), ").",
           call. = FALSE)
    }
  }
  nm
}


#' Coerce and validate a variant-availability matrix
#'
#' @param variant_available `NULL`, or a logical (or 0/1) matrix / data frame
#'   representing `H`, with `m = n_snp` rows and one column per ancestry.  Column
#'   names, when present, are
#'   used to reorder the matrix to `ancestry_names`; positional interpretation is
#'   only used (with a warning) when column names are absent.
#' @param n_snp Number of variants in the union, `m`.
#' @param ancestry_names Character vector of ancestry names, in model order.
#' @param variant_ids Optional variant identifiers used for row names and error
#'   messages.
#' @return The logical `m x K` availability matrix `H`, whose columns are in
#'   `ancestry_names` order and whose entry `(j,k)` is `H_{jk}`.
#' @export
normalize_variant_available <- function(variant_available, n_snp, ancestry_names,
                                        variant_ids = NULL) {
  if (length(n_snp) != 1L || !is.numeric(n_snp) || !is.finite(n_snp) ||
      n_snp < 0 || n_snp != as.integer(n_snp)) {
    stop("n_snp must be one finite, non-negative integer.", call. = FALSE)
  }
  n_snp <- as.integer(n_snp)

  ancestry_names <- as.character(ancestry_names)
  if (!length(ancestry_names) || anyNA(ancestry_names) ||
      any(!nzchar(trimws(ancestry_names))) || anyDuplicated(ancestry_names)) {
    stop("ancestry_names must contain unique, non-missing, non-empty names.",
         call. = FALSE)
  }
  if (any(ancestry_names != trimws(ancestry_names))) {
    stop("ancestry_names must not contain surrounding whitespace.", call. = FALSE)
  }
  K <- length(ancestry_names)

  if (!is.null(variant_ids)) {
    variant_ids <- as.character(variant_ids)
    if (length(variant_ids) != n_snp) {
      stop(sprintf("variant_ids has length %d but n_snp is %d.",
                   length(variant_ids), n_snp), call. = FALSE)
    }
    if (anyNA(variant_ids) || any(!nzchar(trimws(variant_ids)))) {
      stop("variant_ids must contain non-missing, non-empty identifiers.",
           call. = FALSE)
    }
    if (anyDuplicated(variant_ids)) {
      stop("variant_ids contains duplicated identifiers; cannot align ",
           "variant_available by name.", call. = FALSE)
    }
  }

  if (is.null(variant_available)) {
    out <- matrix(TRUE, nrow = n_snp, ncol = K,
                  dimnames = list(variant_ids, ancestry_names))
    return(out)
  }

  if (is.data.frame(variant_available)) variant_available <- as.matrix(variant_available)
  if (!is.matrix(variant_available)) {
    stop("variant_available must be a matrix or data frame.", call. = FALSE)
  }
  if (nrow(variant_available) != n_snp) {
    stop(sprintf("variant_available has %d rows but the model has %d variants.",
                 nrow(variant_available), n_snp), call. = FALSE)
  }
  if (ncol(variant_available) != K) {
    stop(sprintf("variant_available has %d columns but there are %d ancestries.",
                 ncol(variant_available), K), call. = FALSE)
  }

  cn <- colnames(variant_available)
  unnamed_columns <- is.null(cn) ||
    (!anyNA(cn) && all(!nzchar(cn)))
  if (unnamed_columns) {
    warning("variant_available has no column names; assuming its columns are in the ",
            "same order as summary_stat_list (", paste(ancestry_names, collapse = ", "),
            "). Supplying column names is strongly recommended.", call. = FALSE)
    colnames(variant_available) <- ancestry_names
  } else {
    if (anyNA(cn) || any(!nzchar(trimws(cn))) || anyDuplicated(cn)) {
      stop("colnames(variant_available) must be unique, non-missing, and non-empty.",
           call. = FALSE)
    }
    if (!setequal(cn, ancestry_names)) {
      stop("colnames(variant_available) = (", paste(cn, collapse = ", "),
           ") do not match the ancestry names (", paste(ancestry_names, collapse = ", "),
           ").", call. = FALSE)
    }
    variant_available <- variant_available[, ancestry_names, drop = FALSE]
  }

  # Availability is a data contract, not a generic truthiness test.
  # Coercing first would silently turn values such as 2 or -1 into TRUE.
  if (!is.logical(variant_available) && !is.numeric(variant_available)) {
    stop("variant_available must contain only logical or numeric 0/1 values.",
         call. = FALSE)
  }
  if (anyNA(variant_available)) {
    stop("variant_available contains NA; availability must be TRUE or FALSE for ",
         "every variant/ancestry pair.", call. = FALSE)
  }
  if (is.numeric(variant_available) && any(!variant_available %in% c(0, 1))) {
    stop("numeric variant_available must contain only 0 or 1.", call. = FALSE)
  }
  storage.mode(variant_available) <- "logical"

  # ---- ROW alignment -------------------------------------------------------
  # Rows must be matched to the model's variant order BY NAME whenever names are
  # available.  Overwriting the row names instead would silently mis-assign the
  # mask for a correctly-named but differently-ordered availability matrix.
  rn <- rownames(variant_available)
  if (!is.null(variant_ids)) {
    if (!is.null(rn) && anyNA(rn)) {
      stop("rownames(variant_available) contains missing identifiers.",
           call. = FALSE)
    }
    if (!is.null(rn) && any(nzchar(rn))) {
      if (anyNA(rn) || any(!nzchar(trimws(rn))) || anyDuplicated(rn)) {
        stop("rownames(variant_available) contains duplicated or empty identifiers; ",
             "row names must be unique, non-missing, and non-empty.",
             call. = FALSE)
      }
      if (!setequal(rn, variant_ids)) {
        extra <- setdiff(rn, variant_ids); missing <- setdiff(variant_ids, rn)
        stop("rownames(variant_available) do not match the model's variant IDs. ",
             if (length(missing)) paste0(length(missing), " variant(s) missing from ",
                                         "variant_available (e.g. ",
                                         paste(utils::head(missing, 3), collapse = ", "), "). "),
             if (length(extra)) paste0(length(extra), " unexpected row name(s) (e.g. ",
                                       paste(utils::head(extra, 3), collapse = ", "), "). "),
             call. = FALSE)
      }
      # reorder rows to the model's variant order
      variant_available <- variant_available[variant_ids, , drop = FALSE]
    } else {
      warning("variant_available has no row names; assuming its rows are in the ",
              "same order as summary_stat_list. Supplying row names is strongly ",
              "recommended.", call. = FALSE)
      rownames(variant_available) <- variant_ids
    }
  }

  variant_available
}


#' Configuration feasibility mask M_{js}
#'
#' `M_{js}` is `TRUE` when variant `j` is available in every ancestry selected by
#' configuration vector `s`.  It is built entirely from ancestry membership, so
#' it is correct for any `K` and any ancestry ordering.
#'
#' @param variant_available Logical `m x K` matrix `H` (see
#'   [normalize_variant_available()]); entry `(j,k)` is `H_{jk}`.
#' @param column_config Configuration list (see [build_column_config()]).
#' @param ancestry_names Ancestry names in model order.
#' @return Logical `m x (2^K - 1)` matrix whose `(j,s)` entry is `M_{js}`, with
#'   configuration labels as column names.
#' @export
build_config_mask <- function(variant_available, column_config, ancestry_names) {
  stopifnot(ncol(variant_available) == length(ancestry_names))
  sets <- config_ancestry_sets(column_config, ancestry_names)
  labs <- build_config_labels(column_config, ancestry_names)

  mask <- vapply(sets, function(anc) {
    # Variant j is feasible for configuration s only if H_{jk}=1 for every
    # ancestry k selected by s.
    sub <- variant_available[, anc, drop = FALSE]
    apply(sub, 1L, all)
  }, logical(nrow(variant_available)))

  mask <- matrix(as.logical(mask), nrow = nrow(variant_available),
                 dimnames = list(rownames(variant_available), labs))
  mask
}


#' Restricted SNP index sets, one per configuration s
#'
#' @param config_mask Output of [build_config_mask()].
#' @return A named list; the element for configuration `s` holds the 1-based row
#'   indices of variants satisfying `M_{js}=1`.
#' @export
build_config_snp_index <- function(config_mask) {
  out <- lapply(seq_len(ncol(config_mask)), function(c) which(config_mask[, c]))
  names(out) <- colnames(config_mask)
  out
}


#' Availability-aware configuration probabilities phi^H_{s,j}
#'
#' For baseline configuration probabilities `phi_s`, this function computes
#' `phi^H_{s,j} = phi_s M_{js} / sum_{s'} phi_{s'} M_{js'}`.  The returned
#' software object is named `Qstar`, so `Qstar[j,s]` stores `phi^H_{s,j}`.
#'
#' @param config_mask Logical `m x (2^K - 1)` matrix storing `M_{js}`.
#' @param ancestry_weight Numeric length-`(2^K - 1)` vector of baseline
#'   configuration probabilities `phi_s`, in the same order as the mask columns.
#' @param tol Tolerance for the row-sum assertion.
#' @return Numeric `m x (2^K - 1)` row-stochastic matrix `Qstar` storing
#'   `phi^H_{s,j}`; entries with `M_{js}=0` are exactly zero.
#' @export
build_config_prior <- function(config_mask, ancestry_weight, tol = 1e-10) {
  if (!is.matrix(config_mask)) {
    stop("config_mask must be a matrix.", call. = FALSE)
  }
  if (!nrow(config_mask) || !ncol(config_mask)) {
    stop("config_mask must have at least one row and one column.", call. = FALSE)
  }
  if (!is.logical(config_mask) && !is.numeric(config_mask)) {
    stop("config_mask must contain only logical or numeric 0/1 values.",
         call. = FALSE)
  }
  if (anyNA(config_mask)) {
    stop("config_mask contains NA; feasibility must be TRUE or FALSE for every cell.",
         call. = FALSE)
  }
  if (is.numeric(config_mask) && any(!config_mask %in% c(0, 1))) {
    stop("numeric config_mask must contain only 0 or 1.", call. = FALSE)
  }
  storage.mode(config_mask) <- "logical"

  if (length(tol) != 1L || !is.numeric(tol) || !is.finite(tol) || tol < 0) {
    stop("tol must be one finite, non-negative numeric value.", call. = FALSE)
  }

  C <- ncol(config_mask)
  if (!is.numeric(ancestry_weight) || !is.null(dim(ancestry_weight))) {
    stop("ancestry_weight must be a numeric vector.", call. = FALSE)
  }
  if (length(ancestry_weight) != C) {
    stop(sprintf("ancestry_weight has length %d but there are %d configurations.",
                 length(ancestry_weight), C), call. = FALSE)
  }
  if (any(!is.finite(ancestry_weight)) || any(ancestry_weight < 0)) {
    stop("ancestry_weight must be finite and non-negative.", call. = FALSE)
  }

  weight_names <- names(ancestry_weight)
  if (!is.null(weight_names)) {
    if (anyNA(weight_names) || any(!nzchar(trimws(weight_names))) ||
        anyDuplicated(weight_names)) {
      stop("names(ancestry_weight) must be unique, non-missing, and non-empty.",
           call. = FALSE)
    }
    config_names <- colnames(config_mask)
    if (is.null(config_names)) {
      stop("A named ancestry_weight requires configuration names on config_mask.",
           call. = FALSE)
    }
    if (!setequal(weight_names, config_names)) {
      stop("names(ancestry_weight) do not match colnames(config_mask).",
           call. = FALSE)
    }
    ancestry_weight <- ancestry_weight[config_names]
  }

  num <- sweep(config_mask * 1.0, 2L, ancestry_weight, `*`)
  den <- rowSums(num)

  bad <- which(den <= 0)
  if (length(bad)) {
    ids <- rownames(config_mask)
    lab <- if (is.null(ids)) as.character(bad) else ids[bad]
    stop(sprintf(
      "%d variant(s) have no feasible ancestry configuration with positive weight ",
      length(bad)),
      "(first few: ", paste(utils::head(lab, 5L), collapse = ", "), "). ",
      "Such variants must be dropped during preprocessing; they must never reach ",
      "the inference code.", call. = FALSE)
  }

  Qstar <- num / den

  # sum_s phi^H_{s,j}=1 by construction -- assert, do not repair.
  rs <- rowSums(Qstar)
  if (max(abs(rs - 1)) > tol) {
    stop(sprintf("Qstar row sums deviate from 1 by up to %.3e; masking logic is wrong.",
                 max(abs(rs - 1))), call. = FALSE)
  }
  # Cells with M_{js}=0 must be *exactly* zero, not merely small.
  if (any(Qstar[!config_mask] != 0)) {
    stop("Qstar has non-zero mass on infeasible variant/configuration pairs.",
         call. = FALSE)
  }

  dimnames(Qstar) <- dimnames(config_mask)
  Qstar
}


#' Combine pi_{lj} with the availability-aware configuration probabilities
#'
#' Computes the joint prior `pi_{lj} phi^H_{s,j}`.  In software notation,
#' `snp_prior[j]` stores `pi_{lj}` and `Qstar[j,s]` stores `phi^H_{s,j}`.
#' This replaces the upstream
#' `kronecker(snp_prior, t(ancestry_weight), FUN = "*")`, which unconditionally
#' spread each SNP's mass over *all* configurations.
#'
#' No global renormalization is applied: when `sum_j pi_{lj}=1`, the result sums
#' to 1 automatically because `sum_s phi^H_{s,j}=1` for every variant `j`.
#'
#' @param snp_prior Numeric length-`m` SNP-level prior (`prior_weights` or the
#'   component-specific functional-annotation prior `pi_{lj}`).
#' @param Qstar Output of [build_config_prior()].
#' @param tol Tolerance for the total-mass assertion.
#' @param check_total If `TRUE` (default) assert that the result sums to
#'   `sum(snp_prior)`.
#' @return Numeric `m x (2^K - 1)` matrix storing the joint prior
#'   `pi_{lj} phi^H_{s,j}`.
#' @export
build_prior_matrix <- function(snp_prior, Qstar, tol = 1e-8, check_total = TRUE) {
  if (!is.matrix(Qstar) || !is.numeric(Qstar)) {
    stop("Qstar must be a numeric matrix.", call. = FALSE)
  }
  if (!nrow(Qstar) || !ncol(Qstar)) {
    stop("Qstar must have at least one row and one column.", call. = FALSE)
  }
  if (any(!is.finite(Qstar)) || any(Qstar < 0)) {
    stop("Qstar must contain only finite, non-negative values.", call. = FALSE)
  }
  if (length(tol) != 1L || !is.numeric(tol) || !is.finite(tol) || tol < 0) {
    stop("tol must be one finite, non-negative numeric value.", call. = FALSE)
  }
  if (length(check_total) != 1L || !is.logical(check_total) || is.na(check_total)) {
    stop("check_total must be TRUE or FALSE.", call. = FALSE)
  }

  qstar_row_sums <- rowSums(Qstar)
  max_row_error <- max(abs(qstar_row_sums - 1))
  if (max_row_error > tol) {
    stop(sprintf("Qstar is not row-stochastic; maximum row-sum error is %.3e.",
                 max_row_error), call. = FALSE)
  }

  if (!is.numeric(snp_prior) || !is.null(dim(snp_prior))) {
    stop("snp_prior must be a numeric vector.", call. = FALSE)
  }
  if (length(snp_prior) != nrow(Qstar)) {
    stop(sprintf("snp_prior has length %d but Qstar has %d rows.",
                 length(snp_prior), nrow(Qstar)), call. = FALSE)
  }
  prior_names <- names(snp_prior)
  if (!is.null(prior_names)) {
    if (anyNA(prior_names) || any(!nzchar(trimws(prior_names))) ||
        anyDuplicated(prior_names)) {
      stop("names(snp_prior) must be unique, non-missing, and non-empty.",
           call. = FALSE)
    }
    qstar_names <- rownames(Qstar)
    if (is.null(qstar_names)) {
      stop("A named snp_prior requires variant row names on Qstar.", call. = FALSE)
    }
    if (!setequal(prior_names, qstar_names)) {
      stop("names(snp_prior) do not match rownames(Qstar).", call. = FALSE)
    }
    snp_prior <- snp_prior[qstar_names]
  }
  snp_prior <- as.numeric(snp_prior)
  if (any(!is.finite(snp_prior)) || any(snp_prior < 0)) {
    stop("snp_prior must be finite and non-negative.", call. = FALSE)
  }
  if (sum(snp_prior) <= 0) {
    stop("snp_prior must contain positive total prior mass.", call. = FALSE)
  }

  out <- Qstar * snp_prior          # recycles snp_prior down the columns

  if (check_total) {
    got <- sum(out)
    want <- sum(snp_prior)
    if (abs(got - want) > tol * max(1, want)) {
      stop(sprintf(paste0("Joint prior mass is %.12f but should equal sum(snp_prior) = %.12f. ",
                          "Qstar is not row-stochastic -- refusing to renormalize, ",
                          "because that would hide the underlying bug."),
                   got, want), call. = FALSE)
    }
  }

  dimnames(out) <- dimnames(Qstar)
  out
}


#' Default INFUSE baseline configuration probabilities phi_s
#'
#' Reproduces the upstream default: single-ancestry configurations receive
#' `base_fac_ratio` times more weight than configurations one size larger.  For
#' `K = 2`, this gives `(phi_10, phi_01, phi_11) = (3/7, 3/7, 1/7)`.
#'
#' @param n_ancestry Total number of ancestries, `K`; `k` indexes ancestry.
#' @param base_fac_ratio Ratio between successive configuration sizes (default 3).
#' @export
default_ancestry_weight <- function(n_ancestry, base_fac_ratio = 3) {
  base_fac <- 1 / Reduce("+", lapply(seq_len(n_ancestry), function(x) {
    choose(n_ancestry, x) * base_fac_ratio^(n_ancestry - x)
  }))
  unlist(lapply(seq_len(n_ancestry), function(x) {
    rep(base_fac * base_fac_ratio^(n_ancestry - x), choose(n_ancestry, x))
  }))
}


#' Assemble the full availability context used throughout the model
#'
#' Convenience wrapper returning every availability-derived object in one place,
#' so that the core routine and the audit output share exactly one source of
#' truth.
#'
#' @param variant_available `NULL` or the `m x K` availability matrix `H`.
#' @param n_snp Number of variants in the union, `m`.
#' @param ancestry_names Ancestry names in model order.
#' @param ancestry_weight Baseline configuration probabilities `phi_s` (length
#'   `2^K - 1`).
#' @param variant_ids Optional variant identifiers.
#' @return A list with `variant_available`, `column_config`, `ancestry_names`,
#'   `config_names`, `config_ancestry`, `config_mask`, `config_snp_index`,
#'   `Qstar` (which stores `phi^H_{s,j}`) and `n_config`.
#' @export
build_availability_context <- function(variant_available, n_snp, ancestry_names,
                                       ancestry_weight = NULL, variant_ids = NULL) {
  K <- length(ancestry_names)
  column_config <- build_column_config(K)
  if (is.null(ancestry_weight)) ancestry_weight <- default_ancestry_weight(K)

  va <- normalize_variant_available(variant_available, n_snp, ancestry_names, variant_ids)

  # Prefer explicit IDs, then any row names carried by the availability matrix,
  # and only fall back to positional indices when neither is available -- error
  # messages naming a row number are far less useful than ones naming a variant.
  if (is.null(variant_ids)) variant_ids <- rownames(va)

  dead <- which(!apply(va, 1L, any))
  if (length(dead)) {
    lab <- if (is.null(variant_ids)) as.character(dead) else variant_ids[dead]
    stop(sprintf("%d variant(s) are unavailable in every ancestry (first few: %s). ",
                 length(dead), paste(utils::head(lab, 5L), collapse = ", ")),
         "Drop them in preprocessing (see prepare_infuse_union()) before fitting.",
         call. = FALSE)
  }

  mask  <- build_config_mask(va, column_config, ancestry_names)
  Qstar <- build_config_prior(mask, ancestry_weight)
  if (!is.null(names(ancestry_weight))) {
    ancestry_weight <- ancestry_weight[colnames(mask)]
  }

  list(
    variant_available = va,
    column_config     = column_config,
    ancestry_names    = ancestry_names,
    ancestry_weight   = ancestry_weight,
    config_names      = build_config_labels(column_config, ancestry_names),
    config_ancestry   = config_ancestry_sets(column_config, ancestry_names),
    config_mask       = mask,
    config_snp_index  = build_config_snp_index(mask),
    Qstar             = Qstar,
    n_config          = length(column_config),
    all_available     = all(va)
  )
}


#' Validate that every ancestry's rows describe the same variants, in the same order
#'
#' `infuse_core_avail()` indexes every ancestry by the same row number, so a
#' permuted or reversed variant order in one ancestry silently pairs that
#' ancestry's summary statistics with the wrong variants.  That corruption does not
#' always move the PIPs -- a strong ancestry can mask it entirely -- but it does
#' damage the prior covariance (observed: `V` collapsing to a zero variance for the
#' misaligned ancestry).  Dimension checks alone cannot catch it, so identity and
#' order are checked explicitly here.
#'
#' @param summary_stat_list,R_mat_list,annot The model inputs.
#' @param ancestry_names Ancestry names in model order.
#' @param n_snp Expected number of variants.
#' @param variant_ids Reference union-variant order, taken from the first
#'   preprocessed summary-statistics element.
#' @return `TRUE` invisibly; errors otherwise.
#' @export
validate_row_alignment <- function(summary_stat_list, R_mat_list, annot,
                                   ancestry_names, n_snp, variant_ids) {
  variant_ids <- as.character(variant_ids)
  if (anyDuplicated(variant_ids)) {
    stop("summary_stat_list[[1]]$SNP contains duplicated variant identifiers.",
         call. = FALSE)
  }

  first_mismatch <- function(a, b) {
    i <- which(a != b)
    if (!length(i)) return("")
    sprintf(" First mismatch at row %d: '%s' vs '%s'.", i[1], a[i[1]], b[i[1]])
  }

  for (k in seq_along(ancestry_names)) {
    nm <- ancestry_names[k]

    ## ---- summary statistics -------------------------------------------------
    ss <- summary_stat_list[[k]]
    if (nrow(ss) != n_snp) {
      stop(sprintf("summary_stat_list[['%s']] has %d rows; expected %d (the variant union).",
                   nm, nrow(ss), n_snp), call. = FALSE)
    }
    if (!"SNP" %in% names(ss)) {
      stop(sprintf("summary_stat_list[['%s']] has no SNP column, so row alignment cannot be verified.",
                   nm), call. = FALSE)
    }
    snp_k <- as.character(ss$SNP)
    if (!identical(snp_k, variant_ids)) {
      extra <- if (setequal(snp_k, variant_ids))
        " The variant sets match, so the rows are merely in a different order -- reorder them to ancestry 1."
      else
        sprintf(" %d variant(s) differ between the two sets.",
                length(union(setdiff(snp_k, variant_ids), setdiff(variant_ids, snp_k))))
      stop(sprintf("summary_stat_list[['%s']]$SNP is not identical to summary_stat_list[['%s']]$SNP.",
                   nm, ancestry_names[1]), extra, first_mismatch(snp_k, variant_ids),
           call. = FALSE)
    }

    ## ---- LD ------------------------------------------------------------------
    Rk <- R_mat_list[[k]]
    if (nrow(Rk) != n_snp || ncol(Rk) != n_snp) {
      stop(sprintf("R_mat_list[['%s']] is %d x %d; expected %d x %d.",
                   nm, nrow(Rk), ncol(Rk), n_snp, n_snp), call. = FALSE)
    }
    for (which_dn in c("colnames", "rownames")) {
      dn <- if (which_dn == "colnames") colnames(Rk) else rownames(Rk)
      if (is.null(dn) || !any(nzchar(dn))) next        # unnamed LD: nothing to check
      dn <- as.character(dn)
      if (!identical(dn, variant_ids)) {
        extra <- if (setequal(dn, variant_ids))
          " The variant sets match, so the LD is merely in a different order -- reorder it to the summary statistics."
        else ""
        stop(sprintf("%s(R_mat_list[['%s']]) do not match summary_stat_list[['%s']]$SNP.",
                     which_dn, nm, ancestry_names[1]), extra,
             first_mismatch(dn, variant_ids), call. = FALSE)
      }
    }
  }

  ## ---- annotations -----------------------------------------------------------
  if (!is.null(annot)) {
    A <- as.matrix(annot)
    if (nrow(A) != n_snp) {
      stop(sprintf(paste0("annot has %d rows but the variant union has %d variants; ",
                          "align annotations to the union in preprocessing ",
                          "(see align_annotations_to_union())."),
                   nrow(A), n_snp), call. = FALSE)
    }
    rn <- rownames(A)
    if (!is.null(rn) && any(nzchar(rn)) && !identical(as.character(rn), variant_ids)) {
      extra <- if (setequal(as.character(rn), variant_ids))
        " The variant sets match, so the annotations are merely in a different order."
      else ""
      stop("rownames(annot) do not match summary_stat_list[[1]]$SNP.", extra,
           first_mismatch(as.character(rn), variant_ids), call. = FALSE)
    }
  }

  invisible(TRUE)
}
