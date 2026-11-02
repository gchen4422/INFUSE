# -*- coding: utf-8 -*-

#' Estimate joint shrinkage parameter for two-ancestry LD-discrepancy test
#'
#' Estimates a single regularization parameter \eqn{s \in [0,1]} by maximizing
#' the joint null log-likelihood across both ancestries simultaneously. The
#' regularized LD matrix is \eqn{(1-s)\mathbf{R} + s\mathbf{I}}, where
#' \eqn{\mathbf{R} = \mathrm{blockdiag}(R_1, R_2)}.
#'
#' @param z1 Numeric vector of z-scores from ancestry 1.
#' @param z2 Numeric vector of z-scores from ancestry 2.
#' @param R1 Symmetric LD correlation matrix for ancestry 1 (\eqn{p_1 \times p_1}).
#' @param R2 Symmetric LD correlation matrix for ancestry 2 (\eqn{p_2 \times p_2}).
#' @param n1 Integer sample size for ancestry 1. Used for finite-sample
#'   z-score shrinkage; ignored if \code{NULL}.
#' @param n2 Integer sample size for ancestry 2. Used for finite-sample
#'   z-score shrinkage; ignored if \code{NULL}.
#' @param r_tol Numeric tolerance for flooring small negative eigenvalues to
#'   enforce positive semi-definiteness (default \code{1e-8}).
#'
#' @return A scalar \eqn{\hat{s} \in [0,1]}: the estimated regularization
#'   parameter.
#'
#' @keywords internal
estimate_s_rss_joint <- function(z1, z2,
                                 R1, R2,
                                 n1 = NULL, n2 = NULL,
                                 r_tol = 1e-8) {
  #-- 1) finite-sample shrinkage on each z
  shrink <- function(z, n) {
    if (!is.null(n)) {
      if (n <= 1) stop("n must be >1")
      sigma2 <- (n-1)/(z^2 + n - 2)
      z <- sqrt(sigma2) * z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)

  #-- 2) stack z’s and build block-diag R
  z_joint <- c(z1, z2)
  # use Matrix::bdiag for efficiency if you like
  R_joint <- as.matrix(Matrix::bdiag(R1, R2))

  #-- 3) eigen-decomp (cached if already present)
  if (is.null(attr(R_joint, "eigen")))
    attr(R_joint, "eigen") <- eigen(R_joint, symmetric = TRUE)
  eigJ <- attr(R_joint, "eigen")
  # floor tiny negatives
  if (any(eigJ$values < -r_tol))
    warning("R_joint not PSD: flooring small negative eigenvalues")
  d <- eigJ$values
  d[d < r_tol] <- 0

  #-- 4) define joint neg-log-likelihood in s
  ztv <- crossprod(z_joint, eigJ$vectors)  # V^T z
  neglog <- function(s) {
    tmp <- (1 - s)*d + s
    0.5 * sum(log(tmp)) +
      0.5 * sum((ztv^2)/tmp)
  }

  #-- 5) optimize on [0,1]
  out <- optim(0.5, fn = neglog, method = "Brent", lower = 0, upper = 1)
  s_hat <- out$par
  s_hat
}


#' Joint LD-discrepancy test for two ancestries
#'
#' @description
#' A novel multi-ancestry extension of the \code{kriging_rss} diagnostic
#' (Zou et al., 2022) that jointly tests for LD mismatches and allele
#' flips across two ancestries simultaneously. For each variant, the
#' conditional distribution of its z-score given all others is computed under
#' the stacked model \eqn{(z_1, z_2) \sim N(0,\, (1-s)R + sI)}, where
#' \eqn{R = \mathrm{blockdiag}(R_1, R_2)} and the regularization parameter
#' \eqn{s} is estimated jointly by maximum likelihood via
#' \code{\link{estimate_s_rss_joint}}. A mixture-of-normals model is then fit
#' to the standardized residuals, and a log-likelihood ratio (logLR) flags
#' variants whose z-scores are inconsistent with LD — indicative of allele
#' flips or other reference-panel mismatches. The joint logLR combines
#' evidence across both ancestries, increasing sensitivity over single-ancestry
#' diagnostics.
#'
#' @param z1 Numeric vector of z-scores from ancestry 1 (\eqn{p} variants).
#' @param z2 Numeric vector of z-scores from ancestry 2 (\eqn{p} variants,
#'   same order as \code{z1}).
#' @param R1 Symmetric \eqn{p \times p} LD correlation matrix for ancestry 1.
#' @param R2 Symmetric \eqn{p \times p} LD correlation matrix for ancestry 2.
#' @param n1 Integer sample size for ancestry 1. Used for finite-sample
#'   z-score shrinkage; ignored if \code{NULL}.
#' @param n2 Integer sample size for ancestry 2. Used for finite-sample
#'   z-score shrinkage; ignored if \code{NULL}.
#' @param r_tol Numeric tolerance for flooring small negative eigenvalues
#'   (default \code{1e-8}).
#' @param s Numeric regularization parameter in \eqn{[0,1]}. If \code{NULL}
#'   (default), estimated jointly from the data via
#'   \code{\link{estimate_s_rss_joint}}.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{plot}}{A \code{ggplot2} scatter of observed vs. expected
#'     z-scores, colored by ancestry. Variants flagged as potential allele
#'     flips (logLR > 2 and |z_std_diff| > 4) are highlighted in red.}
#'   \item{\code{conditional_dist}}{A data frame with one row per variant
#'     (stacked ancestry 1 then ancestry 2) and columns:
#'     \code{ancestry}, \code{z}, \code{condmean}, \code{condvar},
#'     \code{z_std_diff}, \code{logLR} (per-variant), \code{logLR_joint}
#'     (sum of logLR across both ancestries for the same variant).}
#' }
#'
#' @references
#' Zou, Y., Carbonetto, P., Wang, G. & Stephens, M. (2022).
#' Fine-mapping from summary data with the "Sum of Single Effects" model.
#' \emph{PLOS Genetics}, 18(7), e1010299.
#' \doi{10.1371/journal.pgen.1010299}
#'
#' @importFrom Matrix bdiag
#' @importFrom stats optim dnorm
#' @importFrom ggplot2 ggplot aes geom_point geom_abline theme_bw labs
#' @importFrom mixsqp mixsqp
kriging_rss_joint <- function(z1, z2, R1, R2,
                                  n1 = NULL, n2 = NULL,
                                  r_tol = 1e-8,
                                  s = NULL) {
  #-- 1) finite‐sample shrinkage & NA‐impute
  shrink <- function(z, n) {
    if (!is.null(n)) {
      if (n <= 1) stop("n must be > 1")
      sigma2 <- (n - 1)/(z^2 + n - 2)
      z <- sqrt(sigma2) * z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)

  #-- 2) stack z and block‐diag LD
  z  <- c(z1, z2)
  R  <- as.matrix(Matrix::bdiag(R1, R2))

  #-- 3) estimate s if missing
  if (is.null(s)) {
    s <- estimate_s_rss_joint(z1, z2, R1, R2, n1, n2, r_tol)
  }
  force(s)
  if (s < 0) stop("s must be non‐negative")
  if (s > 1) {
    warning("s > 1; clamping to 0.8")
    s <- 0.8
  }

  #-- 4) eigen‐decompose & enforce PSD
  if (is.null(attr(R, "eigen")))
    attr(R, "eigen") <- eigen(R, symmetric = TRUE)
  eig <- attr(R, "eigen")
  d   <- eig$values
  if (any(d < -r_tol))
    warning("R not PSD: flooring small negative eigenvalues")
  d[d < r_tol] <- 0

  #-- 5) build precision matrix Ω = [(1-s)R + sI]^{-1}
  dinv      <- 1/((1 - s)*d + s)
  dinv[is.infinite(dinv)] <- 0 ## Here we are using the psedo-inverse (if d_i = 0, then 1/d_i = 0, if d_i > 0, then 1/d_i)
  precision <- eig$vectors %*% (t(eig$vectors) * dinv)

  #-- 6) compute conditional means & variances
  p_total <- length(z)
  condmean <- numeric(p_total)
  condvar  <- numeric(p_total)
  for (i in seq_len(p_total)) {
    condmean[i] <- -(1/precision[i,i]) * precision[i,-i] %*% z[-i]
    condvar[i]  <- 1/precision[i,i]
  }
  z_std_diff <- (z - condmean)/sqrt(condvar)

  #-- 7) fit mixture‐of‐normals on residuals
  a_min <- 0.8
  a_max <- if (max(z_std_diff^2) < 1) 2 else 2*sqrt(max(z_std_diff^2))
  npoint <- ceiling(log2(a_max/a_min)/log2(1.05))
  a_grid <- 1.05^(seq(-npoint, 0)) * a_max

  sd_mtx      <- outer(sqrt(condvar), a_grid)
  matrix_llik <- stats::dnorm(z - condmean, sd = sd_mtx, log = TRUE)
  lfactors    <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  w           <- mixsqp::mixsqp(matrix_llik, log = TRUE,
                                control = list(verbose = FALSE))$x

  #-- 8) compute log‐likelihood ratios for allele‐flip
  logl0mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  # rebuild for flipped sign
  matrix_llik <- stats::dnorm(z + condmean, sd = sd_mtx, log = TRUE)
  lfactors    <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  logl1mix    <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  logLRmix    <- logl1mix - logl0mix


  # 7) split into LR1, LR2, form joint
  p1 <- length(z1)
  logLR1       <- logLRmix[1:p1]
  logLR2       <- logLRmix[(p1+1):length(z)]
  logLR_joint  <- logLR1 + logLR2


  #-- 9) assemble results & plot
  ancestry <- rep(c("anc1","anc2"), times = c(length(z1), length(z2)))
  res      <- data.frame(
    ancestry   = ancestry,
    z          = z,
    condmean   = condmean,
    condvar    = condvar,
    z_std_diff = z_std_diff,
    logLR      = logLRmix,
    logLR_joint= logLR_joint
  )

  p <- ggplot2::ggplot(res, ggplot2::aes(x = condmean, y = z, color = ancestry)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Expected z", y = "Observed z")

  # — now flag per‐ancestry
  p1 <- length(z1)
  # ancestry1 positions 1…p1
  idx_flag1 <- which(logLRmix[1:p1] > 2 & abs(z_std_diff[1:(p1)])>4)

  # ancestry2 positions (p1+1)…end
  idx_flag2  <- which(logLRmix[(p1+1):length(z)] > 2 &
                        abs(z_std_diff[(p1+1):length(z)]) > 4)
  idx_flag_joint  <- which(logLR_joint > 2 & abs(z_std_diff) > 4)


  # highlight in red
  if (length(idx_flag1) > 0) {
    p <- p + ggplot2::geom_point(
      data = res[idx_flag1, ],
      ggplot2::aes(x = condmean, y = z),
      color = "red"
    )
  }
  if (length(idx_flag2) > 0) {
    p <- p + ggplot2::geom_point(
      data = res[p1+idx_flag2, ],
      ggplot2::aes(x = condmean, y = z),
      color = "red"
    )
  }
  if (length(idx_flag_joint) > 0) {
    p <- p + ggplot2::geom_point(
      data = res[idx_flag_joint, ],
      ggplot2::aes(x = condmean, y = z),
      color = "red"
    )
  }
  list(plot = p, conditional_dist = res)
}









kriging_rss_joint_ori <- function(z1, z2, R1, R2,
                              n1 = NULL, n2 = NULL,
                              r_tol = 1e-8,
                              s = NULL) {
  #-- 1) finite‐sample shrinkage & NA‐impute
  shrink <- function(z, n) {
    if (!is.null(n)) {
      if (n <= 1) stop("n must be > 1")
      sigma2 <- (n - 1)/(z^2 + n - 2)
      z <- sqrt(sigma2) * z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)

  #-- 2) stack z and block‐diag LD
  z  <- c(z1, z2)
  R  <- as.matrix(Matrix::bdiag(R1, R2))

  #-- 3) estimate s if missing
  if (is.null(s)) {
    s <- estimate_s_rss_joint(z1, z2, R1, R2, n1, n2, r_tol)
  }
  force(s)
  if (s < 0) stop("s must be non‐negative")
  if (s > 1) {
    warning("s > 1; clamping to 0.8")
    s <- 0.8
  }

  #-- 4) eigen‐decompose & enforce PSD
  if (is.null(attr(R, "eigen")))
    attr(R, "eigen") <- eigen(R, symmetric = TRUE)
  eig <- attr(R, "eigen")
  d   <- eig$values
  if (any(d < -r_tol))
    warning("R not PSD: flooring small negative eigenvalues")
  d[d < r_tol] <- 0

  #-- 5) build precision matrix Ω = [(1-s)R + sI]^{-1}
  dinv      <- 1/((1 - s)*d + s)
  dinv[is.infinite(dinv)] <- 0 ## Here we are using the psedo-inverse (if d_i = 0, then 1/d_i = 0, if d_i > 0, then 1/d_i)
  precision <- eig$vectors %*% (t(eig$vectors) * dinv)

  #-- 6) compute conditional means & variances
  p_total <- length(z)
  condmean <- numeric(p_total)
  condvar  <- numeric(p_total)
  for (i in seq_len(p_total)) {
    condmean[i] <- -(1/precision[i,i]) * precision[i,-i] %*% z[-i]
    condvar[i]  <- 1/precision[i,i]
  }
  z_std_diff <- (z - condmean)/sqrt(condvar)

  #-- 7) fit mixture‐of‐normals on residuals
  a_min <- 0.8
  a_max <- if (max(z_std_diff^2) < 1) 2 else 2*sqrt(max(z_std_diff^2))
  npoint <- ceiling(log2(a_max/a_min)/log2(1.05))
  a_grid <- 1.05^(seq(-npoint, 0)) * a_max

  sd_mtx      <- outer(sqrt(condvar), a_grid)
  matrix_llik <- stats::dnorm(z - condmean, sd = sd_mtx, log = TRUE)
  lfactors    <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  w           <- mixsqp::mixsqp(matrix_llik, log = TRUE,
                                control = list(verbose = FALSE))$x

  #-- 8) compute log‐likelihood ratios for allele‐flip
  logl0mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  # rebuild for flipped sign
  matrix_llik <- stats::dnorm(z + condmean, sd = sd_mtx, log = TRUE)
  lfactors    <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  logl1mix    <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  logLRmix    <- logl1mix - logl0mix

  #-- 9) assemble results & plot
  ancestry <- rep(c("anc1","anc2"), times = c(length(z1), length(z2)))
  res      <- data.frame(
    ancestry   = ancestry,
    z          = z,
    condmean   = condmean,
    condvar    = condvar,
    z_std_diff = z_std_diff,
    logLR      = logLRmix
  )

  p <- ggplot2::ggplot(res, ggplot2::aes(x = condmean, y = z, color = ancestry)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Expected z", y = "Observed z")

  # — now flag per‐ancestry
  p1 <- length(z1)
  # ancestry1 positions 1…p1
  idx_flag1 <- which(logLRmix[1:p1] > 2 & abs(z[1:p1]) > 2)
  # ancestry2 positions (p1+1)…end
  idx_flag2  <- which(logLRmix[(p1+1):length(z)] > 2 &
                       abs(z[(p1+1):length(z)]) > 2)

  # highlight in red
  if (length(idx_flag1) > 0) {
    p <- p + ggplot2::geom_point(
      data = res[idx_flag1, ],
      ggplot2::aes(x = condmean, y = z),
      color = "red"
    )
  }
  if (length(idx_flag2) > 0) {
    p <- p + ggplot2::geom_point(
      data = res[p1+idx_flag2, ],
      ggplot2::aes(x = condmean, y = z),
      color = "red"
    )
  }
  list(plot = p, conditional_dist = res)
}

# Benchmarking only: ancestry-specific shrinkage for the joint diagnostic.
# This implementation is the revision sensitivity-analysis comparator.
# Primary INFUSE analyses use kriging_rss_joint().

estimate_s_rss_joint_separate <- function(z1, z2,
                                          R1, R2,
                                          n1 = NULL, n2 = NULL,
                                          r_tol = 1e-8) {
  #-- 1) finite-sample shrinkage on each z (same as original)
  shrink <- function(z, n) {
    if (!is.null(n)) {
      if (n <= 1) stop("n must be >1")
      sigma2 <- (n - 1) / (z^2 + n - 2)
      z <- sqrt(sigma2) * z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)

  #-- 2) fit the same null-MLE objective separately in each ancestry
  estimate_one <- function(z, R) {
    if (is.null(attr(R, "eigen")))
      attr(R, "eigen") <- eigen(R, symmetric = TRUE)
    eig <- attr(R, "eigen")
    if (any(eig$values < -r_tol))
      warning("R not PSD: flooring small negative eigenvalues")
    d <- eig$values
    d[d < r_tol] <- 0
    ztv <- crossprod(z, eig$vectors)

    neglog <- function(s) {
      tmp <- (1 - s) * d + s
      0.5 * sum(log(tmp)) + 0.5 * sum((ztv^2) / tmp)
    }
    stats::optim(0.5, fn = neglog, method = "Brent",
                 lower = 0, upper = 1)$par
  }

  c(s1 = estimate_one(z1, R1), s2 = estimate_one(z2, R2))
}


#' Joint diagnostic with ancestry-specific shrinkage, for benchmarking only
#'
#' Reproduces the separate-s sensitivity analysis used in the revision.
#' Routine INFUSE analyses use kriging_rss_joint(). This benchmark retains
#' joint likelihood-ratio scoring and is distinct from the manuscript's
#' Separate Test, which runs single-ancestry diagnostics independently.
#' @noRd
kriging_rss_joint_separate_s <- function(z1, z2, R1, R2,
                                         n1 = NULL, n2 = NULL,
                                         r_tol = 1e-8,
                                         s1 = NULL, s2 = NULL) {
  if (length(z1) != length(z2))
    stop("z1 and z2 must contain the same ordered variants")
  if (nrow(R1) != length(z1) || ncol(R1) != length(z1) ||
      nrow(R2) != length(z2) || ncol(R2) != length(z2))
    stop("LD dimensions must match their ancestry-specific z-score vectors")

  #-- 1) finite-sample shrinkage & NA-impute (same as original)
  shrink <- function(z, n) {
    if (!is.null(n)) {
      if (n <= 1) stop("n must be > 1")
      sigma2 <- (n - 1) / (z^2 + n - 2)
      z <- sqrt(sigma2) * z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)

  #-- 2) stack z (same as original)
  z <- c(z1, z2)

  #-- 3) estimate s1 and s2 if missing. z has already received the same
  #      finite-sample shrinkage, so pass n=NULL to avoid shrinking twice.
  if (is.null(s1) || is.null(s2)) {
    s_hat <- estimate_s_rss_joint_separate(z1, z2, R1, R2,
                                            n1 = NULL, n2 = NULL,
                                            r_tol = r_tol)
    if (is.null(s1)) s1 <- unname(s_hat["s1"])
    if (is.null(s2)) s2 <- unname(s_hat["s2"])
  }
  if (s1 < 0 || s2 < 0) stop("s1 and s2 must be non-negative")
  if (s1 > 1 || s2 > 1) stop("s1 and s2 must not exceed 1")

  #-- 4/5) the only substantive model change: apply s1 and s2 to the
  #         eigenvalues of their corresponding unregularized LD blocks. This
  #         exactly preserves the original order of operations: floor the LD
  #         eigenvalues first, then add s. In particular, when 0 < s < r_tol,
  #         a null LD direction has precision 1/s; it must not be discarded by
  #         flooring the already-regularized covariance eigenvalue.
  precision_one <- function(R, s) {
    if (is.null(attr(R, "eigen")))
      attr(R, "eigen") <- eigen(R, symmetric = TRUE)
    eig <- attr(R, "eigen")
    d <- eig$values
    if (any(d < -r_tol))
      warning("R not PSD: flooring small negative eigenvalues")
    d[d < r_tol] <- 0
    dinv <- 1 / ((1 - s) * d + s)
    dinv[is.infinite(dinv)] <- 0
    eig$vectors %*% (t(eig$vectors) * dinv)
  }
  precision <- as.matrix(Matrix::bdiag(precision_one(R1, s1),
                                       precision_one(R2, s2)))

  #-- 6) conditional means and variances (same loop as original)
  p_total <- length(z)
  condmean <- numeric(p_total)
  condvar <- numeric(p_total)
  for (i in seq_len(p_total)) {
    condmean[i] <- -(1 / precision[i, i]) * precision[i, -i] %*% z[-i]
    condvar[i] <- 1 / precision[i, i]
  }
  z_std_diff <- (z - condmean) / sqrt(condvar)

  #-- 7) fit the same mixture of normal residual scales
  a_min <- 0.8
  a_max <- if (max(z_std_diff^2) < 1) 2 else 2 * sqrt(max(z_std_diff^2))
  npoint <- ceiling(log2(a_max / a_min) / log2(1.05))
  a_grid <- 1.05^(seq(-npoint, 0)) * a_max

  sd_mtx <- outer(sqrt(condvar), a_grid)
  matrix_llik <- stats::dnorm(z - condmean, sd = sd_mtx, log = TRUE)
  lfactors <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  w <- mixsqp::mixsqp(matrix_llik, log = TRUE,
                      control = list(verbose = FALSE))$x

  #-- 8) compute the same allele-flip log likelihood ratios
  logl0mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  matrix_llik <- stats::dnorm(z + condmean, sd = sd_mtx, log = TRUE)
  lfactors <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  logl1mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  logLRmix <- logl1mix - logl0mix

  #-- 9) split into LR1/LR2 and form the same summed joint LR
  p1 <- length(z1)
  logLR1 <- logLRmix[seq_len(p1)]
  logLR2 <- logLRmix[p1 + seq_len(p1)]
  logLR_joint <- logLR1 + logLR2

  ancestry <- rep(c("anc1", "anc2"), times = c(length(z1), length(z2)))
  res <- data.frame(
    ancestry = ancestry,
    z = z,
    condmean = condmean,
    condvar = condvar,
    z_std_diff = z_std_diff,
    logLR = logLRmix,
    logLR_joint = rep(logLR_joint, 2)
  )

  p <- ggplot2::ggplot(res, ggplot2::aes(x = condmean, y = z,
                                          color = ancestry)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Expected z", y = "Observed z")

  #-- 10) retain the original thresholds
  idx_flag1 <- which(logLR1 > 2 & abs(z_std_diff[seq_len(p1)]) > 4)
  idx_flag2 <- which(logLR2 > 2 &
                       abs(z_std_diff[p1 + seq_len(p1)]) > 4)
  # The original code implicitly recycles the length-p joint LR over the stacked
  # length-2p residual vector. Make the identical rule explicit.
  idx_flag_joint_stacked <- which(rep(logLR_joint, 2) > 2 &
                                    abs(z_std_diff) > 4)
  idx_flag_joint <- sort(unique(((idx_flag_joint_stacked - 1L) %% p1) + 1L))
  joint_flags <- sort(unique(c(idx_flag1, idx_flag2, idx_flag_joint)))

  if (length(idx_flag1) > 0) {
    p <- p + ggplot2::geom_point(data = res[idx_flag1, ],
                                 ggplot2::aes(x = condmean, y = z), color = "red")
  }
  if (length(idx_flag2) > 0) {
    p <- p + ggplot2::geom_point(data = res[p1 + idx_flag2, ],
                                 ggplot2::aes(x = condmean, y = z), color = "red")
  }
  if (length(idx_flag_joint_stacked) > 0) {
    p <- p + ggplot2::geom_point(data = res[idx_flag_joint_stacked, ],
                                 ggplot2::aes(x = condmean, y = z), color = "red")
  }

  list(
    plot = p,
    conditional_dist = res,
    s = c(s1 = s1, s2 = s2),
    flags = list(EUR = idx_flag1, AFR = idx_flag2,
                 joint_only = idx_flag_joint, joint = joint_flags)
  )
}
