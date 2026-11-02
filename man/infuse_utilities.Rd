\name{infuse_utilities}
\alias{compute_b1b2}
\alias{compute_softmax}
\alias{compute_softmax_avail}
\alias{estimate_prior_annot}
\alias{get_purity}
\alias{in_CS}
\alias{in_CS_x}
\alias{kriging_rss_joint_ori}
\alias{estimate_s_rss_joint_separate}
\alias{kriging_rss_joint_separate_s}
\alias{llk_prior_causal_cpp}
\alias{loglik_cpp_R6}
\alias{meSuSie_get_cs}
\alias{meSusie_get_pip}
\alias{meSusie_get_pip_config}
\alias{meSusie_get_pip_either}
\alias{meSusie_plot_pip}
\alias{meSuSieData}
\alias{meSuSieObject}
\alias{multivariate_regression}
\alias{mvlmm_reg}
\alias{n_in_CS_x}
\alias{optimize_llk}
\alias{optimize_prior_variance_avail}
\alias{pre_optim}
\alias{ser_from_V_avail}
\alias{single_effect_regression_annot_avail}
\alias{single_effect_regression_avail}
\alias{test_ELBO}
\alias{test_run_loglik_cpp}
\alias{test_run_loglik_cpp_avail}
\alias{test_run_mvlmm_reg}
\alias{uni_reg}
\alias{vec_to_cov}
\title{Supporting INFUSE model and diagnostic utilities}
\description{Low-level covariance, regression, posterior, credible-set, and LD diagnostic utilities retained for compatibility. For routine analyses, use \code{infuse_core}, \code{prepare_infuse_union}, or \code{kriging_rss_joint}.}
\details{
\code{estimate_s_rss_joint_separate} and
\code{kriging_rss_joint_separate_s} are for benchmarking only. They reproduce
the revision sensitivity analysis with ancestry-specific regularization
parameters while retaining joint likelihood-ratio scoring and flagging rules.
Use \code{\link{kriging_rss_joint}} for routine INFUSE LD-discrepancy detection.
This separate-s joint benchmark is distinct from the manuscript's Separate Test,
which applies \code{susieR::kriging_rss()} independently to each ancestry.
}
\keyword{internal}
