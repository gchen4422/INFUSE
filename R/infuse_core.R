require(progress)
require(R6)
require(nloptr)
require(elasticnet)

#' Multi-ancestry fine-mapping with functional annotations
#'
#' Fits common or union variant sets using a single availability-aware engine.
#' With `variant_available = NULL`, all input variants are treated as available
#' in every ancestry and the model reduces to the original INFUSE formulation.
#' Summary-statistic and LD lists must have matching ancestry names.
#'
#' @param R_mat_list Named list of aligned LD correlation matrices.
#' @param summary_stat_list Named list of aligned summary-statistic data frames
#'   containing SNP, Beta, Se, Z and N.
#' @param L Maximum number of single effects.
#' @param residual_variance Optional residual variances.
#' @param prior_weights Optional SNP prior weights.
#' @param ancestry_weight Prior weights over nonempty ancestry configurations.
#' @param optim_method Prior-variance optimizer.
#' @param estimate_residual_variance Whether to estimate residual variance.
#' @param max_iter Maximum number of iterations.
#' @param cor_method Credible-set purity criterion.
#' @param cor_threshold Minimum credible-set purity.
#' @param annot Optional SNP-by-feature annotation matrix.
#' @param annot_method Annotation-prior estimation method.
#' @param est_annot_prior Annotation-prior update mode.
#' @param variant_available Optional logical variant-by-ancestry matrix.
#'   NULL means all variants are available in every ancestry.
#' @param verbose Whether to print progress.
#' @return An R6 object with PIPs, credible sets, posterior effect moments,
#'   ancestry configurations and availability metadata.
#' @export
infuse_core <- function(R_mat_list, summary_stat_list, L,
                        residual_variance = NULL, prior_weights = NULL,
                        ancestry_weight = NULL, optim_method = "optim",
                        estimate_residual_variance = FALSE, max_iter = 100,
                        cor_method = "min.abs.corr", cor_threshold = 0.5,
                        annot = NULL, annot_method = NULL,
                        est_annot_prior = "fixed", variant_available = NULL,
                        verbose = TRUE) {

  # Reject an invalid annotation-prior mode before any fitting.  Without
  # this guard, the annotation loop can reuse a stale SER_res when the
  # value is neither "fixed" nor "optim", silently returning a wrong fit.
  if(length(est_annot_prior)!=1L || !is.character(est_annot_prior) ||
     is.na(est_annot_prior) || !est_annot_prior %in% c("fixed","optim")){
    stop("est_annot_prior must be either 'fixed' or 'optim'.",call. = FALSE)
  }

  if(verbose){                                                                  # verbose guard
  cat("*************************************************************\n",
      "INFUSE: multi-ancestry fine-mapping\n",
      "Guangxin Chen and Boran Gao\n",
      "Purdue University\n",
      "https://github.com/gchen4422/INFUSE\n",
      "GNU General Public License\n",
      "*************************************************************", sep = "")
  }

  time_start<-Sys.time()
  if(verbose) cat("\n# Start data processing for sufficient statistics \n")
  meSuSieData_obj<-meSuSieData$new(R_mat_list,summary_stat_list)

  n_snp = nrow(summary_stat_list[[1]])
  n_ancestry = length(summary_stat_list)


  ancestry_names = resolve_ancestry_names(summary_stat_list,R_mat_list)
  variant_ids = as.character(summary_stat_list[[1]]$SNP)
  # Row alignment is validated, not assumed.  Everything downstream --
  # XtX.diag, Xty, betahat, the mask -- indexes all ancestries by the SAME
  # row number, so a permuted or reversed SNP order in one ancestry
  # silently pairs that ancestry's statistics with the wrong variants.  It
  # does not necessarily move the PIPs (a strong ancestry can mask it) but
  # it does corrupt the prior covariance, so it must be an error.
  validate_row_alignment(summary_stat_list,R_mat_list,annot,ancestry_names,n_snp,variant_ids)


  if(is.null(prior_weights)){
    prior_weights = rep(1/n_snp,n_snp)
  }
  if(is.null(ancestry_weight)){
    base_fac_ratio = 3
    base_fac = 1/Reduce("+",lapply(seq(1,n_ancestry),function(x)choose(n_ancestry,x)*base_fac_ratio^(n_ancestry-x)))
    ancestry_weight = unlist(lapply(seq(1,n_ancestry),function(x)rep(base_fac*base_fac_ratio^(n_ancestry-x),choose(n_ancestry,x))))
  }


  #          configuration probabilities phi^H_{s,j} (stored as Qstar),
  #          and the per-configuration feasible SNP index sets
  avail = build_availability_context(variant_available = variant_available,
                                     n_snp = n_snp,
                                     ancestry_names = ancestry_names,
                                     ancestry_weight = ancestry_weight,
                                     variant_ids = variant_ids)
  if(verbose){
    cat(sprintf("# %d variants in the union; %d available in all %d ancestries\n",
                n_snp,sum(rowSums(avail$variant_available)==n_ancestry),n_ancestry))
    cat("# feasible variants per configuration: ",
        paste(sprintf("%s=%d",avail$config_names,
                      vapply(avail$config_snp_index,length,integer(1))),collapse=", "),"\n",sep="")
  }


  # was: used_weights = kronecker(prior_weights, t(ancestry_weight), FUN = "*")
  # now: used_weights[j,s] = prior_weights[j] * Qstar[j,s], representing
  #      pi_{lj} phi^H_{s,j}. No global renormalization: for every j,
  #      sum_s phi^H_{s,j}=1 by construction, and that is asserted.
  used_weights = build_prior_matrix(prior_weights, avail$Qstar)

  if(is.null(residual_variance)){
    residual_variance = rep(1,n_ancestry)
  }

  if(verbose) cat("# Create INFUSE object \n")
  meSuSieObject_obj<-meSuSieObject$new(n_snp,L,n_ancestry,residual_variance,used_weights,optim_method,estimate_residual_variance,max_iter)
  meSuSieObject_obj$avail <- avail                                               # the SER functions read this

  if(verbose) cat("# Start data analysis \n")
  pb = progress::progress_bar$new(format =paste(" :elapsed"),clear = TRUE,total = max_iter,show_after = 0)

  n_iter = 0
  for (iter in 1:max_iter) {
    comp_residual<-meSuSieObject_obj$compute_residual(meSuSieData_obj,meSuSieObject_obj)

    for (l_index in seq(1,L,1)) {
      comp_residual = comp_residual + meSuSieObject_obj$Xr[[l_index]]
      SER_res<-single_effect_regression_avail(comp_residual,meSuSieData_obj$XtX.diag, meSuSieObject_obj,l_index)   # _avail
      meSuSieObject_obj$par_update(SER_res,l_index)
      meSuSieObject_obj$compute_KL(SER_res,meSuSieObject_obj$compute_SER_posterior_loglik(meSuSieData_obj,comp_residual,SER_res$b1b2),l_index)
      meSuSieObject_obj$compute_Xr(meSuSieData_obj,SER_res$b1b2$EB1,l_index)
      comp_residual = comp_residual - meSuSieObject_obj$Xr[[l_index]]
    }

    if(verbose) pb$tick(tokens = list(iteration = iter))                         # verbose guard

    updated_sigma2 = meSuSieObject_obj$update_residual_variance(meSuSieData_obj,iter)
    if((meSuSieObject_obj$ELBO[iter+1] - meSuSieObject_obj$ELBO[iter])<0.001){
      break
    }
    if(meSuSieObject_obj$estimate_residual_variance==TRUE){
      meSuSieObject_obj$sigma2 =  updated_sigma2
    }
    n_iter = n_iter + 1
  }

  prior_ELBO <- meSuSieObject_obj$ELBO[iter]
  meSuSieObject_obj$ELBO<-rep(NA,max_iter)
  meSuSieObject_obj$ELBO[1] = prior_ELBO

  # =====================================================================
  # MASTER ANNOTATION WRAPPER
  # Safely bypasses all annotation logic if no annotations are provided
  # =====================================================================
  if (!is.null(annot) && !is.null(annot_method)) {

    # annotations must already be aligned to the union (see
    # prepare_infuse_union / align_annotations_to_union)
    if(nrow(as.matrix(annot))!=n_snp){
      stop(sprintf("annot has %d rows but the variant union has %d variants; align annotations to the union in preprocessing.",
                   nrow(as.matrix(annot)),n_snp),call. = FALSE)
    }

    if (annot_method == "mlk"){
      # AUTO-SPCA BYPASS: Only run SPCA if there are 15 or more columns
      if(ncol(as.matrix(annot)) >= 15) {
        sds <- apply(as.matrix(annot), 2, sd)
        annotation_nonconst <- as.matrix(annot)[, sds > 0]
        annotation_nonconst_scaled <- scale(annotation_nonconst)
        scaling_scales <- attr(annotation_nonconst_scaled, "scaled:scale")

        k_values <- c(5,10,15,20)
        nonzero_values <- c(5,10)
        results_comb <- expand.grid(k = k_values, q = nonzero_values)
        loadings_list <- list()

        X <- as.matrix(annotation_nonconst_scaled)
        storage.mode(X) <- "double"
        G <- crossprod(X) / (nrow(X) - 1)

        for (i in 1:nrow(results_comb)) {
          k_val <- results_comb$k[i]
          q_val <- results_comb$q[i]
          spca_result <- elasticnet::spca(G, K = k_val, type = "Gram", sparse = "varnum", para = rep(q_val, k_val))
          loadings_list[[paste0("k", k_val, "_q", q_val)]] <- spca_result$loadings
          if(verbose) cat("Finished combination: k =", k_val, "and q =", q_val, "\n")
        }
        results_comb$BIC <- NA
        weights_mat = matrix(NA,nrow = nrow(spca_result$loadings) ,ncol = nrow(results_comb))
        rownames(weights_mat) = colnames(annotation_nonconst_scaled)
      }
    }

    updated_prior_list<-list()
    annot_weights_list<-list() # Initialize list to store the annotation weights
    annot_snp_prior<-list()                                                      # keep u[,l] for the audit output

    for (l_index in seq(1,L,1)) {
      input.response<-rowSums(meSuSieObject_obj$alpha[[l_index]])

      if(max(diag(meSuSieObject_obj$V[[l_index]]))<1e-9|max(input.response)>0.999){
        updated_prior <- rep(1/nrow(annot),nrow(annot))
        annot_weights_list[[l_index]] <- rep(0, ncol(as.matrix(annot)))
      }else{
        if(annot_method == "glmnet"){
          response.matrix<-matrix(c(1-input.response, input.response),length(input.response),2)
          try.index<-try(cv.logistic<-glmnet::cv.glmnet(x=as.matrix(annot),y=response.matrix, family='multinomial',alpha=0.5,type.measure = 'deviance'))
          if(class(try.index)[1]!='try-error'){
            cv.index<-which(cv.logistic$lambda==cv.logistic$lambda.min)
            glm.beta<-as.matrix(c(cv.logistic$glmnet.fit$a0[cv.index],cv.logistic$glmnet.fit$beta[[2]][,cv.index]))
            updated_prior<-exp(as.matrix(cbind(1,annot))%*%glm.beta)/sum(exp(as.matrix(cbind(1,annot))%*%glm.beta))

            # Save glmnet weights
            annot_weights_list[[l_index]] <- cv.logistic$glmnet.fit$beta[[2]][,cv.index]
          }else{
            updated_prior <- rep(1/nrow(annot),nrow(annot))
            annot_weights_list[[l_index]] <- rep(0, ncol(as.matrix(annot)))
          }
        }else if (annot_method == "mlk"){

          # DIRECT OPTIMIZATION (For 1-14 annotations)
          if(ncol(as.matrix(annot)) < 15) {
            inter_weight <- optimize_llk(as.matrix(annot), meSuSieObject_obj$alpha[[l_index]], rep(0, ncol(as.matrix(annot))))

            # Save Direct Weights
            names(inter_weight) <- colnames(as.matrix(annot))
            annot_weights_list[[l_index]] <- inter_weight

            z <- as.matrix(annot) %*% inter_weight
            z_stable <- z - max(z)
            updated_prior <- exp(z_stable) / sum(exp(z_stable))

          } else {
            # SPCA OPTIMIZATION (For 15+ annotations)
            for(i in 1:nrow(results_comb)) {
              k <- results_comb$k[i]
              nonzero <- results_comb$q[i]
              current_loadings <- loadings_list[[paste0("k", k, "_q", nonzero)]]

              annotation_spca <- annotation_nonconst_scaled %*% current_loadings
              inter_weight <- optimize_llk(as.matrix(annotation_spca), meSuSieObject_obj$alpha[[l_index]], rep(0, ncol(annotation_spca)))
              gamma <- inter_weight
              V_k <- current_loadings
              beta_std <- V_k %*% gamma
              beta_orig <- beta_std / scaling_scales
              weights_mat[,i] = beta_orig

              n_beta <- length(beta_orig)
              beta_RSS <- sum(beta_orig^2)
              results_comb$BIC[i] <- n_beta * log(beta_RSS / n_beta) + log(n_beta) * (k * nonzero)
            }
            optimal_params <- results_comb[which.min(results_comb$BIC), ]
            k_optim = optimal_params$k
            nonzero_optim = optimal_params$q
            beta_orig = weights_mat[,which.min(results_comb$BIC)]
            idx = match(rownames(weights_mat),colnames(as.matrix(annot)))

            beta_annot = rep(0,ncol(as.matrix(annot)))
            beta_annot[idx] = beta_orig

            # Save SPCA Weights
            names(beta_annot) <- colnames(as.matrix(annot))
            annot_weights_list[[l_index]] <- beta_annot

            z <- as.matrix(annot) %*% beta_annot
            z_stable <- z - max(z)
            updated_prior <- exp(z_stable) / sum(exp(z_stable))

            # PRINT TOP WEIGHTS FOR CAUSAL EFFECT 1 TO CONSOLE
            if (l_index == 1 && verbose) {
              cat("\n--- Top Annotations for Effect L=1 ---\n")
              sorted_weights <- sort(abs(beta_annot), decreasing = TRUE)
              print(sorted_weights[1:min(10, length(sorted_weights))])
            }
          }

        }else if (annot_method == "irwu"){
          # upstream refers to undefined `annot_file_subset`; use the
          # validated, row-aligned annotation matrix supplied to this fit.
          updated_prior<-drop(estimate_prior_annot(as.matrix(annot), rowSums(meSuSieObject_obj$alpha[[l_index]])))
          annot_weights_list[[l_index]] <- rep(0, ncol(as.matrix(annot))) # Placeholder if IRWU used
        }else{
          stop("Unknown annot_method: ",annot_method,call. = FALSE)               # fail loudly instead of leaving updated_prior undefined
        }
      }
      # was: updated_prior_list[[l_index]] = kronecker(updated_prior, t(ancestry_weight), FUN = "*")
      # now: updated_prior_list[[l]][j,s] = pi_{lj} phi^H_{s,j}, where
      #      updated_prior[j] stores the original INFUSE pi_{lj}
      annot_snp_prior[[l_index]] = as.numeric(updated_prior)
      updated_prior_list[[l_index]] = build_prior_matrix(as.numeric(updated_prior), avail$Qstar)
    }

    # SECOND LOOP: Update fine-mapping with the new functional priors
    n_iter = 0
    for (iter in 1:max_iter) {
      comp_residual<-meSuSieObject_obj$compute_residual(meSuSieData_obj,meSuSieObject_obj)
      for (l_index in seq(1,L,1)) {
        comp_residual = comp_residual + meSuSieObject_obj$Xr[[l_index]]
        if(est_annot_prior == "optim"){
          SER_res<-single_effect_regression_annot_avail(comp_residual,meSuSieData_obj$XtX.diag, meSuSieObject_obj,l_index,updated_prior_list[[l_index]],est_annot_prior)   # _avail
        }else if(est_annot_prior == "fixed"){
          V_input<-as.matrix(meSuSieObject_obj$V[[l_index]],ncol=n_ancestry,nrow=n_ancestry)
          SER_res<-single_effect_regression_annot_avail(comp_residual,meSuSieData_obj$XtX.diag, meSuSieObject_obj,l_index,updated_prior_list[[l_index]],est_annot_prior,prior_var = V_input) # _avail
        }
        meSuSieObject_obj$par_update(SER_res,l_index)
        meSuSieObject_obj$compute_KL(SER_res,meSuSieObject_obj$compute_SER_posterior_loglik(meSuSieData_obj,comp_residual,SER_res$b1b2),l_index)
        meSuSieObject_obj$compute_Xr(meSuSieData_obj,SER_res$b1b2$EB1,l_index)
        comp_residual = comp_residual - meSuSieObject_obj$Xr[[l_index]]
      }
      updated_sigma2 = meSuSieObject_obj$update_residual_variance(meSuSieData_obj,iter)
      if((meSuSieObject_obj$ELBO[iter+1] - meSuSieObject_obj$ELBO[iter])<0.001){
        break
      }
      if(meSuSieObject_obj$estimate_residual_variance==TRUE){
        meSuSieObject_obj$sigma2 =  updated_sigma2
      }
      n_iter = n_iter + 1
    }

    # Attach the final weights to the output object
    meSuSieObject_obj$annot_weights <- annot_weights_list
    avail$annot_snp_prior <- annot_snp_prior                                     # u[,l]
    avail$joint_prior     <- updated_prior_list                                  # pi[,,l]

  } # <-- END OF MASTER ANNOTATION WRAPPER

  if(verbose) cat("\n# Data analysis is done, and now generates result \n\n")
  ###Use function in Utility to output result
  meSuSieObject_obj$get_result(meSuSie_get_cs(meSuSieObject_obj,R_mat_list,cor_method=cor_method,cor_threshold=cor_threshold),meSusie_get_pip_either(meSuSieObject_obj),meSusie_get_pip_config(meSuSieObject_obj))


  avail$variant_ids  <- variant_ids
  avail$prior_matrix <- used_weights
  meSuSieObject_obj$availability <- avail
  meSuSieObject_obj$avail        <- avail
  if (!is.null(variant_available)) {
    names(meSuSieObject_obj$pip) <- variant_ids
    colnames(meSuSieObject_obj$pip_config) <- avail$config_names
  } else {
    # Preserve the original alpha naming for common-variant callers.
    rownames(meSuSieObject_obj$pip_config) <- NULL
    meSuSieObject_obj$alpha <- lapply(meSuSieObject_obj$alpha, function(a) {
      dimnames(a) <- list(NULL, c("init", rep("", ncol(a) - 1L)))
      a
    })
  }


  if(verbose) meSuSieObject_obj$mesusie_summary(meSuSieData_obj)

  time_end<-Sys.time()
  if(verbose) cat(c("\n# Total time used for the analysis:",paste0(round(as.numeric(difftime(time_end,time_start,units=c("mins"))),2)," mins\n")))
  return(meSuSieObject_obj)
}


meSuSieObject <- R6Class("meSuSieObject",public = list(
  initialize = function(p,L,N_ancestry, residual_variance,prior_weights,estimate_prior_method,estimate_residual_variance,max_iter){

    self$name_config = Reduce(append , lapply(seq(1,N_ancestry),function(x){
      poss_config = combn(1:N_ancestry,x)
      apply(poss_config,2,function(x)paste(x, collapse = '_'))
    }))


    self$column_config = Reduce(append,lapply(seq(1,N_ancestry),function(x){
      poss_config = combn(1:N_ancestry,x)
      lapply(1:ncol(poss_config), function(y)return(matrix(poss_config[,y])))
    }))
    self$alpha =rep(list(matrix(0,nrow = p,ncol = N_ancestry)),L)
    self$mu1 = rep(list(lapply(self$column_config,function(x){matrix(0,nrow = p,ncol = length(x))})),L)
    self$mu2 = rep(list(lapply(self$column_config,function(x){matrix(0,nrow = p,ncol = length(x))})),L)
    self$EB1 =rep(list(matrix(0,nrow = p,ncol = N_ancestry)),L)
    self$EB2 =rep(list(matrix(0,nrow = p,ncol = N_ancestry)),L)

    self$Xr     = rep(list(matrix(0,nrow = p,ncol=N_ancestry)),L)
    self$KL     = rep(as.numeric(NA),L)
    self$lbf    = rep(as.numeric(NA),L)
    self$lbf_variable = vector("list",L)
    self$ELBO = rep(NA,max_iter)
    self$ELBO[1] = -Inf
    self$sigma2 = residual_variance


    self$V      = rep(list(matrix(0,ncol = N_ancestry,nrow=N_ancestry)),L)
    self$pi     = prior_weights

    self$estimate_prior_method = estimate_prior_method
    self$estimate_residual_variance = estimate_residual_variance

    self$L = L
    self$nSNP = p
    self$nancestry = N_ancestry

    self$cs = list()
    self$pip = rep(as.numeric(NA),p)
    self$pip_config = matrix(as.numeric(NA),p,2^N_ancestry-1)
    self$avail = NULL                                                            # the only field added
  },

  compute_Xb = function (meSuSie_Data,b){
    lapply(1:length(meSuSie_Data$XtX_list),function(x){
      meSuSie_Data$XtX_list[[x]]%*%b[x,]
    })},

  compute_residual = function (meSuSie_Data,meSuSie_Obj) {
    residual<-Reduce(cbind,meSuSie_Data$Xty_list)-Reduce("+",meSuSie_Obj$Xr)
    return(residual)
  },
  compute_Xr = function(meSuSie_Data,b1b2,l_index){
    self$Xr[[l_index]] = Reduce(cbind,lapply(1:self$nancestry,function(ancestry_index)meSuSie_Data$XtX_list[[ancestry_index]]%*%b1b2[,ancestry_index]))
  },

  par_update = function(SER_res,l_index){
    self$alpha[[l_index]]<-SER_res$alpha
    self$mu1[[l_index]]<-SER_res$mu1_multi
    self$mu2[[l_index]]<-SER_res$mu2_multi
    self$lbf_variable[[l_index]]<-SER_res$lbf_multi
    self$V[[l_index]]<-SER_res$V
    self$EB1[[l_index]] =SER_res$b1b2$EB1
    self$EB2[[l_index]] =SER_res$b1b2$EB2

  },

  compute_SER_posterior_loglik = function(meSuSie_Data,comp_residual,b1b2){

    return(Reduce("+",lapply(1:self$nancestry,function(x){
      sum(-0.5/self$sigma2[[x]]*(-2*comp_residual[,x]*b1b2$EB1[,x]+meSuSie_Data$XtX.diag[[x]]*b1b2$EB2[,x]))
    })))

  },

  compute_KL = function(SER_res,value,l_index){

    self$KL[l_index] = -SER_res$loglik+value

  },

  update_residual_variance = function(meSuSie_Data,niter){


    B_1_ancestry = lapply(1:self$nancestry,function(y)Reduce(cbind,lapply(self$EB1,function(x)x[,y])))


    BXXB = lapply(1:self$nancestry,function(x){
      sum((t(B_1_ancestry[[x]])%*% meSuSie_Data$XtX_list[[x]])*t(B_1_ancestry[[x]]))
    })

    betabar = Reduce("+",self$EB1)

    BbarXXBbar = lapply(1:self$nancestry,function(x){
      sum(betabar[,x]*(meSuSie_Data$XtX_list[[x]]%*%betabar[,x]))
    })

    BbarXty = lapply(1:self$nancestry,function(x){
      2*sum(betabar[,x]*meSuSie_Data$Xty_list[[x]])
    })



    B_2_ancestry = lapply(1:self$nancestry,function(y)Reduce(cbind,lapply(self$EB2,function(x)x[,y])))

    XXB2 = lapply(1:self$nancestry,function(x){
      sum(meSuSie_Data$XtX.diag[[x]]*B_2_ancestry[[x]])
    })

    updated_sigma = lapply(1:self$nancestry,function(x){

      ((meSuSie_Data$yty_list[[x]]-BbarXty[[x]]+BbarXXBbar[[x]])+(XXB2[[x]]-BXXB[[x]]))/meSuSie_Data$N_list[[x]]
    })


    self$ELBO[niter+1] = Reduce(sum,lapply(1:self$nancestry,function(x){
      -meSuSie_Data$N_list[[x]]/2*log(2*pi*self$sigma2[x])-(1/(2*self$sigma2[x]))*((meSuSie_Data$yty_list[[x]]-BbarXty[[x]]+BbarXXBbar[[x]])+(XXB2[[x]]-BXXB[[x]]))
    }))-sum(self$KL)

    return( unlist(updated_sigma))
  },

  get_result = function(cs, pip, pip_config){
    self$cs = cs
    self$pip = pip
    self$pip_config = pip_config
  },
  mesusie_summary = function(meSuSie_Data){

    cat(c(paste0("Potential causal SNPs with PIP > 0.5: "),meSuSie_Data$Summary_Stat[[1]]$SNP[which(self$pip>0.5)],"\n\n"))
    cat("Credible sets for effects: \n")
    print(self$cs)
    cat("\n Use meSusie_plot_pip() for Mahattan and PIP Plot")

  }

),lock_objects = F
)

##############################################################
#  Assume var(y) = 1 and causal snps contribute negligible variance
#  therefore se(\beta) = var(y)*diag(xtx) => diag(xtx) = var(y)/se(beta)
#                          or
#  R^2 = z^2/(z^2+N-2) => sigma^2 = var(y)*(N-1)/(z^2+N-2)=>diag(xtx)=sigma^2/se(beta)^2
#
#  xtx = sqrt(diag(xtx))%*%R%*%sqrt(diag(xtx))
#  xty = diag(xtx)\beta
#
#############################################################
meSuSieData <- R6Class("meSuSieData",public = list(
  initialize = function(X,Y,var_y = 1){
    self$R<- X
    self$Summary_Stat <-Y
    self$var_y = var_y


    self$Name_list<-as.list(names(self$R))
    names(self$Name_list)<-names(self$R)
    self$N_ancestry<-length(X)

    self$XtX.diag<-self$XtX_diag(self$Summary_Stat,self$Name_list)
    self$XtX_list<-self$XtX_pro(self$R,self$XtX.diag,self$Name_list)
    self$Xty_list<-self$Xty_pro(self$Summary_Stat,self$XtX.diag,self$Name_list)

    self$N_list<-lapply(self$Summary_Stat,function(x)median(x$N))
    self$yty_list<-lapply(self$N_list,function(x)return(self$var_y*(x-1)))

    return(self)},
  ##First compute XtX.diag
  XtX_diag = function(Summary_Stat,Name_list){
    return( lapply(Name_list,function(x){
      R2 = (Summary_Stat[[x]]$Z^2)/(Summary_Stat[[x]]$Z^2+Summary_Stat[[x]]$N-2)
      sigma2 = self$var_y*(1-R2)*(Summary_Stat[[x]]$N-1)/(Summary_Stat[[x]]$N-2)
      return(sigma2/(Summary_Stat[[x]]$Se)^2)
    }))
  },


  ###Process XtX
  XtX_pro = function(R,XtX.diag,Name_list){
    return(lapply(Name_list,function(x){
      return(diag(sqrt(XtX.diag[[x]]))%*%R[[x]]%*%diag(sqrt(XtX.diag[[x]])))
    }))
  },

  ###Process XtY
  Xty_pro = function(Summary_Stat,XtX.diag,Name_list){
    return(lapply(Name_list,function(x){
      XtX.diag[[x]]*Summary_Stat[[x]]$Beta
    }))

  }
),
lock_objects = F)


# ======================================================================
#  single_effect_regression -> single_effect_regression_avail
#
#  Two changes, both required:
#   (a) the prior-variance optimiser calls test_run_loglik_cpp_avail and
#       passes config_snp_index, so infeasible cells cannot influence V;
#   (b) for each configuration s, the regression runs only on
#       J_s={j:M_{js}=1}, then pads back to m x |s| with
#       log BF^H_{s,j}=-Inf and
#       exactly zero posterior moments -- unavailable variant/ancestry
#       entries are never passed to uni_reg() or mvlmm_reg().
#  The softmax then uses compute_softmax_avail, which takes its log-sum-exp
#  stabiliser over feasible cells only and forces exact zero support.
# ======================================================================

single_effect_regression_avail<-function(XtR,XtX.diag, meSuSieObject_obj,l_index){
  column_config = meSuSieObject_obj$column_config
  N_ancestry = meSuSieObject_obj$nancestry
  avail = meSuSieObject_obj$avail
  Xty_standardized =Reduce(cbind, lapply(1:N_ancestry,function(x){
    XtR[,x]/meSuSieObject_obj$sigma2[x]
  }))

  shat2 = Reduce(cbind, lapply(1:N_ancestry,function(x){
    meSuSieObject_obj$sigma2[x]/XtX.diag[[x]]
  }))

  betahat = shat2 * Xty_standardized

  if(meSuSieObject_obj$estimate_prior_method =="optim"){
    V_mat = optimize_prior_variance_avail(betahat,shat2,meSuSieObject_obj$pi,N_ancestry,
                                          column_config,avail$config_snp_index)
  }

  return(ser_from_V_avail(betahat,shat2,V_mat,meSuSieObject_obj$pi,column_config,
                          avail$config_snp_index,avail$config_mask))
}


single_effect_regression_annot_avail<-function(XtR,XtX.diag, meSuSieObject_obj,l_index,annot_prior,estimate_prior_var_method,prior_var = NULL){
  column_config = meSuSieObject_obj$column_config
  N_ancestry = meSuSieObject_obj$nancestry
  avail = meSuSieObject_obj$avail
  Xty_standardized =Reduce(cbind, lapply(1:N_ancestry,function(x){
    XtR[,x]/meSuSieObject_obj$sigma2[x]
  }))

  shat2 = Reduce(cbind, lapply(1:N_ancestry,function(x){
    meSuSieObject_obj$sigma2[x]/XtX.diag[[x]]
  }))

  betahat = shat2 * Xty_standardized

  if(estimate_prior_var_method =="optim"){
    V_mat = optimize_prior_variance_avail(betahat,shat2,annot_prior,N_ancestry,
                                          column_config,avail$config_snp_index)
  }else if(estimate_prior_var_method =="fixed"){
    V_mat = prior_var
  }else{
    stop("estimate_prior_var_method must be 'fixed' or 'optim'.",call. = FALSE)
  }

  return(ser_from_V_avail(betahat,shat2,V_mat,annot_prior,column_config,
                          avail$config_snp_index,avail$config_mask))
}


# Effect-size covariance optimisation, restricted to J_s for each s.
# Identical structure to upstream (optim/L-BFGS-B for K == 2, the two-stage
# nloptr DIRECT_L + BOBYQA path otherwise); the only change is that the
# objective is test_run_loglik_cpp_avail and takes config_snp_index.
optimize_prior_variance_avail<-function(betahat,shat2,prior_weight,N_ancestry,column_config,config_snp_index){
  opt_par<-pre_optim(N_ancestry,-30,10)
  cfg_idx<-lapply(config_snp_index,as.numeric)

  if(N_ancestry==2){
    update_V<-optim(opt_par$inital_par,fn = test_run_loglik_cpp_avail,gr=NULL,betahat=betahat,shat2=shat2,prior_weight=prior_weight,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list =column_config,config_snp_index = cfg_idx,method = "L-BFGS-B",lower=opt_par$lower_bound,upper=opt_par$upper_bound)
    V_mat = vec_to_cov(update_V$par,opt_par$diag_index, opt_par$nancestry)
  }else{
    intermediate_V<-nloptr::nloptr(opt_par$inital_par, eval_f=test_run_loglik_cpp_avail,eval_grad_f = NULL,lb = opt_par$lower_bound,ub = opt_par$upper_bound,betahat=betahat,shat2=shat2,prior_weight=prior_weight,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list =column_config,config_snp_index = cfg_idx,opts=list("algorithm"= "NLOPT_GN_DIRECT_L","xtol_rel"=1.0e-10))
    update_V<-nloptr::nloptr(intermediate_V$solution, eval_f=test_run_loglik_cpp_avail,eval_grad_f = NULL,lb = opt_par$lower_bound,ub = opt_par$upper_bound,betahat=betahat,shat2=shat2,prior_weight=prior_weight,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list=column_config,config_snp_index = cfg_idx,opts=list("algorithm"= "NLOPT_LN_BOBYQA","xtol_rel"=1.0e-10))
    V_mat = vec_to_cov(update_V$solution,opt_par$diag_index, opt_par$nancestry)
  }
  return(V_mat)
}


# The shared tail of both SER functions.  Upstream inlines this; it is
# factored out here only so the two callers cannot drift apart.
ser_from_V_avail<-function(betahat,shat2,V_mat,prior_weight,column_config,config_snp_index,config_mask){
  p = nrow(betahat)

  reg_out<-lapply(seq_along(column_config),function(c){
    reg_one_config_avail(betahat,shat2,V_mat,
                         cols = as.integer(column_config[[c]]),
                         idx  = as.integer(config_snp_index[[c]]),
                         p    = p)
  })

  lbf<-Reduce(cbind,lapply(reg_out,function(x)x$lbf))
  lbf<-matrix(lbf,nrow = p)
  lbf[!config_mask]<- -Inf                    # kind of redundant, can replace with assertion

  softmax_out<-compute_softmax_avail(lbf,prior_weight,feasible = config_mask)

  mu1_multi = lapply(reg_out,function(x)x$post_mean)
  mu2_multi = lapply(reg_out,function(x)x$post_mean2)
  b1b2 = compute_b1b2(softmax_out$alpha_wmulti,mu1_multi,mu2_multi,column_config,ncol(betahat),p)
  return(list(alpha = softmax_out$alpha_wmulti,mu1_multi = mu1_multi,mu2_multi = mu2_multi,lbf_multi = lbf,V = V_mat,loglik =softmax_out$loglik,b1b2 = b1b2))
}


#' Compatibility entry point for availability-aware INFUSE
#'
#' Uses the same implementation as `infuse_core()`, retaining the argument order
#' of the development entry point. There is no second inference engine.
#' @inheritParams infuse_core
#' @return The fitted object returned by `infuse_core()`.
#' @export
infuse_core_avail <- function(R_mat_list, summary_stat_list, L,
                              variant_available = NULL, residual_variance = NULL,
                              prior_weights = NULL, ancestry_weight = NULL,
                              optim_method = "optim",
                              estimate_residual_variance = FALSE, max_iter = 100,
                              cor_method = "min.abs.corr", cor_threshold = 0.5,
                              annot = NULL, annot_method = NULL,
                              est_annot_prior = "fixed", verbose = TRUE) {
  infuse_core(R_mat_list = R_mat_list, summary_stat_list = summary_stat_list,
              L = L, residual_variance = residual_variance,
              prior_weights = prior_weights, ancestry_weight = ancestry_weight,
              optim_method = optim_method,
              estimate_residual_variance = estimate_residual_variance,
              max_iter = max_iter, cor_method = cor_method,
              cor_threshold = cor_threshold, annot = annot,
              annot_method = annot_method, est_annot_prior = est_annot_prior,
              variant_available = variant_available, verbose = verbose)
}
