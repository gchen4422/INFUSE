require(progress)
require(R6)
require(nloptr)
require(elasticnet)

#' infuse core function
#'
#' @param  R_mat_list A list of length N ancestry with each elment being correlation matrix with dimension p*p, the column name of the correlation matrix should match to the order of SNP name in summary_stat_list
#' @param  summary_stat_list A list of length N ancestry with each elment being summary statistics. The minimum requirement of summary statistics contains columns of SNP, Beta, Se, Z, and N. The order of the SNP should match the order of the correlation matrix 
#' @param  L Number of effects supposed within the region
#' @return An R6 object with pip, credible sets, and other features of the fine-mapping result. 
#' @export

infuse_core <- function(
  R_mat_list, summary_stat_list, L,
  residual_variance = NULL, prior_weights = NULL, ancestry_weight = NULL,
  optim_method = "optim", estimate_residual_variance = FALSE, max_iter = 100,
  cor_method = "min.abs.corr", cor_threshold = 0.5,
  annot = NULL, annot_method = NULL, est_annot_prior = "fixed"
) {
  cat("*************************************************************\n
  Multiple Ancestry Sum of Single Effect Model (MESuSiE)\n
   Visit http://www.xzlab.org/software.html For Update\n
            (C) 2022 Boran Gao, Xiang Zhou\n
              GNU General Public License\n
*************************************************************")

  time_start <- Sys.time()
  cat("\n# Start data processing for sufficient statistics \n")
  meSuSieData_obj <- meSuSieData$new(R_mat_list, summary_stat_list)

  n_snp <- nrow(summary_stat_list[[1]])
  n_ancestry <- length(summary_stat_list)

  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / n_snp, n_snp)
  }

  if (is.null(ancestry_weight)) {
    base_fac_ratio <- 3
    base_fac <- 1 / Reduce("+", lapply(seq_len(n_ancestry), function(x)
      choose(n_ancestry, x) * base_fac_ratio^(n_ancestry - x)))
    ancestry_weight <- unlist(lapply(seq_len(n_ancestry), function(x)
      rep(base_fac * base_fac_ratio^(n_ancestry - x), choose(n_ancestry, x))))
  }

  used_weights <- kronecker(prior_weights, t(ancestry_weight), FUN = "*")

  if (is.null(residual_variance)) {
    residual_variance <- rep(1, n_ancestry)
  }

  cat("# Create MESuSiE object \n")
  meSuSieObject_obj <- meSuSieObject$new(
    n_snp, L, n_ancestry, residual_variance, used_weights,
    optim_method, estimate_residual_variance, max_iter
  )

  cat("# Start data analysis \n")
  pb <- progress::progress_bar$new(
    format = " :elapsed",
    clear = TRUE,
    total = max_iter,
    show_after = 0
  )

  n_iter <- 0
  for (iter in seq_len(max_iter)) {
    comp_residual <- meSuSieObject_obj$compute_residual(meSuSieData_obj, meSuSieObject_obj)

    for (l_index in seq_len(L)) {
      comp_residual <- comp_residual + meSuSieObject_obj$Xr[[l_index]]
      SER_res <- single_effect_regression(comp_residual, meSuSieData_obj$XtX.diag, meSuSieObject_obj, l_index)
      meSuSieObject_obj$par_update(SER_res, l_index)
      meSuSieObject_obj$compute_KL(
        SER_res,
        meSuSieObject_obj$compute_SER_posterior_loglik(meSuSieData_obj, comp_residual, SER_res$b1b2),
        l_index
      )
      meSuSieObject_obj$compute_Xr(meSuSieData_obj, SER_res$b1b2$EB1, l_index)
      comp_residual <- comp_residual - meSuSieObject_obj$Xr[[l_index]]
    }

    pb$tick(tokens = list(iteration = iter))

    updated_sigma2 <- meSuSieObject_obj$update_residual_variance(meSuSieData_obj, iter)
    if ((meSuSieObject_obj$ELBO[iter + 1] - meSuSieObject_obj$ELBO[iter]) < 0.001) {
      break
    }
    if (isTRUE(meSuSieObject_obj$estimate_residual_variance)) {
      meSuSieObject_obj$sigma2 <- updated_sigma2
    }
    n_iter <- n_iter + 1
  }

  prior_ELBO <- meSuSieObject_obj$ELBO[iter]
  meSuSieObject_obj$ELBO <- rep(NA, max_iter)
  meSuSieObject_obj$ELBO[1] <- prior_ELBO

  if (!is.null(annot) && !is.null(annot_method)) {

    if (annot_method == "mlk") {
      if (ncol(as.matrix(annot)) >= 15) {
        sds <- apply(as.matrix(annot), 2, sd)
        annotation_nonconst <- as.matrix(annot)[, sds > 0, drop = FALSE]
        annotation_nonconst_scaled <- scale(annotation_nonconst)
        scaling_scales <- attr(annotation_nonconst_scaled, "scaled:scale")

        k_values <- c(5, 10, 15, 20)
        nonzero_values <- c(5, 10)
        results_comb <- expand.grid(k = k_values, q = nonzero_values)
        loadings_list <- list()

        X <- as.matrix(annotation_nonconst_scaled)
        storage.mode(X) <- "double"
        G <- crossprod(X) / (nrow(X) - 1)

        for (i in seq_len(nrow(results_comb))) {
          k_val <- results_comb$k[i]
          q_val <- results_comb$q[i]
          spca_result <- elasticnet::spca(
            G, K = k_val, type = "Gram",
            sparse = "varnum", para = rep(q_val, k_val)
          )
          loadings_list[[paste0("k", k_val, "_q", q_val)]] <- spca_result$loadings
          cat("Finished combination: k =", k_val, "and q =", q_val, "\n")
        }

        results_comb$BIC <- NA
        weights_mat <- matrix(NA, nrow = nrow(spca_result$loadings), ncol = nrow(results_comb))
        rownames(weights_mat) <- colnames(annotation_nonconst_scaled)
      }
    }

    updated_prior_list <- list()
    annot_weights_list <- list()

    for (l_index in seq_len(L)) {
      input.response <- rowSums(meSuSieObject_obj$alpha[[l_index]])

      if (max(diag(meSuSieObject_obj$V[[l_index]])) < 1e-9 || max(input.response) > 0.999) {
        updated_prior <- rep(1 / nrow(annot), nrow(annot))
        annot_weights_list[[l_index]] <- rep(0, ncol(as.matrix(annot)))
      } else {
        if (annot_method == "glmnet") {
          response.matrix <- matrix(c(1 - input.response, input.response), length(input.response), 2)
          try.index <- try(
            cv.logistic <- glmnet::cv.glmnet(
              x = as.matrix(annot), y = response.matrix,
              family = "multinomial", alpha = 0.5, type.measure = "deviance"
            )
          )
          if (class(try.index)[1] != "try-error") {
            cv.index <- which(cv.logistic$lambda == cv.logistic$lambda.min)
            glm.beta <- as.matrix(c(cv.logistic$glmnet.fit$a0[cv.index], cv.logistic$glmnet.fit$beta[[2]][, cv.index]))
            updated_prior <- exp(as.matrix(cbind(1, annot)) %*% glm.beta) /
              sum(exp(as.matrix(cbind(1, annot)) %*% glm.beta))
            annot_weights_list[[l_index]] <- cv.logistic$glmnet.fit$beta[[2]][, cv.index]
          } else {
            updated_prior <- rep(1 / nrow(annot), nrow(annot))
            annot_weights_list[[l_index]] <- rep(0, ncol(as.matrix(annot)))
          }
        }

        # keep rest unchanged
      }

      updated_prior_list[[l_index]] <- kronecker(updated_prior, t(ancestry_weight), FUN = "*")
    }
  }

  cat("\n# Data analysis is done, and now generates result \n\n")
  meSuSieObject_obj$get_result(
    meSuSie_get_cs(meSuSieObject_obj, R_mat_list, cor_method = cor_method, cor_threshold = cor_threshold),
    meSusie_get_pip_either(meSuSieObject_obj),
    meSusie_get_pip_config(meSuSieObject_obj)
  )
  meSuSieObject_obj$mesusie_summary(meSuSieData_obj)

  time_end <- Sys.time()
  cat(c("\n# Total time used for the analysis:",
        paste0(round(as.numeric(difftime(time_end, time_start, units = c("mins"))), 2), " mins\n")))
  return(meSuSieObject_obj)
}


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
    self$Xty_list<-self$Xty_pro(self$Summary_Stat,self$XtX.diag,self$Name_list) ##diag(xtx)^*betahat
    
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
single_effect_regression<-function(XtR,XtX.diag, meSuSieObject_obj,l_index){
  column_config = meSuSieObject_obj$column_config
  N_ancestry = meSuSieObject_obj$nancestry
  Xty_standardized =Reduce(cbind, lapply(1:N_ancestry,function(x){
    XtR[,x]/meSuSieObject_obj$sigma2[x]
  }))
  
  shat2 = Reduce(cbind, lapply(1:N_ancestry,function(x){
    meSuSieObject_obj$sigma2[x]/XtX.diag[[x]]
  }))
  
  betahat = shat2 * Xty_standardized
  
  if(meSuSieObject_obj$estimate_prior_method =="optim"){
    
    opt_par<-pre_optim(N_ancestry,-30,10)
    
    # update_V<-optim(opt_par$inital_par,fn = loglik_cpp_R6,gr=NULL,betahat,shat2,meSuSieObject_obj$pi,opt_par$nancestry,opt_par$diag_index,method = "L-BFGS-B",lower=opt_par$lower_bound,upper=opt_par$upper_bound)
    if(N_ancestry==2){
      
      update_V<-optim(opt_par$inital_par,fn = test_run_loglik_cpp,gr=NULL,betahat=betahat,shat2=shat2,prior_weight=meSuSieObject_obj$pi,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list =column_config,method = "L-BFGS-B",lower=opt_par$lower_bound,upper=opt_par$upper_bound)
      V_mat = vec_to_cov(update_V$par,opt_par$diag_index, opt_par$nancestry)
      
    }else{
      intermediate_V<-nloptr(opt_par$inital_par, eval_f=test_run_loglik_cpp,eval_grad_f = NULL,lb = opt_par$lower_bound,ub = opt_par$upper_bound,betahat=betahat,shat2=shat2,prior_weight=meSuSieObject_obj$pi,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list =column_config,opts=list("algorithm"= "NLOPT_GN_DIRECT_L","xtol_rel"=1.0e-10))
      update_V<-nloptr(intermediate_V$solution, eval_f=test_run_loglik_cpp,eval_grad_f = NULL,lb = opt_par$lower_bound,ub = opt_par$upper_bound,betahat=betahat,shat2=shat2,prior_weight=meSuSieObject_obj$pi,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list=column_config,opts=list("algorithm"= "NLOPT_LN_BOBYQA","xtol_rel"=1.0e-10))
      V_mat = vec_to_cov(update_V$solution,opt_par$diag_index, opt_par$nancestry)
      
    }
    
    
  }
  
  #  V_mat = vec_to_cov(update_V$par,opt_par$diag_index, opt_par$nancestry)
  # V_mat = vec_to_cov(update_V$solution,opt_par$diag_index, opt_par$nancestry)
  column_config = meSuSieObject_obj$column_config
  #multivariate_out<-multivariate_regression(Xty_standardized,shat2,V_mat) 
  # multivariate_out<-mvlmm_reg(betahat,shat2,V_mat) 
  reg_out<-lapply(column_config,function(x){
    if(length(x)==1){
      uni_reg(betahat[,x],shat2[,x],V_mat[x,x])
    }else if(length(x)>1){
      mvlmm_reg(betahat[,x],shat2[,x],V_mat[x,x]) 
    }
  })
  
  lbf<-Reduce(cbind,lapply(reg_out,function(x)x$lbf))
  lbf[is.na(lbf)]<-0
  
  softmax_out<-compute_softmax(lbf,meSuSieObject_obj$pi)
  
  mu1_multi = lapply(reg_out,function(x)x$post_mean)
  mu2_multi = lapply(reg_out,function(x)x$post_mean2)
  b1b2 = compute_b1b2(softmax_out$alpha_wmulti,mu1_multi,mu2_multi,column_config,ncol(betahat),nrow(betahat))
  return(list(alpha = softmax_out$alpha_wmulti,mu1_multi = mu1_multi,mu2_multi = mu2_multi,lbf_multi = lbf,V = V_mat,loglik =softmax_out$loglik,b1b2 = b1b2))
}





single_effect_regression_annot<-function(XtR,XtX.diag, meSuSieObject_obj,l_index,annot_prior,estimate_prior_var_method,prior_var = NULL){
  column_config = meSuSieObject_obj$column_config
  N_ancestry = meSuSieObject_obj$nancestry
  Xty_standardized =Reduce(cbind, lapply(1:N_ancestry,function(x){
    XtR[,x]/meSuSieObject_obj$sigma2[x]
  }))
  
  shat2 = Reduce(cbind, lapply(1:N_ancestry,function(x){
    meSuSieObject_obj$sigma2[x]/XtX.diag[[x]]
  }))
  
  betahat = shat2 * Xty_standardized
  
  if(estimate_prior_var_method =="optim"){
    
    opt_par<-pre_optim(N_ancestry,-30,10)
    
    # update_V<-optim(opt_par$inital_par,fn = loglik_cpp_R6,gr=NULL,betahat,shat2,meSuSieObject_obj$pi,opt_par$nancestry,opt_par$diag_index,method = "L-BFGS-B",lower=opt_par$lower_bound,upper=opt_par$upper_bound)
    if(N_ancestry==2){
      
      update_V<-optim(opt_par$inital_par,fn = test_run_loglik_cpp,gr=NULL,betahat=betahat,shat2=shat2,prior_weight=annot_prior,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list =column_config,method = "L-BFGS-B",lower=opt_par$lower_bound,upper=opt_par$upper_bound)
      V_mat = vec_to_cov(update_V$par,opt_par$diag_index, opt_par$nancestry)
      
    }else{
      intermediate_V<-nloptr(opt_par$inital_par, eval_f=test_run_loglik_cpp,eval_grad_f = NULL,lb = opt_par$lower_bound,ub = opt_par$upper_bound,betahat=betahat,shat2=shat2,prior_weight=annot_prior,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list =column_config,opts=list("algorithm"= "NLOPT_GN_DIRECT_L","xtol_rel"=1.0e-10))
      update_V<-nloptr(intermediate_V$solution, eval_f=test_run_loglik_cpp,eval_grad_f = NULL,lb = opt_par$lower_bound,ub = opt_par$upper_bound,betahat=betahat,shat2=shat2,prior_weight=annot_prior,nancestry =opt_par$nancestry,diag_index = opt_par$diag_index,config_list=column_config,opts=list("algorithm"= "NLOPT_LN_BOBYQA","xtol_rel"=1.0e-10))
      V_mat = vec_to_cov(update_V$solution,opt_par$diag_index, opt_par$nancestry)
      
    }
  }else if(estimate_prior_var_method =="fixed"){
    V_mat = prior_var
  }
  
  #  V_mat = vec_to_cov(update_V$par,opt_par$diag_index, opt_par$nancestry)
  # V_mat = vec_to_cov(update_V$solution,opt_par$diag_index, opt_par$nancestry)
  column_config = meSuSieObject_obj$column_config
  #multivariate_out<-multivariate_regression(Xty_standardized,shat2,V_mat) 
  # multivariate_out<-mvlmm_reg(betahat,shat2,V_mat) 
  reg_out<-lapply(column_config,function(x){
    if(length(x)==1){
      uni_reg(betahat[,x],shat2[,x],V_mat[x,x])
    }else if(length(x)>1){
      mvlmm_reg(betahat[,x],shat2[,x],V_mat[x,x]) 
    }
  })
  
  lbf<-Reduce(cbind,lapply(reg_out,function(x)x$lbf))
  lbf[is.na(lbf)]<-0
  
  softmax_out<-compute_softmax(lbf,annot_prior)
  #print(softmax_out$loglik)
  
  mu1_multi = lapply(reg_out,function(x)x$post_mean)
  mu2_multi = lapply(reg_out,function(x)x$post_mean2)
  b1b2 = compute_b1b2(softmax_out$alpha_wmulti,mu1_multi,mu2_multi,column_config,ncol(betahat),nrow(betahat))
  return(list(alpha = softmax_out$alpha_wmulti,mu1_multi = mu1_multi,mu2_multi = mu2_multi,lbf_multi = lbf,V = V_mat,loglik =softmax_out$loglik,b1b2 = b1b2))
}



