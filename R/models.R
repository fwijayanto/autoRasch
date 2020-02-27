#' Estimation of The Partial Credit Model (PCM)
#'
#' This function computes the parameter estimates of a partial credit model for dichotomous and polytomous responses
#' by using penalized JML estimation.
#'
#' @param X Input dataset as matrix or data frame with ordinal responses (starting from 0);
#' rows represent individuals, columns represent items.
#' @param isHessian a logical parameter setting whether or not the Hessian matrix is needed to be computed.
#' The Hessian matrix is needed to compute the separation reliability.
#'
#' @return
#' \item{X}{   The dataset that is used for estimation.}
#' \item{mt_vek}{   A vector of the highest response category as many as the number of items.}
#' \item{real_vek}{   A vector of the presence of thresholds in the items-thresholds order: \code{1} if the threshold is present and \code{NA} if unavailable.}
#' \item{itemName}{   The vector of names of items (columns) in the dataset.}
#' \item{loglik}{   The log likelihood of the estimation.}
#' \item{hessian}{   The hessian matrix. Only when the \code{isHessian = TRUE}.}
#' \item{beta}{   A vector of the difficulty parameter of each categories of items (thresholds).}
#' \item{theta}{   A vector of the ability parameters of each individuals.}
#'
#' @examples
#' res <- pcm(poly_inh_dset)
#' res
#' summary(res)
#'
#' @export
pcm <- function(X, isHessian = TRUE){

  result <- pjmle(X = X, fixed_par = c("gamma","deltabeta"), isPenalized_gamma = FALSE, isPenalized_deltabeta = FALSE, isHessian = isHessian)
  class(result) <- c("armodels","pcm","autoRasch")
  return(result)

}

#' Estimation of The Generalized Partial Credit Model
#'
#' This function computes the parameter estimates of a generalized partial credit model for polytomous responses
#' by using penalized JML estimation.
#'
#' @inheritParams pcm
#'
#' @return
#' \item{X}{   The dataset that is used for estimation.}
#' \item{mt_vek}{   A vector of the highest response category as many as the number of items.}
#' \item{real_vek}{   A vector of the presence of thresholds in the items-thresholds order: \code{1} if the threshold is present and \code{NA} if unavailable.}
#' \item{itemName}{   The vector of names of items (columns) in the dataset.}
#' \item{loglik}{   The log likelihood of the estimation.}
#' \item{hessian}{   The hessian matrix. Only when the \code{isHessian = TRUE}.}
#' \item{gamma}{   A vector of the natural logarithm of discrimination parameters of each items.}
#' \item{beta}{   A vector of the difficulty parameter of each items' categories (thresholds).}
#' \item{theta}{   A vector of the ability parameters of each individuals.}
#'
#' @details
#' In the discrimination parameters estimation, instead of estimating the discrimination parameters,
#' we are estimating the natural logarithm of the parameters to avoid negative values, \eqn{\alpha = exp(\gamma)}.
#'
#' @examples
#' res <- gpcm(poly_inh_dset)
#' res
#' summary(res)
#'
#' @export
gpcm <- function(X, isHessian = TRUE){

  result <- pjmle(X = X, fixed_par = c("deltabeta"), isPenalized_deltabeta = FALSE, isHessian = isHessian)
  class(result) <- c("armodels","gpcm","autoRasch")
  return(result)
}



#' Estimation of The Partial Credit Model with DIF
#'
#' This function computes the parameter estimates of a partial credit model with DIF for dichotomous and polytomous responses
#' by using penalized JML estimation.
#'
#' @inheritParams pcm
#' @param dif_map Binary matrix. Respondents membership to DIF groups; rows represent individuals, column represent group partitions.
#'
#' @return
#' \item{X}{   The dataset that is used for estimation.}
#' \item{mt_vek}{   A vector of the highest response category as many as the number of items.}
#' \item{real_vek}{   A vector of the presence of thresholds in the items-thresholds order:
#' \code{1} if the threshold is present and \code{NA} if unavailable.}
#' \item{itemName}{   The vector of names of items (columns) in the dataset.}
#' \item{loglik}{   The log likelihood of the estimation.}
#' \item{hessian}{   The hessian matrix. Only when the \code{isHessian = TRUE}.}
#' \item{beta}{   A vector of the difficulty parameter of each categories of items (thresholds).}
#' \item{theta}{   A vector of the ability parameters of each individuals.}
#'
#' @examples
#' res <- pcm_dif(polydif_inh_dset, groups_map = c(rep(1,245),rep(0,245)))
#' res
#' summary(res)
#'
#' @export
pcm_dif <- function(X, groups_map, isHessian = TRUE){

  result <- pjmle(X = X, fixed_par = c("gamma"), isPenalized_gamma = FALSE, groups_map = groups_map, isHessian = isHessian)
  class(result) <- c("armodels","pcmdif","autoRasch")
  return(result)
}


#' Estimation of The Generalized Partial Credit Model with DIF
#'
#' This function computes the parameter estimates of a generalized partial credit model with DIF for polytomous responses
#' by using penalized JML estimation.
#'
#' @inheritParams pcm_dif
#'
#' @return
#' \item{X}{   The dataset that is used for estimation.}
#' \item{mt_vek}{   A vector of the highest response category as many as the number of items.}
#' \item{real_vek}{   A vector of the presence of thresholds in the items-thresholds order: \code{1} if the threshold is present and \code{NA} if unavailable.}
#' \item{itemName}{   The vector of names of items (columns) in the dataset.}
#' \item{loglik}{   The log likelihood of the estimation.}
#' \item{hessian}{   The hessian matrix. Only when the \code{isHessian = TRUE}.}
#' \item{deltabeta}{   A vector of the DIF parameters of each items on each groups.}
#' \item{gamma}{   A vector of the natural logarithm of discrimination parameters of each items.}
#' \item{beta}{   A vector of the difficulty parameter of each items' categories (thresholds).}
#' \item{theta}{   A vector of the ability parameters of each individuals.}
#'
#' @details
#' In the discrimination parameters estimation, instead of estimating the discrimination parameters,
#' we are estimating the natural logarithm of the parameters to avoid negative values, \eqn{\alpha = exp(\gamma)}.
#'
#' @examples
#' res <- gpcm_dif(polydif_inh_dset, groups_map = c(rep(1,245),rep(0,245)))
#' res
#' summary(res)
#'
#'
#' @export
gpcm_dif <- function(X, groups_map, isHessian = TRUE){

  result <- pjmle(X = X, groups_map = groups_map, isHessian = isHessian)
  class(result) <- c("armodels","gpcmdif","autoRasch")
  return(result)

}

#' Estimation of the generic form of the models
#'
#' This function computes the parameter estimates of the generic form of the models by using penalized JML estimation. It allows users to adjust the default settings of the estimation.
#'
#' @param X Input dataset as matrix or data frame with ordinal responses (starting from 0); rows represent individuals, column represent items.
#' @param ... Parameter settings which are listed in \code{\link[autoRaschOptions]{autoRaschOptions()}}.
#'
#' @return
#' \item{X}{   The dataset that is used for estimation.}
#' \item{name}{   The name of each items in the dataset.}
#' \item{mt_vek}{   }
#' \item{loglik}{   The log likelihood of the estimation.}
#' \item{objtype}{   Type of the model that is used.}
#' \item{deltabeta}{   A vector of the DIF parameters of each items on each groups.}
#' \item{gamma}{   A vector of the natural logarithm of discrimination parameters of each items.}
#' \item{beta}{   A vector of the difficulty parameter of each items' categories (thresholds).}
#' \item{theta}{   A vector of the ability parameters of each individuals.}
#'
#' @details
#' In the discrimination parameters estimation, instead of estimating the discrimination parameters,
#' we are estimating the natural logarithm of the parameters to avoid negative values, \eqn{\alpha = exp(\gamma)}.
#'
#'
#' @export
generic_model <- function(X, ...){

  dotdotdot <- list(...)
  opts_default <- autoRaschOptions()

  opt_names <- names(dotdotdot)[which(names(dotdotdot) %in% names(opts_default))]
  opts <- opts_default
  for (i in opt_names) {
    opts[[i]] <- dotdotdot[[i]]
  }

  return(pjmle(X = X, fixed_par = opts$fixed_par, fixed_theta = opts$fixed_theta, fixed_beta = opts$fixed_beta, fixed_gamma = opts$fixed_gamma,
               fixed_deltabeta = opts$fixed_deltabeta, isPenalized_theta = opts$isPenalized_theta, isPenalized_gamma = opts$isPenalized_gamma,
               isPenalized_deltabeta = opts$isPenalized_deltabeta, groups_map = opts$groups_map, optz_tuner = opts$optz_tuner,
               lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in, lambda_out = opts$lambda_out, lambda_deltabeta = opts$lambda_deltabeta,
               isHessian = opts$isHessian, isTracked = opts$isTracked))

}


#' @param obj The object of class \code{'pcm'}.
#' @param par The parameter that are wanted to be summarized.
#'
#' @rdname pcm
#' @export
summary.pcm <- function(obj, par = c()){

  if(is.null(par) | "theta" %in% par){
    cat("\n\n")
    cat("The estimated ability scores:")
    cat("\n")
    print(obj$theta)
    cat("\n")
    cat("The highest ability score: ",round(max(obj$theta,na.rm = TRUE),4))
    cat("\n")
    cat("The lowest ability score: ",round(min(obj$theta,na.rm = TRUE),4))
  }

  if(is.null(par) | "beta" %in% par){
    cat("\n\n")
    cat("The estimated difficulty scores:")
    cat("\n")
    reported_beta <- obj$beta * obj$real_vek
    beta_mat <- matrix(reported_beta, nrow = length(obj$mt_vek), byrow = TRUE)
    beta_mat <- as.data.frame(round(beta_mat,4), row.names = obj$itemName)
    colnames(beta_mat) <- paste("Th_",c(1:max(obj$mt_vek)),sep = "")
    beta_mat[["Item Loc."]] <- round(apply(beta_mat,1,mean,na.rm=TRUE),4)
    beta_mat$` ` <- apply(beta_mat[,1:max(obj$mt_vek)],1,function(x){if(is.unsorted(na.omit(x))){return("*")}else{return("")}})
    print(beta_mat, quote = FALSE)
    cat("\n")
    cat("The most difficult item: ",obj$itemName[which(beta_mat[,5] == max(beta_mat[,5],na.rm = TRUE))])
    cat("\n")
    cat("The easiest item: ",obj$itemName[which(beta_mat[,5] == min(beta_mat[,5],na.rm = TRUE))])
    cat("\n")
    ntd_items <- length(which(beta_mat[,ncol(beta_mat)] == "*"))
    cat("There are",ntd_items,"items which have disordered thresholds.")
    cat("\n")
    cat("'*' Item has disordered thresholds.")
  }
}

#' @inheritParams summary.pcm
#' @rdname gpcm
#' @export
summary.gpcm <- function(obj, par = c()){

  if(is.null(par) | "theta" %in% par){
    cat("\n\n")
    cat("The estimated ability scores:")
    cat("\n")
    print(obj$theta)
    cat("\n")
    cat("The highest ability score: ",round(max(obj$theta,na.rm = TRUE),4))
    cat("\n")
    cat("The lowest ability score: ",round(min(obj$theta,na.rm = TRUE),4))
  }

  if(is.null(par) | "beta" %in% par){
    cat("\n\n")
    cat("The estimated difficulty scores:")
    cat("\n")
    reported_beta <- obj$beta * obj$real_vek
    beta_mat <- matrix(reported_beta, nrow = length(obj$mt_vek), byrow = TRUE)
    beta_mat <- as.data.frame(round(beta_mat,4), row.names = obj$itemName)
    colnames(beta_mat) <- paste("Th_",c(1:max(obj$mt_vek)),sep = "")
    beta_mat[["Item Loc."]] <- round(apply(beta_mat,1,mean,na.rm=TRUE),4)
    beta_mat$` ` <- apply(beta_mat[,1:max(obj$mt_vek)],1,function(x){if(is.unsorted(na.omit(x))){return("*")}else{return("")}})
    print(beta_mat, quote = FALSE)
    cat("\n")
    cat("The most difficult item: ",obj$itemName[which(beta_mat[,5] == max(beta_mat[,5],na.rm = TRUE))])
    cat("\n")
    cat("The easiest item: ",obj$itemName[which(beta_mat[,5] == min(beta_mat[,5],na.rm = TRUE))])
    cat("\n")
    ntd_items <- length(which(beta_mat[,ncol(beta_mat)] == "*"))
    cat("There are",ntd_items,"items which have disordered thresholds.")
    cat("\n")
    cat("'*' Item has disordered thresholds.")
  }

  if(is.null(par) | "gamma" %in% par){
    cat("\n\n")
    cat("The estimated discrimination parameters:")
    cat("\n")
    alpha_mat <- matrix(exp(obj$gamma), ncol = 1, dimnames = list(c(obj$itemName),c("alpha")))
    print(alpha_mat, quote = FALSE)
  }
}

#' @inheritParams summary.pcm
#' @param th_dif The threshold of a recognizable estimated DIF parameters.
#' @rdname pcm_dif
#' @export
summary.pcm_dif <- function(obj, par = c(), th_dif = NULL){

  if(is.null(par) | "theta" %in% par){
    cat("\n\n")
    cat("The estimated ability scores:")
    cat("\n")
    print(obj$theta)
    cat("\n")
    cat("The highest ability score: ",round(max(obj$theta,na.rm = TRUE),4))
    cat("\n")
    cat("The lowest ability score: ",round(min(obj$theta,na.rm = TRUE),4))
  }

  if(is.null(par) | "beta" %in% par){
    cat("\n\n")
    cat("The estimated difficulty scores:")
    cat("\n")
    reported_beta <- obj$beta * obj$real_vek
    beta_mat <- matrix(reported_beta, nrow = length(obj$mt_vek), byrow = TRUE)
    beta_mat <- as.data.frame(round(beta_mat,4), row.names = obj$itemName)
    colnames(beta_mat) <- paste("Th_",c(1:max(obj$mt_vek)),sep = "")
    beta_mat[["Item Loc."]] <- round(apply(beta_mat,1,mean,na.rm=TRUE),4)
    beta_mat$` ` <- apply(beta_mat[,1:max(obj$mt_vek)],1,function(x){if(is.unsorted(na.omit(x))){return("*")}else{return("")}})
    print(beta_mat, quote = FALSE)
    cat("\n")
    cat("The most difficult item: ",obj$itemName[which(beta_mat[,5] == max(beta_mat[,5],na.rm = TRUE))])
    cat("\n")
    cat("The easiest item: ",obj$itemName[which(beta_mat[,5] == min(beta_mat[,5],na.rm = TRUE))])
    cat("\n")
    ntd_items <- length(which(beta_mat[,ncol(beta_mat)] == "*"))
    cat("There are",ntd_items,"items which have disordered thresholds.")
    cat("\n")
    cat("'*' Item has disordered thresholds.")
  }

  if(is.null(par) | "deltabeta" %in% par){
    cat("\n\n")
    cat("The estimated DIF effects (gap of difficulties):")
    cat("\n")
    if(is.null(th_dif)){
      th_dif <- 1e-5
    }
    delta_mat <- matrix(round(obj$deltabeta,5), ncol = ncol(obj$groups_map), dimnames = list(c(obj$itemName),c(paste("Group",c(1:ncol(obj$groups_map)),sep = ""))))
    delta_mat[which(delta_mat < th_dif)] <- ""
    delta_mat <- as.data.frame(delta_mat)
    print(delta_mat, quote = FALSE)
    cat("\n")
    cat("DIF effect threshold =",th_dif)
  }
}

#' @inheritParams summary.gpcm
#' @inheritParams summary.pcm_dif
#' @rdname gpcm_dif
#' @export
summary.gpcm_dif <- function(obj, par = c(), th_dif = NULL){

  if(is.null(par) | "theta" %in% par){
    cat("\n\n")
    cat("The estimated ability scores:")
    cat("\n")
    print(obj$theta)
    cat("\n")
    cat("The highest ability score: ",round(max(obj$theta,na.rm = TRUE),4))
    cat("\n")
    cat("The lowest ability score: ",round(min(obj$theta,na.rm = TRUE),4))
  }

  if(is.null(par) | "beta" %in% par){
    cat("\n\n")
    cat("The estimated difficulty scores:")
    cat("\n")
    reported_beta <- obj$beta * obj$real_vek
    beta_mat <- matrix(reported_beta, nrow = length(obj$mt_vek), byrow = TRUE)
    beta_mat <- as.data.frame(round(beta_mat,4), row.names = obj$itemName)
    colnames(beta_mat) <- paste("Th_",c(1:max(obj$mt_vek)),sep = "")
    beta_mat[["Item Loc."]] <- round(apply(beta_mat,1,mean,na.rm=TRUE),4)
    beta_mat$` ` <- apply(beta_mat[,1:max(obj$mt_vek)],1,function(x){if(is.unsorted(na.omit(x))){return("*")}else{return("")}})
    print(beta_mat, quote = FALSE)
    cat("\n")
    cat("The most difficult item: ",obj$itemName[which(beta_mat[,5] == max(beta_mat[,5],na.rm = TRUE))])
    cat("\n")
    cat("The easiest item: ",obj$itemName[which(beta_mat[,5] == min(beta_mat[,5],na.rm = TRUE))])
    cat("\n")
    ntd_items <- length(which(beta_mat[,ncol(beta_mat)] == "*"))
    cat("There are",ntd_items,"items which have disordered thresholds.")
    cat("\n")
    cat("'*' Item has disordered thresholds.")
  }

  if(is.null(par) | "gamma" %in% par){
    cat("\n\n")
    cat("The estimated discrimination parameters:")
    cat("\n")
    alpha_mat <- matrix(exp(obj$gamma), ncol = 1, dimnames = list(c(obj$itemName),c("alpha")))
    print(alpha_mat, quote = FALSE)
  }

  if(is.null(par) | "deltabeta" %in% par){
    cat("\n\n")
    cat("The estimated DIF effects (gap of difficulties):")
    cat("\n")
    if(is.null(th_dif)){
      th_dif <- 1e-5
    }
    delta_mat <- matrix(round(obj$deltabeta,5), ncol = ncol(obj$groups_map), dimnames = list(c(obj$itemName),c(paste("Group",c(1:ncol(obj$groups_map)),sep = ""))))
    delta_mat[which(delta_mat < th_dif)] <- ""
    delta_mat <- as.data.frame(delta_mat)
    print(delta_mat, quote = FALSE)
    cat("\n")
    cat("DIF effect threshold =",th_dif)
  }
}

#' @inheritParams summary.pcm
#' @rdname pcm
#' @export
print.pcm <- function(obj, par = c()){
  cls <- class(obj)
  class(obj) <- "list"

  if(is.null(par) | "theta" %in% par){
    cat("\n")
    cat("$theta")
    cat("\n")
    print(obj$theta)
    cat("\n")
  }

  if(is.null(par) | "beta" %in% par){
    cat("\n")
    cat("$beta")
    cat("\n")
    print(obj$beta)
    cat("\n")
  }

}

#' @inheritParams print.pcm
#' @rdname gpcm
#' @export
print.gpcm <- function(obj, par = c()){
  cls <- class(obj)
  class(obj) <- "list"

  if(is.null(par) | "theta" %in% par){
    cat("\n")
    cat("$theta")
    cat("\n")
    print(obj$theta)
    cat("\n")
  }

  if(is.null(par) | "beta" %in% par){
    cat("\n")
    cat("$beta")
    cat("\n")
    print(obj$beta)
    cat("\n")
  }

  if(is.null(par) | "gamma" %in% par){
    cat("\n")
    cat("$gamma")
    cat("\n")
    print(obj$gamma)
    cat("\n")
  }

}

#' @inheritParams print.pcm
#' @rdname pcm_dif
#' @export
print.pcm_dif <- function(obj, par = c()){
  cls <- class(obj)
  class(obj) <- "list"

  if(is.null(par) | "theta" %in% par){
    cat("\n")
    cat("$theta")
    cat("\n")
    print(obj$theta)
    cat("\n")
  }

  if(is.null(par) | "beta" %in% par){
    cat("\n")
    cat("$beta")
    cat("\n")
    print(obj$beta)
    cat("\n")
  }

  if(is.null(par) | "gamma" %in% par){
    cat("\n")
    cat("$gamma")
    cat("\n")
    print(obj$gamma)
    cat("\n")
  }

  if(is.null(par) | "deltabeta" %in% par){
    cat("\n")
    cat("$deltabeta")
    cat("\n")
    print(obj$deltabeta)
    cat("\n")
  }

}

#' @inheritParams print.pcm
#' @rdname gpcm_dif
#' @export
print.gpcm_dif <- function(obj, par = c()){
  cls <- class(obj)
  class(obj) <- "list"

  if(is.null(par) | "theta" %in% par){
    cat("\n")
    cat("$theta")
    cat("\n")
    print(obj$theta)
    cat("\n")
  }

  if(is.null(par) | "beta" %in% par){
    cat("\n")
    cat("$beta")
    cat("\n")
    print(obj$beta)
    cat("\n")
  }

  if(is.null(par) | "deltabeta" %in% par){
    cat("\n")
    cat("$deltabeta")
    cat("\n")
    print(obj$deltabeta)
    cat("\n")
  }

}

