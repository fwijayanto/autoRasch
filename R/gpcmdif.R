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
#' @seealso \code{\link{pcm}}, \code{\link{pcm_dif}}, \code{\link{gpcm}}, \code{\link{gpcm_dif}}
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


#' @param obj The object of class \code{'gpcm_dif'}.
#' @param par The parameter that are wanted to be summarized.
#' @param th_dif The threshold of a recognizable estimated DIF parameters.
#'
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
