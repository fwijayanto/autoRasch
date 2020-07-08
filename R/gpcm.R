#' Estimation of The Generalized Partial Credit Model
#'
#' \code{gpcm()} computes the parameter estimates of a generalized partial credit model for polytomous responses
#' by using penalized JML estimation.
#'
#' @inheritParams pcm
#'
#' @return
#' \item{X}{   The dataset that is used for estimation.}
#' \item{mt_vek}{   A vector of the highest response category as many as the number of items.}
#' \item{real_vek}{   A vector of the presence of thresholds in the items-thresholds order: \code{1}
#' if the threshold is present and \code{NA} if unavailable.}
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
#' @seealso \code{\link{pcm}}, \code{\link{gpcm}}
#'
#' @examples
#' res <- gpcm(pcm_data)
#' res
#' summary(res)
#'
#' @export
gpcm <- function(X, isHessian = TRUE){

  result <- pjmle(X = X, fixed_par = c("deltabeta"), isPenalized_deltabeta = FALSE, isHessian = isHessian)
  class(result) <- c("gpcm","armodels","autoRasch",class(result))
  return(result)
}


#' @param obj The object of class \code{'gpcm'}.
#' @param par The parameter that are wanted to be summarized.
#'
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
