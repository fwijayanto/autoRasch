#' Compute Reliability and Standard Error
#'
#' This function computes the reliability index, separation and the standard error of the models estimation.
#'
#' @param obj Object that resulted from any models estimation, e.g., \code{PCM}, \code{GPCM}, \code{PCMDIF}, and \code{GPCMDIF}.
#'
#' @return
#' A list of two objects, the reliability and the standard error.
#'
#' \emph{reliability}
#' \itemize{
#'    \item{PRI}{   Person reliability index.}
#'    \item{PSR}{   Person separation reliability.}
#'    \item{IRI}{   Item reliability index.}
#'    \item{ISR}{   Item separation reliability.}
#' }
#' \emph{stdError}
#' \itemize{
#'    \item{var_err_pers}{   A matrix of variance error of the estimation.}
#'    \item{std_err_pers}{   A matrix of standard error of the estimation.}
#'    \item{rmsse_pers}{   Root mean square of the standard error per person.}
#'    \item{var_err_item}{   A matrix of variance error of the estimation.}
#'    \item{std_err_item}{   A matrix of standard error of the estimation.}
#'    \item{rmsse_item}{   Root mean square of the standard error per person.}
#'    \item{hessian_theta}{   Hessian matrix of \code{theta} parameter.}
#'    \item{hessian_beta}{   Hessian matrix of \code{beta} parameter.}
#' }
#'
#' @details
#' Person reliability index
#'
#' @examples
#' pcmObject <- pcm(poly_inh_dset)
#' rel <- checkRel(pcmObject)
#' summary(rel)
#'
#'
#' @export
checkRel <- function(obj){

  if(length(which(is.na(obj$real_vek))) != 0){
    rem.idx <- length(obj$theta)+which(is.na(obj$real_vek))

    obj$hessian <- obj$hessian[-c(rem.idx),-c(rem.idx)]
  }

  if(is.null(obj$hessian)){
    stop("autoRasch ERROR: the separation reliability and standard error can not be computed without Hessian matrix.")
  }

  rmseroor <- stdError(obj)

  rmse <- rmseroor$rmsse_pers
  p_var <- var(obj$theta)
  true_pvar <- p_var - (rmse^2)
  true_psd <- sqrt(true_pvar)
  p_sep_coeff <- true_psd/rmse
  p_rel_idx <- (p_sep_coeff^2)/(1+(p_sep_coeff^2))

  rmse_item <- rmseroor$rmsse_item
  i_var <- var(obj$beta)
  true_ivar <- i_var - (rmse_item^2)
  true_isd <- sqrt(true_ivar)
  i_sep_coeff <- true_isd/rmse_item
  i_rel_idx <- (i_sep_coeff^2)/(1+(i_sep_coeff^2))

  result <- list("reliability" = list("PRI" = p_sep_coeff, "PSR" = p_rel_idx, "IRI" = i_sep_coeff, "ISR" = i_rel_idx), "stdError" = rmseroor)
  class(result) <- c("seprel","autoRasch",class(result))
  return(result)
}

stdError <- function(obj){

  hess_theta <- obj$hessian[1:length(obj$theta),1:length(obj$theta)]
  varerr_p <- (diag(solve(hess_theta)))
  stderr_p <- sqrt(varerr_p)
  rmse_p <- sqrt(mean(varerr_p))
  hess_beta <- obj$hessian[(length(obj$theta)+1):(length(obj$theta)+length(which(!is.na(obj$real_vek)))),(length(obj$theta)+1):(length(obj$theta)+length(which(!is.na(obj$real_vek))))]
  varerr_i <- (diag(solve(hess_beta)))
  stderr_i <- sqrt(varerr_i)
  rmse_i <- sqrt(mean(varerr_i))

  return(list("var_err_pers" = varerr_p, "std_err_pers" = stderr_p, "rmsse_pers" = rmse_p, "var_err_item" = varerr_i,
              "std_err_item" = stderr_i, "rmsse_item" = rmse_i, "hessian_theta" = hess_theta, "hessian_beta" = hess_beta))
}

summary.seprel <- function(obj,...){

  res_table <- cbind(c(obj$sepRel$PRI,obj$sepRel$PSR,obj$stdErr$rmse_pers),c(obj$sepRel$IRI,obj$sepRel$ISR,obj$stdErr$rmse_item))
  dimnames(res_table) <- list(c("Reliability Index","Separation Reliability","RMSSE"),c("Person","Item"))
  cat("\n")
  print(res_table)
  cat("\n")

}

