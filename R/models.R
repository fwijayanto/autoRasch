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
#' \item{mt_vek}{   A vector of the highest response category as many as the number of items.}
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



