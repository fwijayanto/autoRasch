#' autoRasch Options
#'
#' Show the default settings used by the functions in \strong{autoRasch} package. These settings can be changed by passing values through the "..." argument.
#'
#' @param x \emph{Character}. A name of single parameter setting that is wanted to be shown.
#'
#' @return
#' \item{fixed_par}{   A vector of parameter types that are set to be fix. It means that these parameters are not estimated.}
#' \item{fixed_theta}{   A vector of \code{theta} values when \code{theta} are set to be fix at the \code{fixed_par}. If it is not set, it will be set to zero.}
#' \item{fixed_beta}{   A vector of \code{theta} values when \code{beta} are set to be fix at the \code{fixed_par}. If it is not set, it will be set to zero.}
#' \item{fixed_gamma}{   A vector of \code{gamma} (natural logarithm of discrimination parameters, \eqn{\alpha = exp(\gamma)}) values when \code{gamma} are set to be fix at the \code{fixed_par}. If it is not set, it will be set to zero.}
#' \item{fixed_deltabeta}{   A vector of \code{deltabeta} values when \code{deltabeta} are set to be fix at the \code{fixed_par}. If it is not set, it will be set to zero.}
#' \item{isPenalized_theta}{   It is a logical parameter whether, in the estimation procedure, \code{theta} is penalized or not.}
#' \item{isPenalized_gamma}{   It is a logical parameter whether, in the estimation procedure, \code{gamma} is penalized or not.}
#' \item{isPenalized_deltabeta}{   It is a logical parameter whether, in the estimation procedure, \code{deltabeta} is penalized or not.}
#' \item{groups_map}{   A matrix \eqn{n x f} to map the subject into DIF groups, where \eqn{n} is number of subjects and \eqn{f} is number of focal groups.}
#' \item{optz_tuner}{   A list of optimization function settings. For complete settings can be seen in \code{\link[stats:optim]{stats::optim()}}.}
#' \item{lambda_theta}{   An integer value to set the regularization parameter to the \code{theta}.}
#' \item{lambda_in}{   An integer value to set the regularization parameter to the \code{gamma} in the included itemset.}
#' \item{lambda_out}{   An integer value to set the regularization parameter to the \code{gamma} in the excluded itemset.}
#' \item{lambda_deltabeta}{   An integer value to set the regularization parameter to the \code{deltabeta}.}
#' \item{isHessian}{   It is a logical parameter whether, in the estimation procedure, need to return the Hessian matrix or not.}
#' \item{isTracked}{   It is a logical parameter whether need to track the process or not.}
#'
#' @examples
#' autoRaschOptions()
#' autoRaschOptions(x = "isTracked")
#'
#' @export
autoRaschOptions <- function(x = NULL){

  aRoptions <- autoRasch_options_default()

  if(!is.null(x)) {
    if(is.character(x)) {
      # lower case only
      #x <- tolower(x)

      # check if x is in names(aRoptions)
      not.ok <- which(!x %in% names(aRoptions))
      if(length(not.ok) > 0L) {
        # only warn if multiple options were requested
        if(length(x) > 1L) {
          warning("autoRasch WARNING: option `", x[not.ok],
                  "' not available")
        }
        x <- x[ -not.ok ]
      }

      # return requested option(s)
      if(length(x) == 0L) {
        return(default)
      } else {
        aRoptions[x]
      }
    } else {
      stop("autoRasch ERROR: `x' must be a character string")
    }
  } else {
    aRoptions
  }
}

autoRasch_options_default <- function(){

  opt <- list(
    fixed_par = c(),
    fixed_theta = c(),
    fixed_beta = c(),
    fixed_gamma = c(),
    fixed_deltabeta = c(),
    isPenalized_gamma = TRUE,
    isPenalized_deltabeta = TRUE,
    isPenalized_theta = TRUE,
    groups_map = c(),
    # resp.th = c()
    # optz_method = c("optim","nlminb"),
    optz_tuner = list(maxit = 2e+4, reltol = 1e-12, fnscale = 10),
    # objtype = "",
    # desc = NULL,
    lambda_theta = 0.05,
    lambda_in = 50,
    lambda_out = 5e-3,
    lambda_deltabeta = 15,
    # lambda_deltagamma = 10000,
    # eps = 0.0,
    # random.init = FALSE,
    # random.init.th = 1e-2,
    isHessian = TRUE,
    isTracked = TRUE
  )

  return(opt)

}


testFUnction <- function(X,...){

  dotdotdot <- list(...)
  temp <- list(x,dotdotdot)
  return(temp)
}
