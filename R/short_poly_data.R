#' A Shorter Inhomogenous Polytomous Dataset
#'
#' The artifical dataset of a polytomous responses (three categories) which contains of three subsets with different discrimination values.
#' To reproduce this dataset: \cr
#' \code{short_poly_data <- generate_data(alpha = c(0.02,0.5,2), nitem = 3, ndim = 3,ncat = 5, theta = c(-6,6), beta = c(-4,4), ntheta = 151)}
#' \cr\cr will lead to similar but not the same dataset, due to the randomization.
#'
#' @docType data
#'
#' @usage data(short_poly_data)
#'
#' @rdname short_poly_data
"short_poly_data"
