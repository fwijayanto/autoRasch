#' The Inhomogenous Polytomous Dataset
#'
#' The artifical dataset of a polytomous responses (five categories) which contains of three subsets with different discrimination values.
#' To reproduce this dataset: \cr
#' \code{poly_inh_dset <- generate_data(responseType = "discriminate", ncat = 5, alpha = c(0.04,0.045,0.05,0.055,0.06,0.065,0.2,0.25,0.3,0.35,0.4,0.45,2.6,2.65,2.7,2.75,2.8,2.85,2.9))}
#' \cr\cr will lead to similar but not the same dataset, due to the randomization.
#'
#' @docType data
#'
#' @usage data(poly_inh_dset)
#'
#' @rdname poly_inh_dset
"poly_inh_dset"
