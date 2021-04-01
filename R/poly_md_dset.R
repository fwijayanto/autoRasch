#' The Uncorrelated Multidimensional Polytomous Dataset
#'
#' The artifical dataset of a polytomous responses (five categories) which contains of three subsets  which represent different dimensions.
#' The three dimensions are uncorrelated. To reproduce this dataset: \cr
#' \code{poly_md_dset <- generate_data(responseType = "multidim.nocorrel", ncat = 5)}
#' \cr\cr will lead to similar but not the same dataset, due to the randomization.
#'
#' @docType data
#'
#' @usage data(poly_md_dset)
#'
#' @rdname poly_md_dset
"poly_md_dset"
