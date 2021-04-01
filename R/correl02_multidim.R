#' The Multidimensional Polytomous Dataset with 0.2 Correlation
#'
#' The artifical dataset of a polytomous responses (five categories) which contains of two subsets  which represent different dimensions.
#' Both dimensions are having 0.2 correlation to each other. To reproduce this dataset: \cr
#' \code{correl02_multidim <- generate_data(responseType = "multidim.withcorrel", corLevel = 0.2)}
#' \cr\cr will lead to similar but not the same dataset, due to the randomization.
#'
#' @docType data
#'
#' @usage data(correl02_multidim)
#'
#' @rdname correl02_multidim
"correl02_multidim"
