#' The Multidimensional Polytomous Dataset with 0.6 Correlation
#'
#' The artifical dataset of a polytomous responses (five categories) which contains of two subsets  which represent different dimensions.
#' Both dimensions are having 0.6 correlation to each other. To reproduce this dataset: \cr
#' \code{correl06_multidim <- generate_data(responseType = "multidim.withcorrel", corLevel = 0.6)}
#' \cr\cr will lead to similar but not the same dataset, due to the randomization.
#'
#' @docType data
#'
#' @usage data(correl06_multidim)
#'
#' @rdname correl06_multidim
"correl06_multidim"
