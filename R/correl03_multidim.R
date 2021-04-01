#' The Multidimensional Polytomous Dataset with 0.3 Correlation
#'
#' The artifical dataset of a polytomous responses (five categories) which contains of two subsets  which represent different dimensions.
#' Both dimensions are having 0.3 correlation to each other. To reproduce this dataset: \cr
#' \code{correl03_multidim <- generate_data(responseType = "multidim.withcorrel", corLevel = 0.3)}
#' \cr\cr will lead to similar but not the same dataset, due to the randomization.
#'
#' @docType data
#'
#' @usage data(correl03_multidim)
#'
#' @rdname correl03_multidim
"correl03_multidim"
