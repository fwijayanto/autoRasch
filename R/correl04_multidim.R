#' The Multidimensional Polytomous Dataset with 0.4 Correlation
#'
#' The artifical dataset of a polytomous responses (five categories) which contains of two subsets  which represent different dimensions.
#' Both dimensions are having 0.4 correlation to each other. To reproduce this dataset: \cr
#' \code{correl04_multidim <- generate_data(responseType = "multidim.withcorrel", corLevel = 0.4)}
#' \cr\cr will lead to similar but not the same dataset, due to the randomization.
#'
#' @docType data
#'
#' @usage data(correl04_multidim)
#'
#' @rdname correl04_multidim
"correl04_multidim"
