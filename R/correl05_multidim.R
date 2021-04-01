#' The Multidimensional Polytomous Dataset with 0.5 Correlation
#'
#' The artifical dataset of a polytomous responses (five categories) which contains of two subsets  which represent different dimensions.
#' Both dimensions are having 0.5 correlation to each other. To reproduce this dataset: \cr
#' \code{correl05_multidim <- generate_data(responseType = "multidim.withcorrel", corLevel = 0.5)}
#' \cr\cr will lead to similar but not the same dataset, due to the randomization.
#'
#' @docType data
#'
#' @usage data(correl05_multidim)
#'
#' @rdname correl05_multidim
"correl05_multidim"
