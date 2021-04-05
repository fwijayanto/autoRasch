#' The Within-item Multidimensional Polytomous Dataset
#'
#' The artifical dataset of a polytomous responses (five categories) which contains of three subsets which represent different dimensions.
#' The three dimensions are uncorrelated, however, some of items relate to more than one dimension. To reproduce this dataset: \cr
#' \code{withinItem_multidim <- generate_data(responseType = "multidim.within", ndim = 3, dim.members = list(c(1:6,13),c(3,7:12),c(5,13:18)))}
#' \cr\cr will lead to similar but not the same dataset, due to the randomization.
#'
#' @docType data
#'
#' @usage data(withinItem_multidim)
#'
#' @rdname withinItem_multidim
"withinItem_multidim"
