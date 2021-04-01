#' The Multi-testlets Polytomous Dataset
#'
#' The artifical dataset of a polytomous responses (five categories) which contains of two subsets which represent different testlets.
#' Both testlets are having difference variance effects. To reproduce this dataset: \cr
#' \code{testlets_dataset <- generate_data(responseType = "testlets", ndim = 2, sdlambda = c(0,4))}
#' \cr\cr will lead to similar but not the same dataset, due to the randomization.
#'
#' @docType data
#'
#' @usage data(testlets_dataset)
#'
#' @rdname testlets_dataset
"testlets_dataset"
