#' Unidimensionality Check
#'
#' This function checks the unidimensionality status using the confirmatory factor analysis.
#'
#' @param x The dataset of responses.
#' @param is.polychor A boolean parameter to set whether the dataset is categorical or not.
#' @param se The standard error type. Please refer to \code{\link[lavaan:cfa]{cfa()}}.
#' @param estimator The type of the estimator. Please refer to \code{\link[lavaan:cfa]{cfa()}}.
#' @param test The test used in the factor analysis. Please refer to \code{\link[lavaan:cfa]{cfa()}}.
#'
#' @return A list of the CFA output and the some of the goodness-of-fit indices (i.e., cfi, tli, rmsea, and srmr)
#'
#' @importFrom lavaan cfa
#' @importFrom lavaan fitMeasures
#'
#' @rdname unidimensionality
#' @export
check.unidim <- function(x, is.polychor = TRUE, se = "robust", estimator = "WLSMV", test="Satorra.Bentler"){

  resid <- x

  if(is.null(dim(resid))){
    resid <- matrix(resid,ncol = 1,dimnames = list(c(1:length(resid)), c("item")))
  }
  var.name <- names(as.data.frame(resid))
  var.names <- substr(paste(paste(var.name,"+",sep = ""),collapse = ""),1,nchar(paste(paste(var.name,"+",sep = ""),collapse = ""))-1)

  model <- paste('dimension  =~',var.names)
  if(is.polychor){
    unidim <- cfa(model, data=resid, ordered = c(var.name), missing = "listwise", se = se, estimator = estimator, test=test,
                  optim.method = "nlminb",optim.force.converged = TRUE, optim.dx.tol = 1e-3,check.gradient = FALSE,
                  bootstrap = 1000, start = "Mplus")
  }else {
    unidim <- cfa(model, data=resid, missing = "pairwise", se = se, estimator = estimator, test=test)
  }

  fit.unidim <- fitMeasures(unidim, c("cfi","tli","rmsea","cfi.scaled","tli.scaled","rmsea.scaled","cfi.robust",
                                      "tli.robust","rmsea.robust","srmr","chisq","pvalue"))

  cfa.cfi <- fit.unidim[1]
  cfa.tli <- fit.unidim[2]
  cfa.rmsea <- fit.unidim[3]
  cfa.cfi.scaled <- fit.unidim[4]
  cfa.tli.scaled <- fit.unidim[5]
  cfa.rmsea.scaled <- fit.unidim[6]
  cfa.cfi.robust <- fit.unidim[7]
  cfa.tli.robust <- fit.unidim[8]
  cfa.rmsea.robust <- fit.unidim[9]
  cfa.srmr <- fit.unidim[10]
  cfa.chisq <- fit.unidim[11]
  cfa.pvalue <- fit.unidim[12]

  MLplot <- list("unidim" = unidim,"cfi" = cfa.cfi,"tli" = cfa.tli,"rmsea" = cfa.rmsea,"cfi.scaled" = cfa.cfi.scaled,
                 "tli.scaled" = cfa.tli.scaled,"rmsea.scaled" = cfa.rmsea.scaled,"cfi.robust" = cfa.cfi.robust,
                 "tli.robust" = cfa.tli.robust,"rmsea.robust" = cfa.rmsea.robust,"srmr" = cfa.srmr,"chisq"=cfa.chisq,"pvalue" = cfa.pvalue)


  return(MLplot)

}
