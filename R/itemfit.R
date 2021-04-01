#' Fit statistics
#'
#' The goodness-of-fit statistics of Rasch analysis for items and persons. It consists of Outfit (Unweighted) Mean Square,
#' Infit (Weighted) Mean Square, Outfit ZSTD (Standardized Unweighted Mean Square), and Outfit ZSTD (Standardized Weighted Mean Square)
#'
#' @param obj The object of class \code{'pcm'} or \code{'pcmdif'}.
#' @param isAlpha Boolean value that indicates whether the discrimination parameters is needed to be estimated or not.
#' The discrimination parameters are estimated using the corresponding models (GPCM or GPCM-DIF).
#' @param isTraced A list of some matrices, i.e., the expected values, the variances, the curtosis, and the standardized residual matrix.
#'
#' @rdname fit
#' @export
fitStats <- function (obj, isAlpha = TRUE, isTraced = FALSE) {
  if(!("pcm" %in% class(obj)) & !("pcmdif" %in% class(obj))){
    stop("autoRasch ERROR: itemfit is only for rasch, pcm and pcmdif object.")
  }
  UseMethod("fitStats", obj)
}


#' @param object The object of class \code{'fit'}.
#' @param ... Further arguments to be passed.
#'
#' @rdname fit
#' @export
summary.fit <- function(object, ...){

  obj <- object

  dotdotdot <- list(...)

  if(!is.null(dotdotdot$type)){
    type <- dotdotdot$type
  } else {
    type <- NULL
  }

  if(is.null(type) | "item" %in% type){
    i.mat <- cbind(obj$i.fit$i.outfitMSQ, obj$i.fit$i.infitMSQ, obj$i.fit$i.outfitZ, obj$i.fit$i.infitZ)
    if(!is.null(obj$alpha)){
      i.mat <- cbind(i.mat, obj$alpha)
      dimnames(i.mat) <- list(c(names(obj$i.fit$i.outfitMSQ)),c("OutfitMSQ","InfitMSQ","OutfitZ","InfitZ","Discrimination"))
    } else {
      dimnames(i.mat) <- list(c(names(obj$i.fit$i.outfitMSQ)),c("OutfitMSQ","InfitMSQ","OutfitZ","InfitZ"))
    }
    cat("\n")
    cat("Item Fit Statistics:")
    cat("\n\n")
    print(round(i.mat,2))
    cat("\n\n")
  }

  if(is.null(type) | "person" %in% type){
    p.mat <- cbind(obj$p.fit$p.outfitMSQ, obj$p.fit$p.infitMSQ, obj$p.fit$p.outfitZ, obj$p.fit$p.infitZ)
    dimnames(p.mat) <- list(c(paste("P",c(1:length(obj$p.fit$p.outfitMSQ)),sep = "")), c("OutfitMSQ","InfitMSQ","OutfitZ","InfitZ"))
    cat("Person Fit Statistics:")
    cat("\n\n")
    print(round(p.mat,2))
    cat("\n")
  }

}

#' @param objFit The object of class \code{'fit'}.
#'
#' @rdname fit
#' @export
itemfit <- function(objFit){
  summary(objFit, type = "item")
}

#' @rdname fit
#' @export
personfit <- function(objFit){
  summary(objFit, type = "person")
}

#' @param toPlot An array with length two \code{c(x,y)}, to choose what to plot. There are five options to plot, which are alpha, outfit, infit, outfitz, and infitz
#' @param  useName A logical statement whether the name of the variable are going to be used in the plot instead of the variable order.
#'
#' @rdname fit
#' @export
plot_fitStats <- function(objFit, toPlot = c("alpha","infit"), useName = FALSE, ...){

  obj <- objFit

  if(!"fit" %in% class(obj) & !"autoRasch" %in% class(obj)){
    stop("The input should be an object of class 'fit'")
  }

  dotdotdot <- list(...)

  if(!is.null(dotdotdot$use.name)){
    use.name <- dotdotdot$use.name
  } else {
    use.name <- FALSE
  }

  if(!is.null(dotdotdot$font.text)){
    font.text <- dotdotdot$font.text
  } else {
    font.text <- 1
  }

  if(length(toPlot) < 2){
    plotx <- "alpha"
    ploty <- "infit"
  } else {
    plotx <- toPlot[1]
    ploty <- toPlot[2]
  }

  # if(!is.null(dotdotdot$main)){
  #   par(mar = c(6.5, 7.5, 2.5, 1), oma = c(0, 0,0, 0))
  # } else {
  #   par(mar = c(6.5, 7.5, 0.5, 1), oma = c(0, 0,0, 0))
  # }

  if(plotx == "outfit"){
    plotx <- obj$i.fit$i.outfitMSQ
    x.lab <- "Outfit"
  } else if(plotx == "infit"){
    plotx <- obj$i.fit$i.infitMSQ
    x.lab <- "Infit"
  } else if(plotx == "outfitz"){
    plotx <- obj$i.fit$i.outfitZ
    x.lab <- "outfitZ"
  } else if(plotx == "infitz"){
    plotx <- obj$i.fit$i.infitZ
    x.lab <- "infitZ"
  } else if(plotx == "alpha"){
    plotx <- obj$alpha
    # x.lab <- expression(paste("alpha (",alpha,")"))
    # x.lab <- expression(paste("estimated ",alpha[i]))
    x.lab <- expression(hat(alpha))
  }
  if(ploty == "outfit"){
    ploty <- obj$i.fit$i.outfitMSQ
    y.lab <- "Outfit"
  } else if(ploty == "infit"){
    ploty <- obj$i.fit$i.infitMSQ
    y.lab <- "Infit"
  } else if(ploty == "outfitz"){
    ploty <- obj$i.fit$i.outfitZ
    y.lab <- "outfitZ"
  } else if(ploty == "infitz"){
    ploty <- obj$i.fit$i.infitZ
    y.lab <- "infitZ"
  } else if(ploty == "alpha"){
    ploty <- obj$alpha
    y.lab <- expression(hat(alpha))
  }


  if(useName){
    text <- names(obj$i.fit$i.outfitMSQ)
  }


  if(!is.null(dotdotdot$xlab)){
    x.lab <- dotdotdot$xlab
  }

  if(!is.null(dotdotdot$ylab)){
    y.lab <- dotdotdot$ylab
  }

  if(!is.null(dotdotdot$xlab) | !is.null(dotdotdot$ylab)){
    plot(plotx, ploty, ... = ...)
  } else {
    suppressWarnings(plot(plotx, ploty, xlab = x.lab, ylab = y.lab, ... = ...))
  }

  if(useName){
    suppressWarnings(text(plotx, ploty, labels = text, font = font.text))
  } else {
    suppressWarnings(text(plotx, ploty, font = font.text))
  }

}


#' Residual Correlation
#'
#' Compute the correlation of the standardized residual to check the local dependency status
#'
#' @param objFit object of class "fit", the output of \code{fitStats()}.
#'
#' @return
#' \item{ld_correl}{  Correlation matrix of the standradized residual.}
#' \item{ld_mean}{  The mean of the correlation.}
#' \item{ld_lowertri}{  The lower triangle of the correlation matrix.}
#'
#' @rdname ld
#'
#' @export
resid_corr <- function(objFit){
  # if(!is.object(objFit$traceMat)){
  #   stop("Please compute Fitness object using isTrace = TRUE!")
  # }
  corLD <- cor(objFit$traceMat$std.res, use = "pairwise.complete.obs")
  ld <- mean((corLD[lower.tri(corLD)]), na.rm = TRUE)
  return(list("ld_correl" = corLD, "ld_mean" = ld, "ld_lowertri" = corLD[lower.tri(corLD)]))
}
