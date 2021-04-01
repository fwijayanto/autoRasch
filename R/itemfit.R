<<<<<<< HEAD
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
=======
#' FIt statistics
#'
#' @param obj The object from class \code{pcm}. The result of the parameter estimation using the PCM.
#' @param isAlpha Boolean value whether discrimination parameter is needed or not.
#'
#' @examples
#' res <- pcm(pcm_data)
#' fit_res <- fitStats(res)
#' summary(fit_res)
#' plot(fit_res, plot.x = "gamma", plot.y = "outfit")
#'
#' @rdname fit
#' @export
fitStats <- function (obj, isAlpha = TRUE) {
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
  if(!("pcm" %in% class(obj)) & !("pcmdif" %in% class(obj))){
    stop("autoRasch ERROR: itemfit is only for rasch and pcm object.")
  }
  UseMethod("fitStats", obj)
}

<<<<<<< HEAD

#' @param object The object of class \code{'fit'}.
#' @param ... Further arguments to be passed.
=======
#' @param object The object from class \code{fit}. The item fit statistics results.
#' @param ... further argument passed or from other method.
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
#'
#' @rdname fit
#' @export
summary.fit <- function(object, ...){

<<<<<<< HEAD
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
=======
  dotdotdot <- list(...)

  i.mat <- cbind(object$i.fit$i.outfitMSQ, object$i.fit$i.infitMSQ, object$i.fit$i.outfitZ, object$i.fit$i.infitZ)
  if(!is.null(object$alpha)){
    i.mat <- cbind(i.mat, object$alpha)
    dimnames(i.mat) <- list(c(names(object$i.fit$i.outfitMSQ)),c("OutfitMSQ","InfitMSQ","OutfitZ","InfitZ","Discrimination"))
  } else {
    dimnames(i.mat) <- list(c(names(object$i.fit$i.outfitMSQ)),c("OutfitMSQ","InfitMSQ","OutfitZ","InfitZ"))
  }
  cat("\n")
  cat("Item Fit Statistics:")
  cat("\n\n")
  print(i.mat, ... = ...)
  cat("\n\n")

  p.mat <- cbind(object$p.fit$p.outfitMSQ, object$p.fit$p.infitMSQ, object$p.fit$p.outfitZ, object$p.fit$p.infitZ)
  dimnames(p.mat) <- list(c(paste("P",c(1:length(object$p.fit$p.outfitMSQ)),sep = "")), c("OutfitMSQ","InfitMSQ","OutfitZ","InfitZ"))
  cat("Person Fit Statistics:")
  cat("\n\n")
  print(p.mat, ... = ...)
  cat("\n")

}

#' @param x The object which is needed to be plot.
#' @param plotx The statistics that is wanted to be plot in the x axis.
#' @param ploty The statistics that is wanted to be plot in the y axis.
#' @param type The type of the plot.
#' @param xlab Thelabel of the x axis.
#' @param ylab Thelabel of the y axis.
#' @param use.name Boolean value whether the plot using the variable names or not.
#' @param fileOutput Boolean value whether the plot is saved to a file or not.
#'
#' @rdname fit
#' @export
plot.fit <- function(x, plotx = "alpha", ploty = "outfit", type = "n", xlab = NULL, ylab = NULL, use.name = FALSE, fileOutput = TRUE, ...){

  dotdotdot <- list(...)
  # if(!is.null(dotdotdot$main)){
  #   par(mar = c(6.5, 7.5, 2.5, 1), oma = c(0, 0,0, 0))
  # } else {
  #   par(mar = c(6.5, 7.5, 0.5, 1), oma = c(0, 0,0, 0))
  # }

>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5

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
    plotx <- x$i.fit$i.outfitMSQ
    x.lab <- "Outfit"
  } else if(plotx == "infit"){
    plotx <- x$i.fit$i.infitMSQ
    x.lab <- "Infit"
  } else if(plotx == "outfitz"){
    plotx <- x$i.fit$i.outfitZ
    x.lab <- "outfitZ"
  } else if(plotx == "infitz"){
    plotx <- x$i.fit$i.infitZ
    x.lab <- "infitZ"
  } else if(plotx == "alpha"){
<<<<<<< HEAD
    plotx <- obj$alpha
=======
    plotx <- x$alpha
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
    # x.lab <- expression(paste("alpha (",alpha,")"))
    # x.lab <- expression(paste("estimated ",alpha[i]))
    x.lab <- expression(hat(alpha))
  }
  if(ploty == "outfit"){
    ploty <- x$i.fit$i.outfitMSQ
    y.lab <- "Outfit"
  } else if(ploty == "infit"){
    ploty <- x$i.fit$i.infitMSQ
    y.lab <- "Infit"
  } else if(ploty == "outfitz"){
    ploty <- x$i.fit$i.outfitZ
    y.lab <- "outfitZ"
  } else if(ploty == "infitz"){
    ploty <- x$i.fit$i.infitZ
    y.lab <- "infitZ"
  } else if(ploty == "alpha"){
<<<<<<< HEAD
    ploty <- obj$alpha
    y.lab <- expression(hat(alpha))
  }


  if(useName){
    text <- names(obj$i.fit$i.outfitMSQ)
=======
    ploty <- x$alpha
    y.lab <- expression(alpha)
  }


  if(use.name){
    text <- names(x$i.fit$i.outfitMSQ)
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
  }


  if(!is.null(dotdotdot$xlab)){
    x.lab <- dotdotdot$xlab
  }

  if(!is.null(dotdotdot$ylab)){
    y.lab <- dotdotdot$ylab
  }

<<<<<<< HEAD
  if(!is.null(dotdotdot$xlab) | !is.null(dotdotdot$ylab)){
    plot(plotx, ploty, ... = ...)
  } else {
    suppressWarnings(plot(plotx, ploty, xlab = x.lab, ylab = y.lab, ... = ...))
  }
=======
  suppressWarnings(plot(plotx, ploty, type = type, xlab = x.lab, ylab = y.lab, ... = ...))
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5

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
