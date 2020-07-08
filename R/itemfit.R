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
  if(!("pcm" %in% class(obj)) & !("pcmdif" %in% class(obj))){
    stop("autoRasch ERROR: itemfit is only for rasch and pcm object.")
  }
  UseMethod("fitStats", obj)
}

#' @param object The object from class \code{fit}. The item fit statistics results.
#' @param ... further argument passed or from other method.
#'
#' @rdname fit
#' @export
summary.fit <- function(object, ...){

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
    plotx <- x$alpha
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
    ploty <- x$alpha
    y.lab <- expression(alpha)
  }


  if(use.name){
    text <- names(x$i.fit$i.outfitMSQ)
  }


  if(!is.null(xlab)){
    x.lab <- xlab
  }

  if(!is.null(ylab)){
    y.lab <- ylab
  }

  suppressWarnings(plot(plotx, ploty, type = type, xlab = x.lab, ylab = y.lab, ... = ...))

  if(use.name){
    suppressWarnings(text(plotx, ploty, labels = text, ... = ...))
  } else {
    suppressWarnings(text(plotx, ploty, ... = ...))
  }

}
