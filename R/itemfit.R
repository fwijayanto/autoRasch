
fitStats <- function (obj, isAlpha = TRUE) {
  if(!("pcm" %in% class(obj))){
    stop("autoRasch ERROR: itemfit is only for pcm.")
  }
  UseMethod("fitStats", obj)
}


#' @rdname fit
#' @export
summary.fit <- function(obj){

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
  print(i.mat)
  cat("\n\n")

  p.mat <- cbind(obj$p.fit$p.outfitMSQ, obj$p.fit$p.infitMSQ, obj$p.fit$p.outfitZ, obj$p.fit$p.infitZ)
  dimnames(p.mat) <- list(c(paste("P",c(1:length(obj$p.fit$p.outfitMSQ)),sep = "")), c("OutfitMSQ","InfitMSQ","OutfitZ","InfitZ"))
  cat("Person Fit Statistics:")
  cat("\n\n")
  print(p.mat)
  cat("\n")

}

#' @rdname fit
#' @export
plot.fit <- function(obj, type = "n", plotx, ploty, xlab = NULL, ylab = NULL, use.name = FALSE, ...){

  dotdotdot <- list(...)
  if(!is.null(dotdotdot$main)){
    par(mar = c(6.5, 7.5, 2.5, 1), oma = c(0, 0,0, 0))
  } else {
    par(mar = c(6.5, 7.5, 0.5, 1), oma = c(0, 0,0, 0))
  }

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
  } else if(plotx == "gamma"){
    plotx <- obj$alpha
    # x.lab <- expression(paste("alpha (",alpha,")"))
    # x.lab <- expression(paste("estimated ",alpha[i]))
    x.lab <- expression(alpha)
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
  } else if(ploty == "gamma"){
    ploty <- obj$alpha
    y.lab <- expression(alpha)
  }


  if(use.name){
    text <- names(obj$i.fit$i.outfitMSQ)
  }


  if(!is.null(xlab)){
    x.lab <- xlab
  }

  if(!is.null(ylab)){
    y.lab <- ylab
  }

  suppressWarnings(plot(plotx, ploty, type = type, xlab = x.lab, ylab = y.lab, mgp =c(5,2,0), ... = ...))

  if(use.name){
    suppressWarnings(text(plotx, ploty, labels = text, ... = ...))
  } else {
    suppressWarnings(text(plotx, ploty, ... = ...))
  }

}
