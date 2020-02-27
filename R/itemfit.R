#' Fit Statistics
#'
#' Compute the fit statistics (Outfit, Infit, etc.) of the PCM dan PCM-DIF model estimation (items and persons)
#'
#' @param obj The object which is resulted from models estimation (PCM or PCM-DIF)
#' @param isAlpha Boolean value. Are the discrimination parameters needed to be estimated using the correspond extension of the 2PL models (GPCM or GPCM-DIF, respectively).
#'
#' @return
#' List of values
#'
#' \item{alpha}{   A vector of estimated discrimination parameters for each items.}
#' \emph{i.fit}{   Item fit statistics.}
#' \itemize{
#'    \item{i.outfitMSQ}{   A vector of Outfit mean square values for each items.}
#'    \item{i.infitMSQ}{   A vector of Infit mean square values for each items.}
#'    \item{i.outfitZ}{   A vector of OutfitZ values for each items.}
#'    \item{i.infitZ}{   A vector of InfitZ values for each items.}
#' }
#' \emph{p.fit}{   Person fit statistics.}
#' \itemize{
#'    \item{p.outfitMSQ}{   A vector of Outfit mean square values for each persons.}
#'    \item{p.infitMSQ}{   A vector of Infit mean square values for each persons.}
#'    \item{p.outfitZ}{   A vector of OutfitZ values for each persons.}
#'    \item{p.infitZ}{   A vector of InfitZ values for each persons.}
#' }
#'
#'
#' @examples
#' pcmdif_res <- pcm_dif(polydif_inh_dset, groups_map = c(rep(1,245),rep(0,245)))
#' fit_res <- fitStats(pcmdif_res)
#' summary(fit_res)
#' plot(fit_res, plotx = "gamma", ploty = "outfit")
#'
#' @rdname fit
#' @export
fitStats <- function (obj, isAlpha = TRUE) {
  if(!("pcm" %in% class(obj)) & !("pcmdif" %in% class(obj))){
    stop("autoRasch ERROR: itemfit is only for rasch, pcm and pcmdif object.")
  }
  UseMethod("fitStats", obj)
}

#' @rdname fit
#' @export
fitStats.pcmdif <- function(obj, isAlpha = TRUE){

  X <- obj$X
  groups_map <- as.matrix(obj$groups_map)

  # Map the parameters
  theta <- obj$theta
  beta <- obj$beta
  deltabeta <- obj$deltabeta
  mt_vek <- obj$mt_vek
  n.th <- max(mt_vek)


  deltabeta.tot <- 0
  for(i in 1:ncol(groups_map)){
    deltabeta.tot <- deltabeta.tot + outer(deltabeta[(((i-1)*ncol(X))+1):(i*ncol(X))],groups_map[,i],"*")
  }
  deltabeta.tot.rep <- rep.int((deltabeta.tot), rep.int(obj$mt_vek,nrow(groups_map)))       #deltabeta.tot.rep is total deltabeta which has been replicated to every categoory

  xna.mat <- matrix(1,nrow = nrow(X), ncol = ncol(X))

  idx <- which(is.na(X))
  xna.mat[idx] <- NA
  XNA <- xna.mat

  t.diff <- outer((-beta), theta, "+")
  t.diff <- t.diff - deltabeta.tot.rep
  disc.diff <- t.diff

  per.cat.list <- matrix(disc.diff, nrow = n.th)

  temp.prob <- as.matrix(per.cat.list[1,])
  temp.l2 <- exp(temp.prob)
  temp.l1 <- exp(temp.l2*0)
  temp.l1 <- cbind(temp.l1,temp.l2)
  for(i in 2:n.th){
    temp.prob <- cbind(temp.prob,(temp.prob[,i-1]+per.cat.list[i,]))
    temp.l1 <- cbind(temp.l1,(exp(temp.prob[,i])))
    temp.l2 <- temp.l2 + (exp(temp.prob[,i]))
  }
  l2 <- (temp.l2+1)

  l1 <- as.vector(t(temp.l1))
  l2 <- rep(l2, each = (n.th+1))

  pmat <- l1/l2

  mt_vek0 <- mt_vek + 1
  mt_seq <- sequence(mt_vek0)-1

  Emat <- pmat * mt_seq
  Emat <- matrix(Emat, nrow = (n.th+1))
  Emat <- colSums(Emat)

  Emat.cat <- rep(Emat, each = (n.th+1))

  Vvect.cat <- ((mt_seq - Emat.cat)^2)*pmat
  Vmat.cat <- matrix(Vvect.cat, nrow = (n.th+1))
  Vmat <- colSums(Vmat.cat)

  Cvect.cat <- ((mt_seq - Emat.cat)^4)*pmat
  Cmat.cat <- matrix(Cvect.cat, nrow = (n.th+1))
  Cmat <- colSums(Cmat.cat)

  Emat <- t(matrix(Emat, nrow = ncol(X)))
  Vmat <- t(matrix(Vmat, nrow = ncol(X)))
  Cmat <- t(matrix(Cmat, nrow = ncol(X)))

  st.res <- (X-Emat)/sqrt(Vmat)
  sq.res <- st.res^2                            #squared standardized residuals
  ifit <- colSums(sq.res, na.rm = TRUE)
  pfit <- rowSums(sq.res, na.rm = TRUE)

  idf <- apply(X, 2, function(x) {length(na.exclude(x))})
  pdf <- apply(X, 1, function(x) {length(na.exclude(x))})

  i.outfitMSQ <- ifit/idf
  p.outfitMSQ <- pfit/pdf

  qsq.outfitMSQ <- (colSums(Cmat/Vmat^2, na.rm=TRUE)/idf^2) - 1/idf
  q.outfitMSQ <- sqrt(qsq.outfitMSQ)
  p.qsq.outfitMSQ <- (rowSums(Cmat/Vmat^2, na.rm=TRUE)/pdf^2) - 1/pdf
  p.q.outfitMSQ <- sqrt(p.qsq.outfitMSQ)

  isumVmat<-colSums(Vmat*XNA, na.rm = TRUE)
  psumVmat<-rowSums(Vmat*XNA, na.rm = TRUE)
  i.infitMSQ <- colSums(sq.res*Vmat, na.rm = TRUE)/isumVmat
  p.infitMSQ <- rowSums(sq.res*Vmat, na.rm = TRUE)/psumVmat

  qsq.infitMSQ <- colSums(Cmat-Vmat^2, na.rm=TRUE)/isumVmat^2
  q.infitMSQ <- sqrt(qsq.infitMSQ)
  p.qsq.infitMSQ <- rowSums(Cmat-Vmat^2, na.rm=TRUE)/psumVmat^2
  p.q.infitMSQ <- sqrt(p.qsq.infitMSQ)

  i.outfitZ <- (i.outfitMSQ^(1/3) - 1)*(3/q.outfitMSQ)+(q.outfitMSQ/3)
  i.infitZ  <- (i.infitMSQ^(1/3)  - 1)*(3/q.infitMSQ) +(q.infitMSQ/3)

  p.outfitZ <- (p.outfitMSQ^(1/3) - 1)*(3/p.q.outfitMSQ)+(p.q.outfitMSQ/3)
  p.infitZ  <- (p.infitMSQ^(1/3)  - 1)*(3/p.q.outfitMSQ) +(p.q.outfitMSQ/3)

  res_fit <- list("i.fit" = list("i.outfitMSQ" = i.outfitMSQ, "i.infitMSQ" = i.infitMSQ, "i.outfitZ" = i.outfitZ, "i.infitZ" = i.infitZ), "p.fit" = list("p.outfitMSQ" = p.outfitMSQ, "p.infitMSQ" = p.infitMSQ, "p.outfitZ" = p.outfitZ, "p.infitZ" = p.infitZ))

  if(isAlpha){
    gpcmdif_res <- gpcm_dif(X = X, groups_map = groups_map, isHessian = FALSE)
    res_fit[["alpha"]] <- exp(gpcmdif_res$gamma)
  }

  class(res_fit) <- c("fit","autoRasch")
  return(res_fit)
}

#' @rdname fit
#' @export
fitStats.pcm <- function(obj, isAlpha = TRUE){

  X <- obj$X

  # Map the parameters
  theta <- obj$theta
  beta <- obj$beta
  mt_vek <- obj$mt_vek
  n.th <- max(mt_vek)

  xna.mat <- matrix(1,nrow = nrow(X), ncol = ncol(X))

  idx <- which(is.na(X))
  xna.mat[idx] <- NA
  XNA <- xna.mat

  t.diff <- outer((-beta), theta, "+")
  disc.diff <- t.diff

  per.cat.list <- matrix(disc.diff, nrow = n.th)

  temp.prob <- as.matrix(per.cat.list[1,])
  temp.l2 <- exp(temp.prob)
  temp.l1 <- exp(temp.l2*0)
  temp.l1 <- cbind(temp.l1,temp.l2)
  for(i in 2:n.th){
    temp.prob <- cbind(temp.prob,(temp.prob[,i-1]+per.cat.list[i,]))
    temp.l1 <- cbind(temp.l1,(exp(temp.prob[,i])))
    temp.l2 <- temp.l2 + (exp(temp.prob[,i]))
  }
  l2 <- (temp.l2+1)

  l1 <- as.vector(t(temp.l1))
  l2 <- rep(l2, each = (n.th+1))

  pmat <- l1/l2

  mt_vek0 <- mt_vek + 1
  mt_seq <- sequence(mt_vek0)-1

  Emat <- pmat * mt_seq
  Emat <- matrix(Emat, nrow = (n.th+1))
  Emat <- colSums(Emat)

  Emat.cat <- rep(Emat, each = (n.th+1))

  Vvect.cat <- ((mt_seq - Emat.cat)^2)*pmat
  Vmat.cat <- matrix(Vvect.cat, nrow = (n.th+1))
  Vmat <- colSums(Vmat.cat)

  Cvect.cat <- ((mt_seq - Emat.cat)^4)*pmat
  Cmat.cat <- matrix(Cvect.cat, nrow = (n.th+1))
  Cmat <- colSums(Cmat.cat)

  Emat <- t(matrix(Emat, nrow = ncol(X)))
  Vmat <- t(matrix(Vmat, nrow = ncol(X)))
  Cmat <- t(matrix(Cmat, nrow = ncol(X)))

  st.res <- (X-Emat)/sqrt(Vmat)
  sq.res <- st.res^2                            #squared standardized residuals
  ifit <- colSums(sq.res, na.rm = TRUE)
  pfit <- rowSums(sq.res, na.rm = TRUE)

  idf <- apply(X, 2, function(x) {length(na.exclude(x))})
  pdf <- apply(X, 1, function(x) {length(na.exclude(x))})

  i.outfitMSQ <- ifit/idf
  p.outfitMSQ <- pfit/pdf

  qsq.outfitMSQ <- (colSums(Cmat/Vmat^2, na.rm=TRUE)/idf^2) - 1/idf
  q.outfitMSQ <- sqrt(qsq.outfitMSQ)
  p.qsq.outfitMSQ <- (rowSums(Cmat/Vmat^2, na.rm=TRUE)/pdf^2) - 1/pdf
  p.q.outfitMSQ <- sqrt(p.qsq.outfitMSQ)

  isumVmat<-colSums(Vmat*XNA, na.rm = TRUE)
  psumVmat<-rowSums(Vmat*XNA, na.rm = TRUE)
  i.infitMSQ <- colSums(sq.res*Vmat, na.rm = TRUE)/isumVmat
  p.infitMSQ <- rowSums(sq.res*Vmat, na.rm = TRUE)/psumVmat

  qsq.infitMSQ <- colSums(Cmat-Vmat^2, na.rm=TRUE)/isumVmat^2
  q.infitMSQ <- sqrt(qsq.infitMSQ)
  p.qsq.infitMSQ <- rowSums(Cmat-Vmat^2, na.rm=TRUE)/psumVmat^2
  p.q.infitMSQ <- sqrt(p.qsq.infitMSQ)

  i.outfitZ <- (i.outfitMSQ^(1/3) - 1)*(3/q.outfitMSQ)+(q.outfitMSQ/3)
  i.infitZ  <- (i.infitMSQ^(1/3)  - 1)*(3/q.infitMSQ) +(q.infitMSQ/3)

  p.outfitZ <- (p.outfitMSQ^(1/3) - 1)*(3/p.q.outfitMSQ)+(p.q.outfitMSQ/3)
  p.infitZ  <- (p.infitMSQ^(1/3)  - 1)*(3/p.q.outfitMSQ) +(p.q.outfitMSQ/3)

  res_fit <- list("i.fit" = list("i.outfitMSQ" = i.outfitMSQ, "i.infitMSQ" = i.infitMSQ, "i.outfitZ" = i.outfitZ, "i.infitZ" = i.infitZ), "p.fit" = list("p.outfitMSQ" = p.outfitMSQ, "p.infitMSQ" = p.infitMSQ, "p.outfitZ" = p.outfitZ, "p.infitZ" = p.infitZ))

  if(isAlpha){
    gpcm_res <- gpcm(X = X, isHessian = FALSE)
    res_fit[["alpha"]] <- exp(gpcm_res$gamma)
  }

  class(res_fit) <- c("fit","autoRasch")
  return(res_fit)
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
