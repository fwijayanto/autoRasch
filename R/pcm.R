#' Estimation of The Partial Credit Model (PCM)
#'
#' This function computes the parameter estimates of a partial credit model for dichotomous and polytomous responses
#' by using penalized JML estimation.
#'
#' @param X Input dataset as matrix or data frame with ordinal responses (starting from 0);
#' rows represent individuals, columns represent items.
#' @param isHessian a logical parameter setting whether or not the Hessian matrix is needed to be computed.
#' The Hessian matrix is needed to compute the separation reliability.
#'
#' @return
#' \strong{\code{pcm()} will return a \code{\link[base:list]{list}} which contains:}
#' \item{X}{   The dataset that is used for estimation.}
#' \item{mt_vek}{   A vector of the highest response category as many as the number of items.}
#' \item{real_vek}{   A vector of the presence of thresholds in the items-thresholds order: \code{1} if the threshold is present and \code{NA} if unavailable.}
#' \item{itemName}{   The vector of names of items (columns) in the dataset.}
#' \item{loglik}{   The log likelihood of the estimation.}
#' \item{hessian}{   The hessian matrix. Only when the \code{isHessian = TRUE}.}
#' \item{beta}{   A vector of the difficulty parameter of each categories of items (thresholds).}
#' \item{theta}{   A vector of the ability parameters of each individuals.}
#'
#' @seealso \code{\link{pcm}}, \code{\link{gpcm}}
#'
#' @examples
#' res <- pcm(poly_inh_dset)
#' res
#' summary(res)
#'
#' @rdname pcm
#' @export
pcm <- function(X, isHessian = TRUE){

  result <- pjmle(X = X, fixed_par = c("gamma"), isPenalized_gamma = FALSE, isHessian = isHessian)
  class(result) <- c("armodels","pcm","autoRasch")
  return(result)

}


#' Fit Statistics PCM
#'
#' \code{fitStats} compute the fit statistics (e.g., Outfit and Infit) of the PCM model estimation (items and persons).
#'
#' @param obj The object of class \code{'pcm'}.
#' @param isAlpha Boolean value that indicates whether the discrimination parameters is needed to be estimated or not.
#' The discrimination parameters are estimated using the corresponding models (GPCM).
#'
#' @return
#' \strong{\code{fitStats()} will return a \code{\link[base:list]{list}} which contains:}
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
#' plot(fit_res, plot.x = "gamma", plot.y = "outfit")
#'
#' @rdname pcm
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


#' @param par The parameter that are wanted to be summarized.
#'
#' @rdname pcm
#' @export
summary.pcm <- function(obj, par = c()){

  if(is.null(par) | "theta" %in% par){
    cat("\n\n")
    cat("The estimated ability scores:")
    cat("\n")
    print(obj$theta)
    cat("\n")
    cat("The highest ability score: ",round(max(obj$theta,na.rm = TRUE),4))
    cat("\n")
    cat("The lowest ability score: ",round(min(obj$theta,na.rm = TRUE),4))
  }

  if(is.null(par) | "beta" %in% par){
    cat("\n\n")
    cat("The estimated difficulty scores:")
    cat("\n")
    reported_beta <- obj$beta * obj$real_vek
    beta_mat <- matrix(reported_beta, nrow = length(obj$mt_vek), byrow = TRUE)
    beta_mat <- as.data.frame(round(beta_mat,4), row.names = obj$itemName)
    colnames(beta_mat) <- paste("Th_",c(1:max(obj$mt_vek)),sep = "")
    beta_mat[["Item Loc."]] <- round(apply(beta_mat,1,mean,na.rm=TRUE),4)
    beta_mat$` ` <- apply(beta_mat[,1:max(obj$mt_vek)],1,function(x){if(is.unsorted(na.omit(x))){return("*")}else{return("")}})
    print(beta_mat, quote = FALSE)
    cat("\n")
    cat("The most difficult item: ",obj$itemName[which(beta_mat[,5] == max(beta_mat[,5],na.rm = TRUE))])
    cat("\n")
    cat("The easiest item: ",obj$itemName[which(beta_mat[,5] == min(beta_mat[,5],na.rm = TRUE))])
    cat("\n")
    ntd_items <- length(which(beta_mat[,ncol(beta_mat)] == "*"))
    cat("There are",ntd_items,"items which have disordered thresholds.")
    cat("\n")
    cat("'*' Item has disordered thresholds.")
  }
}

#' @rdname pcm
#' @export
print.pcm <- function(obj, par = c()){
  cls <- class(obj)
  class(obj) <- "list"

  if(is.null(par) | "theta" %in% par){
    cat("\n")
    cat("$theta")
    cat("\n")
    print(obj$theta)
    cat("\n")
  }

  if(is.null(par) | "beta" %in% par){
    cat("\n")
    cat("$beta")
    cat("\n")
    print(obj$beta)
    cat("\n")
  }

}
