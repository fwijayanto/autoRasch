#' Compute the In-plus-out-of-questionnaire log likelihood (IPOQ-LL)
#'
#' \code{compute_score} computes the the IPOQ-LL score of an instrument (included set) of the given initial survey.
#' While \code{compute_scores} computes the IPOQ-LL score of many (more than one) instruments (included sets) of
#' the given initial survey simultanously.
#'
#' @param X A matrix or data.frame of the observed responses (ordinal or binary response).
#' @param incl_set A vector of the items (columns) number in the data.frame X that are included in the included set.
#' @param type The type of the score. \code{ipoqll} if we ignore the presence of the DIF and \code{ipoqlldif} if we want to consider the DIF effect.
#' @param groups_map Matrix to map the respondents to the DIF groups.
#' @param init_par_iq Initial values of the parameters in the included set before the estimation begin.
#' @param init_par_oq Initial values of the parameters in the excluded set before the estimation begin.
#' @param optz_tuner_iq The optimisation setting of the included set. See \code{\link[stats:optim]{stats::optim()}} \code{control} parameter.
#' @param optz_tuner_oq The optimisation setting of the excluded set. See \code{\link[stats:optim]{stats::optim()}} \code{control} parameter.
#' @param isTracked A logical value whether the computation will be tracked or not.
#'
#' @return
#' \code{compute_score} will return a vector which contains in-questionnaire log likelihood (iqll), out-of-questionnaire log likelihood(oqll),
#' ipoqll, included set's items' number of the given initial survey, the estimated parameters of the included set,
#' and the estimated parameters of the excluded set, respectively.
#'
#' @examples
#' ipoqll_score <- compute_score(pcm_data, incl_set = c(4,5,6))
#' summary(ipoqll_score)
#'
#' @rdname compute_score
#' @export
compute_score <- function(X, incl_set, type = c("ipoqll"), groups_map = c(),
                          init_par_iq = c(), init_par_oq = c(), optz_tuner_iq = c(), optz_tuner_oq = c(),
                          isTracked = FALSE){

  if(is.null(type)){
    type <- "ipoqll"
  }

  if(type[1] == "ipoqll"){
    fixed_par <- c("deltabeta")
    isPenalized_deltabeta <- FALSE
  }

  dset <- as.matrix(X)
  incl_set <- incl_set[!is.na(incl_set)]
  incl_resp <- dset[,c(incl_set)]
  if(length(incl_set) != ncol(dset)){
    excl_resp <- dset[,-c(incl_set)]
  }

  iqll <- pjmle(incl_resp, init_par = init_par_iq, fixed_par = fixed_par,
                optz_tuner = optz_tuner_iq, groups_map = groups_map, isHessian = FALSE, isTracked = isTracked)

  if(ncol(dset) == length(incl_set)){
    loglik_oqll <- NA
  } else {
    oqll <- pjmle(excl_resp, init_par = init_par_oq, fixed_par = c("theta",fixed_par), fixed_theta = iqll$theta,
                  isPenalized_theta = FALSE, optz_tuner = optz_tuner_oq,
                  groups_map = groups_map, isHessian = FALSE, isTracked = isTracked)
    loglik_oqll <- oqll$loglik
  }

  ipoqll <- sum(c(iqll$loglik,loglik_oqll),na.rm = TRUE)
  res <- c(iqll$loglik, loglik_oqll, ipoqll)


  n_par <- sum(nrow(dset),((1+(max(dset,na.rm = TRUE)-min(dset,na.rm = TRUE)))*ncol(dset)))
  iqll_params <- c(iqll$theta, iqll$beta, iqll$gamma)
  if(ncol(dset) == length(incl_set)){
    oqll_params <- NA
  } else {
    oqll_params <- c(oqll$beta, oqll$gamma)
  }
  length(iqll_params) <- n_par
  length(oqll_params) <- n_par - nrow(dset)


  length(incl_set) <- ncol(dset)
  res <- c(res, incl_set, iqll_params, oqll_params)
  class(res) <- c("score",type[1])
  return(res)
}

compute_scores_unparalleled <- function(X, itemsets, type = c("ipoqll"), step_direct = c("fixed","forward","backward"), groups_map = c(),
                           init_par_iq = c(), init_par_oq = c(), optz_tuner_iq = c(), optz_tuner_oq = c()){


  dset <- as.matrix(X)
  incl_sets <- itemsets
  if(is.vector(incl_sets) & length(incl_sets) > ncol(dset)){
    stop("autoRasch ERROR: the number of items in the incl_set can not exceed the initial items.")
  }

  if(is.matrix(incl_sets) | is.null(step_direct)){
    step_direct <- "fixed"
  }

  if(step_direct == "forward"){
    excl_set <- c(1:ncol(dset))[-c(incl_sets)]
    add_items <- t(combn(excl_set,1))
    rep_itemsets <- matrix(rep.int(incl_sets,length(add_items)), nrow = length(add_items), byrow = TRUE)
    incl_sets <- cbind(rep_itemsets,add_items)
  } else if(step_direct == "backward"){
    incl_sets <- t(combn(incl_sets,(length(incl_sets)-1)))
  }

  incl_sets <- as.matrix(incl_sets)

# print("IN")

  scoreList <- foreach(i=1:nrow(incl_sets), .combine = rbind, .errorhandling = "remove", .packages = c("autoRasch"), .export = c()) %dopar% {

    incl_set <- incl_sets[i,]
    incl_set <- incl_set[!is.na(incl_set)]
    incl_set <- sort(incl_set,decreasing = FALSE)

    if(!is.null(init_par_iq) & !is.null(init_par_oq)){
      init_iq <- iqll_init(dset = dset, prev_incl_set = itemsets, prev_par_iq = init_par_iq, prev_par_oq = init_par_oq, incl_set = incl_set, direction = step_direct, type = type[1], groups_map = groups_map, iq_noise = 1e-3)
      init_oq <- oqll_init(dset = dset, prev_incl_set = itemsets, prev_par_iq = init_par_iq, prev_par_oq = init_par_oq, incl_set = incl_set, direction = step_direct, type = type[1], groups_map = groups_map, oq_noise = 1e-3)
    } else {
      init_iq <- c()
      init_oq <- c()
    }

    # cat("incl_set : ", incl_set)
    # cat("\n length.init_iq : ",length(init_iq))
    # cat("\n length.init_oq : ",length(init_oq),"\n")

    score_res <- compute_score(dset, incl_set = incl_set, type = type, groups_map = groups_map,
                              init_par_iq = init_iq, init_par_oq = init_oq,
                              optz_tuner_iq = optz_tuner_iq, optz_tuner_oq = optz_tuner_oq)

    length(incl_set) <- ncol(dset)
    res <- c(score_res)
    return(res)
  }

  res <- scoreList
  # class(res) <- c(paste(type,"s",sep = ""))
  return(res)
}

#' @param itemsets A matrix as a results of a \code{rbind} of \code{incl_set}.
#' @param step_direct The direction of the computation; \code{fixed}, \code{forward}, or \code{backward} .
#' @param cores Number of cores that is used in the paralellization.
#'
#' @return
#' \code{compute_scores} will return a matrix as a result of the \code{rbind} operation of the \code{compute_score}'s result.
#'
#' @rdname compute_score
#' @export
compute_scores <- function(X, itemsets, type = c("ipoqll"), step_direct = c("fixed","forward","backward"), groups_map = c(),
                           init_par_iq = c(), init_par_oq = c(), optz_tuner_iq = c(), optz_tuner_oq = c(), cores = NULL){



  itemsets <- as.matrix(itemsets)

  if(is.null(cores)){
    cores <- nrow(itemsets)
    if(cores > 2){
      cores <- 2
    }
  } else {
    if(cores > detectCores()){
      cores <- detectCores()
    }
  }


  cl <- makeCluster(cores)
  registerDoParallel(cl, cores = cores)

  scoreList <- compute_scores_unparalleled(X = X, itemsets = itemsets, type = type, step_direct = step_direct, groups_map = groups_map,
                                          init_par_iq = init_par_iq, init_par_oq = init_par_oq, optz_tuner_iq = optz_tuner_iq,
                                          optz_tuner_oq = optz_tuner_oq)

  stopCluster(cl)
  registerDoSEQ()

  res <- scoreList
  class(res) <- c(paste(type,"s",sep = ""))
  return(res)
}


iqll_init <- function(dset, prev_incl_set, prev_par_iq, prev_par_oq, incl_set, direction, type, groups_map, iq_noise){

  old.set.iq <- prev_incl_set
  new.set.iq <- incl_set
  old.dataset.iq <- as.matrix(dset[,old.set.iq])
  nn.th <- max(dset,na.rm = TRUE)-min(dset,na.rm = TRUE)
  if(!is.null(groups_map)){
    n.resp.th <- ncol(as.matrix(groups_map))
  }
  nlmPar.old.iq <- prev_par_iq
  nlmPar.new.iq <- c()

  old.set.oq <- c(1:ncol(dset))[-c(prev_incl_set)]
  new.set.oq <- c(1:ncol(dset))[-c(incl_set)]
  old.dataset.oq <- as.matrix(dset[,old.set.oq])
  nlmPar.old.oq <- prev_par_oq

  if(direction == "backward"){
    item.no.iq <- which(!(old.set.iq %in% new.set.iq))
    idx.remv.iq.beta <- c((nrow(old.dataset.iq) + ((nn.th*item.no.iq)-(nn.th-1))):(nrow(old.dataset.iq)+(nn.th*item.no.iq)))
    idx.remv.iq.gamma <- c(nrow(old.dataset.iq)+(ncol(old.dataset.iq)*nn.th)+item.no.iq)
    if(type[1] == "ipoqlldif"){
      idx.remv.iq.deltabeta <- c()
      for(i in 1: n.resp.th){
        idx.remv.iq.deltabeta <- c(idx.remv.iq.deltabeta,(nrow(old.dataset.iq)+((ncol(old.dataset.iq)*nn.th)+ncol(old.dataset.iq))+(ncol(old.dataset.iq)*(i-1))+item.no.iq))
      }
      idx.remv.iq <- c(idx.remv.iq.beta,idx.remv.iq.gamma,idx.remv.iq.deltabeta)
    } else {
      idx.remv.iq <- c(idx.remv.iq.beta,idx.remv.iq.gamma)
    }
    nlmPar.new.iq <- nlmPar.old.iq[-c(idx.remv.iq)]
    nlmPar.new.iq <- nlmPar.new.iq + (runif(length(nlmPar.new.iq),-1,1)*iq_noise)
  } else if(direction == "forward"){
    item.no.iq <- which(!(new.set.iq %in% old.set.iq))
    idx.remv.iq.beta <- c((nrow(old.dataset.iq) + ((nn.th*item.no.iq)-(nn.th-1))):(nrow(old.dataset.iq)+(nn.th*item.no.iq)))
    idx.remv.iq.gamma <- c(length(idx.remv.iq.beta)+nrow(old.dataset.iq)+(ncol(old.dataset.iq)*nn.th)+item.no.iq)
    if(type[1] == "ipoqlldif"){
      idx.remv.iq.deltabeta <- c()
      for(i in 1: n.resp.th){
        idx.remv.iq.deltabeta <- c(idx.remv.iq.deltabeta,(nrow(old.dataset.iq)+((ncol(old.dataset.iq)*nn.th)+ncol(old.dataset.iq))+(ncol(old.dataset.iq)*(i-1))+item.no.iq+length(idx.remv.iq.beta)+length(idx.remv.iq.gamma+(i-1))))
      }
      idx.remv.iq <- c(idx.remv.iq.beta,idx.remv.iq.gamma,idx.remv.iq.deltabeta)
    } else {
      idx.remv.iq <- c(idx.remv.iq.beta,idx.remv.iq.gamma)
    }
    nlmPar.new.iq <- insert.at(nlmPar.old.iq,idx.remv.iq,(length(nlmPar.old.iq)+length(idx.remv.iq)))

    item.no.oq <- which(!(old.set.oq %in% new.set.oq))
    idx.remv.oq.beta <- c((((nn.th*item.no.oq)-(nn.th-1))):((nn.th*item.no.oq)))
    idx.remv.oq.gamma <- c((ncol(old.dataset.oq)*nn.th)+item.no.oq)
    if(type[1] == "ipoqlldif"){
      idx.remv.oq.deltabeta <- c()
      for(i in 1: n.resp.th){
        idx.remv.oq.deltabeta <- c(idx.remv.oq.deltabeta,(((ncol(old.dataset.oq)*nn.th)+ncol(old.dataset.oq))+(ncol(old.dataset.oq)*(i-1))+item.no.oq))
      }
      idx.remv.oq <- c(idx.remv.oq.beta,idx.remv.oq.gamma,idx.remv.oq.deltabeta)
    } else {
      idx.remv.oq <- c(idx.remv.oq.beta,idx.remv.oq.gamma)
    }
    nlmPar.new.iq[idx.remv.iq] <- nlmPar.old.oq[idx.remv.oq]
    nlmPar.new.iq <- nlmPar.new.iq + (runif(length(nlmPar.new.iq),-1,1)*iq_noise)
  }

  return(nlmPar.new.iq)
}

oqll_init <- function(dset, prev_incl_set, prev_par_iq, prev_par_oq, incl_set, direction, type, groups_map, oq_noise){
  old.set.iq <- prev_incl_set
  new.set.iq <- incl_set
  old.dataset.iq <- as.matrix(dset[,old.set.iq])
  nn.th <- max(dset,na.rm = TRUE)-min(dset,na.rm = TRUE)
  if(!is.null(groups_map)){
    n.resp.th <- ncol(as.matrix(groups_map))
  }
  nlmPar.old.iq <- prev_par_iq

  old.set.oq <- c(1:ncol(dset))[-c(prev_incl_set)]
  new.set.oq <- c(1:ncol(dset))[-c(incl_set)]
  old.dataset.oq <- as.matrix(dset[,old.set.oq])
  nlmPar.old.oq <- prev_par_oq
  nlmPar.new.oq <- c()

  if(direction == "forward"){
    item.no.oq <- which(!(old.set.oq %in% new.set.oq))
    idx.remv.oq.beta <- c((((nn.th*item.no.oq)-(nn.th-1))):((nn.th*item.no.oq)))
    idx.remv.oq.gamma <- c((ncol(old.dataset.oq)*nn.th)+item.no.oq)
    if(type[1] == "ipoqlldif"){
      idx.remv.oq.deltabeta <- c()
      for(i in 1: n.resp.th){
        idx.remv.oq.deltabeta <- c(idx.remv.oq.deltabeta,(((ncol(old.dataset.oq)*nn.th)+ncol(old.dataset.oq))+(ncol(old.dataset.oq)*(i-1))+item.no.oq))
      }
      idx.remv.oq <- c(idx.remv.oq.beta,idx.remv.oq.gamma,idx.remv.oq.deltabeta)
    } else {
      idx.remv.oq <- c(idx.remv.oq.beta,idx.remv.oq.gamma)
    }
    nlmPar.new.oq <- nlmPar.old.oq[-c(idx.remv.oq)]
    nlmPar.new.oq <- nlmPar.new.oq + (runif(length(nlmPar.new.oq),-1,1)*oq_noise)
  } else if(direction == "backward"){
    item.no.oq <- which(!(new.set.oq %in% old.set.oq))
    idx.remv.oq.beta <- c((((nn.th*item.no.oq)-(nn.th-1))):((nn.th*item.no.oq)))
    idx.remv.oq.gamma <- c(length(idx.remv.oq.beta)+(ncol(old.dataset.oq)*nn.th)+item.no.oq)
    if(type[1] == "ipoqlldif"){
      idx.remv.oq.deltabeta <- c()
      for(i in 1: n.resp.th){
        idx.remv.oq.deltabeta <- c(idx.remv.oq.deltabeta,(((ncol(old.dataset.oq)*nn.th)+ncol(old.dataset.oq))+(ncol(old.dataset.oq)*(i-1))+item.no.oq+length(idx.remv.oq.beta)+length(idx.remv.oq.gamma)+(i-1)))
      }
      idx.remv.oq <- c(idx.remv.oq.beta,idx.remv.oq.gamma,idx.remv.oq.deltabeta)
    } else {
      idx.remv.oq <- c(idx.remv.oq.beta,idx.remv.oq.gamma)
    }
    nlmPar.new.oq <- insert.at(nlmPar.old.oq,idx.remv.oq,(length(nlmPar.old.oq)+length(idx.remv.oq)))


    item.no.iq <- which(!(old.set.iq %in% new.set.iq))
    idx.remv.iq.beta <- c((nrow(old.dataset.iq) + ((nn.th*item.no.iq)-(nn.th-1))):(nrow(old.dataset.iq)+(nn.th*item.no.iq)))
    idx.remv.iq.gamma <- c(nrow(old.dataset.iq)+(ncol(old.dataset.iq)*nn.th)+item.no.iq)
    if(type[1] == "ipoqlldif"){
      idx.remv.iq.deltabeta <- c()
      for(i in 1: n.resp.th){
        idx.remv.iq.deltabeta <- c(idx.remv.iq.deltabeta,(nrow(old.dataset.iq)+((ncol(old.dataset.iq)*nn.th)+ncol(old.dataset.iq))+(ncol(old.dataset.iq)*(i-1))+item.no.iq))
      }
      idx.remv.iq <- c(idx.remv.iq.beta,idx.remv.iq.gamma,idx.remv.iq.deltabeta)
    } else {
      idx.remv.iq <- c(idx.remv.iq.beta,idx.remv.iq.gamma)
    }
    nlmPar.new.oq[idx.remv.oq] <- nlmPar.old.iq[c(idx.remv.iq)]
    nlmPar.new.oq <- nlmPar.new.oq + (runif(length(nlmPar.new.oq),-1,1)*oq_noise)
  }

  return(nlmPar.new.oq)
}

insert.at <- function(a, pos, max.nlmpar){
  addLast <- FALSE
  pos <- (c(pos)-c(1:length(pos)))

  if(pos[length(pos)] == (max.nlmpar-length(pos))){
    pos <- pos[-c(length(pos))]
    a <- c(a,0)
    addLast <- TRUE
  }
  if(!identical(pos, integer(0))){
    if(pos[1] == 0){
      length.begin <- length(which(pos == 0))
      pos <- pos[-c(1:length.begin)]
      pos <- pos + length.begin
      a <- c(rep.int(0,length.begin),a)
    }
  }
  if(!identical(pos, integer(0)) & !identical(pos, numeric(0))){
    pos.idx <- split(pos,pos)
    dots <- rapply(pos.idx,function(x) x*0, how = "replace")
    pos <- unique(pos)
    stopifnot(length(dots)==length(pos))
    result <- vector("list",2*length(pos)+1)
    result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
    result[c(FALSE,TRUE)] <- dots

    result <- unlist(result)
  } else {
    result <- a
  }
  return(result)
}

#' @param object The object from the class \code{score}. The result of the score computation.
#' @param ... further argument passed or from other method.
#'
#' @rdname compute_score
#' @export
summary.score <- function(object, ...){
  dotdotdot <- list(...)
  cat("\n")
  cat("Score of the itemsets: ")
  cat("\n\n")
  cat("IQ-LL: ", object[1])
  cat("\nOQ-LL: ", object[2])
  cat("\nIPOQ-LL: ", object[3])
  cat("\n\n")
}
