#' Compute the In-plus-out-of-questionnaire log likelihood (with DIF) (IPOQ-LL(-DIF))
#'
#' \code{compute_score} computes the the IPOQ-LL/IPOQ-LL-DIF score of an instrument (included set) of the given initial survey.
#' While \code{compute_scores} computes the IPOQ-LL/IPOQ-LL-DIF score of many (more than one) instruments (included sets) of
#' the given initial survey simultanously.
#'
#' @param X A matrix or data.frame of the observed responses (ordinal or binary response).
#' @param incl_set A vector of the items (columns) number in the data.frame X that are included in the included set.
#' @param type The type of the score. \code{ipoqll} if we ignore the presence of the DIF and \code{ipoqlldif} if we want to consider the DIF effect.
#' @param groups_map Matrix to map the respondents to the DIF groups.
#' @param init_par_iq Initial values of the parameters in the included set before the estimation begin.
#' @param init_par_oq Initial values of the parameters in the excluded set before the estimation begin.
#' @param optim_control_iq The optimisation setting of the included set. See \code{\link[stats:optim]{stats::optim()}} \code{control} parameter.
#' @param optim_control_oq The optimisation setting of the excluded set. See \code{\link[stats:optim]{stats::optim()}} \code{control} parameter.
#' @param setting_par_iq The coordinate descent optimisation setting of the included set. See \code{\link[autoRasch:autoRaschOptions]{autoRasch::autoRaschOptions()}} \code{cd_control} parameter.
#' @param setting_par_oq The coordinate descent optimisation setting of the excluded set. See \code{\link[autoRasch:autoRaschOptions]{autoRasch::autoRaschOptions()}} \code{cd_control} parameter.
#'
#' @return
#' \code{compute_score} will return a vector which contains in-questionnaire log likelihood (IQ-LL(-DIF)), out-of-questionnaire log likelihood(OQ-LL(-DIF)),
#' IPOQ-LL(-DIF), included set's items' number in the given initial survey, the estimated theta parameters, the estimated items' parameters in the included set,
#' and the estimated items' parameters in the excluded set, sequentially.
#'
#' @examples
#' ipoqll_score <- compute_score(short_poly_data,incl_set = c(1:3),type = "ipoqll")
#'
#' \dontrun{
#' ipoqll_scores <- compute_scores(short_poly_data,incl_set = rbind(c(1:3),c(4:6)),
#'                                 type = "ipoqll", cores = 2)
#' }
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#'
#' @rdname compute_score
#' @export
compute_score <- function(X, incl_set, type = c("ipoqll","ipoqlldif"), groups_map = c(),
                          init_par_iq = c(), init_par_oq = c(),
                          optim_control_iq = c(), optim_control_oq = c(),
                          setting_par_iq = c(), setting_par_oq = c(), method = c("fast","novel")){

  if(is.null(type)){
    type <- "ipoqll"
  }

  if(type[1] == "ipoqll"){
    fixed_par <- c("delta")
    isPenalized_delta <- FALSE
    groups_map <- NULL
    scoreName <- "IPOQ-LL"
  } else if(type[1] == "ipoqlldif"){
    fixed_par <- c()
    isPenalized_delta <- TRUE
    if(is.null(groups_map)){
      stop("autoRasch ERROR: to use the `ipoqlldif`, `groups_map` must be provided.")
    }
    groups_map <- as.matrix(groups_map)
    scoreName <- "IPOQ-LL-DIF"
  }

  dset <- as.matrix(X)
  incl_set <- incl_set[!is.na(incl_set)]
  incl_resp <- dset[,incl_set]

  if(length(incl_set) != ncol(dset)){
    excl_resp <- dset[,-c(incl_set)]
  }

  # if(!is.null(init_par_iq) ){
  #
  #   if((minCat <- min(dset,na.rm = TRUE)) != 0){  ### makes sure the response is started at 0
  #     dset <- dset - minCat
  #   }
  #
  #   mt_vek_ori <- apply(dset, 2L, max, na.rm = TRUE)
  #   mt_vek_ori <- rep(max(mt_vek_ori),length(mt_vek_ori))
  #
  #   thetaidx_iq <- c((1):(nrow(dset)))
  #   betaidx_iq <- c((nrow(dset)+1):(nrow(dset)+(max(mt_vek_ori)*length(incl_set))))
  #   gammaidx_iq <- c((nrow(dset)+(max(mt_vek_ori)*length(incl_set))+1):(nrow(dset)+(max(mt_vek_ori)*length(incl_set))+ncol(incl_resp)))
  #
  #   mt_vek_incl <- apply(matrix(dset[,c(incl_set)],ncol = length(c(1:ncol(dset))[c(incl_set)])), 2L, max, na.rm = TRUE)
  #   mt_vek_incl <- rep(max(mt_vek_incl),length(mt_vek_incl))
  #
  #   mt_idx_ori_incl <- rep(c(1:length(mt_vek_incl)),each = max(mt_vek_ori))
  #
  #
  #   betalist_incl <- as.vector(unlist(tapply(init_par_iq[betaidx_iq], mt_idx_ori_incl, function(x){
  #     if(max(mt_vek_ori) > max(mt_vek_incl)){
  #       temp <- x[-c(max(mt_idx_ori_incl))]
  #     } else {
  #       temp <- x
  #     }
  #     return(temp)
  #   })))
  #
  #   if(type[1] == "ipoqll"){
  #     init_par_iq <- c(init_par_iq[thetaidx_iq], betalist_incl, init_par_iq[gammaidx_iq])
  #   } else if(type[1] == "ipoqlldif"){
  #     deltaidx_iq <- c((nrow(dset)+(max(mt_vek_ori)*length(incl_set))+ncol(incl_resp)+1):length(init_par_iq))
  #     init_par_iq <- c(init_par_iq[thetaidx_iq], betalist_incl, init_par_iq[gammaidx_iq], init_par_iq[deltaidx_iq])
  #   }
  #
  # }
  # #####
  # if(!is.null(init_par_oq)){
  #
  #   if((minCat <- min(dset,na.rm = TRUE)) != 0){  ### makes sure the response is started at 0
  #     dset <- dset - minCat
  #   }
  #
  #   mt_vek_ori <- apply(dset, 2L, max, na.rm = TRUE)
  #   mt_vek_ori <- rep(max(mt_vek_ori),length(mt_vek_ori))
  #
  #   betaidx_oq <- c((1):(max(mt_vek_ori)*ncol(excl_resp)))
  #   gammaidx_oq <- c(((max(mt_vek_ori)*ncol(excl_resp))+1):((max(mt_vek_ori)*ncol(excl_resp))+ncol(excl_resp)))
  #
  #   mt_vek_excl <- apply(matrix(dset[,-c(incl_set)],ncol = length(c(1:ncol(dset))[-c(incl_set)])), 2L, max, na.rm = TRUE)
  #   mt_vek_excl <- rep(max(mt_vek_excl),length(mt_vek_excl))
  #
  #   mt_idx_ori_excl <- rep(c(1:length(mt_vek_excl)),each = max(mt_vek_ori))
  #
  #   betalist_excl <- as.vector(unlist(tapply(init_par_oq[betaidx_oq], mt_idx_ori_excl, function(x){
  #     if(max(mt_vek_ori) > max(mt_vek_excl)){
  #       temp <- x[-c(max(mt_vek_ori))]
  #     } else {
  #       temp <- x
  #     }
  #     return(temp)
  #   })))
  #
  #   if(type[1] == "ipoqll"){
  #     init_par_oq <- c(betalist_excl, init_par_oq[gammaidx_oq])
  #   } else if(type[1] == "ipoqlldif"){
  #     deltaidx_oq <- c(((max(mt_vek_ori)*ncol(excl_resp))+ncol(excl_resp)+1):length(init_par_oq))
  #     init_par_oq <- c(betalist_excl, init_par_oq[gammaidx_oq], init_par_oq[deltaidx_oq])
  #   }
  #
  # }

  if(is.null(setting_par_iq)){
    setting_par_iq <- autoRaschOptions()
  } else {
    if("aR_opt" %in% class(setting_par_iq)){
    } else {
      stop("The setting used should be a class of aR_opt!")
    }
  }

  if(type[1] == "ipoqlldif"){
    setting_par_iq$optz_method <- "mixed"
  } else {
    setting_par_iq$optz_method <- "optim"
  }
  setting_par_iq$isHessian <- FALSE
  setting_par_iq$fixed_par <- fixed_par
  setting_par_iq$isPenalized_delta <- isPenalized_delta
  setting_par_iq$groups_map <- groups_map
  setting_par_iq$randomized <- TRUE

  iqll <- pjmle(incl_resp, init_par = init_par_iq, setting = setting_par_iq, method = method)

  if(ncol(dset) == length(incl_set)){
    loglik_oqll <- NA
  } else {
    if(is.null(setting_par_oq)){
      setting_par_oq <- autoRaschOptions()
    } else {
      if("aR_opt" %in% class(setting_par_oq)){
      } else {
        stop("The setting used should be a class of aR_opt!")
      }
    }

    if(type[1] == "ipoqlldif"){
      setting_par_oq$optz_method <- "mixed"
    } else {
      setting_par_oq$optz_method <- "optim"
    }
    setting_par_oq$isHessian <- FALSE
    setting_par_oq$fixed_par <- c("theta",fixed_par)
    setting_par_oq$fixed_theta <- iqll$theta
    setting_par_oq$isPenalized_delta <- isPenalized_delta
    setting_par_oq$isPenalized_theta <- FALSE
    setting_par_oq$groups_map <- groups_map
    setting_par_oq$randomized <- TRUE

    oqll <- pjmle(excl_resp, init_par = init_par_oq, setting = setting_par_oq, method = method)

    loglik_oqll <- oqll$loglik
  }


  ipoqll <- sum(c(iqll$loglik,loglik_oqll),na.rm = TRUE)
  res <- c(iqll$loglik, loglik_oqll, ipoqll)

  ### Parse the parameters estimate
  if(type[1] == "ipoqlldif"){

    ### parse for iq-ll estimate
    n_par <- sum(nrow(dset),((1+ncol(groups_map)+(max(dset,na.rm = TRUE)-min(dset,na.rm = TRUE)))*ncol(dset)))

    if((minCat <- min(dset,na.rm = TRUE)) != 0){  ### makes sure the response is started at 0
      dset <- dset - minCat
    }

    if(method[1] == "novel"){

      mt_vek_ori <- apply(dset, 2L, max, na.rm = TRUE)
      mt_vek_ori <- max(mt_vek_ori)

      mt_vek_incl <- apply(matrix(dset[,c(incl_set)],ncol = length(c(1:ncol(dset))[c(incl_set)])), 2L, max, na.rm = TRUE)

      mt_idx_incl <- rep(c(1:length(mt_vek_incl)),mt_vek_incl)
      betalist_incl <- as.vector(unlist(tapply(iqll$beta, mt_idx_incl, function(x){
        temp <- c(x,rep(0,(mt_vek_ori-length(x))))
        return(temp)
      })))

      betalength <- sum(apply(dset,2,function(x){
        temp <- max(x,na.rm = TRUE)-min(x,na.rm = TRUE)
        return(temp)
      }))


    } else {

      betalength <- ncol(dset)*max(iqll$mt_vek)
      betalist_incl <- iqll$beta

    }

    length(betalist_incl) <- betalength
    gamma.ret <- iqll$gamma
    length(gamma.ret) <- ncol(dset)
    delta.ret <- iqll$delta
    length(delta.ret) <- ncol(groups_map)*ncol(dset)

    iqll_params <- c(iqll$theta, betalist_incl, gamma.ret, delta.ret)

    ### parse for oq-ll estimate
    if(ncol(dset) == length(incl_set)){
      oqll_params <- NA
    } else {
      if(method[1] == "novel"){

        mt_vek_excl <- apply(matrix(dset[,-c(incl_set)],ncol = length(c(1:ncol(dset))[-c(incl_set)])), 2L, max, na.rm = TRUE)

        mt_idx_excl <- rep(c(1:length(mt_vek_excl)),mt_vek_excl)

        betalist_excl <- as.vector(unlist(tapply(oqll$beta, mt_idx_excl, function(x){
          temp <- c(x,rep(0,(mt_vek_ori-length(x))))
          return(temp)
        })))

      } else {

        betalist_excl <- oqll$beta

      }

      length(betalist_excl) <- betalength
      gamma.ret <- oqll$gamma
      length(gamma.ret) <- ncol(dset)
      delta.ret <- oqll$delta
      length(delta.ret) <- ncol(groups_map)*ncol(dset)

      oqll_params <- c(betalist_excl, gamma.ret, delta.ret)

    }
    length(iqll_params) <- n_par
    length(oqll_params) <- n_par - nrow(dset)
  } else {
    n_par <- sum(nrow(dset),((1+(max(dset,na.rm = TRUE)-min(dset,na.rm = TRUE)))*ncol(dset)))

    if((minCat <- min(dset,na.rm = TRUE)) != 0){  ### makes sure the response is started at 0
      dset <- dset - minCat
    }

    if(method[1] == "novel"){

      mt_vek_ori <- apply(dset, 2L, max, na.rm = TRUE)
      mt_vek_ori <- max(mt_vek_ori)

      mt_vek_incl <- apply(matrix(dset[,c(incl_set)],ncol = length(c(1:ncol(dset))[c(incl_set)])), 2L, max, na.rm = TRUE)

      mt_idx_incl <- rep(c(1:length(mt_vek_incl)),mt_vek_incl)

      betalist_incl <- as.vector(unlist(tapply(iqll$beta, mt_idx_incl, function(x){
        temp <- c(x,rep(0,(mt_vek_ori-length(x))))
        return(temp)
      })))

      betalength <- sum(apply(dset,2,function(x){
        temp <- max(x,na.rm = TRUE)-min(x,na.rm = TRUE)
        return(temp)
      }))

    } else {

      betalength <- ncol(dset)*max(iqll$mt_vek)
      betalist_incl <- iqll$beta

    }

    length(betalist_incl) <- betalength
    gamma.ret <- iqll$gamma
    length(gamma.ret) <- ncol(dset)

    iqll_params <- c(iqll$theta, betalist_incl, gamma.ret)


    if(ncol(dset) == length(incl_set)){
      oqll_params <- NA
    } else {

      if(method[1] == "novel"){

        mt_vek_excl <- apply(matrix(dset[,-c(incl_set)],ncol = length(c(1:ncol(dset))[-c(incl_set)])), 2L, max, na.rm = TRUE)

        mt_idx_excl <- rep(c(1:length(mt_vek_excl)),mt_vek_excl)

        betalist_excl <- as.vector(unlist(tapply(oqll$beta, mt_idx_excl, function(x){
          temp <- c(x,rep(0,(mt_vek_ori-length(x))))
          return(temp)
        })))

      } else {

        betalist_excl <- oqll$beta

      }

      length(betalist_excl) <- betalength
      gamma.ret <- oqll$gamma
      length(gamma.ret) <- ncol(dset)

      oqll_params <- c(betalist_excl, gamma.ret)

    }
    length(iqll_params) <- n_par
    length(oqll_params) <- n_par - nrow(dset)

  }

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
  names(res) <- c("IQ-LL","OQ-LL",scoreName,rep("item no.",length(incl_set)),rep("iq-ll par.",length(iqll_params)),rep("oq-ll par.",length(oqll_params)))
  class(res) <- c(class(res),"score",type[1])
  return(res)

}

#' @rdname compute_score
#' @export
# compute_score_fast <- function(X, incl_set, type = c("ipoqll","ipoqlldif"), groups_map = c(),
#                           init_par_iq = c(), init_par_oq = c(),
#                           optim_control_iq = c(), optim_control_oq = c(),
#                           setting_par_iq = c(), setting_par_oq = c()){
#
#   if(is.null(type)){
#     type <- "ipoqll"
#   }
#
#   if(type[1] == "ipoqll"){
#     fixed_par <- c("delta")
#     isPenalized_delta <- FALSE
#     groups_map <- NULL
#     scoreName <- "IPOQ-LL"
#   } else if(type[1] == "ipoqlldif"){
#     fixed_par <- c()
#     isPenalized_delta <- TRUE
#     if(is.null(groups_map)){
#       stop("autoRasch ERROR: to use the `ipoqlldif`, `groups_map` must be provided.")
#     }
#     groups_map <- as.matrix(groups_map)
#     scoreName <- "IPOQ-LL-DIF"
#   }
#
#   dset <- as.matrix(X)
#   incl_set <- incl_set[!is.na(incl_set)]
#   incl_resp <- dset[,incl_set]
#
#   if(length(incl_set) != ncol(dset)){
#     excl_resp <- dset[,-c(incl_set)]
#   }
#
#   # if(!is.null(init_par_iq) ){
#   #
#   #   if((minCat <- min(dset,na.rm = TRUE)) != 0){  ### makes sure the response is started at 0
#   #     dset <- dset - minCat
#   #   }
#   #
#   #   mt_vek_ori <- apply(dset, 2L, max, na.rm = TRUE)
#   #   mt_vek_ori <- rep(max(mt_vek_ori),length(mt_vek_ori))
#   #
#   #   thetaidx_iq <- c((1):(nrow(dset)))
#   #   betaidx_iq <- c((nrow(dset)+1):(nrow(dset)+(max(mt_vek_ori)*length(incl_set))))
#   #   gammaidx_iq <- c((nrow(dset)+(max(mt_vek_ori)*length(incl_set))+1):(nrow(dset)+(max(mt_vek_ori)*length(incl_set))+ncol(incl_resp)))
#   #
#   #   mt_vek_incl <- apply(matrix(dset[,c(incl_set)],ncol = length(c(1:ncol(dset))[c(incl_set)])), 2L, max, na.rm = TRUE)
#   #   mt_vek_incl <- rep(max(mt_vek_incl),length(mt_vek_incl))
#   #
#   #   mt_idx_ori_incl <- rep(c(1:length(mt_vek_incl)),each = max(mt_vek_ori))
#   #
#   #
#   #   betalist_incl <- as.vector(unlist(tapply(init_par_iq[betaidx_iq], mt_idx_ori_incl, function(x){
#   #     if(max(mt_vek_ori) > max(mt_vek_incl)){
#   #       temp <- x[-c(max(mt_idx_ori_incl))]
#   #     } else {
#   #       temp <- x
#   #     }
#   #     return(temp)
#   #   })))
#   #
#   #   if(type[1] == "ipoqll"){
#   #     init_par_iq <- c(init_par_iq[thetaidx_iq], betalist_incl, init_par_iq[gammaidx_iq])
#   #   } else if(type[1] == "ipoqlldif"){
#   #     deltaidx_iq <- c((nrow(dset)+(max(mt_vek_ori)*length(incl_set))+ncol(incl_resp)+1):length(init_par_iq))
#   #     init_par_iq <- c(init_par_iq[thetaidx_iq], betalist_incl, init_par_iq[gammaidx_iq], init_par_iq[deltaidx_iq])
#   #   }
#   #
#   # }
#   #   ######
#   # if(!is.null(init_par_oq)){
#   #
#   #   if((minCat <- min(dset,na.rm = TRUE)) != 0){  ### makes sure the response is started at 0
#   #     dset <- dset - minCat
#   #   }
#   #
#   #   mt_vek_ori <- apply(dset, 2L, max, na.rm = TRUE)
#   #   mt_vek_ori <- rep(max(mt_vek_ori),length(mt_vek_ori))
#   #
#   #   betaidx_oq <- c((1):(max(mt_vek_ori)*ncol(excl_resp)))
#   #   gammaidx_oq <- c(((max(mt_vek_ori)*ncol(excl_resp))+1):((max(mt_vek_ori)*ncol(excl_resp))+ncol(excl_resp)))
#   #
#   #   mt_vek_excl <- apply(matrix(dset[,-c(incl_set)],ncol = length(c(1:ncol(dset))[-c(incl_set)])), 2L, max, na.rm = TRUE)
#   #   mt_vek_excl <- rep(max(mt_vek_excl),length(mt_vek_excl))
#   #
#   #   mt_idx_ori_excl <- rep(c(1:length(mt_vek_excl)),each = max(mt_vek_ori))
#   #
#   #   betalist_excl <- as.vector(unlist(tapply(init_par_oq[betaidx_oq], mt_idx_ori_excl, function(x){
#   #     if(max(mt_vek_ori) > max(mt_vek_excl)){
#   #       temp <- x[-c(max(mt_vek_ori))]
#   #     } else {
#   #       temp <- x
#   #     }
#   #     return(temp)
#   #   })))
#   #
#   #   if(type[1] == "ipoqll"){
#   #     init_par_oq <- c(betalist_excl, init_par_oq[gammaidx_oq])
#   #   } else if(type[1] == "ipoqlldif"){
#   #     deltaidx_oq <- c(((max(mt_vek_ori)*ncol(excl_resp))+ncol(excl_resp)+1):length(init_par_oq))
#   #     init_par_oq <- c(betalist_excl, init_par_oq[gammaidx_oq], init_par_oq[deltaidx_oq])
#   #   }
#   #
#   # }
#
#   ### setting for iq-ll estimation
#   if(is.null(setting_par_iq)){
#     setting_par_iq <- autoRaschOptions()
#   } else {
#     if("aR_opt" %in% class(setting_par_iq)){
#     } else {
#       stop("The setting used should be a class of aR_opt!")
#     }
#   }
#
#   if(type[1] == "ipoqlldif"){
#     setting_par_iq$optz_method <- "mixed"
#   } else {
#     setting_par_iq$optz_method <- "optim"
#   }
#   setting_par_iq$isHessian <- FALSE
#   setting_par_iq$fixed_par <- fixed_par
#   setting_par_iq$isPenalized_delta <- isPenalized_delta
#   setting_par_iq$groups_map <- groups_map
#   setting_par_iq$randomized <- TRUE
#
#   ### compute the iq-ll/-dif
#   iqll <- pjmle(incl_resp, init_par = init_par_iq, setting = setting_par_iq)
#
#   ### setting for oq-ll estimation
#   if(ncol(dset) == length(incl_set)){
#     loglik_oqll <- NA
#   } else {
#     if(is.null(setting_par_oq)){
#       setting_par_oq <- autoRaschOptions()
#     } else {
#       if("aR_opt" %in% class(setting_par_oq)){
#       } else {
#         stop("The setting used should be a class of aR_opt!")
#       }
#     }
#
#     if(type[1] == "ipoqlldif"){
#       setting_par_oq$optz_method <- "mixed"
#     } else {
#       setting_par_oq$optz_method <- "optim"
#     }
#     setting_par_oq$isHessian <- FALSE
#     setting_par_oq$fixed_par <- c("theta",fixed_par)
#     setting_par_oq$fixed_theta <- iqll$theta
#     setting_par_oq$isPenalized_delta <- isPenalized_delta
#     setting_par_oq$isPenalized_theta <- FALSE
#     setting_par_oq$groups_map <- groups_map
#     setting_par_oq$randomized <- TRUE
#
#     ### compute the oq-ll/-dif
#     oqll <- pjmle(excl_resp, init_par = init_par_oq, setting = setting_par_oq)
#
#     loglik_oqll <- oqll$loglik
#   }
#
#
#   ### compute the ipoq-ll/-dif
#   ipoqll <- sum(c(iqll$loglik,loglik_oqll),na.rm = TRUE)
#   res <- c(iqll$loglik, loglik_oqll, ipoqll)
#
#
#   if(type[1] == "ipoqlldif"){
#     n_par <- sum(nrow(dset),((1+ncol(groups_map)+(max(dset,na.rm = TRUE)-min(dset,na.rm = TRUE)))*ncol(dset)))
#
#     if((minCat <- min(dset,na.rm = TRUE)) != 0){  ### makes sure the response is started at 0
#       dset <- dset - minCat
#     }
#
#     # mt_vek_ori <- apply(dset, 2L, max, na.rm = TRUE)
#     # mt_vek_ori <- max(mt_vek_ori)
#     #
#     # print(mt_vek_ori)
#
#     # mt_vek_incl <- apply(matrix(dset[,c(incl_set)],ncol = length(c(1:ncol(dset))[c(incl_set)])), 2L, max, na.rm = TRUE)
#     # mt_vek_incl <- rep(mt_vek_ori, length(incl_set))
#     # mt_vek_incl <- iqll$mt_vek
#
#     # print(mt_vek_incl)
#
#     # mt_idx_incl <- rep(c(1:length(mt_vek_incl)),mt_vek_incl)
#     # betalist_incl <- as.vector(unlist(tapply(iqll$beta, mt_idx_incl, function(x){
#     #   temp <- c(x,rep(0,(mt_vek_ori-length(x))))
#     #   return(temp)
#     # })))
#     #
#     # betalength <- sum(apply(dset,2,function(x){
#     #   temp <- max(x,na.rm = TRUE)-min(x,na.rm = TRUE)
#     #   return(temp)
#     # }))
#     # betalength <- sum(mt_vek_incl)
#     # length(betalist_incl) <- betalength
#     betalength <- ncol(dset)*max(iqll$mt_vek)
#     beta.ret <- iqll$beta
#     length(beta.ret) <- betalength
#     gamma.ret <- iqll$gamma
#     length(gamma.ret) <- ncol(dset)
#     delta.ret <- iqll$delta
#     length(delta.ret) <- ncol(groups_map)*ncol(dset)
#
#     # iqll_params <- c(iqll$theta, betalist_incl, iqll$gamma, iqll$delta)
#     # iqll_params <- c(iqll$theta, betalist_incl, gamma.ret, delta.ret)
#     iqll_params <- c(iqll$theta, beta.ret, gamma.ret, delta.ret)
#
#
#     if(ncol(dset) == length(incl_set)){
#       oqll_params <- NA
#     } else {
#       # mt_vek_excl <- apply(matrix(dset[,-c(incl_set)],ncol = length(c(1:ncol(dset))[-c(incl_set)])), 2L, max, na.rm = TRUE)
#       # mt_vek_excl <- rep(max(dset[,-c(incl_set)],na.rm = TRUE), length(c(1:ncol(dset))[-c(incl_set)]))
#
#
#       # mt_idx_excl <- rep(c(1:length(mt_vek_excl)),mt_vek_excl)
#       #
#       # betalist_excl <- as.vector(unlist(tapply(oqll$beta, mt_idx_excl, function(x){
#       #   temp <- c(x,rep(0,(mt_vek_ori-length(x))))
#       #   return(temp)
#       # })))
#
#       # length(betalist_excl) <- sum(betalength)
#       beta.ret <- oqll$beta
#       length(beta.ret) <- betalength
#       gamma.ret <- oqll$gamma
#       length(gamma.ret) <- ncol(dset)
#       delta.ret <- oqll$delta
#       length(delta.ret) <- ncol(groups_map)*ncol(dset)
#
#       # oqll_params <- c(betalist_excl, oqll$gamma, oqll$delta)
#       # oqll_params <- c(betalist_excl, gamma.ret, delta.ret)
#       oqll_params <- c(beta.ret, gamma.ret, delta.ret)
#     }
#     length(iqll_params) <- n_par
#     length(oqll_params) <- n_par - nrow(dset)
#   } else {
#     n_par <- sum(nrow(dset),((1+(max(dset,na.rm = TRUE)-min(dset,na.rm = TRUE)))*ncol(dset)))
#
#     if((minCat <- min(dset,na.rm = TRUE)) != 0){  ### makes sure the response is started at 0
#       dset <- dset - minCat
#     }
#
#     # mt_vek_ori <- apply(dset, 2L, max, na.rm = TRUE)
#     # mt_vek_ori <- max(mt_vek_ori)
#     #
#     # # mt_vek_incl <- apply(matrix(dset[,c(incl_set)],ncol = length(c(1:ncol(dset))[c(incl_set)])), 2L, max, na.rm = TRUE)
#     # # mt_vek_incl <- apply(matrix(dset[,c(incl_set)],ncol = length(c(1:ncol(dset))[c(incl_set)])), 2L, max, na.rm = TRUE)
#     # mt_vek_incl <- rep(mt_vek_ori, length(incl_set))
#
#     # mt_idx_incl <- rep(c(1:length(mt_vek_incl)),mt_vek_incl)
#     #
#     # betalist_incl <- as.vector(unlist(tapply(iqll$beta, mt_idx_incl, function(x){
#     #   temp <- c(x,rep(0,(mt_vek_ori-length(x))))
#     #   return(temp)
#     # })))
#     #
#     # betalength <- sum(apply(dset,2,function(x){
#     #   temp <- max(x,na.rm = TRUE)-min(x,na.rm = TRUE)
#     #   return(temp)
#     # }))
#     # length(betalist_incl) <- betalength
#
#     betalength <- ncol(dset)*max(iqll$mt_vek)
#     beta.ret <- iqll$beta
#     length(beta.ret) <- betalength
#     gamma.ret <- iqll$gamma
#     length(gamma.ret) <- ncol(dset)
#
#     # iqll_params <- c(iqll$theta, betalist_incl, iqll$gamma)
#     # iqll_params <- c(iqll$theta, betalist_incl, gamma.ret)
#     iqll_params <- c(iqll$theta, beta.ret, gamma.ret)
#
#     if(ncol(dset) == length(incl_set)){
#       oqll_params <- NA
#     } else {
#       # mt_vek_excl <- apply(matrix(dset[,-c(incl_set)],ncol = length(c(1:ncol(dset))[-c(incl_set)])), 2L, max, na.rm = TRUE)
#       # mt_vek_excl <- apply(matrix(dset[,-c(incl_set)],ncol = length(c(1:ncol(dset))[-c(incl_set)])), 2L, max, na.rm = TRUE)
#       # mt_vek_excl <- rep(max(dset[,-c(incl_set)],na.rm = TRUE), length(c(1:ncol(dset))[-c(incl_set)]))
#       #
#       # mt_idx_excl <- rep(c(1:length(mt_vek_excl)),mt_vek_excl)
#       #
#       # betalist_excl <- as.vector(unlist(tapply(oqll$beta, mt_idx_excl, function(x){
#       #   temp <- c(x,rep(0,(mt_vek_ori-length(x))))
#       #   return(temp)
#       # })))
#
#       # length(betalist_excl) <- sum(betalength)
#       beta.ret <- oqll$beta
#       length(beta.ret) <- betalength
#       gamma.ret <- oqll$gamma
#       length(gamma.ret) <- ncol(dset)
#
#       # oqll_params <- c(betalist_excl, oqll$gamma)
#       # oqll_params <- c(betalist_excl, gamma.ret)
#       oqll_params <- c(beta.ret, gamma.ret)
#     }
#     length(iqll_params) <- n_par
#     length(oqll_params) <- n_par - nrow(dset)
#
#   }
#
#   n_par <- sum(nrow(dset),((1+(max(dset,na.rm = TRUE)-min(dset,na.rm = TRUE)))*ncol(dset)))
#   iqll_params <- c(iqll$theta, iqll$beta, iqll$gamma)
#   if(ncol(dset) == length(incl_set)){
#     oqll_params <- NA
#   } else {
#     oqll_params <- c(oqll$beta, oqll$gamma)
#   }
#   length(iqll_params) <- n_par
#   length(oqll_params) <- n_par - nrow(dset)
#
#
#   length(incl_set) <- ncol(dset)
#   res <- c(res, incl_set, iqll_params, oqll_params)
#   names(res) <- c("IQ-LL","OQ-LL",scoreName,rep("item no.",length(incl_set)),rep("iq-ll par.",length(iqll_params)),rep("oq-ll par.",length(oqll_params)))
#   class(res) <- c(class(res),"score",type[1])
#   return(res)
# }

compute_scores_unparalleled <- function(X, incl_sets, type = c("ipoqll","ipoqlldif"),
                                        step_direct = c("fixed","forward","backward"), groups_map = c(),
                                        init_par_iq = c(), init_par_oq = c(), optim_control_iq = c(), optim_control_oq = c(),
                                        setting_par_iq = c(), setting_par_oq = c(), method = c("fast","novel")){


  dset <- as.matrix(X)
  # incl_sets <- itemsets

  if(is.vector(incl_sets) & length(incl_sets) > ncol(dset)){
    stop("autoRasch ERROR: the number of items in the incl_set can not exceed the initial items.")
  }

  if(type[1] == "ipoqlldif"){
    if(is.null(groups_map)){
      stop("autoRasch ERROR: to use the `ipoqlldif`, `groups_map` must be provided.")
    }
  } else {
    type <- "ipoqll"
  }

  if(is.matrix(incl_sets) | is.null(step_direct)){
    step_direct <- "fixed"
  }

  if(step_direct[1] == "forward"){
    excl_set <- c(1:ncol(dset))[-c(incl_sets)]
    add_items <- t(combn(excl_set,1))
    rep_itemsets <- matrix(rep.int(incl_sets,length(add_items)), nrow = length(add_items), byrow = TRUE)
    incl_sets <- cbind(rep_itemsets,add_items)
  } else if(step_direct[1] == "backward"){
    incl_sets <- t(combn(incl_sets,(length(incl_sets)-1)))
  } else if(step_direct[1] == "fixed"){
    if((any(class(incl_sets) == "matrix")) & (dim(incl_sets)[2] == 1 | dim(incl_sets)[1] == 1)){
      incl_sets <- matrix(as.vector(incl_sets), ncol = 1)
    }
  }

  incl_sets <- as.matrix(incl_sets)

  i <- NULL

  scoreList <- foreach(i=1:nrow(incl_sets), .combine = rbind, .errorhandling = "stop") %dopar% {

    incl_set <- incl_sets[i,]
    incl_set <- incl_set[!is.na(incl_set)]
    incl_set <- sort(incl_set,decreasing = FALSE)

    # if(!is.null(init_par_iq) & !is.null(init_par_oq)){
    #
    #
    #   init_iq <- iqll_init(dset = dset, prev_incl_set = incl_sets, prev_par_iq = init_par_iq, prev_par_oq = init_par_oq,
    #                        incl_set = incl_set, direction = step_direct, type = type[1], groups_map = groups_map,
    #                        iq_noise = 1e-3)
    #
    #
    #
    #   init_oq <- oqll_init(dset = dset, prev_incl_set = incl_sets, prev_par_iq = init_par_iq, prev_par_oq = init_par_oq,
    #                        incl_set = incl_set, direction = step_direct, type = type[1], groups_map = groups_map,
    #                        oq_noise = 1e-3)
    # } else {
      init_iq <- c()
      init_oq <- c()
    # }

    # cat("incl_set : ", incl_set)
    # cat("\n length.init_iq : ",length(init_iq))
    # cat("\n length.init_oq : ",length(init_oq),"\n")

    # if(method[1] == "novel"){
        score_res <- compute_score(dset, incl_set = incl_set, type = type, groups_map = groups_map,
                              init_par_iq = init_iq, init_par_oq = init_oq,
                              optim_control_iq = optim_control_iq, optim_control_oq = optim_control_oq,
                              setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq,
                              method = method)
    # } else {
    #   score_res <- compute_score_fast(dset, incl_set = incl_set, type = type, groups_map = groups_map,
    #                              init_par_iq = init_iq, init_par_oq = init_oq,
    #                              optim_control_iq = optim_control_iq, optim_control_oq = optim_control_oq,
    #                              setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq)
    # }

    length(incl_set) <- ncol(dset)
    res <- c(score_res)
    return(res)
  }


  res <- scoreList

  return(res)
}

#' @param incl_sets A matrix as a results of a \code{rbind} of \code{incl_set}.
#' @param cores Number of cores that is used in the paralellization.
#' @param step_direct How will you compute the criterion score. \code{fixed} for the given itemset,
#' \code{forward} computes all the scores of the possible combination of items if an item is added to the given set,
#' \code{backward}  computes all the scores of the possible combination of items if an item is removed to the given set.
#'
#' @return
#' \code{compute_scores} will return a matrix as a result of the \code{rbind} operation of the \code{compute_score}'s result.
#'
#' @import doParallel
#' @import foreach
#'
#' @rdname compute_score
#' @export
compute_scores <- function(X, incl_sets, type = c("ipoqll","ipoqlldif"),
                           step_direct = c("fixed","forward","backward"), groups_map = c(),
                           init_par_iq = c(), init_par_oq = c(), optim_control_iq = c(), optim_control_oq = c(),
                           setting_par_iq = c(), setting_par_oq = c(),
                           cores = NULL, method = c("fast","novel")){


  incl_sets <- as.matrix(incl_sets)

  if(is.null(cores)){
    cores <- nrow(incl_sets)
    if(cores > 2){
      cores <- 2
    }
  } else {
    if(cores > detectCores()){
      cores <- detectCores()
    }
  }


  cl <- parallel::makeCluster(cores)
  # oFuture::registerDoFuture()
  # future::plan(future::cluster, workers = cl)
  doParallel::registerDoParallel(cl=cl, cores = cores)

  scoreList <- compute_scores_unparalleled(X = X, incl_sets = incl_sets, type = type,
                                           step_direct = step_direct, groups_map = groups_map,
                                          init_par_iq = init_par_iq, init_par_oq = init_par_oq,
                                          optim_control_iq = optim_control_iq, optim_control_oq = optim_control_oq,
                                          setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq, method = method)

  parallel::stopCluster(cl)
  foreach::registerDoSEQ()

  res <- scoreList

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
      idx.remv.iq.delta <- c()
      for(i in 1: n.resp.th){
        idx.remv.iq.delta <- c(idx.remv.iq.delta,(nrow(old.dataset.iq)+((ncol(old.dataset.iq)*nn.th)+ncol(old.dataset.iq))+(ncol(old.dataset.iq)*(i-1))+item.no.iq))
      }
      idx.remv.iq <- c(idx.remv.iq.beta,idx.remv.iq.gamma,idx.remv.iq.delta)
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
      idx.remv.iq.delta <- c()
      for(i in 1: n.resp.th){
        idx.remv.iq.delta <- c(idx.remv.iq.delta,(nrow(old.dataset.iq)+((ncol(old.dataset.iq)*nn.th)+ncol(old.dataset.iq))+(ncol(old.dataset.iq)*(i-1))+item.no.iq+length(idx.remv.iq.beta)+length(idx.remv.iq.gamma+(i-1))))
      }
      idx.remv.iq <- c(idx.remv.iq.beta,idx.remv.iq.gamma,idx.remv.iq.delta)
    } else {
      idx.remv.iq <- c(idx.remv.iq.beta,idx.remv.iq.gamma)
    }
    nlmPar.new.iq <- insert.at(nlmPar.old.iq,idx.remv.iq,(length(nlmPar.old.iq)+length(idx.remv.iq)))

    item.no.oq <- which(!(old.set.oq %in% new.set.oq))
    idx.remv.oq.beta <- c((((nn.th*item.no.oq)-(nn.th-1))):((nn.th*item.no.oq)))
    idx.remv.oq.gamma <- c((ncol(old.dataset.oq)*nn.th)+item.no.oq)
    if(type[1] == "ipoqlldif"){
      idx.remv.oq.delta <- c()
      for(i in 1: n.resp.th){
        idx.remv.oq.delta <- c(idx.remv.oq.delta,(((ncol(old.dataset.oq)*nn.th)+ncol(old.dataset.oq))+(ncol(old.dataset.oq)*(i-1))+item.no.oq))
      }
      idx.remv.oq <- c(idx.remv.oq.beta,idx.remv.oq.gamma,idx.remv.oq.delta)
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
      idx.remv.oq.delta <- c()
      for(i in 1: n.resp.th){
        idx.remv.oq.delta <- c(idx.remv.oq.delta,(((ncol(old.dataset.oq)*nn.th)+ncol(old.dataset.oq))+(ncol(old.dataset.oq)*(i-1))+item.no.oq))
      }
      idx.remv.oq <- c(idx.remv.oq.beta,idx.remv.oq.gamma,idx.remv.oq.delta)
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
      idx.remv.oq.delta <- c()
      for(i in 1: n.resp.th){
        idx.remv.oq.delta <- c(idx.remv.oq.delta,(((ncol(old.dataset.oq)*nn.th)+ncol(old.dataset.oq))+(ncol(old.dataset.oq)*(i-1))+item.no.oq+length(idx.remv.oq.beta)+length(idx.remv.oq.gamma)+(i-1)))
      }
      idx.remv.oq <- c(idx.remv.oq.beta,idx.remv.oq.gamma,idx.remv.oq.delta)
    } else {
      idx.remv.oq <- c(idx.remv.oq.beta,idx.remv.oq.gamma)
    }
    nlmPar.new.oq <- insert.at(nlmPar.old.oq,idx.remv.oq,(length(nlmPar.old.oq)+length(idx.remv.oq)))


    item.no.iq <- which(!(old.set.iq %in% new.set.iq))
    idx.remv.iq.beta <- c((nrow(old.dataset.iq) + ((nn.th*item.no.iq)-(nn.th-1))):(nrow(old.dataset.iq)+(nn.th*item.no.iq)))
    idx.remv.iq.gamma <- c(nrow(old.dataset.iq)+(ncol(old.dataset.iq)*nn.th)+item.no.iq)
    if(type[1] == "ipoqlldif"){
      idx.remv.iq.delta <- c()
      for(i in 1: n.resp.th){
        idx.remv.iq.delta <- c(idx.remv.iq.delta,(nrow(old.dataset.iq)+((ncol(old.dataset.iq)*nn.th)+ncol(old.dataset.iq))+(ncol(old.dataset.iq)*(i-1))+item.no.iq))
      }
      idx.remv.iq <- c(idx.remv.iq.beta,idx.remv.iq.gamma,idx.remv.iq.delta)
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
