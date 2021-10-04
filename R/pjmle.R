# Penalized Joint Maximum likelihood Estimation
#
# Estimate the theta, beta and gamma parameter of the full dataset
# X is the dataset which is want to be estimated
# init_par is a vector contains the initial values of the estimated parameters
# setting contains the parameter setting used  in the estimation. See autoRaschOptions().

pjmle <- function(X, init_par = c(), setting = c()){

  dset <- as.data.frame(X)                          ### makes sure that the dataset has a matrix format

  if(is.null(setting)){
    opts <- autoRaschOptions()
  } else {
    opts <- setting
  }

  ### preprocessing the data
  dataPrep <- data_prep(dset = dset, fixed_par = opts$fixed_par, groups_map = opts$groups_map)

  ### Intializing the parameters ###
  if(opts$randomized){
    theta <- runif(nrow(dataPrep$dset),-1,1)*opts$random.init.th
    beta <- runif(dataPrep$allcat,-1,1)*opts$random.init.th
    gamma <- runif(ncol(dataPrep$dset),-1,1)*opts$random.init.th
    delta <- rep(0,(ncol(dataPrep$dset)*ncol(dataPrep$groups_map)))
  } else {
    theta <- rep(0,nrow(dataPrep$dset))
    beta <- rep(0,dataPrep$allcat)
    gamma <- rep(0,ncol(dataPrep$dset)) #gpcm uses different gamma for each item
    delta <- rep(0,(ncol(dataPrep$dset)*ncol(dataPrep$groups_map)))
  }

  ### setting the optimized parameter
  nlmPar <- c(theta,beta,gamma,delta)
  length_theta <- length(theta)
  length_beta <- length(beta)
  length_gamma <- length(gamma)
  length_delta <- length(delta)
  length_array <- c(length_theta,length_beta,length_gamma,length_delta)
  fullPar_arr <- c("theta","beta","gamma","delta")

  if(is.null(opts$fixed_par)){
    estPar_arr <- fullPar_arr
  } else {
    estPar_arr <- fullPar_arr[-c(which(fullPar_arr %in% opts$fixed_par))]
  }
  fixLength_arr <- length_array[c(which(fullPar_arr %in% opts$fixed_par))]

  fixValue <- c()
  estLength_array <- length_array

  if(!identical(grep("delta",opts$fixed_par), integer(0))){
    nlmPar <- nlmPar[-c((sum(length_array[1:3])+1):(sum(length_array[1:4])))]
    estLength_array <- estLength_array[-c(4)]
    if(!is.null(opts$fixed_delta)){
      fixValue <- c(opts$fixed_delta,fixValue)
    } else {
      fixValue <- c(rep(0,(ncol(dset)*ncol(dataPrep$groups_map))),fixValue)
    }
  }
  if(!identical(grep("^gamma",opts$fixed_par), integer(0))){
    if(length(nlmPar) == length(length_array)){
      nlmPar <- nlmPar[-c((sum(length_array[1:2])+1):(sum(length_array[1:3])),(sum(length_array[1:4])+1):(sum(length_array[1:5])))]
      estLength_array <- estLength_array[-c(3,5)]
      if(!is.null(opts$fixed_gamma)){
        fixValue <- c(opts$fixed_gamma, rep(0,(ncol(dset)*ncol(dataPrep$groups_map))), fixValue)
      } else {
        fixValue <- c(rep(0,ncol(dset)), rep(0,(ncol(dset)*ncol(dataPrep$groups_map))), fixValue)
      }
    } else {
      nlmPar <- nlmPar[-c((sum(length_array[1:2])+1):(sum(length_array[1:3])))]
      estLength_array <- estLength_array[-c(3)]
      if(!is.null(opts$fixed_gamma)){
        fixValue <- c(opts$fixed_gamma, fixValue)
        # print(opts$fixed_gamma)
      } else {
        fixValue <- c(rep(0,ncol(dset)), fixValue)
      }
    }
  }
  if(!identical(grep("^beta",opts$fixed_par), integer(0))){
    nlmPar <- nlmPar[-c((sum(length_array[1])+1):(sum(length_array[1:2])))]
    estLength_array <- estLength_array[-c(2)]
    fixValue <- c(opts$fixed_beta, fixValue)
  }
  if(!identical(grep("theta",opts$fixed_par), integer(0))){
    nlmPar <- nlmPar[-c((1):(sum(length_array[1])))]
    estLength_array <- estLength_array[-c(1)]
    fixValue <- c(opts$fixed_theta, fixValue)
  }


  if(!is.null(init_par)) {
    nlmPar <- init_par
  }

  ### Maximizing the loglikelihood function ###
  nameCol <- colnames(as.data.frame(X))
  output <- list("X" = X, "mt_vek" = dataPrep$mt_vek, "real_vek" = dataPrep$na_catVec, "itemName" = nameCol,
                 "penalty.coeff" = list("lambda_theta" = opts$lambda_theta,"lambda_in" = opts$lambda_in,
                                        "lambda_out" = opts$lambda_out,"lambda_delta" = opts$lambda_delta))

  if(!is.null(opts$groups_map)){
    output[["groups_map"]] <- dataPrep$groups_map
  }

  if(opts$optz_method == "optim"){
    (minRes <- optim(nlmPar, loglik_fun, gr = grad_fun, hessian = opts$isHessian, dset = dataPrep$dset,
                     lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = opts$eps,
                     lambda_delta = opts$lambda_delta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                     fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                     groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_idx, fixed_par = opts$fixed_par, fixValue = fixValue,
                     isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                     isPenalized_delta = opts$isPenalized_delta, method = "BFGS", control = opts$optim_control))

    est <- minRes$par
    obj <- minRes$value
    iter <- minRes$iterations
    conv <- minRes$convergence
    counts <- minRes$counts
    output[["loglik"]] <- -obj
    if(opts$isHessian){
      output[["hessian"]] <- minRes$hessian
    }
  } else if(opts$optz_method == "mixed"){
    (minRes <- mixed.min(nlmPar, dset = dset, dataPrep = dataPrep, opts = opts, estPar_arr = estPar_arr,
                         estLength_array = estLength_array, eps = opts$eps,fixLength_arr = fixLength_arr,
                         length_array = length_array, fixValue = fixValue,
                         maxit.cd.lower = opts$cd_control$maxit.cd.lower, init.step = opts$cd_control$init.step,
                         abs.tol = opts$cd_control$abs.tol, maxit.cd.higher = opts$cd_control$maxit.cd.higher,
                         maxit.optim = opts$cd_control$maxit.optim, scale.down = opts$cd_control$scale.down,
                         checkNonZero = opts$cd_control$checkNonZero,
                         max.diff.par = opts$cd_control$max.diff.par))

    est <- minRes$par
    obj <- minRes$value
    iter <- minRes$iterations
    conv <- minRes$convergence
    counts <- minRes$counts
    output[["loglik"]] <- -obj
    if(opts$isHessian){
      output[["hessian"]] <- minRes$hessian
    }
  }
  print("...done!")

  if(!identical(grep("delta",estPar_arr), integer(0))){
    parNo <- grep("delta",estPar_arr)
    if(parNo == 1){
      output[["delta"]] <- est[c((1):(sum(estLength_array[1])))]
    } else {
      output[["delta"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  }
  if(!identical(grep("^gamma",estPar_arr), integer(0))){
    parNo <- grep("^gamma",estPar_arr)
    if(parNo == 1){
      output[["gamma"]] <- est[c((1):(sum(estLength_array[1])))]
    } else {
      output[["gamma"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  }
  if(!identical(grep("^beta",estPar_arr), integer(0))){
    parNo <- grep("^beta",estPar_arr)
    if(parNo == 1){
      output[["beta"]] <- est[c((1):(sum(estLength_array[1])))]#*dataPrep$na_catVec
    } else {
      output[["beta"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]#*dataPrep$na_catVec
    }
  }
  if(!identical(grep("theta",estPar_arr), integer(0))){
    parNo <- grep("theta",estPar_arr)
    if(parNo == 1){
      output[["theta"]] <- est[c((1):(sum(estLength_array[1])))]
    } else {
      output[["theta"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  }

  return(output)

}

data_prep <- function(dset, fixed_par, groups_map){

  if(!is.null(groups_map)){
    groups_map <- as.matrix(groups_map)
    if(nrow(groups_map) != nrow(dset)){
      stop("autoRasch ERROR: `groups_map` must be has the same number of rows with observation.")
    }
    if(max(groups_map, na.rm = TRUE) > 1 | min(groups_map, na.rm = TRUE) < 0 | length(levels(factor(unlist(groups_map)))) > 2){
      stop("autoRasch ERROR: `groups_map` must be binary matrix with only 1 and 0 values.")
    }
  } else {
    groups_map <- matrix(0,nrow = nrow(dset),ncol = 1)
  }

  ### Manages catagory's vector
  mt_vek <- apply(dset, 2L, max, na.rm = TRUE)      ### create vector of max categories for each item
  n_th <- max(mt_vek)                               ### number of thresholds
  mt_vek_min <- apply(dset, 2L, min, na.rm = TRUE)  ### create vector of min categories for each item
  true_catVec <- (mt_vek-n_th)*(-1)
  allcat <- sum(mt_vek)                             ### number of items * categories (assumption : item has the same number of categories)
  if(n_th != 1){
    na_catVec <- as.vector(apply(rbind(mt_vek_min,true_catVec),2,function(x){     ### To remove the unavailable categories if items contain different number of categories (used for output filter)
      temp <- c(rep(NA,x[1]),rep(1,(n_th-sum(x))),rep(NA,x[2]))
      return(temp)
    }))
  } else {
    na_catVec <- mt_vek
  }

  mt_idx <- rep(c(1:length(mt_vek)),mt_vek)
  dimResp <- dim(dset)

  ### Create the xn.mat (matrix of x_ni) and the xna.mat (matrix of missing value) to map the responses and the missing values

  ### XNA is infused to XN
  xn.vec <- xna.vec <- c()
  for(rIdx in 1:dimResp[1]){
    for(cIdx in 1:dimResp[2]){
      if(is.na(dset[rIdx,cIdx])){
        xn.cell <- xna.cell <- rep(NA,mt_vek[cIdx])
      } else {
        xn.cell <- c(rep(1,dset[rIdx,cIdx]),rep(0,(mt_vek[cIdx]-dset[rIdx,cIdx])))
        xna.cell <- rep(1,(mt_vek[cIdx]))
      }
      xn.vec <- c(xn.vec,xn.cell)
      xna.vec <- c(xna.vec,xna.cell)
    }
  }
  xn.mat <- matrix(xn.vec, nrow = dimResp[1], byrow = TRUE)
  xna.mat <- matrix(xna.vec, nrow = dimResp[1], byrow = TRUE)

  XN <- as.vector(t(xn.mat))
  XNA <- as.vector(t(xna.mat))

  ret <- list("dset" = dset, "mt_vek" = mt_vek, "mt_idx" = mt_idx, "dimResp" = dimResp, "groups_map" = groups_map, "n_th" = n_th, "na_catVec" = na_catVec, "allcat" = allcat, "XN" = XN, "XNA" = XNA)#,"XREAL" = XREAL)

  return(ret)

}

coord.descent <- function(nlmPar, dset, dataPrep, opts, fixed_par = c(), step.vec = NULL,
                          estPar_arr = c(), estLength_array = c(), eps = c(),fixLength_arr = c(), fixValue = c(),
                          maxit.cd.lower = NULL, init.step = NULL, scale.down = NULL, checkNonZero = NULL){

  nlmpar.new <- nlmPar


  if(checkNonZero){
    n.par.idx <- c(1:length(nlmPar))
  } else {
    n.par.idx <- which(abs(nlmPar) >= 1e-2)
    if(length(n.par.idx)==0){
      n.par.idx <- c(1:length(nlmPar))
    } else {
      nlmpar.new[-c(n.par.idx)] <- 0
    }
  }

  if(is.null(step.vec)){
    step.vec <- rep(init.step,length(nlmPar))
  }

  opts$isPenalized_delta <- TRUE

  delta.vector <- c()
  for(n.par in n.par.idx){
    stepsize <- step.vec[n.par]
    ll.val.old <- loglik_fun(nlmpar.new, dset = dataPrep$dset, lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = eps,
                             lambda_delta = opts$lambda_delta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                             fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                             groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_id, fixed_par = fixed_par, fixValue = fixValue,
                             isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                             isPenalized_delta = opts$isPenalized_delta)
    ll.val.baseline <- ll.val.old
    nlmpar.old <- nlmpar.newminus <- nlmpar.newplus <- nlmpar.new

    nlmpar.newminus[n.par] <- nlmpar.old[n.par] - stepsize
    nlmpar.newplus[n.par] <- nlmpar.old[n.par] - (-stepsize)

    ll.val.newminus <- loglik_fun(nlmpar.newminus, dset = dataPrep$dset, lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = eps,
                                  lambda_delta = opts$lambda_delta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                                  fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                                  groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_id, fixed_par = fixed_par, fixValue = fixValue,
                                  isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                                  isPenalized_delta = opts$isPenalized_delta)
    ll.val.newplus <- loglik_fun(nlmpar.newplus, dset = dataPrep$dset, lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = eps,
                                 lambda_delta = opts$lambda_delta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                                 fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                                 groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_id, fixed_par = fixed_par, fixValue = fixValue,
                                 isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                                 isPenalized_delta = opts$isPenalized_delta)

    if(ll.val.newminus < ll.val.old){
      stepsize <- stepsize
      ll.val.new <- ll.val.newminus
      nlmpar.new[n.par] <- nlmpar.newminus[n.par]
    } else if(ll.val.newplus < ll.val.old){
      stepsize <- -stepsize
      ll.val.new <- ll.val.newplus
      nlmpar.new[n.par] <- nlmpar.newplus[n.par]
    } else {
      stepsize <- 0
      step.vec[n.par] <- step.vec[n.par]*scale.down
      ll.val.new <- ll.val.old
    }



    i <- 0
    if(stepsize == 0){
    } else {
      i <- 1
      while((diff <- abs(abs(ll.val.new) - abs(ll.val.old))) > 1e-12 & ll.val.new < ll.val.old & (i < maxit.cd.lower)){
        nlmpar.old<- nlmpar.new
        ll.val.old <- ll.val.new
        nlmpar.new[n.par] <- nlmpar.old[n.par] - stepsize

        ll.val.new <- loglik_fun(nlmpar.new, dset = dataPrep$dset, lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = eps,
                                 lambda_delta = opts$lambda_delta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                                 fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                                 groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_id, fixed_par = fixed_par, fixValue = fixValue,
                                 isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                                 isPenalized_delta = opts$isPenalized_delta)

        diff <- abs(ll.val.new) - abs(ll.val.old)

        if(diff > 0){
          nlmpar.new[n.par] <- nlmpar.new[n.par] + stepsize
          stepsize <- 0
        }


        i <- i+1
      }

    }
  }

  return(list("value" = ll.val.new, "par" = nlmpar.new, iterations = i, step.vec = step.vec))

}

mixed.min <- function(nlmPar, dset, dataPrep, opts,
                      estPar_arr = NULL, estLength_array = NULL, eps = NULL,fixLength_arr = NULL, fixValue = c(),
                      maxit.cd.lower = NULL, maxit.cd.higher = NULL, maxit.optim = NULL, init.step = NULL, abs.tol = NULL, scale.down = NULL,
                      checkNonZero = TRUE, length_array = NULL, max.diff.par = NULL){

  ll.val.old <- 1e+10
  ll.value <- c()

  par.diff <- c()
  par.max <- c()

  par.arr <- c("theta", "beta","gamma","delta")
  par.fix.idx <-  which(par.arr %in% opts$fixed_par)

  fixed_par_gpcm <- unique(c(par.arr[par.fix.idx],"delta"))
  estPar_arr_gpcm <- par.arr[-c(which(par.arr %in% fixed_par_gpcm))]

  fixLength_arr_gpcm <- length_array[c(which(par.arr %in% fixed_par_gpcm),5)]
  estLength_array_gpcm <- length_array[which(par.arr %in% estPar_arr_gpcm)]

  fixed_par_gpcmdif <- c("theta","beta","gamma")
  estPar_arr_gpcmdif <- c("delta")
  fixLength_arr_gpcmdif <- length_array[c(1:3,5)]
  estLength_array_gpcmdif <- length_array[4]

  par.old <- c(nlmPar)
  fixed_delta <- nlmPar[(sum(estLength_array_gpcm)+1):(sum(estLength_array_gpcm)+estLength_array_gpcmdif)]

  nlmPar <- nlmPar[1:sum(estLength_array_gpcm)]

  fixValue <- c(fixed_delta)
  if(length(which("gamma" %in% fixed_par_gpcm)) != 0){
    fixValue <- c(rep(0,(dataPrep$dimResp[2])),fixValue)
  }
  if(!is.null(opts$fixed_beta)){
    fixValue <- c(opts$fixed_beta,fixValue)
  }
  if(!is.null(opts$fixed_theta)){
    fixValue <- c(opts$fixed_theta,fixValue)
  }

  (gpcm <- optim(nlmPar, fn = loglik_fun, gr = grad_fun, hessian = opts$isHessian, dset = dataPrep$dset,
                 fixValue = fixValue, fixed_par = fixed_par_gpcm, estPar_arr = estPar_arr_gpcm, fixLength_arr = fixLength_arr_gpcm, estLength_array = estLength_array_gpcm,
                 lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = 0, lambda_delta = opts$lambda_delta,
                 isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta, isPenalized_delta = FALSE,
                 allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                 groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_idx,
                 method = "BFGS", control = list(maxit = maxit.optim,reltol=1e-12,fnscale = 10)))

  if(length(which("gamma" %in% fixed_par_gpcm)) != 0){
    fixValue <- c(opts$fixed_theta,gpcm$par,rep(0,(dataPrep$dimResp[2])))
  } else {
    fixValue <- c(opts$fixed_theta,gpcm$par)
  }

  nlmPar <- fixed_delta

  (gpcmdif <- coord.descent(nlmPar, dset = dset, dataPrep = dataPrep, opts = opts, fixed_par = fixed_par_gpcmdif,
                            fixValue = fixValue, estPar_arr = estPar_arr_gpcmdif, estLength_array = estLength_array_gpcmdif, fixLength_arr = fixLength_arr_gpcmdif,
                            maxit.cd.lower = maxit.cd.lower, init.step = init.step, eps = eps,
                            scale.down = scale.down, checkNonZero = TRUE))


  ll.val.new <- gpcmdif$value
  par.new <- c(gpcm$par,gpcmdif$par)

  i <- 0

  while((diff <- (abs(abs(ll.val.new) - abs(ll.val.old)))) > abs.tol & i < maxit.cd.higher &
        ((diff.par <- (max(abs(abs(par.new) - abs(par.old))))) > max.diff.par)){
    ll.val.old <- ll.val.new
    par.old <- par.new


    nlmPar <- gpcm$par
    if(is.null(opts$fixed_theta)){
      if(length(which("gamma" %in% fixed_par_gpcm)) != 0){
        fixValue <- c(rep(0,(dataPrep$dimResp[2])),gpcmdif$par)
      } else {
        fixValue <- c(gpcmdif$par)
      }
    } else {
      if(length(which("gamma" %in% fixed_par_gpcm)) != 0){
        fixValue <- c(opts$fixed_theta,rep(0,(dataPrep$dimResp[2])),gpcmdif$par)
      } else {
        fixValue <- c(opts$fixed_theta,gpcmdif$par)
      }
    }
    (gpcm <- optim(nlmPar, loglik_fun, gr = grad_fun, hessian = opts$isHessian, dset = dataPrep$dset,
                   fixValue = fixValue, fixed_par = fixed_par_gpcm, estPar_arr = estPar_arr_gpcm, fixLength_arr = fixLength_arr_gpcm, estLength_array = estLength_array_gpcm,
                   lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = 0, lambda_delta = opts$lambda_delta,
                   isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta, isPenalized_delta = FALSE,
                   allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                   groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_idx,
                   method = "BFGS", control = list(maxit = 1e+4,reltol=1e-12,fnscale = 10)))

    if(length(which("gamma" %in% fixed_par_gpcm)) != 0){
      fixValue <- c(opts$fixed_theta,gpcm$par,rep(0,(dataPrep$dimResp[2])))
    } else {
      fixValue <- c(opts$fixed_theta,gpcm$par)
    }
    nlmPar <- gpcmdif$par

    if(i < 20){
      checkNonZero <- TRUE
    } else {
      checkNonZero <- FALSE
    }


    (gpcmdif <- coord.descent(nlmPar, dset = dset, dataPrep = dataPrep, opts = opts, fixed_par = fixed_par_gpcmdif,
                              fixValue = fixValue, estPar_arr = estPar_arr_gpcmdif, estLength_array = estLength_array_gpcmdif, fixLength_arr = fixLength_arr_gpcmdif,
                              maxit.cd.lower = maxit.cd.lower, init.step = init.step, step.vec = gpcmdif$step.vec, eps = eps,
                              scale.down = scale.down, checkNonZero = checkNonZero))



    ll.val.new <- gpcmdif$value
    par.new <- c(gpcm$par,gpcmdif$par)
    par.diff <- c(par.diff, mean(abs(par.new - par.old)))
    par.max <- c(par.max, max(abs(abs(par.new) - abs(par.old))))
    ll.value <- c(ll.value,ll.val.new)

  }

  pardif <- gpcmdif$par
  pardif_idx <- which(abs(pardif) < 1e-10)
  pardif[pardif_idx] <- 0
  par.new <- c(gpcm$par,pardif)

  if(ll.val.new > ll.val.old){
    ll.val.new <- ll.val.old
    par.new <- par.old
  }

  return(list("value" = ll.val.new, "par" = par.new, iterations = i))

}

pjmle_fast <- function(X, init_par = c(), setting = c()){

  dset <- as.data.frame(X)                          ### makes sure that the dataset has a matrix format

  if(is.null(setting)){
    opts <- autoRaschOptions()
  } else {
    opts <- setting
  }

  ### preprocessing the data
  dataPrep <- data_prep_fast(dset = dset, fixed_par = opts$fixed_par, groups_map = opts$groups_map)

  ### Intializing the parameters ###
  if(opts$randomized){
    theta <- runif(nrow(dataPrep$dset),-1,1)*opts$random.init.th
    beta <- runif(dataPrep$allcat,-1,1)*opts$random.init.th
    gamma <- runif(ncol(dataPrep$dset),-1,1)*opts$random.init.th
    delta <- rep(0,(ncol(dataPrep$dset)*ncol(dataPrep$groups_map)))
  } else {
    theta <- rep(0,nrow(dataPrep$dset))
    beta <- rep(0,dataPrep$allcat)
    gamma <- rep(0,ncol(dataPrep$dset)) #gpcm uses different gamma for each item
    delta <- rep(0,(ncol(dataPrep$dset)*ncol(dataPrep$groups_map)))
  }

  ### setting the optimized parameter
  nlmPar <- c(theta,beta,gamma,delta)
  length_theta <- length(theta)
  length_beta <- length(beta)
  length_gamma <- length(gamma)
  length_delta <- length(delta)
  length_array <- c(length_theta,length_beta,length_gamma,length_delta)
  fullPar_arr <- c("theta","beta","gamma","delta")

  if(is.null(opts$fixed_par)){
    estPar_arr <- fullPar_arr
  } else {
    estPar_arr <- fullPar_arr[-c(which(fullPar_arr %in% opts$fixed_par))]
  }
  fixLength_arr <- length_array[c(which(fullPar_arr %in% opts$fixed_par))]

  fixValue <- c()
  estLength_array <- length_array

  if(!identical(grep("delta",opts$fixed_par), integer(0))){
    nlmPar <- nlmPar[-c((sum(length_array[1:3])+1):(sum(length_array[1:4])))]
    estLength_array <- estLength_array[-c(4)]
    if(!is.null(opts$fixed_delta)){
      fixValue <- c(opts$fixed_delta,fixValue)
    } else {
      fixValue <- c(rep(0,(ncol(dset)*ncol(dataPrep$groups_map))),fixValue)
    }
  }
  if(!identical(grep("^gamma",opts$fixed_par), integer(0))){
    if(length(nlmPar) == length(length_array)){
      nlmPar <- nlmPar[-c((sum(length_array[1:2])+1):(sum(length_array[1:3])),(sum(length_array[1:4])+1):(sum(length_array[1:5])))]
      estLength_array <- estLength_array[-c(3,5)]
      if(!is.null(opts$fixed_gamma)){
        fixValue <- c(opts$fixed_gamma, rep(0,(ncol(dset)*ncol(dataPrep$groups_map))), fixValue)
      } else {
        fixValue <- c(rep(0,ncol(dset)), rep(0,(ncol(dset)*ncol(dataPrep$groups_map))), fixValue)
      }
    } else {
      nlmPar <- nlmPar[-c((sum(length_array[1:2])+1):(sum(length_array[1:3])))]
      estLength_array <- estLength_array[-c(3)]
      if(!is.null(opts$fixed_gamma)){
        fixValue <- c(opts$fixed_gamma, fixValue)
        # print(opts$fixed_gamma)
      } else {
        fixValue <- c(rep(0,ncol(dset)), fixValue)
      }
    }
  }
  if(!identical(grep("^beta",opts$fixed_par), integer(0))){
    nlmPar <- nlmPar[-c((sum(length_array[1])+1):(sum(length_array[1:2])))]
    estLength_array <- estLength_array[-c(2)]
    fixValue <- c(opts$fixed_beta, fixValue)
  }
  if(!identical(grep("theta",opts$fixed_par), integer(0))){
    nlmPar <- nlmPar[-c((1):(sum(length_array[1])))]
    estLength_array <- estLength_array[-c(1)]
    fixValue <- c(opts$fixed_theta, fixValue)
  }


  if(!is.null(init_par)) {
    nlmPar <- init_par
  }

  ### Maximizing the loglikelihood function ###
  nameCol <- colnames(as.data.frame(X))
  output <- list("X" = X, "mt_vek" = dataPrep$mt_vek, "real_vek" = dataPrep$na_catVec, "itemName" = nameCol,
                 "penalty.coeff" = list("lambda_theta" = opts$lambda_theta,"lambda_in" = opts$lambda_in,
                                        "lambda_out" = opts$lambda_out,"lambda_delta" = opts$lambda_delta))

  if(!is.null(opts$groups_map)){
    output[["groups_map"]] <- dataPrep$groups_map
  }

  if(opts$optz_method == "optim"){
    (minRes <- optim(nlmPar, loglik_fun_fast, gr = grad_fun_fast, hessian = opts$isHessian, dset = dataPrep$dset,
                                          lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = opts$eps,
                                          lambda_delta = opts$lambda_delta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                                          fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                                          groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_idx, fixed_par = opts$fixed_par, fixValue = fixValue,
                                          isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                                          isPenalized_delta = opts$isPenalized_delta, method = "BFGS", control = opts$optim_control))

    est <- minRes$par
    obj <- minRes$value
    iter <- minRes$iterations
    conv <- minRes$convergence
    counts <- minRes$counts
    output[["loglik"]] <- -obj
    if(opts$isHessian){
      output[["hessian"]] <- minRes$hessian
    }
  } else if(opts$optz_method == "mixed"){
    (minRes <- mixed.min.fast(nlmPar, dset = dset, dataPrep = dataPrep, opts = opts, estPar_arr = estPar_arr,
                         estLength_array = estLength_array, eps = opts$eps,fixLength_arr = fixLength_arr,
                         length_array = length_array, fixValue = fixValue,
                         maxit.cd.lower = opts$cd_control$maxit.cd.lower, init.step = opts$cd_control$init.step,
                         abs.tol = opts$cd_control$abs.tol, maxit.cd.higher = opts$cd_control$maxit.cd.higher,
                         maxit.optim = opts$cd_control$maxit.optim, scale.down = opts$cd_control$scale.down,
                         checkNonZero = opts$cd_control$checkNonZero,
                         max.diff.par = opts$cd_control$max.diff.par))

    est <- minRes$par
    obj <- minRes$value
    iter <- minRes$iterations
    conv <- minRes$convergence
    counts <- minRes$counts
    output[["loglik"]] <- -obj
    if(opts$isHessian){
      output[["hessian"]] <- minRes$hessian
    }
  }
  print("...done!")

  if(!identical(grep("delta",estPar_arr), integer(0))){
    parNo <- grep("delta",estPar_arr)
    if(parNo == 1){
      output[["delta"]] <- est[c((1):(sum(estLength_array[1])))]
    } else {
      output[["delta"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  }
  if(!identical(grep("^gamma",estPar_arr), integer(0))){
    parNo <- grep("^gamma",estPar_arr)
    if(parNo == 1){
      output[["gamma"]] <- est[c((1):(sum(estLength_array[1])))]
    } else {
      output[["gamma"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  }
  if(!identical(grep("^beta",estPar_arr), integer(0))){
    parNo <- grep("^beta",estPar_arr)
    if(parNo == 1){
      output[["beta"]] <- est[c((1):(sum(estLength_array[1])))]*dataPrep$na_catVec
    } else {
      output[["beta"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]*dataPrep$na_catVec
    }
  }
  if(!identical(grep("theta",estPar_arr), integer(0))){
    parNo <- grep("theta",estPar_arr)
    if(parNo == 1){
      output[["theta"]] <- est[c((1):(sum(estLength_array[1])))]
    } else {
      output[["theta"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  }

  return(output)

}

data_prep_fast <- function(dset, fixed_par, groups_map){

  if(!is.null(groups_map)){
      groups_map <- as.matrix(groups_map)
    if(nrow(groups_map) != nrow(dset)){
      stop("autoRasch ERROR: `groups_map` must be has the same number of rows with observation.")
    }
    if(max(groups_map, na.rm = TRUE) > 1 | min(groups_map, na.rm = TRUE) < 0 | length(levels(factor(unlist(groups_map)))) > 2){
      stop("autoRasch ERROR: `groups_map` must be binary matrix with only 1 and 0 values.")
    }
  } else {
    groups_map <- matrix(0,nrow = nrow(dset),ncol = 1)
  }

  ### Manages catagory's vector
  mt_vek <- apply(dset, 2L, max, na.rm = TRUE)      ### create vector of max categories for each item
  n_th <- max(mt_vek)                               ### number of thresholds
  mt_vek_min <- apply(dset, 2L, min, na.rm = TRUE)  ### create vector of min categories for each item
  true_catVec <- (mt_vek-n_th)*(-1)
  if(n_th != 1){
    na_catVec <- as.vector(apply(rbind(mt_vek_min,true_catVec),2,function(x){     ### To remove the unavailable categories if items contain different number of categories (used for output filter)
      temp <- c(rep(NA,x[1]),rep(1,(n_th-sum(x))),rep(NA,x[2]))
      return(temp)
    }))
  } else {
    na_catVec <- mt_vek
  }

  mt_vek <- rep(max(dset,na.rm = TRUE),ncol(dset))
  allcat <- sum(mt_vek)                             ### number of items * categories (assumption : item has the same number of categories)

  mt_idx <- rep(c(1:length(mt_vek)),mt_vek)
  dimResp <- dim(dset)

  ### Create the xn.mat (matrix of x_ni) and the xna.mat (matrix of missing value) to map the responses and the missing values

  ### XNA is infused to XN
  xn.vec <- xna.vec <- c()
  for(rIdx in 1:dimResp[1]){
    for(cIdx in 1:dimResp[2]){
      if(is.na(dset[rIdx,cIdx])){
        xn.cell <- xna.cell <- rep(NA,mt_vek[cIdx])
      } else {
        xn.cell <- c(rep(1,dset[rIdx,cIdx]),rep(0,(mt_vek[cIdx]-dset[rIdx,cIdx])))
        xna.cell <- rep(1,(mt_vek[cIdx]))
      }
      xn.vec <- c(xn.vec,xn.cell)
      xna.vec <- c(xna.vec,xna.cell)
    }
  }
  xn.mat <- matrix(xn.vec, nrow = dimResp[1], byrow = TRUE)
  xna.mat <- matrix(xna.vec, nrow = dimResp[1], byrow = TRUE)

  XN <- as.vector(t(xn.mat))
  XNA <- as.vector(t(xna.mat))

  ret <- list("dset" = dset, "mt_vek" = mt_vek, "mt_idx" = mt_idx, "dimResp" = dimResp, "groups_map" = groups_map, "n_th" = n_th, "na_catVec" = na_catVec, "allcat" = allcat, "XN" = XN, "XNA" = XNA)#,"XREAL" = XREAL)

  return(ret)

}

coord.descent.fast <- function(nlmPar, dset, dataPrep, opts, fixed_par = c(), step.vec = NULL,
                          estPar_arr = c(), estLength_array = c(), eps = c(),fixLength_arr = c(), fixValue = c(),
                          maxit.cd.lower = NULL, init.step = NULL, scale.down = NULL, checkNonZero = NULL){

  nlmpar.new <- nlmPar


  if(checkNonZero){
    n.par.idx <- c(1:length(nlmPar))
  } else {
    n.par.idx <- which(abs(nlmPar) >= 1e-2)
    if(length(n.par.idx)==0){
      n.par.idx <- c(1:length(nlmPar))
    } else {
      nlmpar.new[-c(n.par.idx)] <- 0
    }
  }

  if(is.null(step.vec)){
    step.vec <- rep(init.step,length(nlmPar))
  }

  opts$isPenalized_delta <- TRUE

  delta.vector <- c()
  for(n.par in n.par.idx){
    stepsize <- step.vec[n.par]
    ll.val.old <- loglik_fun_fast(nlmpar.new, dset = dataPrep$dset, lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = eps,
                             lambda_delta = opts$lambda_delta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                             fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                             groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_id, fixed_par = fixed_par, fixValue = fixValue,
                             isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                             isPenalized_delta = opts$isPenalized_delta)
    ll.val.baseline <- ll.val.old
    nlmpar.old <- nlmpar.newminus <- nlmpar.newplus <- nlmpar.new

    nlmpar.newminus[n.par] <- nlmpar.old[n.par] - stepsize
    nlmpar.newplus[n.par] <- nlmpar.old[n.par] - (-stepsize)

    ll.val.newminus <- loglik_fun_fast(nlmpar.newminus, dset = dataPrep$dset, lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = eps,
                             lambda_delta = opts$lambda_delta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                             fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                             groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_id, fixed_par = fixed_par, fixValue = fixValue,
                             isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                             isPenalized_delta = opts$isPenalized_delta)
    ll.val.newplus <- loglik_fun_fast(nlmpar.newplus, dset = dataPrep$dset, lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = eps,
                                  lambda_delta = opts$lambda_delta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                                  fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                                  groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_id, fixed_par = fixed_par, fixValue = fixValue,
                                  isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                                  isPenalized_delta = opts$isPenalized_delta)

    if(ll.val.newminus < ll.val.old){
      stepsize <- stepsize
      ll.val.new <- ll.val.newminus
      nlmpar.new[n.par] <- nlmpar.newminus[n.par]
    } else if(ll.val.newplus < ll.val.old){
      stepsize <- -stepsize
      ll.val.new <- ll.val.newplus
      nlmpar.new[n.par] <- nlmpar.newplus[n.par]
    } else {
      stepsize <- 0
      step.vec[n.par] <- step.vec[n.par]*scale.down
      ll.val.new <- ll.val.old
    }



    i <- 0
    if(stepsize == 0){
    } else {
      i <- 1
      while((diff <- abs(abs(ll.val.new) - abs(ll.val.old))) > 1e-12 & ll.val.new < ll.val.old & (i < maxit.cd.lower)){
        nlmpar.old<- nlmpar.new
        ll.val.old <- ll.val.new
        nlmpar.new[n.par] <- nlmpar.old[n.par] - stepsize

        ll.val.new <- loglik_fun_fast(nlmpar.new, dset = dataPrep$dset, lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = eps,
                                 lambda_delta = opts$lambda_delta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                                 fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                                 groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_id, fixed_par = fixed_par, fixValue = fixValue,
                                 isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                                 isPenalized_delta = opts$isPenalized_delta)

        diff <- abs(ll.val.new) - abs(ll.val.old)

        if(diff > 0){
          nlmpar.new[n.par] <- nlmpar.new[n.par] + stepsize
          stepsize <- 0
        }


        i <- i+1
      }

    }
  }

  return(list("value" = ll.val.new, "par" = nlmpar.new, iterations = i, step.vec = step.vec))

}

mixed.min.fast <- function(nlmPar, dset, dataPrep, opts,
                      estPar_arr = NULL, estLength_array = NULL, eps = NULL,fixLength_arr = NULL, fixValue = c(),
                      maxit.cd.lower = NULL, maxit.cd.higher = NULL, maxit.optim = NULL, init.step = NULL, abs.tol = NULL, scale.down = NULL,
                      checkNonZero = TRUE, length_array = NULL, max.diff.par = NULL){

  ll.val.old <- 1e+10
  ll.value <- c()

  par.diff <- c()
  par.max <- c()

  par.arr <- c("theta", "beta","gamma","delta")
  par.fix.idx <-  which(par.arr %in% opts$fixed_par)

  fixed_par_gpcm <- unique(c(par.arr[par.fix.idx],"delta"))
  estPar_arr_gpcm <- par.arr[-c(which(par.arr %in% fixed_par_gpcm))]

  fixLength_arr_gpcm <- length_array[c(which(par.arr %in% fixed_par_gpcm),5)]
  estLength_array_gpcm <- length_array[which(par.arr %in% estPar_arr_gpcm)]

  fixed_par_gpcmdif <- c("theta","beta","gamma")
  estPar_arr_gpcmdif <- c("delta")
  fixLength_arr_gpcmdif <- length_array[c(1:3,5)]
  estLength_array_gpcmdif <- length_array[4]

  par.old <- c(nlmPar)
  fixed_delta <- nlmPar[(sum(estLength_array_gpcm)+1):(sum(estLength_array_gpcm)+estLength_array_gpcmdif)]

  nlmPar <- nlmPar[1:sum(estLength_array_gpcm)]

  fixValue <- c(fixed_delta)
  if(length(which("gamma" %in% fixed_par_gpcm)) != 0){
    fixValue <- c(rep(0,(dataPrep$dimResp[2])),fixValue)
  }
  if(!is.null(opts$fixed_beta)){
    fixValue <- c(opts$fixed_beta,fixValue)
  }
  if(!is.null(opts$fixed_theta)){
    fixValue <- c(opts$fixed_theta,fixValue)
  }

  (gpcm <- optim(nlmPar, fn = loglik_fun_fast, gr = grad_fun_fast, hessian = opts$isHessian, dset = dataPrep$dset,
                                  fixValue = fixValue, fixed_par = fixed_par_gpcm, estPar_arr = estPar_arr_gpcm, fixLength_arr = fixLength_arr_gpcm, estLength_array = estLength_array_gpcm,
                                  lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = 0, lambda_delta = opts$lambda_delta,
                                  isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta, isPenalized_delta = FALSE,
                                  allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                                  groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_idx,
                                  method = "BFGS", control = list(maxit = maxit.optim,reltol=1e-12,fnscale = 10)))

  if(length(which("gamma" %in% fixed_par_gpcm)) != 0){
    fixValue <- c(opts$fixed_theta,gpcm$par,rep(0,(dataPrep$dimResp[2])))
  } else {
    fixValue <- c(opts$fixed_theta,gpcm$par)
  }

  nlmPar <- fixed_delta

  (gpcmdif <- coord.descent.fast(nlmPar, dset = dset, dataPrep = dataPrep, opts = opts, fixed_par = fixed_par_gpcmdif,
                              fixValue = fixValue, estPar_arr = estPar_arr_gpcmdif, estLength_array = estLength_array_gpcmdif, fixLength_arr = fixLength_arr_gpcmdif,
                              maxit.cd.lower = maxit.cd.lower, init.step = init.step, eps = eps,
                              scale.down = scale.down, checkNonZero = TRUE))


  ll.val.new <- gpcmdif$value
  par.new <- c(gpcm$par,gpcmdif$par)

  i <- 0

  while((diff <- (abs(abs(ll.val.new) - abs(ll.val.old)))) > abs.tol & i < maxit.cd.higher &
        ((diff.par <- (max(abs(abs(par.new) - abs(par.old))))) > max.diff.par)){
      ll.val.old <- ll.val.new
      par.old <- par.new


      nlmPar <- gpcm$par
      if(is.null(opts$fixed_theta)){
        if(length(which("gamma" %in% fixed_par_gpcm)) != 0){
          fixValue <- c(rep(0,(dataPrep$dimResp[2])),gpcmdif$par)
        } else {
          fixValue <- c(gpcmdif$par)
        }
      } else {
        if(length(which("gamma" %in% fixed_par_gpcm)) != 0){
          fixValue <- c(opts$fixed_theta,rep(0,(dataPrep$dimResp[2])),gpcmdif$par)
        } else {
          fixValue <- c(opts$fixed_theta,gpcmdif$par)
        }
      }
      (gpcm <- optim(nlmPar, loglik_fun_fast, gr = grad_fun_fast, hessian = opts$isHessian, dset = dataPrep$dset,
                     fixValue = fixValue, fixed_par = fixed_par_gpcm, estPar_arr = estPar_arr_gpcm, fixLength_arr = fixLength_arr_gpcm, estLength_array = estLength_array_gpcm,
                     lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out, eps = 0, lambda_delta = opts$lambda_delta,
                     isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta, isPenalized_delta = FALSE,
                     allcat = dataPrep$allcat, dimResp = dataPrep$dimResp, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA, #XREAL = dataPrep$XREAL,
                     groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, mt_idx = dataPrep$mt_idx,
                     method = "BFGS", control = list(maxit = 1e+4,reltol=1e-12,fnscale = 10)))

      if(length(which("gamma" %in% fixed_par_gpcm)) != 0){
        fixValue <- c(opts$fixed_theta,gpcm$par,rep(0,(dataPrep$dimResp[2])))
      } else {
        fixValue <- c(opts$fixed_theta,gpcm$par)
      }
      nlmPar <- gpcmdif$par

      if(i < 20){
        checkNonZero <- TRUE
      } else {
        checkNonZero <- FALSE
      }


      (gpcmdif <- coord.descent.fast(nlmPar, dset = dset, dataPrep = dataPrep, opts = opts, fixed_par = fixed_par_gpcmdif,
                                fixValue = fixValue, estPar_arr = estPar_arr_gpcmdif, estLength_array = estLength_array_gpcmdif, fixLength_arr = fixLength_arr_gpcmdif,
                                maxit.cd.lower = maxit.cd.lower, init.step = init.step, step.vec = gpcmdif$step.vec, eps = eps,
                                scale.down = scale.down, checkNonZero = checkNonZero))



      ll.val.new <- gpcmdif$value
      par.new <- c(gpcm$par,gpcmdif$par)
      par.diff <- c(par.diff, mean(abs(par.new - par.old)))
      par.max <- c(par.max, max(abs(abs(par.new) - abs(par.old))))
      ll.value <- c(ll.value,ll.val.new)

  }

  pardif <- gpcmdif$par
  pardif_idx <- which(abs(pardif) < 1e-10)
  pardif[pardif_idx] <- 0
  par.new <- c(gpcm$par,pardif)

  if(ll.val.new > ll.val.old){
    ll.val.new <- ll.val.old
    par.new <- par.old
  }

  return(list("value" = ll.val.new, "par" = par.new, iterations = i))

}
