# Joint (Full) Maximum Loglikelihood Estimation for PCM
#
# Estimate the theta, beta and gamma parameter of the full dataset
# X is the dataset which is want to be estimated
# opts$fixed_par is the fixed parameter if any whether it is theta, beta, gamma or the combination of them.
# fixVal is the fixed values of the parameters which are fixed.
# gamma_penalized signs whether the gamma/gamma parameters penalized on the likelihood or not.

# jml.gpcm.dif3.est <- function(X, opts$fixed_par = c(), opts$fixed_theta = c(), opts$fixed_beta = c(), opts$fixed_gamma = c(), opts$fixed_deltabeta = c(), gamma_penalized = TRUE, deltabeta_penalized = TRUE, theta_penalized = TRUE, resp.info = c(), resp.th = c(0), opt.method = c("nlminb"), plot.ll = FALSE, opt.tuner = list(iter.max = 20000, eval.max = 30000, rel.tol = 1e-10, step.max = 0.000001), opt.plot = list(), psi = 0.0078, max.iter = 150, objtype = "", desc = NULL,THETA.PCOEFF = 0.05, GAMMA.PCOEFF = 50, SMALLGAMMA.COEFF = 0.000005, DELTABETA.PCOEFF = 10000, DELTAGAMMA.PCOEFF = 10000, eps = 0.0, random.init = FALSE, random.init.th = 1e-2, init.par = c(), hessian = FALSE,tracked = TRUE){

pjmle <- function(X, init_par = c(), ...){

  dset <- as.data.frame(X)                          ### makes sure that the dataset has a matrix format

  dotdotdot <- list(...)

  # print(length(dotdotdot))
  # print(length(dotdotdot[[1]]))
  #
  # if((length(dotdotdot) == 1 & length(dotdotdot[[1]]) == 1) | length(dotdotdot) > 1){
  #   dotdotdot <- dotdotdot
  # } else {
  #   dotdotdot <- dotdotdot[[1]]
  # }

  opts_default <- autoRaschOptions()

  opt_names <- names(dotdotdot)[which(names(dotdotdot) %in% names(opts_default))]
  opts <- opts_default
  for (i in opt_names) {
    if(!is.null(dotdotdot[[i]])){
      opts[[i]] <- dotdotdot[[i]]
      # print(i)
    }
  }
  opts[["fixed_par"]] <- c(opts[["fixed_par"]], "deltagamma")

  # THETA.PCOEFF <- THETA.PCOEFF
  # GAMMA.PCOEFF <- GAMMA.PCOEFF
  # SMALLGAMMA.COEFF <- SMALLGAMMA.COEFF
  # DELTABETA.PCOEFF <- DELTABETA.PCOEFF
  # DELTAGAMMA.PCOEFF <- DELTAGAMMA.PCOEFF
  # eps <- eps

  dataPrep <- data_prep(dset = dset, fixed_par = opts$fixed_par, groups_map = opts$groups_map)

  ### Intializing the parameters ###
  # if(random.init){
  #   theta <- runif(nrow(dataPrep$dset),-1,1)*random.init.th
  #   beta <- runif(allcat,-1,1)*random.init.th
  #   gamma <- runif(ncol(dataPrep$dset),-1,1)*random.init.th
  #   deltabeta <- runif((ncol(dataPrep$dset)*ncol(dataPrep$groups_map)),-1,1)*random.init.th
  #   deltagamma <- runif((ncol(dataPrep$dset)*ncol(dataPrep$groups_map)),-1,1)*random.init.th
  # } else {
    theta <- rep(0,nrow(dataPrep$dset))
    beta <- rep(0,dataPrep$allcat)
    gamma <- rep(0,ncol(dataPrep$dset)) #gpcm uses different gamma for each item
    deltabeta <- rep(0,(ncol(dataPrep$dset)*ncol(dataPrep$groups_map)))
    deltagamma <- rep(0,(ncol(dataPrep$dset)*ncol(dataPrep$groups_map)))
  # }


  ### setting the optimized parameter
  nlmPar <- c(theta,beta,gamma,deltabeta,deltagamma)
  length_theta <- length(theta)
  length_beta <- length(beta)
  length_gamma <- length(gamma)
  length_deltabeta <- length(deltabeta)
  length_deltagamma <- length(deltagamma)
  length_array <- c(length_theta,length_beta,length_gamma,length_deltabeta,length_deltagamma)
  fullPar_arr <- c("theta","beta","gamma","deltabeta","deltagamma")
  if(is.null(opts$fixed_par)){
    estPar_arr <- fullPar_arr
  } else {
    estPar_arr <- fullPar_arr[-c(which(fullPar_arr %in% opts$fixed_par))]
  }
  fixLength_arr <- length_array[c(which(fullPar_arr %in% opts$fixed_par))]

  fixValue <- c()
  estLength_array <- length_array

  if(!identical(grep("deltagamma",opts$fixed_par), integer(0))){
    nlmPar <- nlmPar[-c((sum(length_array[1:4])+1):(sum(length_array[1:5])))]
    estLength_array <- estLength_array[-c(5)]
    fixValue <- c(rep(0,(ncol(dset)*ncol(dataPrep$groups_map))),fixValue)
  }
  if(!identical(grep("deltabeta",opts$fixed_par), integer(0))){
    nlmPar <- nlmPar[-c((sum(length_array[1:3])+1):(sum(length_array[1:4])))]
    estLength_array <- estLength_array[-c(4)]
    if(!is.null(opts$fixed_deltabeta)){
      fixValue <- c(opts$fixed_deltabeta,fixValue)
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

  # fixValue <- fixValue

  ### Maximizing the loglikelihood function ###
  print("Starting to estimate...")
  nameCol <- colnames(as.data.frame(X))
  output <- list("X" = X, "mt_vek" = dataPrep$mt_vek, "real_vek" = dataPrep$na_catVec, "itemName" = nameCol)
  if(!is.null(dotdotdot$groups_map)){
    output[["groups_map"]] <- dataPrep$groups_map
  }
  # if(opt.method == "nlminb"){
  #   time.v <- system.time(minRes <- nlminb(nlmPar, rasch.gpcm.dif3.ll, gradient = rasch.gpcm.dif3.grad, dset = dset, THETA.PCOEFF = THETA.PCOEFF, GAMMA.PCOEFF = GAMMA.PCOEFF, SMALLGAMMA.COEFF = SMALLGAMMA.COEFF, DELTABETA.PCOEFF = DELTABETA.PCOEFF, DELTAGAMMA.PCOEFF = DELTAGAMMA.PCOEFF, eps = eps, estPar_arr = estPar_arr, estLength_array = estLength_array, fixLength_arr = fixLength_arr, allcat = allcat, n_th = n_th, idx.1st.cat = idx.1st.cat, XN = XN, XNA = XNA, groupMap = groupMap, mt_vek = mt_vek, opts$fixed_par = opts$fixed_par, fixVal = fixValue, gamma_penalized = gamma_penalized, theta_penalized = theta_penalized, control = opt.tuner))
  #   est <- minRes$par
  #   obj <- minRes$objective
  #   iter <- minRes$iterations
  #   conv <- 0
  #   counts <- 0
  #   output[["loglik"]] <- -obj
  # } else if(opt.method == "optim"){
    time.v <- system.time(minRes <- optim(nlmPar, loglik_fun, gr = grad_fun, hessian = opts$isHessian, dset = dataPrep$dset,
                                          lambda_theta = opts$lambda_theta, lambda_in = opts$lambda_in,lambda_out = opts$lambda_out,
                                          lambda_deltabeta = opts$lambda_deltabeta, estPar_arr = estPar_arr, estLength_array = estLength_array,
                                          fixLength_arr = fixLength_arr, allcat = dataPrep$allcat, n_th = dataPrep$n_th, XN = dataPrep$XN, XNA = dataPrep$XNA,
                                          groups_map = dataPrep$groups_map, mt_vek = dataPrep$mt_vek, fixed_par = opts$fixed_par, fixValue = fixValue,
                                          isPenalized_gamma = opts$isPenalized_gamma, isPenalized_theta = opts$isPenalized_theta,
                                          isPenalized_deltabeta = opts$isPenalized_deltabeta, method = "BFGS", control = opts$optz_tuner, tracked =  opts$isTracked))

    est <- minRes$par
    obj <- minRes$value
    iter <- minRes$iterations
    conv <- minRes$convergence
    counts <- minRes$counts
    output[["loglik"]] <- -obj
    if(opts$isHessian){
      output[["hessian"]] <- minRes$hessian
    }
  # }
  print("...done!")

  ### Mapping for the output
  if(!identical(grep("deltagamma",estPar_arr), integer(0))){
    parNo <- grep("deltagamma",estPar_arr)
    if(parNo == 1){
      output[["deltagamma"]] <- est[c((1):(sum(estLength_array[1])))]
    } else {
      output[["deltagamma"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  }
  if(!identical(grep("deltabeta",estPar_arr), integer(0))){
    parNo <- grep("deltabeta",estPar_arr)
    if(parNo == 1){
      output[["deltabeta"]] <- est[c((1):(sum(estLength_array[1])))]
    } else {
      output[["deltabeta"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
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

  if((minCat <- min(dset,na.rm = TRUE)) != 0){  ### makes sure the response is started at 0
    dset <- dset - minCat
  }

  ### Checking the DIF group setting
  if(length(which(c("deltabeta") %in% fixed_par)) == 0){
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
  mt_vek <- rep(max(mt_vek),length(mt_vek))
  allcat <- sum(mt_vek)                             ### number of items * categories (assumption : item has the same number of categories)
  if(n_th != 1){
    na_catVec <- as.vector(apply(rbind(mt_vek_min,true_catVec),2,function(x){     ### To remove the unavailable categories if items contain different number of categories (used for output filter)
      temp <- c(rep(NA,x[1]),rep(1,(n_th-sum(x))),rep(NA,x[2]))
      return(temp)
    }))
  } else {
    na_catVec <- mt_vek
  }

  ### Create the xn.mat (matrix of x_ni) and the xna.mat (matrix of missing value) to map the responses and the missing values
  xn.mat <- matrix(0,nrow = nrow(dset), ncol = allcat) ## response position
  xna.mat <- matrix(1,nrow = nrow(dset), ncol = allcat)## NA position

  for(i in 1:n_th){
    idx <- which(dset==i)
    new.idx <- (nrow(dset)*n_th*(ceiling(idx/nrow(dset))-1))+(idx%%nrow(dset))+(ifelse(idx%%nrow(dset) == 0,1,0)*nrow(dset)) ### locate the categories through vector multiplication
    full.idx <- c()
    for(j in 1:i){
      next.idx <- new.idx+(nrow(dset)*(j-1))
      full.idx <- c(full.idx,next.idx)
    }
    xn.mat[full.idx] <- 1
  }

  idx <- which(is.na(dset))
  new.idx <- (nrow(dset)*n_th*(ceiling(idx/nrow(dset))-1))+(idx%%nrow(dset))+(ifelse(idx%%nrow(dset) == 0,1,0)*nrow(dset))
  full.idx <- c()
  for(j in 1:mt_vek[1]){
    next.idx <- new.idx+(nrow(dset)*(j-1))
    full.idx <- c(full.idx,next.idx)
  }
  xn.mat[full.idx] <- NA
  xna.mat[full.idx] <- NA
  XN <- as.vector(t(xn.mat))
  XNA <- as.vector(t(xna.mat))
  ### end of xn.mat and xna.mat creation

  ret <- list("dset" = dset, "mt_vek" = mt_vek, "groups_map" = groups_map, "n_th" = n_th, "na_catVec" = na_catVec, "allcat" = allcat, "XN" = XN, "XNA" = XNA)

  return(ret)

}

# pjmle <- function(X, init_par = c(), isHessian = TRUE,
#                   lambda_theta = 0.05, lambda_in = 50,lambda_out = 5e-3,lambda_deltagamma = 1e+6, eps = 0,
#                   lambda_deltabeta = 15, fixed_par = c(),
#                   isPenalized_gamma = TRUE, isPenalized_theta = TRUE, groups_map = c(),
#                   isPenalized_deltabeta = TRUE, optz_tuner = list(maxit = 1e+7, reltol = 1e-12, fnscale = 10), isTracked = TRUE,
#                   fixed_theta = c(), fixed_beta = c(), fixed_gamma = c(), fixed_deltabeta = c()){
#
#   dset <- as.matrix(X)                          ### makes sure that the dataset has a matrix format
#
#
#   # THETA.PCOEFF <- THETA.PCOEFF
#   # GAMMA.PCOEFF <- GAMMA.PCOEFF
#   # SMALLGAMMA.COEFF <- SMALLGAMMA.COEFF
#   # DELTABETA.PCOEFF <- DELTABETA.PCOEFF
#   # DELTAGAMMA.PCOEFF <- DELTAGAMMA.PCOEFF
#   # eps <- eps
#
#   if((minCat <- min(dset,na.rm = TRUE)) != 0){  ### makes sure the response is started at 0
#     dset <- dset - minCat
#   }
#
#   ### Checking the DIF group setting
#   if(length(which(c("deltabeta") %in% fixed_par)) == 0){
#     groups_map <- as.matrix(groups_map)
#     if(nrow(groups_map) != nrow(dset)){
#       stop("autoRasch ERROR: `groups_map` must be has the same number of rows with observation.")
#     }
#     if(max(groups_map, na.rm = TRUE) > 1 | min(groups_map, na.rm = TRUE) < 0 | length(levels(factor(unlist(groups_map)))) > 2){
#       stop("autoRasch ERROR: `groups_map` must be binary matrix with only 1 and 0 values.")
#     }
#   } else {
#     groups_map <- matrix(0,nrow = nrow(dset),ncol = 1)
#   }
#
#   ### Manages catagory's vector
#   mt_vek <- apply(dset, 2L, max, na.rm = TRUE)      ### create vector of max categories for each item
#   n_th <- max(mt_vek)                               ### number of thresholds
#   mt_vek_min <- apply(dset, 2L, min, na.rm = TRUE)  ### create vector of min categories for each item
#   true_catVec <- (mt_vek-n_th)*(-1)
#   mt_vek <- rep(max(mt_vek),length(mt_vek))
#   allcat <- sum(mt_vek)                             ### number of items * categories (assumption : item has the same number of categories)
#   if(n_th != 1){
#     na_catVec <- as.vector(apply(rbind(mt_vek_min,true_catVec),2,function(x){     ### To remove the unavailable categories if items contain different number of categories (used for output filter)
#       temp <- c(rep(NA,x[1]),rep(1,(n_th-sum(x))),rep(NA,x[2]))
#       return(temp)
#     }))
#   } else {
#     na_catVec <- mt_vek
#   }
#
#   ### Create the xn.mat (matrix of x_ni) and the xna.mat (matrix of missing value) to map the responses and the missing values
#   xn.mat <- matrix(0,nrow = nrow(dset), ncol = allcat) ## response position
#   xna.mat <- matrix(1,nrow = nrow(dset), ncol = allcat)## NA position
#
#   for(i in 1:n_th){
#     idx <- which(dset==i)
#     new.idx <- (nrow(dset)*n_th*(ceiling(idx/nrow(dset))-1))+(idx%%nrow(dset))+(ifelse(idx%%nrow(dset) == 0,1,0)*nrow(dset)) ### locate the categories through vector multiplication
#     full.idx <- c()
#     for(j in 1:i){
#       next.idx <- new.idx+(nrow(dset)*(j-1))
#       full.idx <- c(full.idx,next.idx)
#     }
#     xn.mat[full.idx] <- 1
#   }
#
#   idx <- which(is.na(dset))
#   new.idx <- (nrow(dset)*n_th*(ceiling(idx/nrow(dset))-1))+(idx%%nrow(dset))+(ifelse(idx%%nrow(dset) == 0,1,0)*nrow(dset))
#   full.idx <- c()
#   for(j in 1:mt_vek[1]){
#     next.idx <- new.idx+(nrow(dset)*(j-1))
#     full.idx <- c(full.idx,next.idx)
#   }
#   xn.mat[full.idx] <- NA
#   xna.mat[full.idx] <- NA
#   XN <- as.vector(t(xn.mat))
#   XNA <- as.vector(t(xna.mat))
#   ### end of xn.mat and xna.mat creation
#
#
#   # dataPrep <- data_prep(dset = dset, fixed_par = fixed_par, groups_map = groups_map)
#
#   ### Intializing the parameters ###
#   # if(random.init){
#   #   theta <- runif(nrow(dataPrep$dset),-1,1)*random.init.th
#   #   beta <- runif(allcat,-1,1)*random.init.th
#   #   gamma <- runif(ncol(dataPrep$dset),-1,1)*random.init.th
#   #   deltabeta <- runif((ncol(dataPrep$dset)*ncol(dataPrep$groups_map)),-1,1)*random.init.th
#   #   deltagamma <- runif((ncol(dataPrep$dset)*ncol(dataPrep$groups_map)),-1,1)*random.init.th
#   # } else {
#   theta <- rep(0,nrow(dset))
#   beta <- rep(0,allcat)
#   gamma <- rep(0,ncol(dset)) #gpcm uses different gamma for each item
#   deltabeta <- rep(0,(ncol(dset)*ncol(groups_map)))
#   deltagamma <- rep(0,(ncol(dset)*ncol(groups_map)))
#   # }
#
#
#   ### setting the optimized parameter
#   nlmPar <- c(theta,beta,gamma,deltabeta,deltagamma)
#   length_theta <- length(theta)
#   length_beta <- length(beta)
#   length_gamma <- length(gamma)
#   length_deltabeta <- length(deltabeta)
#   length_deltagamma <- length(deltagamma)
#   length_array <- c(length_theta,length_beta,length_gamma,length_deltabeta,length_deltagamma)
#   fullPar_arr <- c("theta","beta","gamma","deltabeta","deltagamma")
#   if(is.null(fixed_par)){
#     estPar_arr <- fullPar_arr
#   } else {
#     estPar_arr <- fullPar_arr[-c(which(fullPar_arr %in% fixed_par))]
#   }
#   fixLength_arr <- length_array[c(which(fullPar_arr %in% fixed_par))]
#
#   fixValue <- c()
#   estLength_array <- length_array
#
#   if(!identical(grep("deltagamma",fixed_par), integer(0))){
#     nlmPar <- nlmPar[-c((sum(length_array[1:4])+1):(sum(length_array[1:5])))]
#     estLength_array <- estLength_array[-c(5)]
#     fixValue <- c(rep(0,(ncol(dset)*ncol(groups_map))),fixValue)
#   }
#   if(!identical(grep("deltabeta",fixed_par), integer(0))){
#     nlmPar <- nlmPar[-c((sum(length_array[1:3])+1):(sum(length_array[1:4])))]
#     estLength_array <- estLength_array[-c(4)]
#     if(!is.null(fixed_deltabeta)){
#       fixValue <- c(fixed_deltabeta,fixValue)
#     } else {
#       fixValue <- c(rep(0,(ncol(dset)*ncol(groups_map))),fixValue)
#     }
#   }
#   if(!identical(grep("^gamma",fixed_par), integer(0))){
#     if(length(nlmPar) == length(length_array)){
#       nlmPar <- nlmPar[-c((sum(length_array[1:2])+1):(sum(length_array[1:3])),(sum(length_array[1:4])+1):(sum(length_array[1:5])))]
#       estLength_array <- estLength_array[-c(3,5)]
#       if(!is.null(fixed_gamma)){
#         fixValue <- c(fixed_gamma, rep(0,(ncol(dset)*ncol(groups_map))), fixValue)
#       } else {
#         fixValue <- c(rep(0,ncol(dset)), rep(0,(ncol(dset)*ncol(groups_map))), fixValue)
#       }
#     } else {
#       nlmPar <- nlmPar[-c((sum(length_array[1:2])+1):(sum(length_array[1:3])))]
#       estLength_array <- estLength_array[-c(3)]
#       if(!is.null(fixed_gamma)){
#         fixValue <- c(fixed_gamma, fixValue)
#       } else {
#         fixValue <- c(rep(0,ncol(dset)), fixValue)
#       }
#     }
#   }
#   if(!identical(grep("^beta",fixed_par), integer(0))){
#     nlmPar <- nlmPar[-c((sum(length_array[1])+1):(sum(length_array[1:2])))]
#     estLength_array <- estLength_array[-c(2)]
#     fixValue <- c(fixed_beta, fixValue)
#   }
#   if(!identical(grep("theta",fixed_par), integer(0))){
#     nlmPar <- nlmPar[-c((1):(sum(length_array[1])))]
#     estLength_array <- estLength_array[-c(1)]
#     fixValue <- c(fixed_theta, fixValue)
#   }
#
#   if(!is.null(init_par)) {
#     nlmPar <- init_par
#   }
#
#   # fixValue <- fixValue
#
#   # print(list(nlmPar, hessian = isHessian, dset = dset,
#   #   THETA.PCOEFF = lambda_theta, GAMMA.PCOEFF = lambda_in, SMALLGAMMA.COEFF = lambda_out,
#   #   DELTABETA.PCOEFF = lambda_deltabeta, DELTAGAMMA.PCOEFF = lambda_deltabeta, eps = eps,
#   #   estPar.arr = estPar_arr, estLength.array = estLength_array, fixLength.arr = fixLength_arr,
#   #   allcat = allcat, n.th = n_th, idx.1st.cat = c(), XN = XN, XNA = XNA, groupMap = groups_map,
#   #   mt_vek = mt_vek, fixPar = fixed_par, fixVal = fixValue, gamma_penalized = isPenalized_gamma,
#   #   theta_penalized = isPenalized_theta, deltabeta_penalized = isPenalized_deltabeta,
#   #   method = "BFGS", control = optz_tuner, tracked = isTracked))
#
#   ### Maximizing the loglikelihood function ###
#   print("Starting to estimate...")
#   nameCol <- colnames(as.data.frame(X))
#   output <- list("X" = X, "mt_vek" = mt_vek, "itemName" = nameCol)
#   # if(opt.method == "nlminb"){
#   #   time.v <- system.time(minRes <- nlminb(nlmPar, rasch.gpcm.dif3.ll, gradient = rasch.gpcm.dif3.grad, dset = dset, THETA.PCOEFF = THETA.PCOEFF, GAMMA.PCOEFF = GAMMA.PCOEFF, SMALLGAMMA.COEFF = SMALLGAMMA.COEFF, DELTABETA.PCOEFF = DELTABETA.PCOEFF, DELTAGAMMA.PCOEFF = DELTAGAMMA.PCOEFF, eps = eps, estPar_arr = estPar_arr, estLength_array = estLength_array, fixLength_arr = fixLength_arr, allcat = allcat, n_th = n_th, idx.1st.cat = idx.1st.cat, XN = XN, XNA = XNA, groupMap = groupMap, mt_vek = mt_vek, fixed_par = fixed_par, fixVal = fixValue, gamma_penalized = gamma_penalized, theta_penalized = theta_penalized, control = opt.tuner))
#   #   est <- minRes$par
#   #   obj <- minRes$objective
#   #   iter <- minRes$iterations
#   #   conv <- 0
#   #   counts <- 0
#   #   output[["loglik"]] <- -obj
#   # } else if(opt.method == "optim"){
#   time.v <- system.time(minRes <- optim(nlmPar, loglik_fun, gr = grad_fun, hessian = isHessian, dset = dset,
#                                         lambda_theta = lambda_theta, lambda_in = lambda_in, lambda_out = lambda_out,
#                                         lambda_deltabeta = lambda_deltabeta, lambda_deltagamma = 1e+6, eps = 0,
#                                         estPar_arr = estPar_arr, estLength_array = estLength_array,
#                                         fixLength_arr = fixLength_arr, allcat = allcat, n_th = n_th, XN = XN, XNA = XNA,
#                                         groups_map = groups_map, mt_vek = mt_vek, fixed_par = fixed_par, fixValue = fixValue,
#                                         isPenalized_gamma = isPenalized_gamma, isPenalized_theta = isPenalized_theta,
#                                         isPenalized_deltabeta = isPenalized_deltabeta, method = "BFGS", control = optz_tuner,
#                                         tracked =  isTracked))
#   # time.v <- system.time(minRes <- optim(nlmPar, rasch.gpcm.dif3.ll, gr = rasch.gpcm.dif3.grad, hessian = isHessian, dset = dset,
#   #                                       THETA.PCOEFF = lambda_theta, GAMMA.PCOEFF = lambda_in, SMALLGAMMA.COEFF = lambda_out,
#   #                                       DELTABETA.PCOEFF = lambda_deltabeta, DELTAGAMMA.PCOEFF = lambda_deltabeta, eps = eps,
#   #                                       estPar.arr = estPar_arr, estLength.array = estLength_array, fixLength.arr = fixLength_arr,
#   #                                       allcat = allcat, n.th = n_th, idx.1st.cat = c(), XN = XN, XNA = XNA, groupMap = groups_map,
#   #                                       mt_vek = mt_vek, fixPar = fixed_par, fixVal = fixValue, gamma_penalized = isPenalized_gamma,
#   #                                       theta_penalized = isPenalized_theta, deltabeta_penalized = isPenalized_deltabeta,
#   #                                       method = "BFGS", control = optz_tuner, tracked = isTracked))#list(maxit = 10000, lmm = 10, reltol = 1e-15, abstol = 1e-9))) #factr = 1e9, pgtol = 0.2
#
#   est <- minRes$par
#   obj <- minRes$value
#   iter <- minRes$iterations
#   conv <- minRes$convergence
#   counts <- minRes$counts
#   output[["loglik"]] <- -obj
#   if(isHessian){
#     output[["hessian"]] <- minRes$hessian
#   }
#   # }
#   print("...done!")
#
#   ### Mapping for the output
#   if(!identical(grep("deltagamma",estPar_arr), integer(0))){
#     parNo <- grep("deltagamma",estPar_arr)
#     if(parNo == 1){
#       output[["deltagamma"]] <- est[c((1):(sum(estLength_array[1])))]
#     } else {
#       output[["deltagamma"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
#     }
#   }
#   if(!identical(grep("deltabeta",estPar_arr), integer(0))){
#     parNo <- grep("deltabeta",estPar_arr)
#     if(parNo == 1){
#       output[["deltabeta"]] <- est[c((1):(sum(estLength_array[1])))]
#     } else {
#       output[["deltabeta"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
#     }
#   }
#   if(!identical(grep("^gamma",estPar_arr), integer(0))){
#     parNo <- grep("^gamma",estPar_arr)
#     if(parNo == 1){
#       output[["gamma"]] <- est[c((1):(sum(estLength_array[1])))]
#     } else {
#       output[["gamma"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
#     }
#   }
#   if(!identical(grep("^beta",estPar_arr), integer(0))){
#     parNo <- grep("^beta",estPar_arr)
#     if(parNo == 1){
#       output[["beta"]] <- est[c((1):(sum(estLength_array[1])))]*na_catVec
#     } else {
#       output[["beta"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]*na_catVec
#     }
#   }
#   if(!identical(grep("theta",estPar_arr), integer(0))){
#     parNo <- grep("theta",estPar_arr)
#     if(parNo == 1){
#       output[["theta"]] <- est[c((1):(sum(estLength_array[1])))]
#     } else {
#       output[["theta"]] <- est[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
#     }
#   }
#
#   return(output)
#
# }
