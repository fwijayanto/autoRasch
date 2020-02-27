##### Complete Log likelihood and Gradient Function of GPCM with DIF implementation ######

##########################################################################################
# Mapping the full parameters to theta, beta, gamma, deltabeta and deltagamma
#
# nlmPar = parameters those are wanted to be estimated on the loglikelihood function. First, they are mapped.
# estPar_arr = array of names of parameters which are being estimated
# estLength_array = array of size of each parameter which is being estimated
# fixed_par = array of names of parameters which are set to be fixed
# fixLength_arr = array of size of each parameter which is set to be fixed
# fixValue = array of values for the fixed parameters
##########################################################################################

par_map <- function(nlmPar, estPar_arr, estLength_array, fixValue, fixed_par, fixLength_arr){
  output <- list()
  if(!identical(grep("deltagamma",estPar_arr), integer(0))){
    parNo <- grep("deltagamma",estPar_arr)
    if(parNo == 1){
      output[["deltagamma"]] <- nlmPar[c((1):(sum(estLength_array[1])))]
    } else {
      output[["deltagamma"]] <- nlmPar[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  } else {
    parNo <- grep("deltagamma",fixed_par)
    if(parNo == 1){
      output[["deltagamma"]] <- fixValue[c((1):(sum(fixLength_arr[1])))]
    } else {
      output[["deltagamma"]] <- fixValue[c((sum(fixLength_arr[1:(parNo-1)])+1):(sum(fixLength_arr[1:parNo])))]
    }
  }
  if(!identical(grep("deltabeta",estPar_arr), integer(0))){
    parNo <- grep("deltabeta",estPar_arr)
    if(parNo == 1){
      output[["deltabeta"]] <- nlmPar[c((1):(sum(estLength_array[1])))]
    } else {
      output[["deltabeta"]] <- nlmPar[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  } else {
    parNo <- grep("deltabeta",fixed_par)
    if(parNo == 1){
      output[["deltabeta"]] <- fixValue[c((1):(sum(fixLength_arr[1])))]
    } else {
      output[["deltabeta"]] <- fixValue[c((sum(fixLength_arr[1:(parNo-1)])+1):(sum(fixLength_arr[1:parNo])))]
    }
  }
  if(!identical(grep("^gamma",estPar_arr), integer(0))){
    parNo <- grep("^gamma",estPar_arr)
    if(parNo == 1){
      output[["gamma"]] <- nlmPar[c((1):(sum(estLength_array[1])))]
    } else {
      output[["gamma"]] <- nlmPar[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  } else {
    parNo <- grep("^gamma",fixed_par)
    if(parNo == 1){
      output[["gamma"]] <- fixValue[c((1):(sum(fixLength_arr[1])))]
    } else {
      output[["gamma"]] <- fixValue[c((sum(fixLength_arr[1:(parNo-1)])+1):(sum(fixLength_arr[1:parNo])))]
    }
  }
  if(!identical(grep("^beta",estPar_arr), integer(0))){
    parNo <- grep("^beta",estPar_arr)
    if(parNo == 1){
      output[["beta"]] <- nlmPar[c((1):(sum(estLength_array[1])))]
    } else {
      output[["beta"]] <- nlmPar[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  } else {
    parNo <- grep("^beta",fixed_par)
    if(parNo == 1){
      output[["beta"]] <- fixValue[c((1):(sum(fixLength_arr[1])))]
    } else {
      output[["beta"]] <- fixValue[c((sum(fixLength_arr[1:(parNo-1)])+1):(sum(fixLength_arr[1:parNo])))]
    }
  }
  if(!identical(grep("theta",estPar_arr), integer(0))){
    parNo <- grep("theta",estPar_arr)
    if(parNo == 1){
      output[["theta"]] <- nlmPar[c((1):(sum(estLength_array[1])))]
    } else {
      output[["theta"]] <- nlmPar[c((sum(estLength_array[1:(parNo-1)])+1):(sum(estLength_array[1:parNo])))]
    }
  } else {
    parNo <- grep("theta",fixed_par)
    if(parNo == 1){
      output[["theta"]] <- fixValue[c((1):(sum(fixLength_arr[1])))]
    } else {
      output[["theta"]] <- fixValue[c((sum(fixLength_arr[1:(parNo-1)])+1):(sum(fixLength_arr[1:parNo])))]
    }
  }

  return(output)
}

##########################################################################################
#THE LOG-LIKELIHOOD FUNCTION
#
# nlmPar = parameters which are going to be estimated
# dset = the dataset
# zeroIdx = the nlmPar index which the delta(beta/gamma) very small and need to be rounded to zero. (in evaluation since the gradient descent failed in operation)
# allcat = the number of beta parameters for all items
# n_th = the number of beta for each items (the assumption for now is that every item has the same number of thresholds)
# XN = matrix for mapping x_vi (true response of person v to item i)
# XNA = matrix for mapping the NA responses
# groups_map = a vector or matrix to map person v
##########################################################################################

loglik_fun <- function(nlmPar, dset, lambda_theta, lambda_in, lambda_out, lambda_deltabeta, lambda_deltagamma = 1e+6, eps = 0,
                       estPar_arr, fixLength_arr, estLength_array, allcat, n_th, XN, XNA, groups_map, mt_vek, fixed_par, fixValue,
                       isPenalized_gamma, isPenalized_theta, isPenalized_deltabeta, tracked){


  ### map the nlmPar
  map_nlmPar <- par_map(nlmPar, estPar_arr = estPar_arr, estLength_array = estLength_array, fixValue = fixValue, fixed_par = fixed_par, fixLength_arr = fixLength_arr)
  theta <- map_nlmPar$theta
  beta <- map_nlmPar$beta
  gamma <- map_nlmPar$gamma
  deltabeta <- map_nlmPar$deltabeta
  deltagamma <- map_nlmPar$deltagamma

  ### get the value of alpha
  exp_gamma <- exp(gamma)
  exp_gamma <- rep.int(exp_gamma, mt_vek)

  groups_map <- as.matrix(groups_map) #groups_map should be formed as matrix

  ### take the total of DIF effect for every group for beta or gamma
  deltabeta_tot <- 0
  deltagamma_tot <- 0
  for(i in 1:ncol(groups_map)){
    deltabeta_tot <- deltabeta_tot + outer(deltabeta[(((i-1)*ncol(dset))+1):(i*ncol(dset))],groups_map[,i],"*")
    deltagamma_tot <- deltagamma_tot + outer(deltagamma[(((i-1)*ncol(dset))+1):(i*ncol(dset))],groups_map[,i],"*")
  }
  deltabeta_tot_rep <- rep.int((deltabeta_tot), rep.int(mt_vek,nrow(groups_map)))       #deltabeta_tot_rep is total deltabeta which has been replicated to every categoory
  exp_deltagamma_tot_rep <- rep.int((exp(deltagamma_tot)), rep.int(mt_vek,nrow(groups_map)))     #deltagamma_tot.rep is total deltagamma which has been replicated to every categoory

  ### take the value of DIF on alpha
  #exp_deltagamma_tot_rep <- exp(deltagamma_tot.rep)

  ### compute the theta - (beta+deltabeta)
  t_diff <- rep(theta,each = allcat) - rep.int(beta,length(theta))
  t_diff <- t_diff - deltabeta_tot_rep

  ### multiplied by (exp(gamma+deltagamma))
  disc_diff <- t_diff * exp_gamma
  disc_diff <- disc_diff * exp_deltagamma_tot_rep

  ### map the corresponding NA value of the dataset to the matrix
  disc_diff <- XNA * disc_diff

  ### compute the first part of the log-likelihood (simple addition part)
  l1 <- sum((XN * disc_diff), na.rm = TRUE)

  ### compute the second part of the log-likelohood (with log)
  ### begin
  per_cat_list <- matrix(disc_diff, nrow = n_th)
  temp_prob <- as.matrix(per_cat_list[1,])
  temp_l2 <- exp(temp_prob)
  if(n_th > 1){
    for(i in 2:n_th){
      temp_prob <- cbind(temp_prob,(temp_prob[,i-1]+per_cat_list[i,]))
      temp_l2 <- temp_l2 + (exp(temp_prob[,i]))
    }
  }

  l2 <- sum(log(temp_l2+1), na.rm = TRUE)

  ### end

  st1st2 <- c(l1,l2)

  lnL <- st1st2[1] - st1st2[2]

  if(isPenalized_theta){
    lnL <- lnL - (lambda_theta*(sum(theta^2)))
  }

  if(isPenalized_gamma & isPenalized_theta){
    lnL <- lnL - (lambda_in*(sum(gamma^2)))
  } else if(isPenalized_gamma){
    lnL <- lnL - (lambda_out*(sum(gamma^2)))
  }

  if(isPenalized_deltabeta){
    lnL <- lnL - (lambda_deltabeta*(sum(abs(deltabeta)^(1+eps))))
  }

  #lnL <- lnL - (lambda_deltagamma*(sum(abs(deltagamma)^(1+eps))))

  # if(isPenalized_gamma & isPenalized_theta & isPenalized_deltabeta){
  #   lnL <- st1st2[1] - st1st2[2] - (lambda_theta*(sum(theta^2))) - (lambda_in*(sum(gamma^2)))- (lambda_deltabeta*(sum(abs(deltabeta)^(1+eps))))- (lambda_deltagamma*(sum(abs(deltagamma)^(1+eps))))
  # } else if(isPenalized_gamma){
  #   lnL <- st1st2[1] - st1st2[2] - (lambda_out*(sum(gamma^2)))- (lambda_deltabeta*(sum(abs(deltabeta)^(1+eps))))- (lambda_deltagamma*(sum(abs(deltagamma)^(1+eps))))
  # } else if(isPenalized_theta){
  #   lnL <- st1st2[1] - st1st2[2] - (lambda_theta*(sum(theta^2)))- (lambda_deltabeta*(sum(abs(deltabeta)^(1+eps))))- (lambda_deltagamma*(sum(abs(deltagamma)^(1+eps))))
  # } else {
  #   lnL <- st1st2[1] - st1st2[2]
  # }

  return(-lnL)
}

############################################################################
#THE GRADIENT FUNCTION
#
#
############################################################################

#rasch.gpcm.dif3.grad <- function(nlmPar,dset = c(),zeroIdx = c(), lambda_theta, lambda_in, lambda_out, lambda_deltabeta, lambda_deltagamma, eps, estPar_arr = c(), fixLength_arr = c(), estLength_array = c(), allcat = c(), n_th = c(), idx.1st.cat = c(), XN = c(), XNA = c(), groups_map = c(), mt_vek = c(),fixed_par = c(),fixValue = c(), isPenalized_gamma = TRUE, isPenalized_theta = TRUE,tracked = TRUE, isPenalized_deltabeta = TRUE){
grad_fun <- function(nlmPar, dset, lambda_theta, lambda_in, lambda_out, lambda_deltabeta, lambda_deltagamma = 1e+6, eps = 0,
                         estPar_arr, fixLength_arr, estLength_array, allcat, n_th, XN, XNA, groups_map, mt_vek, fixed_par, fixValue,
                         isPenalized_gamma, isPenalized_theta, isPenalized_deltabeta, tracked){

  lambda_theta <- lambda_theta*2
  lambda_in <- lambda_in*2
  lambda_out <- lambda_out*2

  map_nlmPar <- par_map(nlmPar, estPar_arr = estPar_arr, estLength_array = estLength_array, fixValue = fixValue, fixed_par = fixed_par, fixLength_arr = fixLength_arr)

  theta <- map_nlmPar$theta
  beta <- map_nlmPar$beta
  gamma <- map_nlmPar$gamma
  deltabeta <- map_nlmPar$deltabeta
  deltagamma <- map_nlmPar$deltagamma

  exp_gamma <- exp(gamma)
  exp_gamma <- rep.int(exp_gamma, mt_vek)
  t_exp_gamma_mat <- rep.int(1,(length(XN)))*exp_gamma

  groups_map <- as.matrix(groups_map)

  deltabeta_tot <- 0
  deltagamma_tot <- 0
  for(i in 1:ncol(groups_map)){
    deltabeta_tot <- deltabeta_tot + outer(deltabeta[(((i-1)*ncol(dset))+1):(i*ncol(dset))],groups_map[,i],"*")
    deltagamma_tot <- deltagamma_tot + outer(deltagamma[(((i-1)*ncol(dset))+1):(i*ncol(dset))],groups_map[,i],"*")
  }
  deltabeta_tot_rep <- rep.int((deltabeta_tot), rep.int(mt_vek,nrow(groups_map)))
  exp_deltagamma_tot_rep <- rep.int((exp(deltagamma_tot)), rep.int(mt_vek,nrow(groups_map)))
  # exp_deltagamma_tot_rep <- exp(deltagamma_tot.rep)
  total_alpha_mat <- exp_deltagamma_tot_rep * exp_gamma

  # t_diff <- outer((-beta), theta, "+")
  # tdiff <- as.vector(t_diff)
  t_diff <- rep(theta,each = allcat) - rep.int(beta,length(theta))

  #t_diff <- t_diff - group.rep.deltabeta  #delta over beta should be sum allover delta on certain group
  t_diff <- t_diff - deltabeta_tot_rep          #delta.tot.rep is total delta which has been replicated to every categoory
  disc_diff <- t_diff * total_alpha_mat
  # disc_diff <- t_diff * exp_gamma
  # disc_diff <- disc_diff * exp_deltagamma_tot_rep
  disc_diff <- XNA * disc_diff

  ##### first part of gradient ###

  ## theta ##
  t1_theta_mat <- (XN) * total_alpha_mat
  t1_theta_mat <- matrix(t1_theta_mat, nrow = allcat)
  t1_theta <- colSums(t1_theta_mat, na.rm = TRUE)

  ## beta ##
  t1_beta_mat <- t1_theta_mat
  t1_beta_mat <- matrix(t1_beta_mat, nrow = allcat)
  t1_beta <- rowSums(t1_beta_mat, na.rm = TRUE)

  ## gamma ##
  t1_gamma_mat <- t_diff * t1_beta_mat
  t1_gamma_mat <- matrix(t1_gamma_mat, nrow = allcat)
  t1_gamma <- rowSums(t1_gamma_mat, na.rm = TRUE)
  t1_gamma <- colSums(matrix(t1_gamma, nrow = n_th),na.rm = TRUE)

  t1_deltabeta <- t1_deltagamma <- c()
  for(i in 1:ncol(groups_map)){
    ## deltabeta ##
    t1_deltabeta_mat <- t(t(t1_theta_mat) * groups_map[,i])
    t1_deltabeta_temp <- rowSums(t1_deltabeta_mat, na.rm = TRUE)
    t1_deltabeta_temp <- colSums(matrix(t1_deltabeta_temp, nrow = n_th),na.rm = TRUE)
    t1_deltabeta <- c(t1_deltabeta, t1_deltabeta_temp)

    ## deltagamma ##
    t1_deltagamma_mat <- t(t(t1_gamma_mat) * groups_map[,i])
    t1_deltagamma_temp <- rowSums(t1_deltagamma_mat, na.rm = TRUE)
    t1_deltagamma_temp <- colSums(matrix(t1_deltagamma_temp, nrow = n_th),na.rm = TRUE)
    t1_deltagamma <- c(t1_deltagamma, t1_deltagamma_temp)
  }


  ##### first part of gradient end ###

  per_cat_list_dif <- matrix(disc_diff,nrow = n_th)
  per_cat_list_gamma <- matrix(t_exp_gamma_mat, nrow = n_th)
  per_cat_list_totalalpha <- matrix(total_alpha_mat, nrow = n_th)
  per_cat_list_grmap <- matrix(rep(groups_map,each = ncol(dset)), ncol = ncol(groups_map))
  per_cat_list_tdif <- matrix(t_diff, nrow = n_th)

  # temp_prob <- as.matrix(per_cat_list_dif[1,])
  # temp_gamma <- as.matrix(per_cat_list_totalalpha[1,])
  # temp_denom <- exp(temp_prob)
  # temp_theta_nom <- (temp_denom*temp_gamma)
  # temp_deltabeta_nom <- (per_cat_list_grmap*as.vector(temp_denom*temp_gamma))
  #
  # if(n_th > 1){
  #   for(i in 2:n_th){
  #     temp_prob <- cbind(temp_prob,(temp_prob[,i-1]+per_cat_list_dif[i,]))
  #     temp_gamma <- cbind(temp_gamma,(temp_gamma[,i-1]+per_cat_list_totalalpha[i,]))
  #     expTempProb <- exp(temp_prob[,i])
  #     temp_denom <- temp_denom + (expTempProb)
  #     temp_theta_nom <- temp_theta_nom + (expTempProb*temp_gamma[,i])
  #     temp_deltabeta_nom <- temp_deltabeta_nom + (per_cat_list_grmap*as.vector(expTempProb*temp_gamma[,i]))
  #   }
  # }

  temp_prob <- (per_cat_list_dif[1,])
  temp_gamma <- (per_cat_list_totalalpha[1,])
  temp_denom <- exp(temp_prob)
  temp_theta_nom <- (temp_denom*temp_gamma)
  temp_deltabeta_nom <- (per_cat_list_grmap*as.vector(temp_denom*temp_gamma))
  temp_gamma_nom <- (temp_prob*exp(temp_prob))
  temp_deltagamma_nom <- per_cat_list_grmap*as.vector(temp_gamma_nom)
  if(n_th > 1){
    for(i in 2:n_th){
      temp_prob <- temp_prob+per_cat_list_dif[i,]
      temp_gamma <- temp_gamma+per_cat_list_totalalpha[i,]
      expTempProb <- exp(temp_prob)
      temp_denom <- temp_denom + (expTempProb)
      temp_theta_nom <- temp_theta_nom + (expTempProb*temp_gamma)
      temp_gamma_nom <- temp_gamma_nom + (expTempProb*temp_prob)
      temp_deltabeta_nom <- temp_deltabeta_nom + (per_cat_list_grmap*as.vector(expTempProb*temp_gamma))
      temp_deltagamma_nom <- temp_deltagamma_nom + (per_cat_list_grmap*as.vector(expTempProb*temp_prob))
    }
  }


  # temp_prob <- as.matrix(colSums(per_cat_list_dif))
  # temp_tot <- as.vector(exp(temp_prob))
  # temp_beta_nom <- temp_tot * per_cat_list_totalalpha[n_th,]
  # if(n_th > 1){
  #   for(i in (n_th-1):1){
  #     temp_prob <- cbind(temp_prob[,1] - per_cat_list_dif[i+1,],temp_prob)
  #     temp_tot <- temp_tot + as.vector(exp(temp_prob[,1]))
  #     temp_beta_nom <- rbind(as.vector(temp_tot) * as.vector(per_cat_list_totalalpha[i,]),temp_beta_nom)
  #   }
  # }

  temp_prob <- as.matrix(colSums(per_cat_list_dif))
  temp_tot <- as.vector(exp(temp_prob))
  temp_beta_nom <- temp_tot * per_cat_list_totalalpha[n_th,]
  if(n_th > 1){
    for(i in (n_th-1):1){
      temp_prob <-temp_prob - per_cat_list_dif[i+1,]
      temp_tot <- temp_tot + as.vector(exp(temp_prob))
      temp_beta_nom <- rbind(as.vector(temp_tot) * as.vector(per_cat_list_totalalpha[i,]),temp_beta_nom)
    }
  }


  # temp_prob <- as.matrix(per_cat_list_dif[1,])
  # temp_gamma_nom <- (temp_prob*exp(temp_prob))
  # temp_deltagamma_nom <- per_cat_list_grmap*as.vector(temp_gamma_nom)
  #
  # if(n_th > 1){
  #   for(i in 2:n_th){
  #     temp_prob <- cbind(temp_prob,(temp_prob[,i-1]+per_cat_list_dif[i,]))
  #     expTempProb <- exp(temp_prob[,i])
  #     temp_gamma_nom <- temp_gamma_nom + (expTempProb*temp_prob[,i])
  #     temp_deltagamma_nom <- temp_deltagamma_nom + (per_cat_list_grmap*as.vector(expTempProb*temp_prob[,i]))
  #   }
  # }

  Nom_theta <- matrix(temp_theta_nom,nrow = ncol(dset))
  Nom_beta <- matrix(temp_beta_nom,nrow = allcat)
  Nom_gamma <- matrix(temp_gamma_nom,nrow = ncol(dset))


  Denom <- matrix(temp_denom,nrow = ncol(dset))+1
  Denom_mat <- matrix(rep(Denom, each = n_th), nrow = allcat)


  t2_theta_mat <- (Nom_theta/Denom)
  t2_theta <- colSums(t2_theta_mat, na.rm = TRUE)

  t2_beta_mat <- (Nom_beta/Denom_mat)
  t2_beta <- rowSums(t2_beta_mat, na.rm = TRUE)

  t2_gamma_mat <- (Nom_gamma/Denom)
  t2_gamma <- rowSums(t2_gamma_mat, na.rm = TRUE)

  t2_deltabeta <- t2_deltagamma <- c()
  for(i in 1:ncol(groups_map)){
    Nom_deltabeta <- matrix(temp_deltabeta_nom[,i], nrow = ncol(dset))
    t2_deltabeta_mat <- (Nom_deltabeta/Denom)
    t2_deltabeta <- c(t2_deltabeta, rowSums(t2_deltabeta_mat, na.rm = TRUE))

    Nom_deltagamma <- matrix(temp_deltagamma_nom[,i], nrow = ncol(dset))
    t2_deltagamma_mat <- (Nom_deltagamma/Denom)
    t2_deltagamma <- c(t2_deltagamma, rowSums(t2_deltagamma_mat, na.rm = TRUE))
  }


  if(isPenalized_theta){
    # grad_theta <- t1_theta - t2_theta - (0.1*theta)
    grad_theta <- t1_theta - t2_theta - (lambda_theta*theta)
  } else {
    grad_theta <- t1_theta - t2_theta
  }

  if(isPenalized_theta){
    grad_beta <- (-t1_beta) + t2_beta# - (50/(length(beta)))
  } else {
    grad_beta <- (-t1_beta) + t2_beta
  }

  if(isPenalized_gamma & isPenalized_theta){
    # grad_gamma <- t1_gamma - t2_gamma - (100*gamma)
    grad_gamma <- t1_gamma - t2_gamma - (lambda_in*(gamma))
  }else if(isPenalized_gamma){
    grad_gamma <- t1_gamma - t2_gamma - (lambda_out*(gamma))
  }else {
    grad_gamma <- t1_gamma - t2_gamma
  }

  if(isPenalized_deltabeta){
    grad_deltabeta <- (-t1_deltabeta) + t2_deltabeta - (lambda_deltabeta*(1+eps)*sign(deltabeta)*(abs(deltabeta)^eps))
  } else {
    grad_deltabeta <- (-t1_deltabeta) + t2_deltabeta
  }

  grad_deltagamma <- (t1_deltagamma) - t2_deltagamma #- (lambda_deltagamma*(1+eps)*sign(deltagamma)*(abs(deltagamma)^eps))

  output <- c()

  if(!identical(grep("deltagamma",estPar_arr), integer(0))){
    output <- c(grad_deltagamma,output)
  }
  if(!identical(grep("deltabeta",estPar_arr), integer(0))){
    output <- c(grad_deltabeta,output)
  }
  if(!identical(grep("^gamma",estPar_arr), integer(0))){
    output <- c(grad_gamma,output)
  }
  if(!identical(grep("^beta",estPar_arr), integer(0))){
    output <- c(grad_beta,output)
  }
  if(!identical(grep("theta",estPar_arr), integer(0))){
    output <- c(grad_theta,output)
  }


  # if(!is.null(zeroIdx)){
  #   grad_tot <- output[-c(zeroIdx)]
  # } else {
    grad_tot <- output
  # }

  return(-grad_tot)

}

GPCMDIF_LL <- function(dset = c(), fixValue.theta, fixValue.beta, fixValue.gamma, fixValue.deltabeta, lambda_theta = 0.05, lambda_in = 50, lambda_out = 5e-3, lambda_deltabeta = 17, lambda_deltagamma = 100000, eps = 0, resp.info = c(), resp.th = c(1,1), isPenalized_gamma = TRUE, isPenalized_theta = TRUE, isPenalized_deltabeta = TRUE){

  lambda_theta <- lambda_theta
  lambda_in <- lambda_in
  lambda_out <- lambda_out
  lambda_deltabeta <- lambda_deltabeta
  lambda_deltagamma <- lambda_deltagamma
  eps <- eps

  groups_map <- as.matrix(resp.info)
  mt_vek <- apply(dset, 2L, max, na.rm = TRUE)   #number of categories - 1 for each item
  mt_vek <- rep(max(mt_vek),length(mt_vek))
  allcat <- sum(mt_vek)         #number of items * categories (assumption : item has the same number of categories)
  n_th <- max(mt_vek)

  xn.mat <- matrix(0,nrow = nrow(dset), ncol = allcat) ## response position
  xna.mat <- matrix(1,nrow = nrow(dset), ncol = allcat)## NA position

  for(i in 1:n_th){
    idx <- which(dset==i)
    new.idx <- (nrow(dset)*n_th*(ceiling(idx/nrow(dset))-1))+(idx%%nrow(dset))+(ifelse(idx%%nrow(dset) == 0,1,0)*nrow(dset))
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
  # xn.mat[full.idx] <- 9
  XN <- t(xn.mat)  #need to be transposed from the original dataset form
  XN <- as.vector(XN)
  XNA <- t(xna.mat)
  XNA <- as.vector(XNA)

  theta <- fixValue.theta
  beta <- fixValue.beta
  gamma <- fixValue.gamma
  deltabeta <- fixValue.deltabeta
  deltagamma <- rep(0,length(deltabeta))

  ### get the value of alpha
  exp_gamma <- exp(gamma)
  exp_gamma <- rep.int(exp_gamma, mt_vek)

  groups_map <- as.matrix(groups_map) #groups_map should be formed as matrix

  ### take the total of DIF effect for every group for beta or gamma
  deltabeta_tot <- 0
  deltagamma_tot <- 0
  for(i in 1:ncol(groups_map)){
    deltabeta_tot <- deltabeta_tot + outer(deltabeta[(((i-1)*ncol(dset))+1):(i*ncol(dset))],groups_map[,i],"*")
    deltagamma_tot <- deltagamma_tot + outer(deltagamma[(((i-1)*ncol(dset))+1):(i*ncol(dset))],groups_map[,i],"*")
  }
  deltabeta_tot_rep <- rep.int((deltabeta_tot), rep.int(mt_vek,nrow(groups_map)))       #deltabeta_tot_rep is total deltabeta which has been replicated to every categoory
  exp_deltagamma_tot_rep <- rep.int((exp(deltagamma_tot)), rep.int(mt_vek,nrow(groups_map)))     #deltagamma_tot.rep is total deltagamma which has been replicated to every categoory

  ### take the value of DIF on alpha
  #exp_deltagamma_tot_rep <- exp(deltagamma_tot.rep)

  ### compute the theta - (beta+deltabeta)
  t_diff <- rep(theta,each = allcat) - rep.int(beta,length(theta))
  t_diff <- t_diff - deltabeta_tot_rep

  ### multiplied by (exp(gamma+deltagamma))
  disc_diff <- t_diff * exp_gamma
  disc_diff <- disc_diff * exp_deltagamma_tot_rep

  ### map the corresponding NA value of the dataset to the matrix
  disc_diff <- XNA * disc_diff

  ### compute the first part of the log-likelihood (simple addition part)
  l1 <- sum((XN * disc_diff), na.rm = TRUE)

  ### compute the second part of the log-likelohood (with log)
  ### begin
  per_cat_list <- matrix(disc_diff, nrow = n_th)
  temp_prob <- as.matrix(per_cat_list[1,])
  temp_l2 <- exp(temp_prob)
  if(n_th > 1){
    for(i in 2:n_th){
      temp_prob <- cbind(temp_prob,(temp_prob[,i-1]+per_cat_list[i,]))
      temp_l2 <- temp_l2 + (exp(temp_prob[,i]))
    }
  }

  l2 <- sum(log(temp_l2+1), na.rm = TRUE)

  ### end

  # print(l1)
  # print(l2)
  # print((lambda_theta*(sum(theta^2))))
  # print((lambda_in*(sum(gamma^2))))
  # print(lambda_deltabeta*(sum(abs(deltabeta)^(1+eps))))
  # print(lambda_deltagamma*(sum(abs(deltagamma)^(1+eps))))
  st1st2 <- c(l1,l2)

  if(isPenalized_gamma & isPenalized_theta & isPenalized_deltabeta){
    lnL <- st1st2[1] - st1st2[2] - (lambda_theta*(sum(theta^2))) - (lambda_in*(sum(gamma^2)))- (lambda_deltabeta*(sum(abs(deltabeta)^(1+eps))))- (lambda_deltagamma*(sum(abs(deltagamma)^(1+eps))))
  } else if(isPenalized_gamma){
    lnL <- st1st2[1] - st1st2[2] - (lambda_out*(sum(gamma^2)))- (lambda_deltabeta*(sum(abs(deltabeta)^(1+eps))))- (lambda_deltagamma*(sum(abs(deltagamma)^(1+eps))))
  } else if(isPenalized_theta){
    lnL <- st1st2[1] - st1st2[2] - (lambda_theta*(sum(theta^2)))- (lambda_deltabeta*(sum(abs(deltabeta)^(1+eps))))- (lambda_deltagamma*(sum(abs(deltagamma)^(1+eps))))
  } else {
    print("in")
    lnL <- st1st2[1] - st1st2[2]
  }

  # loglik <- c(get("loglik"),-lnL)
  # iter <- get("iter") + 1

  # assign("loglik",loglik, envir = .GlobalEnv)
  # assign("iter",iter, envir = .GlobalEnv)

  return(lnL)
}
