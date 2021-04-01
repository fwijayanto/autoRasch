#' Stepwise Selection Search
#'
#' To search itemset that give maximum value of the criterion
#'
#' @param X A matrix or data.frame of the observed responses (ordinal or binary response).
#' @param criterion The criterion that should be used. The default is ipoqll.
#' @param incl_set A vector of initial items in the included set to start the search. The default is to start with full items.
#' @param groups_map A matrix or vector to map the subject to the DIFs groups.
#' @param cores An integer value of number of cores should be used for computation. The default is 2.
#' @param isLegacy A logical value whether using the previous step estimated parameters as the initial parameter values for the upcoming search step.
#' @param isTracked A logical value whether the progress need to be tracked or not.
#' @param isContinued A logical value whether this search is continuing another unfinished search.
#' @param prevData The filename of the temporary .RData file of the unfinished search.
#' @param fileOutput The filename if it is wished to save the output results in file (.RData and .csv) and FALSE if not.
#' @param tempFile The filename of the temporary file to track the search progress. The default is \code{"temp_stepSearch.RData"} which also automatically produces \code{"temp_stepSearch.csv"}.
#' @param optim_control_iq The optimisation setting of the included set. See \code{\link[stats:optim]{stats::optim()}} \code{control} parameter.
#' @param optim_control_oq The optimisation setting of the excluded set. See \code{\link[stats:optim]{stats::optim()}} \code{control} parameter.
#' @param isConvert A logical value whether it is wanted to recompute the score of the search results using IPOQ-LL-DIF criterion.
#' @param setting_par_iq a list of the optimization control setting parameters for the included set. See \code{setting} parameter in \code{\link[autoRasch:autoRaschOptions]{autoRaschOptions()}}.
#' @param setting_par_oq a list of the optimization control setting parameters for the included set. See \code{setting} parameter in \code{\link[autoRasch:autoRaschOptions]{autoRaschOptions()}}.
#'
#' @return
#' A matrix of the itemsets that obtain the highest scores for each number of items in the included set and it scores (IQ-LL,OQ-LL, and IPOQ-LL).
#'
#'
#' @details
#' To search the itemset that give the maximum score.
#'
<<<<<<< HEAD
#' @import doParallel
#' @import foreach
#'
#' @rdname search
#' @export
stepwise_search <- function(X, criterion = c("ipoqll","ipoqlldif") , incl_set = c(), groups_map = c(), cores = NULL,
                            optim_control_iq = c(), optim_control_oq = c(), isTracked = TRUE, isContinued = FALSE,
                            prevData = c(), isLegacy = FALSE, fileOutput = FALSE, tempFile = "temp_stepSearch.RData",
                            isConvert = FALSE, setting_par_iq = c(), setting_par_oq = c()){

# stepwise_search <- function(X, criterion = c("ipoqll","ipoqlldif") , incl_set = c(), groups_map = c(), cores = NULL,
#                             optim_control_iq = c(), optim_control_oq = c(), isTracked = TRUE, isContinued = FALSE,
#                             prevData = c(), isLegacy = FALSE, fileOutput = FALSE, tempFile = "temp_stepSearch.RData",
#                             isConvert = FALSE, lambda_delta = 10,lambda_delta_out = 10,
#                             random.init = TRUE, random.init.th = 1e-2, eps = 0, lambda_out =1, opt.method = "optim",
#                             max.iter = 2, abs.tol = 1e-12, scale.down = 0.3, psi.mult.up = 0.3,
#                             checkNonZero = TRUE, step = 0.075, maxit = 2,max.diff.par = 1e-8){
=======
#' @examples
#' #pcmdata_search <- stepwise_search(X = pcm_data, incl_set = c(1:ncol(pcm_data)), cores = 2)
#' #plot(pcmdata_search)
#'
#' #To search only using backward search
#' #pcmdata_search <- backward_search(X = pcm_data, incl_set = c(1:ncol(pcm_data)))
#'
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import utils
#'
#' @rdname search
#' @export
stepwise_search <- function(X, criterion = c("ipoqll") , incl_set = c(), groups_map = c(), cores = 20,
                            optz_tuner_iq = list(maxit=5e+2, reltol=1e-21, fnscale = 10),
                            optz_tuner_oq = list(maxit=5e+2, reltol=1e-21, fnscale = 10), isTracked = TRUE,
                            isContinued = FALSE, prevData = c(), isLegacy = TRUE, fileOutput = FALSE,
                            tempFile = "temp_stepSearch.RData", isConvert = FALSE){
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5

  namecsv <- paste(paste(strsplit(tempFile, "(\\.)")[[1]][1:(length(strsplit(tempFile, "(\\.)")[[1]])-1)],collapse = "."),".csv",sep = "")
  fullitem <- c(1:ncol(X))

  if(is.null(incl_set)){
    incl_set <- c(1:ncol(X))
  }

  if(isLegacy){
    optim_control_iq <- list(maxit=5e+2, reltol=1e-21, fnscale = 10)
    optim_control_oq <- list(maxit=5e+2, reltol=1e-21, fnscale = 10)
  } else {
    optim_control_iq <- list(maxit=1e+4, reltol=1e-12, fnscale = 10)
    optim_control_oq <- list(maxit=2e+4, reltol=1e-12, fnscale = 10)
  }

  if(!is.null(optim_control_iq) & !is.null(optim_control_oq)){
    optim_control_iq <- optim_control_iq
    optim_control_oq <- optim_control_oq
  }

<<<<<<< HEAD
  if(criterion[1] == "ipoqlldif"){
    n_par <- sum(nrow(X),((1+ncol(groups_map)+(max(X,na.rm = TRUE)-min(X,na.rm = TRUE)))*ncol(X)))
  } else {
    n_par <- sum(nrow(X),((1+(max(X,na.rm = TRUE)-min(X,na.rm = TRUE)))*ncol(X)))
    # print(n_par)
  }
=======
  n_par <- sum(nrow(X),((1+(max(X,na.rm = TRUE)-min(X,na.rm = TRUE)))*ncol(X)))

>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5

  trace <- list()
  trace[["isLegacy"]] <- isLegacy

  # To handle mechanism of continuing unfinished search
  if(isContinued){

    load(prevData)
    scoreMat <- trace_data$scoreMat
    incl_set <- trace_data$traceStatus$current_set
    isLegacy <- trace_data$traceStatus$isLegacy
    i <- length(incl_set)

    if(trace_data$traceStatus$next_step == "forward"){

      if(is.null(cores)){
        cores <- 2
      } else {
        if(cores > parallel::detectCores()){
          cores <- parallel::detectCores()
        }
      }

      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl = cl, cores = cores)

      trace[["next_step"]] <- "forward"


<<<<<<< HEAD
      if(isLegacy & (i+1) < length(fullitem)){
=======
      if(isLegacy & (i+3) < length(fullitem)){
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
        init_par_iq <- c(na.omit(scoreMat[i,c((3+ncol(X)+1):(3+ncol(X)+n_par))]))
        init_par_oq <- c(na.omit(scoreMat[i,c((3+ncol(X)+n_par+1):(3+ncol(X)+n_par+(n_par-nrow(X))))]))
      } else {
        init_par_iq <- init_par_oq <- c()
      }

      if(isTracked){
        # cat(":::",format(Sys.time(),"%d/%m/%Y %H:%M:%S"),":::","do forward selection...")
        cat("do forward...")
        cat("\n")
      }

      # time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, itemsets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
      #                                                            init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
      #                                                            optim_control_oq = optim_control_oq, lambda_delta = lambda_delta, lambda_delta_out = lambda_delta_out,
      #                                                            random.init = random.init, random.init.th = random.init.th, eps = eps, lambda_out = lambda_out, opt.method = opt.method, max.iter = max.iter, abs.tol = abs.tol, scale.down = scale.down, psi.mult.up = psi.mult.up,                 checkNonZero = checkNonZero, step = step, maxit = maxit,max.diff.par = max.diff.par))

      time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, incl_sets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
                                                                 init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
                                                                 optim_control_oq = optim_control_oq, setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq))

      best.fwd <- res.fwd[order(res.fwd[,3],decreasing = TRUE),][1,]


      # Condition required to continue forward search
      while(((i+1) < (ncol(X)-1) & round(best.fwd[3],4) > round(scoreMat[i+1,3],4) &
             !identical(as.vector(na.omit(scoreMat[i+1,(4:(3+ncol(X)))])),(as.vector(na.omit(best.fwd[4:(3+ncol(X))]))))) |
            (is.na(scoreMat[i+1,3]))){

        i <- i + 1

        scoreMat[i,] <- best.fwd
        new.incl_set <- best.fwd[4:(3+i)]
        incl_set <- new.incl_set
        excl_set <- fullitem[-incl_set]

        if(isTracked){
          cat((i),":",paste(incl_set,collapse = ","))
          cat("\n")
        }

        trace[["current"]] <- "forward"
        trace[["current_set"]] <- incl_set
        trace_data <- list("scoreMat" = scoreMat, "traceStatus" = trace)
        save(trace_data, file = tempFile)
        write.csv(scoreMat, file = namecsv)


<<<<<<< HEAD
        if(isLegacy & (i+1) < length(fullitem)){
=======
        if(isLegacy & (i+3) < length(fullitem)){
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
          init_par_iq <- c(na.omit(scoreMat[i,c((3+ncol(X)+1):(3+ncol(X)+n_par))]))
          init_par_oq <- c(na.omit(scoreMat[i,c((3+ncol(X)+n_par+1):(3+ncol(X)+n_par+(n_par-nrow(X))))]))
        } else {
          init_par_iq <- init_par_oq <- c()
        }

        if(isTracked){
          # cat(":::",format(Sys.time(),"%d/%m/%Y %H:%M:%S"),":::","do forward selection...")
          cat("do forward...")
          cat("\n")
        }

        # time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, incl_sets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
        #                                                            init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
        #                                                            optim_control_oq = optim_control_oq, lambda_delta = lambda_delta, lambda_delta_out = lambda_delta_out,
        #                                                            random.init = random.init, random.init.th = random.init.th, eps = eps, lambda_out = lambda_out, opt.method = opt.method, max.iter = max.iter, abs.tol = abs.tol, scale.down = scale.down, psi.mult.up = psi.mult.up,                 checkNonZero = checkNonZero, step = step, maxit = maxit,max.diff.par = max.diff.par))

        time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, incl_sets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
                                                                   init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
                                                                   optim_control_oq = optim_control_oq, setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq))

        best.fwd <- res.fwd[order(res.fwd[,3],decreasing = TRUE),][1,]

        trace[["next_step"]] <- "forward"

      }

      trace[["next_step"]] <- "backward"
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()

    }


  } else {

    # If not continuing the unfinished search then it will start with the specified incl_set
    if(isTracked){
      # cat(":::",format(Sys.time(),"%d/%m/%Y %H:%M:%S"),":::","do forward selection...")
      cat("do full items estimation...")
      cat("\n")
    }
    # score <- compute_score(X = X, incl_set = incl_set, type = criterion , groups_map = groups_map, optim_control_iq = c(), optim_control_oq = c(),
    #                        lambda_delta = lambda_delta, lambda_delta_out = lambda_delta_out,
    #                        random.init = random.init, random.init.th = random.init.th, eps = eps, lambda_out = lambda_out,
    #                        opt.method = opt.method, max.iter = max.iter, abs.tol = abs.tol, scale.down = scale.down,
    #                        psi.mult.up = psi.mult.up, checkNonZero = checkNonZero, step = step, maxit = maxit,max.diff.par = max.diff.par)

    score <- compute_score(X = X, incl_set = incl_set, type = criterion , groups_map = groups_map,
                           optim_control_iq = optim_control_iq, optim_control_oq = optim_control_oq,
                           setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq)

    scoreMat <- matrix(NA,nrow = ncol(X), ncol = length(score))
    scoreMat[ncol(X),] <- score

    if(isTracked){
      cat((length(incl_set)),":",paste(incl_set,collapse = ","))
      cat("\n")
    }

  }

  excl_set <- fullitem[-c(incl_set)]
  # main iteration begins with the backward search
  i <- (length(incl_set)-1)

  # setting up the parallelization
  if(is.null(cores)){
    cores <- 2
  } else {
    if(cores > parallel::detectCores()){
      cores <- parallel::detectCores()
    }
  }
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl=cl, cores = cores)

  # iteration will stop if the number of items in the incl_set reaches zero
  while(i >= 1){

    #### Begin backward ####

<<<<<<< HEAD
    # if(isLegacy & (i+1) < length(fullitem)){
    if(isLegacy & ((i+1) < length(fullitem) & i > 1)){
=======
    if(isLegacy & ((i+1) < length(fullitem)) & ((i-1) > 1)){
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
      init_par_iq <- c(na.omit(scoreMat[i+1,c((3+ncol(X)+1):(3+ncol(X)+n_par))]))
      # init_par_oq <- c(na.omit(scoreMat[i+1,c((3+ncol(X)+n_par+1):(3+ncol(X)+n_par+(n_par-nrow(X))))]))
      init_par_oq <- c(na.omit(scoreMat[i+1,c((3+ncol(X)+n_par+1):(ncol(scoreMat)))]))
    } else {
      init_par_iq <- init_par_oq <- c()
    }

    # assign("incl_set_ori", incl_set, envir = .GlobalEnv)
    # assign("n_par", n_par, envir = .GlobalEnv)
    # assign("init_par_iq", init_par_iq, envir = .GlobalEnv)
    # assign("init_par_oq", init_par_oq, envir = .GlobalEnv)

    if(isTracked){
      # cat(":::",format(Sys.time(),"%d/%m/%Y %H:%M:%S"),":::","do backward elimination...")
      cat("do backward...")
      cat("\n")
    }

    # Compute all of the ipoq-ll possible for onestep backward elimination
    # time <- system.time(res.bck <- compute_scores_unparalleled(X = X, itemsets = incl_set, type = criterion, step_direct = c("backward"), groups_map = groups_map,
    #                                                            init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
    #                                                            optim_control_oq = optim_control_oq, lambda_delta = lambda_delta, lambda_delta_out = lambda_delta_out,
    #                                                            random.init = random.init, random.init.th = random.init.th, eps = eps, lambda_out = lambda_out, opt.method = opt.method, max.iter = max.iter, abs.tol = abs.tol, scale.down = scale.down, psi.mult.up = psi.mult.up,                 checkNonZero = checkNonZero, step = step, maxit = maxit,max.diff.par = max.diff.par))

    time <- system.time(res.bck <- compute_scores_unparalleled(X = X, incl_sets = incl_set, type = criterion, step_direct = c("backward"), groups_map = groups_map,
                                                               init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
                                                               optim_control_oq = optim_control_oq, setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq))

    best.bck <- res.bck[order(res.bck[,3],decreasing = TRUE),][1,] # keep the combination with the highest ipoq-ll
    new.incl_set <- best.bck[4:(3+i)] # update the incl_set by removing the item
    incl_set <- new.incl_set
    excl_set <- fullitem[-incl_set] # update the excl_set by adding the item that is removed from incl_set


    if(best.bck[3] > scoreMat[i,3] | is.na(scoreMat[i,3])){ ### the condition that is needed to keep the highest backward score

      # Stores the results to the matrix
      scoreMat[i,] <- best.bck

      trace[["current"]] <- "backward"
      trace[["current_set"]] <- incl_set
      trace_data <- list("scoreMat" = scoreMat, "traceStatus" = trace)
      save(trace_data, file = tempFile)
      write.csv(scoreMat, file = namecsv)

      # print the information that is being tracked
      if(isTracked){
        cat(i,":",paste(incl_set,collapse = ","))
        cat("\n")
      }

      # Do forward search only if there are more than two items in the excluded set.
      if(i < (ncol(X)-1)){

        #### Begin forward ####

        trace[["next_step"]] <- "forward"

<<<<<<< HEAD
        if(isLegacy & (i+2) < length(fullitem)){
=======
        if(isLegacy & ((i+1) < length(fullitem)) & ((i-1) > 1)){
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
          init_par_iq <- c(na.omit(scoreMat[i,c((3+ncol(X)+1):(3+ncol(X)+n_par))]))
          init_par_oq <- c(na.omit(scoreMat[i,c((3+ncol(X)+n_par+1):(3+ncol(X)+n_par+(n_par-nrow(X))))]))
        } else {
          init_par_iq <- init_par_oq <- c()
        }

        if(isTracked){
          cat("do forward...")
          cat("\n")
        }

        # time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, itemsets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
        #                                                            init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
        #                                                            optim_control_oq = optim_control_oq, lambda_delta = lambda_delta, lambda_delta_out = lambda_delta_out,
        #                                                            random.init = random.init, random.init.th = random.init.th, eps = eps, lambda_out = lambda_out, opt.method = opt.method, max.iter = max.iter, abs.tol = abs.tol, scale.down = scale.down, psi.mult.up = psi.mult.up,                 checkNonZero = checkNonZero, step = step, maxit = maxit,max.diff.par = max.diff.par))

        time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, incl_sets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
                                                                   init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
                                                                   optim_control_oq = optim_control_oq, setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq))

        best.fwd <- res.fwd[order(res.fwd[,3],decreasing = TRUE),][1,]


        # Do forward loop until can not find higher score anymore
        while(((i+1) < (ncol(X)-1) & round(best.fwd[3],4) > round(scoreMat[i+1,3],4) &
               !identical(as.vector(na.omit(scoreMat[i+1,(4:(3+ncol(X)))])),(as.vector(na.omit(best.fwd[4:(3+ncol(X))]))))) |
              (is.na(scoreMat[i+1,3]))){

          i <- i + 1

          scoreMat[i,] <- best.fwd
          new.incl_set <- best.fwd[4:(3+i)]
          incl_set <- new.incl_set
          excl_set <- fullitem[-incl_set]

          trace[["current"]] <- "forward"
          trace[["current_set"]] <- incl_set
          trace_data <- list("scoreMat" = scoreMat, "traceStatus" = trace)
          save(trace_data, file = tempFile)
          write.csv(scoreMat, file = namecsv)

          if(isTracked){
            cat((i),":",paste(incl_set,collapse = ","))
            cat("\n")
          }

<<<<<<< HEAD
          if(isLegacy & (i+2) < length(fullitem)){
=======
          if(isLegacy & ((i+1) < length(fullitem)) & ((i-1) > 1)){
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
            init_par_iq <- c(na.omit(scoreMat[i,c((3+ncol(X)+1):(3+ncol(X)+n_par))]))
            init_par_oq <- c(na.omit(scoreMat[i,c((3+ncol(X)+n_par+1):(3+ncol(X)+n_par+(n_par-nrow(X))))]))
          } else {
            init_par_iq <- init_par_oq <- c()
          }

          if(isTracked){
            cat("do forward..")
            cat("\n")
          }

          # time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, itemsets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
          #                                                            init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
          #                                                            optim_control_oq = optim_control_oq, lambda_delta = lambda_delta, lambda_delta_out = lambda_delta_out,
          #                                                            random.init = random.init, random.init.th = random.init.th, eps = eps, lambda_out = lambda_out, opt.method = opt.method, max.iter = max.iter, abs.tol = abs.tol, scale.down = scale.down, psi.mult.up = psi.mult.up,                 checkNonZero = checkNonZero, step = step, maxit = maxit,max.diff.par = max.diff.par))

          time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, incl_sets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
                                                                     init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
                                                                     optim_control_oq = optim_control_oq, setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq))

          best.fwd <- res.fwd[order(res.fwd[,3],decreasing = TRUE),][1,]

          trace[["next_step"]] <- "forward"
        }
      }
      #### End forward ####

    } else if(best.bck[3] < scoreMat[i,3] & !is.na(scoreMat[i,3])){
      incl_set <- scoreMat[i,4:(3+i)]
      trace[["current"]] <- "backward"
      trace[["current_set"]] <- incl_set
      excl_set <- fullitem[-incl_set]
    }
    trace[["next_step"]] <- "backward"
    i <- i - 1

  }

  #### End backward ####

  parallel::stopCluster(cl)
  foreach::registerDoSEQ()

  trace_data <- list("scoreMat" = scoreMat, "traceStatus" = trace)
  save(trace_data, file = tempFile)
  write.csv(scoreMat, file = namecsv)

  res_search <- scoreMat[,1:(3+ncol(X))]

  if(isTracked){
    cat("::: End of search :::")
    cat("\n")
  }


  if(isLegacy & !isConvert){
    if(isTracked){
      cat("::: Recompute score... :::")
      cat("\n")
    }
    res_search <- compute_scores(X, incl_sets = res_search[,4:(3+ncol(X))], type = criterion[1], step_direct = c("fixed"), groups_map = groups_map,
                                 cores = cores, init_par_iq = c(), init_par_oq = c(), optim_control_iq = optim_control_iq,
                                 optim_control_oq = optim_control_oq, setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq)
    res_search <- res_search[,1:(3+ncol(X))]
  }

<<<<<<< HEAD
  if(isConvert){
    if(criterion[1] == "ipoqll"){
      if(!is.null(groups_map)){
        if(isTracked){
          cat("::: Start converting.... :::")
          cat("\n")
        }
        res_search <- compute_scores(X, incl_sets = res_search[,4:(3+ncol(X))], type = criterion[1], step_direct = c("fixed"), groups_map = groups_map,
                                     cores = cores, init_par_iq = c(), init_par_oq = c(), optim_control_iq = optim_control_iq,
                                     optim_control_oq = optim_control_oq, setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq)
        res_search <- res_search[,1:(3+ncol(X))]
      } else {
        "Can not convert the score, groups_map is not provided."
      }
    } else {
      cat("The criterion is already IPOQ-LL-DIF.")
      cat("\n")
    }
  }

  if(criterion[1] == "ipoqlldif"){
    class(res_search) <- c(class(res_search),"search", "ipoqlldif", "autoRasch")
  } else {
    class(res_search) <- c(class(res_search),"search", "ipoqll", "autoRasch")
  }

=======
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
  # If saved in file
  if(fileOutput != FALSE & !is.null(fileOutput)){
    save(res_search, file = tempFile)
    namecsv <- paste(paste(strsplit(tempFile, "(\\.)")[[1]][1:(length(strsplit(tempFile, "(\\.)")[[1]])-1)],collapse = "."),".csv",sep = "")
    write.csv(res_search, file = namecsv)
  }

<<<<<<< HEAD
  # if("ipoqll" %in% class(res_search)){
  #   dimnames(res_search) <- list(c(1:ncol(res_search)),c("IQ-LL","OQ-LL","IPOQ-LL",paste("V",c(1:ncol(res_search)),sep = '')))
  # } else if("ipoqlldif" %in% class(res_search)){
  #   dimnames(res_search) <- list(c(1:ncol(res_search)),c("IQ-LL-DIF","OQ-LL-DIF","IPOQ-LL-DIF",paste("V",c(1:ncol(res_search)),sep = '')))
  # }

=======
  colnames(res_search) <- c("IQ-LL","OQ-LL","IPOQ-LL",paste("Item",c(1:nrow(res_search)),sep = ""))

  class(res_search) <- c("search", class(res_search))
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5

  return(res_search)

}

<<<<<<< HEAD
#' @examples
#' #search_res <- backward_search(short_poly_data,criterion = "ipoqll", incl_set = c(1:9), cores = 2)
#' #plot_search(search_res,type="l")
#'
=======

>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
#' @rdname search
#'
#' @import doParallel
#' @import foreach
#'
#' @export
<<<<<<< HEAD
backward_search <- function(X, criterion = c("ipoqll","ipoqlldif") , incl_set = c(), groups_map = c(), cores = NULL,
                            optim_control_iq = c(), optim_control_oq = c(), isTracked = TRUE, isContinued = FALSE,
                            prevData = c(), isLegacy = FALSE, fileOutput = FALSE, tempFile = "temp_backSearch.RData",
                            isConvert = FALSE, setting_par_iq = c(), setting_par_oq = c()){
# backward_search <- function(X, criterion = c("ipoqll","ipoqlldif") , incl_set = c(), groups_map = c(), cores = NULL,
#                             optim_control_iq = c(), optim_control_oq = c(), isTracked = TRUE, isContinued = FALSE,
#                             prevData = c(), isLegacy = FALSE, fileOutput = FALSE, tempFile = "temp_backSearch.RData",
#                             isConvert = FALSE, lambda_delta = 7, random.init = FALSE, random.init.th = 1e-2, eps = 0){
=======
backward_search <- function(X, criterion = c("ipoqll") , incl_set = c(), groups_map = c(), cores = 20,
                            optz_tuner_iq = list(maxit=5e+2, reltol=1e-21, fnscale = 10),
                            optz_tuner_oq = list(maxit=5e+2, reltol=1e-21, fnscale = 10), isTracked = TRUE,
                            isContinued = FALSE, prevData = c(), isLegacy = TRUE, fileOutput = TRUE,
                            tempFile = "temp_backSearch.RData", isConvert = FALSE){
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5


  namecsv <- paste(paste(strsplit(tempFile, "(\\.)")[[1]][1:(length(strsplit(tempFile, "(\\.)")[[1]])-1)],collapse = "."),".csv",sep = "")
  fullitem <- c(1:ncol(X))

  if(isLegacy){
    optim_control_iq <- list(maxit=5e+2, reltol=1e-21, fnscale = 10)
    optim_control_oq <- list(maxit=5e+2, reltol=1e-21, fnscale = 10)
  } else {
    optim_control_iq <- list(maxit=1e+4, reltol=1e-12, fnscale = 10)
    optim_control_oq <- list(maxit=2e+4, reltol=1e-12, fnscale = 10)
  }

  if(!is.null(optim_control_iq) & !is.null(optim_control_oq)){
    optim_control_iq <- optim_control_iq
    optim_control_oq <- optim_control_oq
  }

  n_par <- sum(nrow(X),((1+(max(X,na.rm = TRUE)-min(X,na.rm = TRUE)))*ncol(X)))

  trace <- list()
  trace[["isLegacy"]] <- isLegacy

  # To handle mechanism of continuing unfinished search
  if(isContinued){

    load(prevData)
    scoreMat <- trace_data$scoreMat
    incl_set <- trace_data$traceStatus$current_set
    isLegacy <- trace_data$traceStatus$isLegacy
    i <- length(incl_set)

  } else {

    # If not continuing the unfinished search then it will start with the specified incl_set
    if(isTracked){
      # cat(":::",format(Sys.time(),"%d/%m/%Y %H:%M:%S"),":::","do forward selection...")
      cat("do full items estimation...")
      cat("\n")
    }
    score <- compute_score(X = X, incl_set = incl_set, type = criterion , groups_map = groups_map,
                           optim_control_iq = optim_control_iq, optim_control_oq = optim_control_oq,
                           setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq)

    scoreMat <- matrix(NA,nrow = ncol(X), ncol = length(score))
    scoreMat[ncol(X),] <- score

    if(isTracked){
      cat((length(incl_set)),":",paste(incl_set,collapse = ","))
      cat("\n")
    }

  }

  excl_set <- fullitem[-c(incl_set)]
  # main iteration begins with the backward search
  i <- (length(incl_set)-1)

  # setting up the parallelization
  if(is.null(cores)){
    cores <- 2
  } else {
    if(cores > parallel::detectCores()){
      cores <- parallel::detectCores()
    }
  }
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl = cl, cores = cores)

  # iteration will stop if the number of items in the incl_set reaches zero
  while(i >= 1){

    #### Begin backward ####

    if(isLegacy & (i+1) < length(fullitem)){
      init_par_iq <- c(na.omit(scoreMat[i+1,c((3+ncol(X)+1):(3+ncol(X)+n_par))]))
      init_par_oq <- c(na.omit(scoreMat[i+1,c((3+ncol(X)+n_par+1):(3+ncol(X)+n_par+(n_par-nrow(X))))]))
    } else {
      init_par_iq <- init_par_oq <- c()
    }

    if(isTracked){
      # cat(":::",format(Sys.time(),"%d/%m/%Y %H:%M:%S"),":::","do backward elimination...")
      cat("do backward...")
      cat("\n")
    }

    # Compute all of the ipoq-ll possible for onestep backward elimination
    time <- system.time(res.bck <- compute_scores_unparalleled(X = X, incl_sets = incl_set, type = criterion, step_direct = c("backward"), groups_map = groups_map,
                                                               init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq,
                                                               optim_control_oq = optim_control_oq, setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq))


    best.bck <- res.bck[order(res.bck[,3],decreasing = TRUE),][1,] # keep the combination with the highest ipoq-ll
    new.incl_set <- best.bck[4:(3+i)] # update the incl_set by removing the item
    incl_set <- new.incl_set
    excl_set <- fullitem[-incl_set] # update the excl_set by adding the item that is removed from incl_set


    if(best.bck[3] > scoreMat[i,3] | is.na(scoreMat[i,3])){ ### the condition that is needed to keep the highest backward score

      # Stores the results to the matrix
      scoreMat[i,] <- best.bck

      trace[["current"]] <- "backward"
      trace[["current_set"]] <- incl_set
      trace_data <- list("scoreMat" = scoreMat, "traceStatus" = trace)
      save(trace_data, file = tempFile)
      write.csv(scoreMat, file = namecsv)

      # print the information that is being tracked
      if(isTracked){
        cat(i,":",paste(incl_set,collapse = ","))
        cat("\n")
      }

    } else if(best.bck[3] < scoreMat[i,3] & !is.na(scoreMat[i,3])){
      incl_set <- scoreMat[i,4:(3+i)]
      trace[["current"]] <- "backward"
      trace[["current_set"]] <- incl_set
      excl_set <- fullitem[-incl_set]
    }
    trace[["next_step"]] <- "backward"
    i <- i - 1

  }

  #### End backward ####

  parallel::stopCluster(cl)
  foreach::registerDoSEQ()

  trace_data <- list("scoreMat" = scoreMat, "traceStatus" = trace)
  save(trace_data, file = tempFile)
  write.csv(scoreMat, file = namecsv)

  res_search <- scoreMat[,1:(3+ncol(X))]

  if(isTracked){
    cat("::: End of search :::")
    cat("\n")
  }

  if(isLegacy & !isConvert){
    if(isTracked){
      cat("::: Recompute score... :::")
      cat("\n")
    }
    res_search <- compute_scores(X, incl_sets = res_search[,4:(3+ncol(X))], type = criterion[1], step_direct = c("fixed"), groups_map = groups_map,
                                 init_par_iq = init_par_iq, init_par_oq = init_par_oq, optim_control_iq = optim_control_iq, optim_control_oq = optim_control_oq, cores = NULL,
                                 setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq)
    res_search <- res_search[,1:(3+ncol(X))]
  }

<<<<<<< HEAD
  if(isConvert){
    if(criterion[1] == "ipoqll"){
      if(!is.null(groups_map)){
        if(isTracked){
          cat("::: Start converting.... :::")
          cat("\n")
        }
        res_search <- compute_scores(X, incl_sets = res_search[,4:(3+ncol(X))], type = c("ipoqlldif"), step_direct = c("fixed"), groups_map = groups_map,
                                     init_par_iq = c(), init_par_oq = c(), optim_control_iq = c(), optim_control_oq = c(), cores = NULL,
                                     setting_par_iq = setting_par_iq, setting_par_oq = setting_par_oq)
        res_search <- res_search[,1:(3+ncol(X))]
      } else {
        "Can not convert the score, groups_map is not provided."
      }
    } else {
      cat("The criterion is already IPOQ-LL-DIF.")
      cat("\n")
    }
  }

  if(criterion[1] == "ipoqlldif"){
    class(res_search) <- c(class(res_search),"search", "ipoqlldif", "autoRasch")
  } else {
    class(res_search) <- c(class(res_search),"search", "ipoqll", "autoRasch")
  }

=======
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
  # If saved in file
  if(fileOutput != FALSE & !is.null(fileOutput)){
    save(res_search, file = tempFile)
    namecsv <- paste(paste(strsplit(tempFile, "(\\.)")[[1]][1:(length(strsplit(tempFile, "(\\.)")[[1]])-1)],collapse = "."),".csv",sep = "")
    write.csv(res_search, file = namecsv)
  }

<<<<<<< HEAD
  # if("ipoqll" %in% class(res_search)){
  #   dimnames(res_search) <- list(c(1:ncol(res_search)),c("IQ-LL","OQ-LL","IPOQ-LL",paste("V",c(1:ncol(res_search)),sep = '')))
  # } else if("ipoqlldif" %in% class(res_search)){
  #   dimnames(res_search) <- list(c(1:ncol(res_search)),c("IQ-LL-DIF","OQ-LL-DIF","IPOQ-LL-DIF",paste("V",c(1:ncol(res_search)),sep = '')))
  # }


=======
  colnames(res_search) <- c("IQ-LL","OQ-LL","IPOQ-LL",paste("Item",c(1:nrow(res_search)),sep = ""))

  class(res_search) <- c("search",class(res_search))
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5

  return(res_search)

}

#' @param object The object of class \code{'search'}.
<<<<<<< HEAD
#' @param ... Further arguments to be passed.
=======
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
#'
#' @rdname search
#' @export
summary.search <- function(object, ...){
<<<<<<< HEAD
  obj <- object

  idx_max <- which(obj[,3] == max(obj[,3], na.rm = TRUE))
=======
  idx_max <- which(object[,3] == max(object[,3], na.rm = TRUE))
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5

  cat("\n")
  cat("Maximum IPOQ-LL score is obtained with ",idx_max," items in the included set.")
  cat("\n")
  cat("Items no.: ",paste(na.omit(object[idx_max,4:ncol(object)]),collapse = ","))
  cat("\n\n")
  print(matrix(object[idx_max,1:3],ncol = 3,byrow = TRUE,dimnames = list(c(""),c("IQ-LL","OQ-LL","IPOQ-LL"))), ... = ...)
  cat("\n\n")
}

<<<<<<< HEAD
#' @param x The object of class \code{'search'}.
#'
#' @rdname search
#' @export
print.search <- function(x,...){
  obj <- x
=======
print.search <- function(obj, ...){
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
  class(obj) <- c("matrix")
  print(obj, ... = ...)
}

<<<<<<< HEAD
#' @param obj An object of class "search".
#' @param remOrdered A logical statement whether show the order of the items removal or not.
#' @param locateMax A logical statement whether the location of the maximum score is needed to be marked or not.
=======
#' @param x The object of class \code{search}.
#' @param use.name Boolean value whether the plot will use variable name or not.
#' @param remOrdered Boolean value whether the plot will show the removal order of the items or not.
#' @param isMax Boolean value whether there is a line mark of the highest IPOQ-LL or not.
#' @param ... further argument passed or from other method.
#' @param type The type of the plot.
#' @param xlab The label of the x axis.
#' @param ylab The label of the y axis.
#' @param xlim The range value of the x axis plot.
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
#'
#' @rdname search
#'
#' @export
<<<<<<< HEAD
plot_search <- function(obj, remOrdered = TRUE, locateMax = TRUE, ...){
=======
plot.search <- function(x, use.name = FALSE, remOrdered = TRUE, isMax = TRUE, ylab = "IPOQ-LL", xlab = expression('|S'['in']*'|'),
                        xlim = c(), type = "l", ...){
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5

  obj <- x
  dotdotdot <- list(...)
<<<<<<< HEAD

  # if(!is.null(dotdotdot$xlab)){
  #   xlab <- dotdotdot$xlab
  # } else {
  #   xlab <- expression('|S'['in']*'|')
  # }
  #
  # if(!is.null(dotdotdot$ylab)){
  #   ylab <- dotdotdot$ylab
  # } else {
  #   if("ipoqll" %in% class(x)){
  #     ylab <- "IPOQ-LL"
  #   } else if("ipoqlldif" %in% class(x)){
  #     ylab <- "IPOQ-LL-DIF"
  #   }
  # }

  # if(!is.null(dotdotdot$remOrdered)){
  #   remOrdered <- dotdotdot$remOrdered
  # } else {
  #   remOrdered <- TRUE
  # }
  #
  # if(!is.null(dotdotdot$isMax)){
  #   isMax <- dotdotdot$isMax
  # } else {
  #   isMax <- TRUE
  # }
  #
  # if(!is.null(dotdotdot$type)){
  #   type <- dotdotdot$type
  # } else {
  #   type <- "l"
  # }

  # if(!is.null(dotdotdot$xlim)){
  #   xlim <- dotdotdot$xlim
  # } else {
  #   xlim <- c()
  # }

  # ygap <- (max(obj[,3],na.rm = TRUE)-min(obj[,3],na.rm = TRUE))
  #
  # if(!is.null(xlim)){
  #   if(xlim[2] >= xlim[1]){
  #     xlim <- xlim[c(2,1)]
  #   } else {
  #     xlim <- xlim
  #   }
  # } else {
  #   xlim <- c(nrow(obj),1)
  # }

  if((is.null(dotdotdot$xlab) |is.null(dotdotdot$ylab)) & (is.null(dotdotdot$xlim))){
    if("ipoqll" %in% class(obj)){
      ylab <- "IPOQ-LL"
    } else if("ipoqlldif" %in% class(obj)){
      ylab <- "IPOQ-LL-DIF"
    }
    suppressWarnings(plot(x = c(1:nrow(obj)), y = obj[,3], ylab = ylab, xlab = expression('|S'['in']*'|'), xlim = c(nrow(obj),1), ... = ...))
  } else if((is.null(dotdotdot$xlab) | is.null(dotdotdot$ylab))){
    if("ipoqll" %in% class(obj)){
      ylab <- "IPOQ-LL"
    } else if("ipoqlldif" %in% class(obj)){
      ylab <- "IPOQ-LL-DIF"
    }
    suppressWarnings(plot(x = c(1:nrow(obj)), y = obj[,3], ylab = ylab, xlab = expression('|S'['in']*'|'), ... = ...))
  } else if((is.null(dotdotdot$xlim))){
    suppressWarnings(plot(x = c(1:nrow(obj)), y = obj[,3], xlim = c(nrow(obj),1), ... = ...))
  }

=======

  ygap <- (max(obj[,3],na.rm = TRUE)-min(obj[,3],na.rm = TRUE))

  if(!is.null(xlim)){
    if(xlim[2] >= xlim[1]){
      xlim <- xlim[c(2,1)]
    } else {
      xlim <- xlim
    }
  } else {
    xlim <- c(nrow(obj),1)
  }

  suppressWarnings(plot(x = c(1:nrow(obj)), y = obj[,3], ylab = ylab, xlab = xlab, xlim = xlim, type = type, ... = ...))

>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
  if(remOrdered){
    for(i in (length(obj[,3])-1):1){
      prev.item <- obj[i+1,4:(length(obj[,3])+3)]
      next.item <- obj[i,4:(length(obj[,3])+3)]
      labels.rem.idx <- which(!(prev.item %in% next.item))
      labels.add.idx <- which(!(next.item %in% prev.item))
      labels.rem <- paste(prev.item[labels.rem.idx],collapse="\n")
      if(length(labels.rem.idx) > 1){
        labels.add <- paste("+",next.item[labels.add.idx],collapse="\n",sep="")
        labels <- paste(labels.rem,"\n",labels.add,sep="")
      } else {
        labels <- labels.rem
      }
<<<<<<< HEAD
      text(i,obj[i,3], pos = 1, labels = labels, cex = ifelse(!is.null(dotdotdot$cex),(0.8*dotdotdot$cex),(0.8)))
=======
      text(i,obj[i,3], pos = 1, labels = labels, ... = ...)
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
    }
  }

  if(locateMax){
    maxPos <- which(obj[,3] == max(obj[,3],na.rm = TRUE))
    lines(c(maxPos,maxPos),c(-100000,obj[maxPos,3]), col = 3, lty = 2, ... = ...)
  }

}
<<<<<<< HEAD
# plot.search <- function(obj, type = "l", xlab = NULL, ylab = NULL, xlim = c(), ylim = c(), use.name = FALSE, remOrdered = TRUE, isMax = TRUE, ...){
#
#   dotdotdot <- list(...)
#   if(!is.null(dotdotdot$main)){
#     par(mar = c(6.5, 7.5, 2.5, 1), oma = c(0, 0,0, 0))
#   } else {
#     par(mar = c(6.5, 7.5, 0.5, 1), oma = c(0, 0,0, 0))
#   }
#
#   if(!is.null(xlab)){
#     x.lab <- xlab
#   } else {
#     x.lab <- "Number of included items"
#   }
#
#   if(!is.null(ylab)){
#     y.lab <- ylab
#   } else {
#     if("ipoqll" %in% class(obj)){
#       y.lab <- "IPOQ-LL"
#     } else {
#       y.lab <- "IPOQ-LL-DIF"
#     }
#   }
#
#   ygap <- (max(obj[,3],na.rm = TRUE)-min(obj[,3],na.rm = TRUE))
#
#   if(!is.null(xlim)){
#     xlim <- xlim
#   } else {
#     xlim <- c(nrow(obj),1)
#   }
#
#   if(!is.null(ylim)){
#     ylim <- ylim
#   } else {
#     ylim <- c(min(obj[,3],na.rm = TRUE)-(ygap/12),max(obj[,3],na.rm = TRUE)+(ygap/12))
#   }
#
#   suppressWarnings(plot(x = c(1:nrow(obj)), y = obj[,3], xlim = xlim, ylim = ylim, type = type, xlab = x.lab, ylab = y.lab, mgp =c(5,2,0), ... = ...))
#
#   if(remOrdered){
#     remList <- list()
#     k <- 1
#     for(i in c((nrow(obj)-1):(1))){
#       j <- i+1
#       remList[[k]] <- c(obj[j,4:(j+3)][-c(which(obj[j,4:(j+3)] %in% obj[i,4:(i+3)]))])
#       k <- k+1
#     }
#
#     for(i in c((nrow(obj)-1):(1))){
#         if(i%%2==0){
#           text(x = i, y = (obj[i,3]+(ygap/20)), adj = c(0.5,0), labels = paste(remList[[nrow(obj)-i]],collapse = ",\n"), col = c(), ... = ...)
#         } else {
#           text(x = i, y = (obj[i,3]-(ygap/20)), adj = c(0.5,1), labels = paste(remList[[nrow(obj)-i]],collapse = ",\n"), col = c(), ... = ...)
#         }
#     }
=======

# plot.search <- function(obj, fileOutput = FALSE, use.name = FALSE, remOrdered = TRUE, isMax = TRUE, ...){
#
#   if(fileOutput){
#     cairo_pdf("test.pdf", width = 2, height = 1.43, bg = "transparent")
#     cex.axis <- 0.5
#     cex.lab <- 0.6
#     cex <- 0.4
#     mar = c(1.4, 1.45, 0.05, 0.05)
#     mgp = c(0.8,0.1,0)
#     par(mar = mar, oma = c(0, 0,0, 0))
#     line.x <- 0.5
#     line.y <- 0.7
#     tcl <- -0.2
#     lwd.line <- 0.5
#     lwd.axis <- 0.4
#     mgp.x <- c(0.5,-0.1,0)
#     mgp.y <- c(0.5,0.1,0)
#   } else {
#     cex.axis <- 1.3
#     cex.lab <- 1.5
#     cex <- 1
#     mar = c(4, 4.5, 0.5, 0.5)
#     mgp = c(1,0.1,0)
#     par(mar = mar, oma = c(0, 0,0, 0))
#     line.x <- 3
#     line.y <- 3
#     tcl <- -0.5
#     lwd.line <- 2
#     lwd.axis <- 2
#     mgp.x <- c(1,1,0)
#     mgp.y <- c(1,1,0)
#   }
#   is_title <- FALSE
#
#   diff_y <- max(obj[,3])-min(obj[,3])
#   ymax <- max(obj[,3]) + (0.1*diff_y)
#   ymin <- min(obj[,3]) - (0.2*diff_y)
#
#   plot(c(1:length(obj[,3])),obj[,3], main = ifelse(is_title,"Stepwise Search",""), xlim = c(length(obj[,3]),1), xlab = "", ylab = "",
#        cex.main = cex, cex = 0.1, cex.axis = cex.axis, cex.lab = cex.lab, mgp =mgp, lwd = lwd.line, tcl = tcl, xaxt = "n", yaxt = "n",
#        type = "l")
#   axis(1, lwd = lwd.axis, tcl = tcl,cex.axis = cex.axis, cex = cex, cex.lab = cex.lab, mgp = mgp.x)
#   axis(2, lwd = lwd.axis, tcl = tcl,cex.axis = cex.axis, cex = cex, cex.lab = cex.lab, mgp = mgp.y)
#   box(lwd=lwd.axis)
#
#   for(i in (length(obj[,3])-1):1){
#     prev.item <- obj[i+1,4:(length(obj[,3])+3)]
#     next.item <- obj[i,4:(length(obj[,3])+3)]
#     labels.rem.idx <- which(!(prev.item %in% next.item))
#     labels.add.idx <- which(!(next.item %in% prev.item))
#     labels.rem <- paste(prev.item[labels.rem.idx],collapse="\n")
#     if(length(labels.rem.idx) > 1){
#       labels.add <- paste("+",next.item[labels.add.idx],collapse="\n",sep="")
#       labels <- paste(labels.rem,"\n",labels.add,sep="")
#     } else {
#       labels <- labels.rem
#     }
#     text(i,obj[i,3], pos = 1, labels = labels, cex = cex)
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
#   }
#
#   if(isMax){
#     maxPos <- which(obj[,3] == max(obj[,3],na.rm = TRUE))
#     lines(c(maxPos,maxPos),c(-100000,obj[maxPos,3]), col = 3, lty = 2, ... = ...)
#   }
#
<<<<<<< HEAD
=======
#   title(xlab = "|Sin|", line = line.x, cex.lab = cex.lab)
#   title(ylab = "IPOQ-LL", line = line.y, cex.lab = cex.lab)
#
#   if(fileOutput){
#     dev.off()
#   }
>>>>>>> 7bda44cf6ff72132fa57077ea53e1ef9d6063ea5
# }
