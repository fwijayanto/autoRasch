#' Stepwise Selection Search
#'
#' To search itemset that give maximum value of the criterion
#'
#' @param X A matrix or data.frame of the observed responses (ordinal or binary response).
#' @param criterion The criterion that should be used, either ipoqll or ipoqlldif. The default is ipoqll.
#' @param incl_set A vector of initial items in the included set to start the search. The default is to start with full items.
#' @param groups_map A matrix or vector to map the subject to the DIFs groups.
#' @param cores An integer value of number of cores should be used for computation. The default is 2.
#' @param isLegacy A logical value whether using the previous step estimated parameters as the initial parameter values for the upcoming search step.
#' @param isTracked A logical value whether the progress need to be tracked or not.
#' @param isContinued A logical value whether this search is continuing another unfinished search.
#' @param prevData The filename of the temporary .RData file of the unfinished search.
#' @param fileOutput The filename if it is wished to save the output results in file (.RData and .csv) and FALSE if not.
#' @param tempFile The filename of the temporary file to track the search progress. The default is \code{"temp_stepSearch.RData"} which also automatically produces \code{"temp_stepSearch.csv"}.
#' @param optz_tuner_iq The optimisation setting of the included set. See \code{\link[stats:optim]{stats::optim()}} \code{control} parameter.
#' @param optz_tuner_oq The optimisation setting of the excluded set. See \code{\link[stats:optim]{stats::optim()}} \code{control} parameter.
#' @param isConvert A logical value whether it is wanted to recompute the score of the search results using IPOQ-LL-DIF criterion.
#'
#' @return
#' A matrix of the itemsets that obtain the highest scores for each number of items in the included set.
#'
#'
#' @details
#' To search the itemset that give the maximum score.
#'
#' @rdname search
#' @export
stepwise_search <- function(X, criterion = c("ipoqll","ipoqlldif") , incl_set = c(), groups_map = c(), cores = NULL,
                            optz_tuner_iq = c(), optz_tuner_oq = c(), isTracked = TRUE, isContinued = FALSE,
                            prevData = c(), isLegacy = FALSE, fileOutput = FALSE, tempFile = "temp_stepSearch.RData",
                            isConvert = FALSE){

  namecsv <- paste(paste(strsplit(tempFile, "(\\.)")[[1]][1:(length(strsplit(tempFile, "(\\.)")[[1]])-1)],collapse = "."),".csv",sep = "")
  fullitem <- c(1:ncol(X))

  if(is.null(incl_set)){
    incl_set <- c(1:ncol(X))
  }

  if(isLegacy){
    optz_tuner_iq <- list(maxit=5e+2, reltol=1e-21, fnscale = 10)
    optz_tuner_oq <- list(maxit=5e+2, reltol=1e-21, fnscale = 10)
  } else {
    optz_tuner_iq <- list(maxit=1e+4, reltol=1e-12, fnscale = 10)
    optz_tuner_oq <- list(maxit=2e+4, reltol=1e-12, fnscale = 10)
  }

  if(!is.null(optz_tuner_iq) & !is.null(optz_tuner_oq)){
    optz_tuner_iq <- optz_tuner_iq
    optz_tuner_oq <- optz_tuner_oq
  }

  if(criterion[1] == "ipoqlldif"){
    n_par <- sum(nrow(X),((1+ncol(groups_map)+(max(X,na.rm = TRUE)-min(X,na.rm = TRUE)))*ncol(X)))
  } else {
    n_par <- sum(nrow(X),((1+(max(X,na.rm = TRUE)-min(X,na.rm = TRUE)))*ncol(X)))
  }

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
      doParallel::registerDoParallel(cl, cores = cores)

      trace[["next_step"]] <- "forward"


      if(isLegacy & (i+1) != length(fullitem)){
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

      time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, itemsets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
                                                                 init_par_iq = init_par_iq, init_par_oq = init_par_oq, optz_tuner_iq = optz_tuner_iq,
                                                                 optz_tuner_oq = optz_tuner_oq))

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


        if(isLegacy & (i+1) != length(fullitem)){
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

        time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, itemsets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
                                                                   init_par_iq = init_par_iq, init_par_oq = init_par_oq, optz_tuner_iq = optz_tuner_iq,
                                                                   optz_tuner_oq = optz_tuner_oq))

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
    score <- compute_score(X = X, incl_set = incl_set, type = criterion , groups_map = groups_map, optz_tuner_iq = c(), optz_tuner_oq = c())

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
  doParallel::registerDoParallel(cl, cores = cores)

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
    time <- system.time(res.bck <- compute_scores_unparalleled(X = X, itemsets = incl_set, type = criterion, step_direct = c("backward"), groups_map = groups_map,
                                                               init_par_iq = init_par_iq, init_par_oq = init_par_oq, optz_tuner_iq = optz_tuner_iq,
                                                               optz_tuner_oq = optz_tuner_oq))


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
      if(i < (ncol(X)-2)){

        #### Begin forward ####

        trace[["next_step"]] <- "forward"

        if(isLegacy & (i+1) < length(fullitem)){
          init_par_iq <- c(na.omit(scoreMat[i,c((3+ncol(X)+1):(3+ncol(X)+n_par))]))
          init_par_oq <- c(na.omit(scoreMat[i,c((3+ncol(X)+n_par+1):(3+ncol(X)+n_par+(n_par-nrow(X))))]))
        } else {
          init_par_iq <- init_par_oq <- c()
        }

        if(isTracked){
          cat("do forward...")
          cat("\n")
        }

        time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, itemsets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
                                                                   init_par_iq = init_par_iq, init_par_oq = init_par_oq, optz_tuner_iq = optz_tuner_iq,
                                                                   optz_tuner_oq = optz_tuner_oq))

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

          if(isLegacy & (i+1) < length(fullitem)){
            init_par_iq <- c(na.omit(scoreMat[i,c((3+ncol(X)+1):(3+ncol(X)+n_par))]))
            init_par_oq <- c(na.omit(scoreMat[i,c((3+ncol(X)+n_par+1):(3+ncol(X)+n_par+(n_par-nrow(X))))]))
          } else {
            init_par_iq <- init_par_oq <- c()
          }

          if(isTracked){
            cat("do forward..")
            cat("\n")
          }

          time <- system.time(res.fwd <- compute_scores_unparalleled(X = X, itemsets = incl_set, type = criterion, step_direct = c("forward"), groups_map = groups_map,
                                                                     init_par_iq = init_par_iq, init_par_oq = init_par_oq, optz_tuner_iq = optz_tuner_iq,
                                                                     optz_tuner_oq = optz_tuner_oq))

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
    res_search <- compute_scores(X, itemsets = res_search[,4:(3+ncol(X))], type = criterion[1], step_direct = c("fixed"), groups_map = groups_map,
                                 init_par_iq = c(), init_par_oq = c(), optz_tuner_iq = c(), optz_tuner_oq = c(), cores = NULL)
    res_search <- res_search[,1:(3+ncol(X))]
  }

  if(isConvert){
    if(criterion[1] == "ipoqll"){
      if(!is.null(groups_map)){
        if(isTracked){
          cat("::: Start converting.... :::")
          cat("\n")
        }
        res_search <- compute_scores(X, itemsets = res_search[,4:(3+ncol(X))], type = c("ipoqlldif"), step_direct = c("fixed"), groups_map = groups_map,
                                   init_par_iq = c(), init_par_oq = c(), optz_tuner_iq = c(), optz_tuner_oq = c(), cores = NULL)
        res_search <- res_search[,1:(3+ncol(X))]
      } else {
        "Can not convert the score, groups_map is not provided."
      }
    } else {
      cat("The criterion is already IPOQ-LL-DIF.")
      cat("\n")
    }
  }

  # If saved in file
  if(fileOutput != FALSE & !is.null(fileOutput)){
    save(res_search, file = fileOutput)
    namecsv <- paste(paste(strsplit(fileOutput, "(\\.)")[[1]][1:(length(strsplit(fileOutput, "(\\.)")[[1]])-1)],collapse = "."),".csv",sep = "")
    write.csv(res_search, file = namecsv)
  }

  if((isConvert & !is.null(groups_map)) | criterion[1] == "ipoqlldif"){
    class(res_search) <- c("search", "ipoqlldif", "autoRasch")
  } else {
    class(res_search) <- c("search", "ipoqll", "autoRasch")
  }

  if("ipoqll" %in% class(obj)){
    dimnames(res_search) <- list(c(1:ncol(res_search)),c("IQ-LL","OQ-LL","IPOQ-LL",paste("V",c(1:ncol(res_search)),sep = '')))
  } else if("ipoqlldif" %in% class(obj)){
    dimnames(res_search) <- list(c(1:ncol(res_search)),c("IQ-LL-DIF","OQ-LL-DIF","IPOQ-LL-DIF",paste("V",c(1:ncol(res_search)),sep = '')))
  }
  return(res_search)

}

#' @examples
#' search_res <- backward_search(poly_inh_dset,criterion = "ipoqll", incl_set = c(1:18), cores = 2, isLegacy = TRUE)
#' @rdname search
#' @export
backward_search <- function(X, criterion = c("ipoqll","ipoqlldif") , incl_set = c(), groups_map = c(), cores = NULL,
                            optz_tuner_iq = c(), optz_tuner_oq = c(), isTracked = TRUE, isContinued = FALSE,
                            prevData = c(), isLegacy = FALSE, fileOutput = FALSE, tempFile = "temp_backSearch.RData",
                            isConvert = FALSE){


  namecsv <- paste(paste(strsplit(tempFile, "(\\.)")[[1]][1:(length(strsplit(tempFile, "(\\.)")[[1]])-1)],collapse = "."),".csv",sep = "")
  fullitem <- c(1:ncol(X))

  if(isLegacy){
    optz_tuner_iq <- list(maxit=5e+2, reltol=1e-21, fnscale = 10)
    optz_tuner_oq <- list(maxit=5e+2, reltol=1e-21, fnscale = 10)
  } else {
    optz_tuner_iq <- list(maxit=1e+4, reltol=1e-12, fnscale = 10)
    optz_tuner_oq <- list(maxit=2e+4, reltol=1e-12, fnscale = 10)
  }

  if(!is.null(optz_tuner_iq) & !is.null(optz_tuner_oq)){
    optz_tuner_iq <- optz_tuner_iq
    optz_tuner_oq <- optz_tuner_oq
  }

  if(criterion[1] == "ipoqlldif"){
    n_par <- sum(nrow(X),((1+ncol(groups_map)+(max(X,na.rm = TRUE)-min(X,na.rm = TRUE)))*ncol(X)))
  } else {
    n_par <- sum(nrow(X),((1+(max(X,na.rm = TRUE)-min(X,na.rm = TRUE)))*ncol(X)))
  }

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
    score <- compute_score(X = X, incl_set = incl_set, type = criterion , groups_map = groups_map, optz_tuner_iq = c(), optz_tuner_oq = c())

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
  doParallel::registerDoParallel(cl, cores = cores)

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
    time <- system.time(res.bck <- compute_scores_unparalleled(X = X, itemsets = incl_set, type = criterion, step_direct = c("backward"), groups_map = groups_map,
                                                               init_par_iq = init_par_iq, init_par_oq = init_par_oq, optz_tuner_iq = optz_tuner_iq,
                                                               optz_tuner_oq = optz_tuner_oq))


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
    res_search <- compute_scores(X, itemsets = res_search[,4:(3+ncol(X))], type = criterion[1], step_direct = c("fixed"), groups_map = groups_map,
                                 init_par_iq = c(), init_par_oq = c(), optz_tuner_iq = c(), optz_tuner_oq = c(), cores = NULL)
    res_search <- res_search[,1:(3+ncol(X))]
  }

  if(isConvert){
    if(criterion[1] == "ipoqll"){
      if(!is.null(groups_map)){
        if(isTracked){
          cat("::: Start converting.... :::")
          cat("\n")
        }
        res_search <- compute_scores(X, itemsets = res_search[,4:(3+ncol(X))], type = c("ipoqlldif"), step_direct = c("fixed"), groups_map = groups_map,
                                     init_par_iq = c(), init_par_oq = c(), optz_tuner_iq = c(), optz_tuner_oq = c(), cores = NULL)
        res_search <- res_search[,1:(3+ncol(X))]
      } else {
        "Can not convert the score, groups_map is not provided."
      }
    } else {
      cat("The criterion is already IPOQ-LL-DIF.")
      cat("\n")
    }
  }

  # If saved in file
  if(fileOutput != FALSE & !is.null(fileOutput)){
    save(res_search, file = fileOutput)
    namecsv <- paste(paste(strsplit(fileOutput, "(\\.)")[[1]][1:(length(strsplit(fileOutput, "(\\.)")[[1]])-1)],collapse = "."),".csv",sep = "")
    write.csv(res_search, file = namecsv)
  }

  if((isConvert & !is.null(groups_map)) | criterion[1] == "ipoqlldif"){
    class(res_search) <- c("search", "ipoqlldif", "autoRasch")
  } else {
    class(res_search) <- c("search", "ipoqll", "autoRasch")
  }

  if("ipoqll" %in% class(obj)){
    dimnames(res_search) <- list(c(1:ncol(res_search)),c("IQ-LL","OQ-LL","IPOQ-LL",paste("V",c(1:ncol(res_search)),sep = '')))
  } else if("ipoqlldif" %in% class(obj)){
    dimnames(res_search) <- list(c(1:ncol(res_search)),c("IQ-LL-DIF","OQ-LL-DIF","IPOQ-LL-DIF",paste("V",c(1:ncol(res_search)),sep = '')))
  }
  return(res_search)

}

#' @rdname search
#' @export
summary.search <- function(obj){
  idx_max <- which(obj[,3] == max(obj[,3], na.rm = TRUE))

  cat("\n")
  if("ipoqll" %in% class(obj)){
    cat("Maximum IPOQ-LL score is obtained with ",idx_max," items in the included set.")
  } else if("ipoqlldif" %in% class(obj)){
    cat("Maximum IPOQ-LL-DIF score is obtained with ",idx_max," items in the included set.")
  }
  cat("\n")
  cat("Items no.: ",paste(na.omit(obj[idx_max,4:ncol(obj)]),collapse = ","))
  cat("\n\n")
  if("ipoqll" %in% class(obj)){
    print(matrix(obj[idx_max,1:3],ncol = 3,byrow = TRUE,dimnames = list(c(""),c("IQ-LL","OQ-LL","IPOQ-LL"))))
  } else if("ipoqlldif" %in% class(obj)){
    print(matrix(obj[idx_max,1:3],ncol = 3,byrow = TRUE,dimnames = list(c(""),c("IQ-LL-DIF","OQ-LL-DIF","IPOQ-LL-DIF"))))
  }
  cat("\n\n")
}

print.search <- function(obj,...){
  class(obj) <- c("matrix")
  print(obj)
}

#' @rdname search
#' @export
plot.search <- function(obj, type = "l", xlab = NULL, ylab = NULL, xlim = c(), ylim = c(), use.name = FALSE, remOrdered = TRUE, isMax = TRUE, ...){

  dotdotdot <- list(...)
  if(!is.null(dotdotdot$main)){
    par(mar = c(6.5, 7.5, 2.5, 1), oma = c(0, 0,0, 0))
  } else {
    par(mar = c(6.5, 7.5, 0.5, 1), oma = c(0, 0,0, 0))
  }

  if(!is.null(xlab)){
    x.lab <- xlab
  } else {
    x.lab <- "Number of included items"
  }

  if(!is.null(ylab)){
    y.lab <- ylab
  } else {
    if("ipoqll" %in% class(obj)){
      y.lab <- "IPOQ-LL"
    } else {
      y.lab <- "IPOQ-LL-DIF"
    }
  }

  ygap <- (max(obj[,3],na.rm = TRUE)-min(obj[,3],na.rm = TRUE))

  if(!is.null(xlim)){
    xlim <- xlim
  } else {
    xlim <- c(nrow(obj),1)
  }

  if(!is.null(ylim)){
    ylim <- ylim
  } else {
    ylim <- c(min(obj[,3],na.rm = TRUE)-(ygap/12),max(obj[,3],na.rm = TRUE)+(ygap/12))
  }

  suppressWarnings(plot(x = c(1:nrow(obj)), y = obj[,3], xlim = xlim, ylim = ylim, type = type, xlab = x.lab, ylab = y.lab, mgp =c(5,2,0), ... = ...))

  if(remOrdered){
    remList <- list()
    k <- 1
    for(i in c((nrow(obj)-1):(1))){
      j <- i+1
      remList[[k]] <- c(obj[j,4:(j+3)][-c(which(obj[j,4:(j+3)] %in% obj[i,4:(i+3)]))])
      k <- k+1
    }

    for(i in c((nrow(obj)-1):(1))){
        if(i%%2==0){
          text(x = i, y = (obj[i,3]+(ygap/20)), adj = c(0.5,0), labels = paste(remList[[nrow(obj)-i]],collapse = ",\n"), col = c(), ... = ...)
        } else {
          text(x = i, y = (obj[i,3]-(ygap/20)), adj = c(0.5,1), labels = paste(remList[[nrow(obj)-i]],collapse = ",\n"), col = c(), ... = ...)
        }
    }
  }

  if(isMax){
    maxPos <- which(obj[,3] == max(obj[,3],na.rm = TRUE))
    lines(c(maxPos,maxPos),c(-100000,obj[maxPos,3]), col = 3, lty = 2, ... = ...)
  }

}
