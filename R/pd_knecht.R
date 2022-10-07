# -------------------------------------------------------------------------#
# knecht ----
# -------------------------------------------------------------------------#

#' Title
#'
#' @param FLAGjobDone 
#' @param FLAGjobPurged 
#' @param FLAGjobRecover 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family knecht
#'
#' @examples
knechtStatusMessage <- function(FLAGjobDone, FLAGjobPurged, FLAGjobRecover) {
  if (FLAGjobPurged)   return("Job is done and purged")
  if (FLAGjobDone)     return("Job is done")
  if (FLAGjobRecover)  return("Job is running")
  if (!FLAGjobRecover) return("Job is not yet started")
}

#' @export
knechts <- c(paste0("knecht",c(1,2,3,5,6)), paste0("ruprecht",c(1,2,3)))


#' Fit model on knecht
#'
#' @param pd 
#' @param .outputFolder 
#' @param n_startsPerNode 
#' @param n_nodes 
#' @param id 
#' @param type 
#'
#' @return Characters
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family knecht
#' @importFrom dMod distributed_computing
#'
#' @examples
pd_knecht_mstrust <- function(pd, .outputFolder, nStartsPerCore = 4, 
                              identifier = "mstrust", FLAGforcePurge = FALSE, 
                              opt.parameter_startpoints = "sample",
                              iterlim = 500,
                              machine) {
  
  curwd <- getwd()
  on.exit(setwd(curwd))
  
  # .. General job handling -----
  jobnm <- paste0("mstrust_", identifier, "_", gsub("-","_",gsub("(S\\d+(-\\d+)?).*", "\\1", basename(.outputFolder))))
  dir.create(jobnm, FALSE)
  setwd(jobnm)
  
  fileJobDone    <- dMod_files(.outputFolder, identifier)[["mstrust"]]
  fileJobPurged  <- file.path(dirname(fileJobDone), paste0(".", jobnm, "jobPurged"))
  fileJobRecover <- paste0(jobnm, "_1.R")
  
  FLAGjobDone    <- file.exists(fileJobDone)
  FLAGjobPurged  <- file.exists(fileJobPurged)
  FLAGjobRecover <- file.exists(fileJobRecover) | FLAGjobDone | FLAGjobPurged
  
  cat(knechtStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  
  # Assign Global variables: Important, in future, this might be a source of bugs, if other knecht-functions are written
  assign("nStartsPerCore",nStartsPerCore,.GlobalEnv)
  assign("opt.parameter_startpoints",opt.parameter_startpoints,.GlobalEnv)
  assign("iterlim",iterlim,.GlobalEnv)
  
  # Start mstrust job
  file.copy(file.path(pd$filenameParts$.currentFolder, pd$filenameParts$.compiledFolder, "/"), ".", recursive = TRUE)
  job <- dMod::runbg(
    {
      loadDLL(pd$obj_data);
      
      seed <- which(c(paste0("knecht",1:6), paste0("ruprecht",1:2)) == system("hostname", intern = TRUE))
      FLAGincludeCurrent <- TRUE
      
      cores <- dMod::detectFreeCores()
      
      if (identical(opt.parameter_startpoints, "sample")){
        center <- pe_sample_parameter_startpoints(pd$pe, n_starts = cores*nStartsPerCore, seed = seed, 
                                                    FLAGincludeCurrent = FLAGincludeCurrent)
      } else {
        center <- opt.parameter_startpoints
      }
      parlower <- petab_getParameterBoundaries(pd$pe, "lower")
      parupper <- petab_getParameterBoundaries(pd$pe, "upper")
      
      center <- center[,setdiff(names(center)     , names(pd$fixed))]
      parlower <- parlower[setdiff(names(parlower), names(pd$fixed))]
      parupper <- parupper[setdiff(names(parupper), names(pd$fixed))]
      
      
      mstrust(objfun = pd$obj, center = center, studyname = paste0("fit", seed),
              fixed = pd$fixed,
              rinit = 0.1, rmax = 10, cores = cores,
              iterlim = iterlim, 
              optmethod = "trust", 
              output = TRUE, cautiousMode = TRUE,
              stats = FALSE, 
              parlower = parlower, parupper = parupper)
    },
    filename = jobnm, 
    machine = machine, 
    recover = FLAGjobRecover,
    compile = FALSE, wait = FALSE
  )
  unlink(list.files(".", "\\.o$|\\.so$|\\.c$|\\.rds$"))
  
  if (!FLAGjobRecover) return("Job submitted")
  if (FLAGforcePurge) {
    # bit ugly code duplication...
    if (readline("Purge job. Are you sure? Type yes: ") == "yes"){
    job$purge()
    setwd(curwd)
    unlink(jobnm, T)
    return("Job was purged")}
  }
  # .. Get results -----
  if (!FLAGjobDone & !FLAGjobPurged) {
    if (job$check()) {
      Sys.sleep(0.1) # avoid being blocked
      job$get()
      f <- list.files(pattern = "_result.RData")
      f <- lapply(f, function(x) {
        load(x)
        .runbgOutput
      })
      f <- f[sapply(f, function(x) !inherits(x, "try-error"))]
      f <- do.call(c, f)
      try(cf_as.parframe(f))
      fits <- f
      setwd(curwd)
      dMod_saveMstrust(fit = fits, path = .outputFolder, 
                                             identifier = identifier, FLAGoverwrite = TRUE)
      savedFits <- readRDS(fileJobDone)
      return("Job done")
    }
  }
  
  FLAGjobDone    <- file.exists(fileJobDone)
  if (FLAGjobDone & !FLAGjobPurged) {
    if (readline("Purge job. Are you sure? Type yes: ") == "yes"){
      job$purge()
      setwd(curwd)
      unlink(jobnm, T)
      writeLines("jobPurged", fileJobPurged)
      return("Job purged\n")
    }
  }
}


#' Run profiles on knecht
#' 
#' Will profile the parameters stored in pd$pars
#'
#' @param pd 
#' @param .outputFolder 
#' @param FLAGforcePurge 
#' @param FLAGfixParsOnBoundary Fix parameters which went to the boundary. Don't fit and don't 
#'   use for reoptimization. A bit heuristic, but [dMod::profile()] does not support boundaries. 
#'   Therefore, if those parameters were not fixed, one potentially not start from the "optimum"
#' @param profpars I suggest to use either 1. names(pd$pars) or 2. hard code the result 
#'   from [pd_profile_getParsNotYetProfiled()]
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family knecht
#' @importFrom dMod profile_pars_per_node
#'
#'
#' @examples
pd_knecht_profile <- function(pd, .outputFolder, FLAGforcePurge = FALSE, FLAGfixParsOnBoundary = TRUE, 
                               profpars = names(pd$pars),
                               passwdEnv = Sys.getenv("hurensohn"), machine = "knecht") {
  # Fix pars which went to boundary
  if (FLAGfixParsOnBoundary){
    fixed_boundary <- pd_pars_getFixedOnBoundary(pd, tol = 1e-2)
    pd$fixed       <- c(pd$fixed, fixed_boundary)
    pd$pars        <- pd$pars[setdiff(names(pd$pars), names(pd$fixed))]
  }
  
  # .. Set up job -----
  cat("* use 5 digit identifier instead of 3\n")
  cat("* use fitrank in identifier")
  jobnm <- paste0("profile_", gsub("(S[0-9-]+-[0-9]+).*", "\\1", basename(.outputFolder)))
  
  var_list <- dMod::profile_pars_per_node(profpars, 16)
  
  fileJobDone   <- dMod_files(.outputFolder, profpars[1])$profile
  fileJobPurged <- file.path(dirname(dMod_files(.outputFolder)$profile), ".jobPurged")
  fileJobRecover <- file.path(paste0(jobnm, "_folder"), paste0(jobnm, ".R"))
  
  FLAGjobDone    <- file.exists(fileJobDone)
  FLAGjobPurged  <- file.exists(fileJobPurged)
  FLAGjobRecover <- file.exists(fileJobRecover) | FLAGjobDone | FLAGjobPurged
  
  cat(knechtStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  
  assign("profpars",profpars,.GlobalEnv)
  assign("var_list",var_list,.GlobalEnv)
  assign("jobnm",jobnm,.GlobalEnv)
  file.copy(file.path(pd$filenameParts$.currentFolder, pd$filenameParts$.compiledFolder, "/"), ".", recursive = TRUE)
  
  job <- distributed_computing(
    {
      loadDLL(pd$obj);
      whichPar <- profpars[(as.numeric(var_1):as.numeric(var_2))]
      cf_profile(pd$obj, pd$pars, whichPar, verbose = TRUE,
                 fixed = pd$fixed,
                 method = "optimize",
                 algoControl = list(gamma = 0.5, reoptimize = TRUE, correction = 0.5),
                 stepControl = list(limit = 100, min = log10(1.005), stepsize = log10(1.005)),
                 optControl = list(iterlim = 20),
                 cautiousMode = TRUE,
                 cores = 16,
                 path = file.path("~", paste0(jobnm, "_folder")))
    },
    jobname = jobnm, 
    partition = "single", cores = 16, nodes = 1, walltime = "12:00:00",
    ssh_passwd = passwdEnv, machine = machine, 
    var_values = var_list, no_rep = NULL, 
    recover = FLAGjobRecover,
    compile = F
  )
  unlink(list.files(".", "\\.o$|\\.so$|\\.c$|\\.rds$"))
  
  if (!FLAGjobRecover) return("Job submitted")
  if (FLAGforcePurge) {
    if (readline("Purge job. Are you sure? Type yes: ") == "yes"){
      # bit ugly code duplication...
      job$purge(purge_local = TRUE)
      return("Job was purged")
    }
  }
  
  # Get results
  if (!FLAGjobDone & !FLAGjobPurged) {
    if (job$check()) {
      job$get()
      # Copy profiles
      prof_files <- list.files(file.path(paste0(jobnm, "_folder"),  "results", "Results", "profile"),
                               full.names = TRUE)
      dir.create(dirname(dMod_files(.outputFolder)$profile), F,T)
      for (pf in prof_files) file.copy(pf, file.path(.outputFolder, "Results", "profile"))
      return("run readPd() again to load the results into the pd")
    }}
  
  # Purge job 
  prof_done <- list.files(file.path(.outputFolder, "Results", "profile"))
  prof_done <- gsub("^profiles-|.rds$","", prof_done)
  
  FLAGjobDone <- length(setdiff(profpars, prof_done)) == 0
  if (FLAGjobDone & !FLAGjobPurged) {
    if (readline("Purge job. Are you sure? Type yes: ") == "yes"){
      job$purge(purge_local = TRUE)
      writeLines("jobPurged", fileJobPurged)
      return("Job purged\n")
    }
  }
}




#' Fit model on knecht
#'
#' @param pd 
#' @param .outputFolder 
#' @param n_startsPerNode 
#' @param n_nodes 
#' @param id 
#' @param type 
#'
#' @return Characters
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family knecht
#' @importFrom dMod distributed_computing profile_pars_per_node
#'
#' @examples
pd_knecht_L1 <- function(pd, .outputFolder, n_nodes = 6, lambdas = 10^(seq(log10(0.0001), log10(100), length.out = n_nodes*16-1)), 
                          identifier = "L1", FLAGforcePurge = FALSE,
                          passwdEnv = Sys.getenv("hurensohn"), machine = "knecht") {
  # .. General job handling -----
  jobnm <- paste0("mstrust_", identifier, "_", gsub("(S\\d+).*", "\\1", basename(.outputFolder)))
  
  fileJobDone    <- dMod_files(.outputFolder, identifier)$L1
  fileJobPurged  <- file.path(dirname(fileJobDone), paste0(".", jobnm, "jobPurged"))
  fileJobRecover <- file.path(paste0(jobnm, "_folder"), paste0(jobnm, ".R"))
  
  FLAGjobDone    <- file.exists(fileJobDone)
  FLAGjobPurged  <- file.exists(fileJobPurged)
  FLAGjobRecover <- file.exists(fileJobRecover) | FLAGjobDone | FLAGjobPurged
  
  cat(knechtStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  
  # Assign Global variables: Important, in future, this might be a source of bugs, if other knecht-functions are written
  assign("lambdas",lambdas,.GlobalEnv)
  
  # Start mstrust job
  file.copy(file.path(pd$filenameParts$.currentFolder, pd$filenameParts$.compiledFolder, "/"), ".", recursive = TRUE)
  var_list <- dMod::profile_pars_per_node(lambdas, 16)
  
  job <- dMod::distributed_computing(
    {
      node <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
      
      loadDLL(pd$objfns$obj_data);
      
      parlower <- petab_getParameterBoundaries(pd$pe, "lower")
      parupper <- petab_getParameterBoundaries(pd$pe, "upper")
      
      ncores <- 16
      parallel::mclapply(X = var_1:var_2, mc.cores = ncores, FUN = function(idx) {
        lambda <- lambdas[idx]
        
        fit <- trustL1(
          objfun = pd$obj, parinit = pd$pars,
          mu = pd$L1$muL1, one.sided = FALSE, lambda = lambda,
          rinit = 0.1, rmax = 10, iterlim = 500,
          parupper = parupper, parlower = parlower
        )
        
        fit <- c(list(lambdaL1 = lambda), fit[c("value", "argument", "iterations", "converged")])
        dput(fit, file = sprintf("L1-%02i-%02i.R",node , idx))
        
        fit
      })
    },
    jobname = jobnm, 
    partition = "single", cores = 16, nodes = 1, walltime = "12:00:00",
    ssh_passwd = passwdEnv, machine = machine, 
    var_values = var_list, no_rep = NULL, 
    recover = FLAGjobRecover,
    compile = F
  )
  unlink(list.files(".", "\\.o$|\\.so$|\\.c$|\\.rds$"))
  
  if (!FLAGjobRecover) return("Job submitted")
  if (FLAGforcePurge) {
    # bit ugly code duplication...
    job$purge(purge_local = TRUE)
    return("Job was purged")
  }
  # .. Get results -----
  if (!FLAGjobDone & !FLAGjobPurged) {
    if (job$check()) {
      Sys.sleep(5) # avoid being blocked
      job$get()
      
      # Gebastelt ...
      fits <- list.files(file.path(paste0(jobnm, "_folder"),"results/"), "^L1.*R$", full.names = T) %>% lapply(source, local = TRUE) %>% lapply(function(x) x$value)
      fits <- lapply(fits, function(f) {data.table(as.data.table(f[setdiff(names(f), "argument")]), as.data.table(as.list(f$argument))) })
      fits <- rbindlist(fits)
      fits <- cf_parframe(fits, metanames = cf_parf_metaNames0$l1)
      dMod_saveL1(L1 = fits, path = .outputFolder, identifier = identifier, FLAGoverwrite = TRUE)
      
      return("Job done. # [ ] TODO You can check out the results by running `readPd` which will load the fit into pd$result$L1. Re-run this function once more to purge the job.")
    }
  }
  
  if (FLAGjobDone & !FLAGjobPurged) {
    if (readline("Purge job. Are you sure? Type yes: ") == "yes"){
      job$purge(purge_local = TRUE)
      writeLines("jobPurged", fileJobPurged)
      return("Job purged\n")
    }
  }
  
}

#' Fit model on knecht
#'
#' @param pd 
#' @param .outputFolder 
#' @param n_startsPerNode 
#' @param n_nodes 
#' @param id 
#' @param type 
#'
#' @return Characters
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family knecht
#' @importFrom dMod distributed_computing profile_pars_per_node
#'
#' @examples
pd_knecht_L1_mstrust <- function(pd, .outputFolder, 
                                  nstarts = 240,
                                  n_nodes = length(pd$pe$meta$L1$parameters_L1reference) * 1.5, 
                                  lambdas = 10^(seq(log10(0.001), log10(1000), length.out = n_nodes)), 
                                  identifier = "L1mstrust", FLAGforcePurge = FALSE,
                                  passwdEnv = Sys.getenv("hurensohn"), machine = "knecht") {
  # .. General job handling -----
  jobnm <- paste0("mstrust_", identifier, "_", gsub("(S\\d+).*", "\\1", basename(.outputFolder)))
  
  fileJobDone    <- dMod_files(.outputFolder, identifier)$L1
  fileJobPurged  <- file.path(dirname(fileJobDone), paste0(".", jobnm, "jobPurged"))
  fileJobRecover <- file.path(paste0(jobnm, "_folder"), paste0(jobnm, ".R"))
  
  FLAGjobDone    <- file.exists(fileJobDone)
  FLAGjobPurged  <- file.exists(fileJobPurged)
  FLAGjobRecover <- file.exists(fileJobRecover) | FLAGjobDone | FLAGjobPurged
  
  cat(knechtStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  
  # Assign Global variables: Important, in future, this might be a source of bugs, if other knecht-functions are written
  assign("lambdas",lambdas,.GlobalEnv)
  assign("nstarts",nstarts,.GlobalEnv)
  
  # Start mstrust job
  file.copy(file.path(pd$filenameParts$.currentFolder, pd$filenameParts$.compiledFolder, "/"), ".", recursive = TRUE)
  
  job <- dMod::distributed_computing(
    {
      node <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
      
      loadDLL(pd$objfns$obj_data);
      
      parlower <- petab_getParameterBoundaries(pd$pe, "lower")
      parupper <- petab_getParameterBoundaries(pd$pe, "upper")
      
      lambda <- lambdas[node]
      
      ncores <- 16
      
      startpars <- pepy_sample_parameter_startpoints(pd$pe, n_starts =  nstarts, FLAGincludeCurrent = TRUE)
      parallel::mclapply(X = seq_len(nstarts), mc.cores = ncores, FUN = function(idx) {
        
        pars <- as.parvec(startpars[idx])
        
        fit <- trustL1(
          objfun = pd$obj, parinit = pars,
          mu = pd$L1$muL1, one.sided = FALSE, lambda = lambda,
          rinit = 0.1, rmax = 10, iterlim = 500,
          parupper = parupper, parlower = parlower
        )
        
        fit <- c(list(lambdaL1 = lambda, index = idx), fit[c("value", "argument", "iterations", "converged")])
        dput(fit, file = sprintf("L1-%02i-%02i.R",node , idx))
        
        fit
      })
    },
    jobname = jobnm, 
    partition = "single", cores = 16, nodes = 1, walltime = "12:00:00",
    ssh_passwd = passwdEnv, machine = machine, 
    var_values = NULL, no_rep = n_nodes, 
    recover = FLAGjobRecover,
    compile = F
  )
  unlink(list.files(".", "\\.o$|\\.so$|\\.c$|\\.rds$"))
  
  if (!FLAGjobRecover) return("Job submitted")
  if (FLAGforcePurge) {
    # bit ugly code duplication...
    job$purge(purge_local = TRUE)
    return("Job was purged")
  }
  # .. Get results -----
  if (!FLAGjobDone & !FLAGjobPurged) {
    if (job$check()) {
      Sys.sleep(0.1) # avoid being blocked
      job$get()
      
      # Gebastelt ...
      fits <- list.files(file.path(paste0(jobnm, "_folder"),"results/"), "^L1.*R$", full.names = T) %>% lapply(source, local = TRUE) %>% lapply(function(x) x$value)
      fits <- lapply(fits, function(f) {data.table(as.data.table(f[setdiff(names(f), "argument")]), as.data.table(as.list(f$argument))) })
      fits <- rbindlist(fits)
      fits <- cf_parframe(fits, metanames = cf_parf_metaNames0$l1)
      dMod_saveL1(L1 = fits, path = .outputFolder, identifier = identifier, FLAGoverwrite = TRUE)
      
      return("Job done. # [ ] TODO You can check out the results by running `readPd` which will load the fit into pd$result$L1. Re-run this function once more to purge the job.")
    }
  }
  
  if (FLAGjobDone & !FLAGjobPurged) {
    if (readline("Purge job. Are you sure? Type yes: ") == "yes"){
      job$purge(purge_local = TRUE)
      writeLines("jobPurged", fileJobPurged)
      return("Job purged\n")
    }
  }
  
}


#' Fit model on knecht
#'
#' @param pd 
#' @param .outputFolder 
#' @param n_startsPerNode 
#' @param n_nodes 
#' @param id 
#' @param type 
#'
#' @return Characters
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family knecht
#' @importFrom dMod distributed_computing
#'
#' @examples
pd_knecht_L1_fitUnbiasedEachMstrust <- function(pd, .outputFolder, n_startsPerNode = 16*3, 
                                                 identifier = "L1UB", FLAGforcePurge = FALSE,
                                                 passwdEnv = Sys.getenv("hurensohn"), machine = "knecht") {
  
  # .. General job handling -----
  jobnm <- paste0("L1UB_", identifier, "_", gsub("(S\\d+).*", "\\1", basename(.outputFolder)))
  
  fileJobDone    <- dMod_files(.outputFolder, paste0(identifier,1))[["mstrust"]]
  fileJobPurged  <- file.path(dirname(fileJobDone), paste0(".", jobnm, "jobPurged"))
  fileJobRecover <- file.path(paste0(jobnm, "_folder"), paste0(jobnm, ".R"))
  
  FLAGjobDone    <- file.exists(fileJobDone)
  FLAGjobPurged  <- file.exists(fileJobPurged)
  FLAGjobRecover <- file.exists(fileJobRecover) | FLAGjobDone | FLAGjobPurged
  
  cat(knechtStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  
  # Hacking var_list: want to have different nodes
  n_nodes <- nrow(L1_getModelCandidates(pd$result$L1))
  var_list <- list(node = 1:n_nodes)
  
  # Assign Global variables: Important, in future, this might be a source of bugs, if other knecht-functions are written
  assign("n_startsPerNode",n_startsPerNode,.GlobalEnv)
  assign("n_nodes"        ,n_nodes        ,.GlobalEnv)
  assign("var_list"       ,var_list       ,.GlobalEnv)
  
  # Start mstrust job
  file.copy(file.path(pd$filenameParts$.currentFolder, pd$filenameParts$.compiledFolder, "/"), ".", recursive = TRUE)
  
  job <- dMod::distributed_computing(
    {
      loadDLL(pd$obj_data);
      
      node <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1
      # node <- var_1
      
      # Determine fixed pars and fix them
      fixed_L1 <- L1_getModelCandidates(pd$result$L1)
      
      parametersFixed <- fixed_L1[node,,drop = TRUE]
      parametersFixed <- names(parametersFixed)[parametersFixed]
      
      fit_par <- pd$pars[setdiff(names(pd$pars), parametersFixed)]
      fit_fix <- c(pd$fixed, pd$pars[parametersFixed] * 0)
      
      pd$pars  <- fit_par
      pd$fixed <- fit_fix
      
      # Sample only narrow
      pd$pe$parameters[,`:=`(initializationPriorParameters = paste0(pd$pars[parameterId] - 1,";", pd$pars[parameterId] + 1))]
      pd$pe$parameters[initializationPriorParameters == "NA;NA",`:=`(initializationPriorParameters = "0;1")] # Petab needs these to be set  but they don't matter
      
      # [ ] "current" from FLAGincludeCurrent should be the L1 fit rather than the global fit...
      center <- pepy_sample_parameter_startpoints(pd$pe, n_starts = n_startsPerNode, seed = node, 
                                                  FLAGincludeCurrent = TRUE)
      parlower <- petab_getParameterBoundaries(pd$pe, "lower")
      parupper <- petab_getParameterBoundaries(pd$pe, "upper")
      
      # only take free paramters
      center <- center[,setdiff(names(center)     , names(pd$fixed))] # redundant, is alredy taken care of by pepy_sample_parameter_startpoints. actually this implementation is not clean, as it does not take pd but pe as input, but uses pd
      parlower <- parlower[setdiff(names(parlower), names(pd$fixed))]
      parupper <- parupper[setdiff(names(parupper), names(pd$fixed))]
      
      fit <- mstrust(objfun = pd$obj, center = center, studyname = paste0("fit", node),
                     fixed = pd$fixed,
                     rinit = 0.1, rmax = 10, cores = 16,
                     iterlim = 500,
                     optmethod = "trust",
                     output = TRUE, cautiousMode = TRUE,
                     stats = FALSE,
                     parlower = parlower, parupper = parupper)
      parf <- try(cf_as.parframe(fit))
      parf <- parframe(cbind(L1modelCandidate = node , parf), parameters = attr(parf, "parameters"))
      parf
      
    },
    jobname = jobnm, 
    partition = "single", cores = 16, nodes = 1, walltime = "12:00:00",
    ssh_passwd = passwdEnv, machine = machine, 
    var_values = var_list, 
    recover = FLAGjobRecover,
    compile = F
  )
  unlink(list.files(".", "\\.o$|\\.so$|\\.c$|\\.rds$"))
  
  if (!FLAGjobRecover) return("Job submitted")
  if (FLAGforcePurge) {
    # bit ugly code duplication...
    job$purge(purge_local = TRUE)
    return("Job was purged")
  }
  # .. Get results -----
  if (!FLAGjobDone & !FLAGjobPurged) {
    if (job$check()) {
      Sys.sleep(5) # avoid being blocked
      job$get()
      parfs <- knecht_result
      lapply(seq_along(parfs), function(idx) dMod_saveMstrust(parfs[[idx]], .outputFolder, paste0(identifier, idx)))
      return("Job done. You can check out the results by running `readPd` which will load the fit into pd$result$fits. Re-run this function once more to purge the job.")
    }
  }
  
  if (FLAGjobDone & !FLAGjobPurged) {
    if (readline("Purge job. Are you sure? Type yes: ") == "yes"){
      job$purge(purge_local = TRUE)
      writeLines("jobPurged", fileJobPurged)
      return("Job purged\n")
    }
  }
}


#' Fit model on knecht
#'
#' @param pd 
#' @param .outputFolder 
#' @param n_startsPerNode 
#' @param n_nodes 
#' @param id 
#' @param type 
#'
#' @return Characters
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family knecht
#' @importFrom dMod distributed_computing
#'
#' @examples
pd_knecht_L1_fitUnbiasedEachOnce <- function(pd, .outputFolder, n_startsPerNode = 16*3, 
                                              identifier = "L1UBSingle", FLAGforcePurge = FALSE,
                                              passwdEnv = Sys.getenv("hurensohn"), machine = "knecht") {
  
  # .. General job handling -----
  jobnm <- paste0("L1UB_", identifier, "_", gsub("(S\\d+).*", "\\1", basename(.outputFolder)))
  
  fileJobDone    <- dMod_files(.outputFolder, identifier)[["mstrust"]]
  fileJobPurged  <- file.path(dirname(fileJobDone), paste0(".", jobnm, "jobPurged"))
  fileJobRecover <- file.path(paste0(jobnm, "_folder"), paste0(jobnm, ".R"))
  
  FLAGjobDone    <- file.exists(fileJobDone)
  FLAGjobPurged  <- file.exists(fileJobPurged)
  FLAGjobRecover <- file.exists(fileJobRecover) | FLAGjobDone | FLAGjobPurged
  
  cat(knechtStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  
  # Hacking var_list: want to have different nodes
  n_nodes <- nrow(L1_getModelCandidates(pd$result$L1))
  var_list <- dMod::profile_pars_per_node(1:n_nodes, 16)
  
  # Assign Global variables: Important, in future, this might be a source of bugs, if other knecht-functions are written
  assign("n_startsPerNode",n_startsPerNode,.GlobalEnv)
  assign("n_nodes"        ,n_nodes        ,.GlobalEnv)
  assign("var_list"       ,var_list       ,.GlobalEnv)
  
  # Start mstrust job
  file.copy(file.path(pd$filenameParts$.currentFolder, pd$filenameParts$.compiledFolder, "/"), ".", recursive = TRUE)
  
  
  stop("update arguments in trust (are mstrust args, not yet fixed)")
  
  job <- dMod::distributed_computing(
    {
      .pd <- copy(pd)
      loadDLL(pd$obj_data);
      
      # node <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1
      parallel::mclapply(var_1:var_2, function(node) {
        pd <- copy(.pd)
        # Determine fixed pars and fix them
        fixed_L1 <- L1_getModelCandidates(pd$result$L1)
        
        parametersFixed <- fixed_L1[node,,drop = TRUE]
        parametersFixed <- names(parametersFixed)[parametersFixed]
        
        fit_par <- pd$pars[setdiff(names(pd$pars), parametersFixed)]
        fit_fix <- c(pd$fixed, pd$pars[parametersFixed] * 0)
        
        pd$pars <- fit_par
        pd$fixed <- fit_fix
        
        # Sample only narrow
        pd$pe$parameters[,`:=`(initializationPriorParameters = paste0(pd$pars[parameterId] - 1,";", pd$pars[parameterId] + 1))]
        pd$pe$parameters[initializationPriorParameters == "NA;NA",`:=`(initializationPriorParameters = "0;1")] # Petab needs these to be set  but they don't matter
        
        parlower <- petab_getParameterBoundaries(pd$pe, "lower")
        parupper <- petab_getParameterBoundaries(pd$pe, "upper")
        
        # only take free paramters
        parlower <- parlower[setdiff(names(parlower), names(pd$fixed))]
        parupper <- parupper[setdiff(names(parupper), names(pd$fixed))]
        
        
        # trust(objfun = pd$obj, center = pd$pars, studyname = paste0("fit", node),
        #       fixed = pd$fixed,
        #       rinit = 0.1, rmax = 10, cores = 16,
        #       iterlim = 500, 
        #       optmethod = "trust", 
        #       output = TRUE, cautiousMode = TRUE,
        #       stats = FALSE, 
        #       parlower = parlower, parupper = parupper)
      })
    },
    jobname = jobnm, 
    partition = "single", cores = 16, nodes = 1, walltime = "12:00:00",
    ssh_passwd = passwdEnv, machine = machine, 
    var_values = var_list, 
    recover = FLAGjobRecover,
    compile = F
  )
  unlink(list.files(".", "\\.o$|\\.so$|\\.c$|\\.rds$"))
  
  if (!FLAGjobRecover) return("Job submitted")
  if (FLAGforcePurge) {
    # bit ugly code duplication...
    job$purge(purge_local = TRUE)
    return("Job was purged")
  }
  # .. Get results -----
  if (!FLAGjobDone & !FLAGjobPurged) {
    if (job$check()) {
      Sys.sleep(5) # avoid being blocked
      job$get()
      fitlist  <- if (exists("knecht_result")) do.call(c, knecht_result) else {NULL
        #cf_dMod_rescueFits()
        #   fitlist <- list.files(file.path(paste0(jobnm, "_folder"), "results","fit"), "\\.Rda$", recursive = TRUE, full.names = TRUE)
        #   fitlist <- lapply(fitlist, function(x) try(local(load(x))))
      }
      fits <- fitlist
      fits <- fits[vapply(fits, is.list, TRUE)]
      class(fits) <- "parlist"
      fits <- cf_as.parframe(fits)
      dMod_saveMstrust(fit = fits, path = .outputFolder, 
                                             identifier = identifier, FLAGoverwrite = TRUE)
      
      return("Job done. You can check out the results by running `readPd` which will load the fit into pd$result$fits. Re-run this function once more to purge the job.")
    }
  }
  
  if (FLAGjobDone & !FLAGjobPurged) {
    if (readline("Purge job. Are you sure? Type yes: ") == "yes"){
      job$purge(purge_local = TRUE)
      writeLines("jobPurged", fileJobPurged)
      return("Job purged\n")
    }
  }
}






