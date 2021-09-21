# -------------------------------------------------------------------------#
# pd-"Class" ----
# -------------------------------------------------------------------------#
# [ ] pd should be its own "class" with clearly defined names
# pd <- function()


# -------------------------------------------------------------------------#
# Files ----
# -------------------------------------------------------------------------#

#' Read a pd and load dlls
#'
#' @param filename
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Cluster
#' @family L1
#' @family pd helpers
#' @family parframe handling
#' @importFrom dMod loadDLL
#' @importFrom conveniencefunctions dMod_readProfiles dMod_readMstrust
#' @examples
readPd <- function(filename) {
  # 0 wd wrangling
  wd <- getwd()
  on.exit({setwd(wd)})
  # 1 Read RDS
  pd <- readRDS(filename)
  # Backwards compatibility
  if (is.null(pd$objfns))
    pd$objfns$obj_data <- pd$obj_data
  # 2 Load DLLs
  setwd(dirname(filename))
  dMod::loadDLL(pd$objfns$obj_data)
  setwd(wd)
  
  # 3 Read potential results if not yet in pd
  # If in pd already, some filenameParts variable might have been derived from them,
  #   so it would be dangerous to reload them if they were postprocessed
  
  # Fits
  path <- file.path(dirname(dirname(filename)), "Results", "mstrust")
  files <- list.files(path, "mstrustList.*\\.rds", full.names = TRUE)
  for (f in files) {
    identifier <- gsub("mstrustList-|\\.rds", "", basename(f))
    if (is.null(pd$result[[identifier]])) pd$result[[identifier]] <- readRDS(f)
  }
  
  # Profiles
  path <- dirname(dirname(filename))
  if (is.null(pd$result$profile) && dir.exists(dirname(conveniencefunctions::dMod_files(path)$profile)))
    pd$result$profiles <- conveniencefunctions::dMod_readProfiles(path)
  
  # L1 - get L1 scan results
  path <- dirname(dirname(filename))
  if (is.null(pd$result$L1) && dir.exists(dirname(conveniencefunctions::dMod_files(path)$L1)))
    pd$result$L1 <- conveniencefunctions::dMod_readL1(path)
  # L1 - get unbiased models
  
  
  # 4 Set Parameters to most relevant fit: mstrust || base_obsParsFitted || base
  pars <- if (!is.null(pd$result$mstrust)) as.parvec(pd$result$mstrust) else if (!is.null(pd$result$base_obsParsFitted)) as.parvec(pd$result$base_obsParsFitted) else pd$pars
  pd$pars <- unclass_parvec(pars)
  
  pd
}



#' Title
#'
#' @param pd
#'
#' @return
#' @export
#'
#' @examples
writePd <- function(pd) {
  rdsfile <- pd_files(pd$filenameParts)$rdsfile
  saveRDS(pd, rdsfile)
}


#' Title
#'
#' @param modelname
#' @param .compiledFolder
#' @param type
#'
#' @return
#' @export
#'
#' @examples
pd_files <- function(filenameParts) {
  if (is.null(filenameParts$type)) filenameParts$type <- "indiv"
  list(rdsfile = file.path(filenameParts$.currentFolder,
                           filenameParts$.compiledFolder,
                           paste0(filenameParts$modelname, "_", filenameParts$type, ".rds")),
       obsParsFitted = file.path(filenameParts$.currentFolder,
                                 filenameParts$.compiledFolder,
                                 paste0(".obsparsfitted")),
       traceFile = file.path(filenameParts$.currentFolder,
                             filenameParts$.compiledFolder,
                             paste0("traceFile.txt"))
  )
}

#' Title
#'
#' @param pd 
#'
#' @return
#' @export
#'
#' @examples
pd_guessPetabYaml <- function(pd) {
  if (!is.null(pd$filenameParts$petabYaml)) {
    pd$filenameParts$petabYaml
  } else { 
    wup <- list.files(file.path(pd$filenameParts$.currentFolder, pd$filenameParts$.compiledFolder, ".."), "yaml", recursive = T, full.names = TRUE) # guessing
    grep("report", wup, value = TRUE, invert = TRUE)
  }
}



# -------------------------------------------------------------------------#
# Manipulating fns ----
# -------------------------------------------------------------------------#

#' Quickly set "controls" amd update all high-level functions
#'
#' @param pd
#' @param rtol
#' @param atol
#' @param maxsteps
#' @param objtimes
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom dMod "controls<-"
#'
#' @examples
#' optionsOde <- list(method = "lsoda", rtol = rtol, atol = atol, maxsteps = maxsteps)
#' optionsSens <- list(method = "lsodes", rtol = rtol, atol = atol, maxsteps = maxsteps)
pdIndiv_updateControls <- function(pd,
                                   optionsOde = NULL,
                                   optionsSens = NULL,
                                   objtimes = NULL) {
  
  # Set integrator controls
  if (!is.null(optionsOde))  dMod::controls(pd$dModAtoms$fns$x, name = "optionsOde") <- optionsOde
  if (!is.null(optionsSens)) dMod::controls(pd$dModAtoms$fns$x, name = "optionsSens") <- optionsSens
  if (!is.null(objtimes))    dMod::controls(pd$objfns$obj_data, "times") <- objtimes
  
  pdIndiv_rebuildPrdObj(pd)
  
}

#' Rebuild the high-level prediction and objective function
#'
#' @param pd pd_indiv
#'
#' @return pd with updated p, prd and obj_data
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom dMod P_indiv PRD_indiv normL2_indiv objtimes controls
#'
#' @examples
pdIndiv_rebuildPrdObj <- function(pd, Nobjtimes = 100) {
  
  # Rebuild p
  p <- dMod::P_indiv((pd$dModAtoms$fns$p1 * pd$dModAtoms$fns$p0), pd$dModAtoms$gridlist$est.grid, pd$dModAtoms$gridlist$fix.grid)
  
  # Rebuild high-level prediction function
  prd0 <- Reduce("*", pd$dModAtoms$fns)
  prd <- PRD_indiv(prd0, pd$dModAtoms$gridlist$est.grid, pd$dModAtoms$gridlist$fix.grid)
  
  # Rebuild obj_data
  times <- c(pd$pe$measurementData$time, pd$pe$meta$presimTimes)
  tobj <- if (!is.null(pd$objfns$obj_data)) dMod::controls(pd$objfns$obj_data, name = "times") else dMod::objtimes(times, Nobjtimes = Nobjtimes)
  obj_data <- normL2_indiv(pd$dModAtoms$data, prd0,
                           pd$dModAtoms$e,
                           est.grid = pd$dModAtoms$gridlist$est.grid,
                           fix.grid = pd$dModAtoms$gridlist$fix.grid,
                           times = tobj)
  obj_prior <- petab_createObjPrior(pd$pe)
  
  # Update p, prd and obj_data obj_prior
  pd$p        <- p
  pd$prd      <- prd
  pd$objfns$obj_data <- obj_data
  pd$objfns$obj_prior <- obj_prior
  
  # Rebuild obj
  pd$obj <- Reduce("+", pd$objfns)
  
  pd
}


#' Title
#'
#' @param pd 
#' @param parameters 
#'
#' @return
#' @export
#'
#' @examples
#' parameters <- c("k_233", "k_234", "k_244")
pd_addObjPrior <- function(pd, parameters, FLAGuseNominalValueAsCenter) {
  
  cat("Todo: Implement this properly\n*write about prior in logfile\n*think about rebuildPrdObj\n*Use obj as default objective function from beginning")
  
  pe_aux <- copy(pd$pe)
  pe_aux$parameters <-pe_aux$parameters[parameterId %in% parameters]
  pd$objfns$obj_prior <- petab_createObjPrior(pe_aux, FLAGuseNominalValueAsCenter = FLAGuseNominalValueAsCenter)
  pd$obj <- Reduce("+", pd$objfns)
  pd
}



# -------------------------------------------------------------------------#
# Helpers ----
# -------------------------------------------------------------------------#

#' Wrapper around predtimes
#'
#' @param pd
#' @param N
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family pd helpers
#' @importFrom dMod predtimes
#'
#' @examples
pd_predtimes <- function(pd, N = 100) {
  datatimes <- unique(sort(pd$pe$measurementData$time))
  eventtimes <- NULL
  dMod::predtimes(datatimes,eventtimes,N)
}

#' Wrapper around predtimes
#'
#' @param pd
#' @param N
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family pd helpers
#' @importFrom dMod objtimes
#'
#' @examples
pd_objtimes <- function(pd, N = 100) {
  datatimes <- unique(sort(pd$pe$measurementData$time))
  eventtimes <- NULL
  dMod::objtimes(datatimes,eventtimes,N)
}



#' Determine which pars were not yet profiled
#'
#' @param pd 
#' @param .outputFolder 
#' @param FLAGreturnVector 
#'
#' @return either vector of p^arameters or output printed to console
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom conveniencefunctions dMod_files
#'
#' @examples
pd_profile_getParsNotYetProfiled <- function(pd, .outputFolder, FLAGreturnVector = FALSE) {
  prof_available <- list.files(dirname(conveniencefunctions::dMod_files(.outputFolder, identifier = "1")$profile))
  prof_available <- gsub("^profiles-","" ,prof_available)
  prof_available <- gsub(".rds$"     ,"" ,prof_available)
  
  profpars <- names(pd$pars)
  profpars <- setdiff(profpars, prof_available)
  if (FLAGreturnVector) {
    return(profpars)
  } else {
    cat("profpars <- "); dput(profpars)
    invisible()
  }
}


# -------------------------------------------------------------------------#
# Parameter handling ----
# -------------------------------------------------------------------------#

#' Update parameters
#'
#' @param pd
#' @param parsEst
#' @param FLAGupdatePE
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family parframe handling
#' @importFrom conveniencefunctions pars2parframe
#'
#' @examples
pd_updateEstPars <- function(pd, parsEst, FLAGupdatePE = TRUE, FLAGsavePd = FALSE) {
  # [ ] don't overwrite, just load from Results
  if(is.parframe(parsEst)) {
    pd$result$mstrust <- parsEst
    parsEst <- as.parvec(parsEst)
    pd$pars[names(parsEst)] <- parsEst
  } else {
    pd$pars[names(parsEst)] <- parsEst
    pd$result$base <- conveniencefunctions::pars2parframe(pd$pars, parameterSetId = "Base", obj = pd$obj)
    }

  if (FLAGupdatePE) {
    cat("pd$pe pars have been *set*")
    petab_setPars_estScale(pd$pe, parsEst)}
  if (FLAGsavePd) writePd(pd)
  pd
}



#' Title
#'
#' @param conveniencefunctions
#' @param pars2parframe
#' @param pd
#'
#' @return [dMod::parframe()]
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family parframe handling
#' @importFrom conveniencefunctions pars2parframe
#'
#' @examples
pd_parf_collectPars <- function(pd, parameterSetId = "base_obsParsFitted") {
  parf <- NULL
  if (!is.null(pd$result$base_obsParsFitted) && parameterSetId == "base_obsParsFitted"){
    parf <- cbind(parameterSetId = "base_obsParsFitted", pd$result$base_obsParsFitted, stringsAsFactors = FALSE)
  } else if (!is.null(pd$result$base) && parameterSetId == "base") {
    parf <- cbind(parameterSetId = "base", pd$result$base, stringsAsFactors = FALSE)
  } else {
    pd <- pd_updateEstPars(pd, pd$pars, FLAGupdatePE = FALSE)
    parf <- pd$result$base
  }
  
  parf <- dMod::parframe(parf,parameters = names(c(pd$pars, pd$fixed)), metanames = c("parameterSetId","fitrank","step","stepsize","index","value","converged","iterations"))
  parf
}


#' Title
#'
#' @param pd
#' @param fitrankRange
#'
#' @return [dMod::parframe()]
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family parframe handling
#' @importFrom conveniencefunctions cf_parf_getStepRepresentatives
#' @importFrom data.table as.data.table
#' @importFrom dMod parframe
#'
#' @examples
pd_parf_collectMstrust <- function(pd, fitrankRange = 1:15, tol = 1) {
  parf0 <- pd$result$mstrust[fitrankRange]
  fitidxs = conveniencefunctions::cf_parf_getStepRepresentatives(parf0, tol = tol)
  pars <- parf0[fitidxs]
  pars <- data.table::as.data.table(pars)
  pars[,`:=`(parameterSetId = paste0("step", step, ",", "rank", fitrank))]
  pars <- dMod::parframe(pars, parameters = names(pd$pars), metanames = setdiff(names(pars), names(pd$pars)))
  pars
}



#' Title
#'
#' @param pd
#' @param rows
#' @param parameters
#'
#' @return [dMod::parframe()]
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family parframe handling
#' @importFrom data.table as.data.table setcolorder rbindlist
#' @importFrom dMod parframe
#'
#' @examples
pd_parf_collectProfile <- function(pd, rows = c("profile_endpoints", "optimum"), parameters = NULL) {
  
  parameternames <- attr(pd$result$profiles, "parameters")
  if (!length(parameters)) parameters <- unique(pd$result$profiles$whichPar)
  
  pars <- copy(pd$result$profiles)
  pars <- data.table::as.data.table(pars)
  pars <- pars[whichPar %in% parameters]
  pars[,`:=`(profileDirection = if(constraint < 0) "left" else if (constraint > 0) "right" else "optimum"), by = 1:nrow(pars)]
  
  parsEnd <- parsOpt <- NULL
  if ("profile_endpoints" %in% rows){
    parsEnd <- pars[profileDirection != "optimum",.SD[which.max(abs(constraint))], by = c("whichPar", "profileDirection")]
    parsEnd[,`:=`(parameterSetId = paste0("profile_endpoints", whichPar))]
    data.table::setcolorder(parsEnd, c(names(pars), "parameterSetId"))
  }
  if ("optimum" %in% rows){
    parsOpt <- pars[which(profileDirection == "optimum")[1]]
    parsOpt[,`:=`(parameterSetId = "optimum")]
    data.table::setcolorder(parsOpt, c(names(pars), "parameterSetId"))
  }
  
  pars <- data.table::rbindlist(list(parsEnd, parsOpt))
  # [ ] Idea: return whichPar as well. Then, when summarizing profiles in pd_predictAndPlot2 by min and max,
  #   one could have a look in which profile the minimum/maximum values occured, i.e. which profiles contribute most to prediction uncertainty
  #   Don't know if this is informative, but for now I'd copy the prediction code which would lead to a lot of code duplication which I don't want
  #   And I shouldn't be distracted from my main goals now
  pars <- pars[,.SD, .SDcols = c("parameterSetId", "profileDirection", "value", parameternames)]
  pars <- dMod::parframe(pars, parameters = parameternames, metanames = setdiff(names(pars), parameternames))
  pars
}


#' Title
#'
#' @param pd
#' @param rows numeric vector of row indices. NULL keeps all rows
#'
#' @return [dMod::parframe()]
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family parframe handling
#'
#' @examples
pd_parf_collectL1 <- function(pd, rows = NULL) {
  if (is.null(rows)) rows <- 1:nrow(pd$result$L1)
  pars <- pd$result$L1[rows]
  pars$parameterSetId <- factor(pars$lambdaL1)
  pars
}


#' Collect parameters for simultaneous prediction or plotting
#'
#' @param pd
#' @param opt.base
#' @param opt.fit
#' @param opt.profile
#'
#' @return pd with pd$parf
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family parframe handling
#' @importFrom conveniencefunctions cf_parf_rbindlist
pd_parf_collect <- function(pd,
                            opt.base =    pd_parf_opt.base(   include = TRUE , parameterSetId = "Base"),
                            opt.mstrust = pd_parf_opt.mstrust(include = TRUE , fitrankRange = 1:20, tol = 1),
                            opt.profile = pd_parf_opt.profile(include = FALSE, rows = "profile_endpoints", parameters = NULL),
                            opt.L1 =      pd_parf_opt.L1(     include = FALSE)
) {
  
  parf_base <- parf_fit <- parf_profile <- parf_L1 <- NULL
  
  if (opt.base$include) {
    args <- c(list(pd = pd), opt.base[setdiff(names(opt.base), "include")])
    parf_base <- do.call(pd_parf_collectPars, args)}
  
  if (opt.mstrust$include && !is.null(pd$result$mstrust)) {
    args <- c(list(pd = pd), opt.mstrust[setdiff(names(opt.mstrust), "include")])
    parf_fit <- do.call(pd_parf_collectMstrust, args)}
  
  if (opt.profile$include && !is.null(pd$result$profiles)) {
    args <- c(list(pd = pd), opt.profile[setdiff(names(opt.profile), "include")])
    parf_profile <- do.call(pd_parf_collectProfile, args)}
  
  if (opt.L1$include && !is.null(pd$result$L1)) {
    args <- c(list(pd = pd), opt.L1[setdiff(names(opt.L1), "include")])
    parf_L1 <- do.call(pd_parf_collectL1, args)}
  
  parf <- conveniencefunctions::cf_parf_rbindlist(list(parf_base, parf_fit, parf_profile, parf_L1))
  
  # debatable
  parf <- parf[order(parf$value)]
  
  parf
}

#' @export
pd_parf_opt.base <- function(include = TRUE, parameterSetId = "base_obsParsFitted") {
  list(include = include, parameterSetId = parameterSetId)
}
#' @export
pd_parf_opt.mstrust <- function(include = TRUE, fitrankRange = 1:20, tol = 1) {
  list(include = include, fitrankRange = fitrankRange, tol = tol)
}

#' @export
pd_parf_opt.profile <- function(include = FALSE, rows = c("profile_endpoints", "optimum"), parameters = NULL) {
  list(include = include, rows = rows, parameters = parameters)
}

#' @export
pd_parf_opt.L1 <- function(include = FALSE, rows = NULL) {
  list(include = include, rows = rows)
}


#' Get parameters which are fixed on the boundary
#'
#' @param pd pd with pd$pars from a fit
#' @param tol tolerance how close each parameter on the est-scale should be to the boundary in order to be fixed
#'
#' @return vector of fixed parameters
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
pd_pars_getFixedOnBoundary <- function(pd, tol = 1e-2) {
  parlower <- petab_getParameterBoundaries(pd$pe, "lower")
  parlower <- parlower[names(pd$pars)]
  parlower <- pd$pars - parlower
  parlower <- parlower[parlower <= tol]
  
  parupper <- petab_getParameterBoundaries(pd$pe, "upper")
  parupper <- parupper[names(pd$pars)]
  parupper <- -(pd$pars - parupper)
  parupper <- parupper[parupper <= tol]
  
  fixedOnBoundary <- pd$pars[c(names(parlower), names(parupper))]
  fixedOnBoundary
}






# -------------------------------------------------------------------------#
# Fitting functions ----
# -------------------------------------------------------------------------#

#' Fit only obsPars
#'
#' @param pd 
#' @param FLAGoverwrite TRUE or FALSE
#'
#' @return pd with pd$result$pd$base_obsParsFitted
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family pd fitting
#' @importFrom conveniencefunctions dMod_files dMod_saveMstrust
#' @importFrom dMod trust
#'
#' @examples
pd_fitObsPars <- function(pd, FLAGoverwrite = FALSE, iterlim = 500) {
  .outputFolder <- dirname(pd$filenameParts$.compiledFolder)
  fit_file <- conveniencefunctions::dMod_files(.outputFolder, "base_obsParsFitted")$mstrust
  if (!FLAGoverwrite && file.exists(fit_file)) {
    cat("FitObsPars: Previous parameters were loaded")
    return(readPd(pd_files(pd$filenameParts)$rdsfile))
  }
  
  obspars <- petab_getParameterType(pd$pe)
  obspars <- obspars[parameterType %in% c("observableParameters", "noiseParameters"), parameterId]
  obspars <- intersect(obspars, pd$pe$parameters[estimate == 1, parameterId])
  
  fit_par <- pd$pars[obspars]
  fit_fix <- pd$pars[setdiff(names(pd$pars), obspars)]
  
  parlower <- petab_getParameterBoundaries(pd$pe, "lower")[obspars]
  parupper <- petab_getParameterBoundaries(pd$pe, "upper")[obspars]
  
  fit <- dMod::trust(pd$obj, fit_par, 1,10, iterlim = iterlim, fixed = fit_fix, parlower = parlower, parupper = parupper, printIter = TRUE)
  fit$argument <- c(fit$argument, fit_fix)
  parf_base_obsParsFitted <- as.parframe(structure(list(fit), class = c("parlist", "list")))
  conveniencefunctions::dMod_saveMstrust(fit = parf_base_obsParsFitted, path = file.path(.outputFolder), identifier = "base_obsParsFitted", FLAGoverwrite = FLAGoverwrite)
  
  readPd(pd_files(pd$filenameParts)$rdsfile)
}


#' Run a fit only for observation parameters
#'
#' @param pd
#' @param NFLAGsavePd 3 as an option doesn't make sense hre
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @family pd
#' @family pd fittting
#'
#' @importFrom dMod trust
#'
#' @examples
pd_fit <- function(pd, NFLAGsavePd = 1, iterlim = 1000, printIter = TRUE, traceFile = NULL) {
  # [ ] dont save pd, save a parframe in Results instead
  
  fit_par <- pd$pars
  fit_fix <- pd$fixed
  
  parlower <- petab_getParameterBoundaries(pd$pe, "lower")
  parupper <- petab_getParameterBoundaries(pd$pe, "upper")
  
  cat("dig into bad fitting of S302 of TGFb...\n")
  
  fit <- dMod::trust(pd$obj, fit_par, 1,10, iterlim = iterlim, fixed = fit_fix,
                     parlower = parlower, parupper = parupper, printIter = printIter, traceFile = traceFile)
  if (!fit$converged) warning("Fit not converged, please try increasing 'iterlim' (was ", iterlim,")")
  
  pd <- pd_updateEstPars(pd, parsEst = fit$argument, FLAGupdatePE = TRUE, FLAGsavePd = NFLAGsavePd > 0)
  pd
}


#' Run mstrust
#'
#' @param pd
#' @param NFLAGsavePd 3 as an option doesn't make sense hre
#'
#' @return
#' @export
#' @author svenja kemmer
#' @md
#'
#' @family pd
#' @family pd fitting
#'
#' @importFrom dMod mstrust
#'
#' @examples
pd_mstrust <- function(pd, NFLAGsavePd = T, iterlim = 1000, nfits = 5, id) {
  
  .outputFolder <- paste0("SelectionProblem/", id)
  
  fit_par <- pd$pars
  fit_fix <- pd$fixed
  
  parlower <- petab_getParameterBoundaries(pd$pe, "lower")
  parupper <- petab_getParameterBoundaries(pd$pe, "upper")
  
  # file.copy(file.path(pd$filenameParts$.currentFolder, pd$filenameParts$.compiledFolder, "/"), ".", recursive = TRUE)
  # loadDLL(pd$obj_data);
  # 
  # center <- pepy_sample_parameter_startpoints(pe = pd$pe, n_starts = nfits, seed = 1, FLAGincludeCurrent = 1)
  center <- fit_par                                             

  out <- dMod::mstrust(objfun = pd$obj,
                       center = center, 
                       studyname = paste0("SelectionProblem/", id, "/Results/trial"),
                       rinit = 1, rmax = 10, 
                       iterlim = iterlim, 
                       fits = nfits,
                       cores = detectFreeCores(),
                       fixed = fit_fix,
                       parlower = parlower, 
                       parupper = parupper)
                           
  myframe <- as.parframe(out)
  conveniencefunctions::dMod_saveMstrust(fit = myframe, path = .outputFolder, 
                                         identifier = paste0(nfits, "fits"), FLAGoverwrite = TRUE)
  bestfit <- as.parvec(myframe, 1)
                       
  if (!myframe$converged[1]) warning("Fit not converged, please try increasing 'iterlim' (was ", iterlim,")")
  
  pd <- pd_updateEstPars(pd, parsEst = myframe, FLAGupdatePE = TRUE, FLAGsavePd = NFLAGsavePd)
  pd
}


# fit_hierarchical -----
# could be implemented as objective function which automatically fits observable parameters
# pd_normHierarchical <- function(pd){
# pd <- petab_exampleRead("01", "pd")
# # normIndiv_hierarchical <- function (pd) {
#   force(pd)
# 
#   # Get Objects from pd
#   est.grid <- data.table(pd$dModAtoms$gridlist$est.grid)
#   fix.grid <- data.table(pd$dModAtoms$gridlist$fix.grid)
#   setkeyv(est.grid, c("ID", "condition"))
#   setkeyv(fix.grid, c("ID", "condition"))
#   xp0      <- (pd$dModAtoms$fns$x * pd$dModAtoms$fns$p0)
#   gxp0     <- (pd$dModAtoms$fns$g * pd$dModAtoms$fns$x * pd$dModAtoms$fns$p0)
#   g0       <- pd$dModAtoms$fns$g
#   data     <- pd$dModAtoms$data
#   errmodel <- pd$dModAtoms$e
# 
#   # Parameter types
#   pepa <- pd$pe$parameters
#   pepa <- pepa[petab_getParameterType(pd$pe), on = "parameterId"]
#   pepa <- pepa[estimate == 1]
#   parametersObsErr <- pepa[parameterType != "other", parameterId]
#   parametersOther  <- pepa[parameterType != "other", parameterId]
#   parsObsErr <- pd$pars[parametersObsErr]
# 
#   # Salat from dMod
#   timesD          <- pd$times
#   x.conditions    <- est.grid$condition
#   data.conditions <- names(data)
#   e.conditions    <- names(attr(errmodel, "mappings"))
#   controls        <- list(times = timesD, attr.name = attr.name, conditions = intersect(x.conditions, data.conditions))
# 
#   myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = controls$conditions,
#                    simcores = 1,
#                    NFLAGbrowser = 0,
#                    FLAGverbose = FALSE,
#                    FLAGNaNInfwarnings = FALSE) {
# 
#     arglist <- list(...)
#     arglist <- arglist[match.fnargs(arglist, "pars")]
#     pars <- arglist[[1]]
# 
#     # Logic
#     # * pars_optObs <- c(parsObsErr); fixed <- c(pars, fixed)
#     # * predictionXP0 <- xp0(times, pars_optObs)
#     # * normObsParsErr <- function(parsObsParsErr, predictionXP0) norm(g(predictionXP0, pars_optObs), data)
#     # * fit <- trust(normObsParsErr, parsObsErr)
#     # * parsObsErr <<- fit$argument
#     # * pars_optOuter <- c(pars, parsObsErr); fixed <- fixed
#     # * norm(g*x*p0, pars_optOuter), data)
# 
# 
#     cn <- (setNames(nm = conditions))[[1]]
#     objlists <- lapply(setNames(nm = conditions), function(cn) {
# 
#       if (FLAGbrowser) browser()
# 
#       ID <- est.grid[condition == cn, ID]
#       if (FLAGverbose) cat(ID, cn, "\n", sep = " ---- ")
#       dummy <- make_pars(pars, fixed, est.grid, fix.grid, ID)
#       pars_ <- dummy$pars
#       fixed_ <- dummy$fixed
# 
#       if (!length(pars_)) return(init_empty_objlist(pars, deriv = deriv, FLAGchisquare = TRUE)) # No pars_ can happen if one fits only condition specific parameters and in this condition there are none
# 
#       prediction <- try(prd0(times = controls$times, pars = pars_, fixed = fixed_, deriv = deriv))
# 
#       if (inherits(prediction, "try-error"))
#         stop("Prediction failed in \n>>>condition = ", cn, "\n>>>ID = ", ID, "\n\nTry iterating p(pars), (x*p)(pars), ... to find the problem.")
# 
#       prediction <- prediction[[1]]
#       prediction <- check_and_sanitize_prediction(prediction, data, cn, FLAGNaNInfwarnings)
# 
#       err <- NULL
#       if (any(is.na(data[[cn]]$sigma))) {
#         err <- errmodel(out = prediction, pars = getParameters(prediction), conditions = cn, deriv=deriv)
#         mywrss <- nll(res(data[[cn]], prediction, err[[1]]), deriv = deriv, pars = pars)
#       } else {
#         mywrss <- nll(res(data[[cn]], prediction), deriv = deriv, pars = pars)
#       }
# 
#       if (deriv) mywrss <- renameDerivParsInObjlist(mywrss, dummy$parnames)
# 
#       mywrss
#     })
# 
#     # Sum all objlists
#     out <- Reduce("+", objlists)
# 
#     # Consider fixed: return only derivs wrt pouter
#     out$gradient <- out$gradient[names(pars)]
#     out$hessian <- out$hessian[names(pars), names(pars)]
# 
#     # Populate attributes
#     attr(out, controls$attr.name) <- out$value
#     ll_conditions <- data.frame(
#       logl = vapply(setNames(objlists, conditions), function(.x) .x$value, 1),
#       chi2 = vapply(setNames(objlists, conditions), function(.x) attr(.x, "chisquare"), 1))
#     ll_sum <- data.frame(logl = sum(ll_conditions$logl),
#                          chi2 = sum(ll_conditions$chi2))
#     attributes(out) <- c(attributes(out), list(ll_cond_df = ll_conditions))
#     attributes(out) <- c(attributes(out), list(ll_sum_df = ll_sum))
#     # attr(out, "AIC") <- out$value + length(pars) * 2
#     # attr(out, "BIC") <- out$value + length(pars) * log(nrow(as.data.frame(data)))
#     return(out)
#   }
# 
#   class(myfn) <- c("objfn", "fn")
#   attr(myfn, "conditions") <- data.conditions
#   attr(myfn, "parameters") <- attr(prd0, "parameters")
#   attr(myfn, "modelname") <- modelname(prd0, errmodel)
#   return(myfn)
# }


# -------------------------------------------------------------------------#
# Cluster ----
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
#' @family Cluster
#'
#' @examples
clusterStatusMessage <- function(FLAGjobDone, FLAGjobPurged, FLAGjobRecover) {
  if (FLAGjobPurged)   return("Job is done and purged")
  if (FLAGjobDone)     return("Job is done")
  if (FLAGjobRecover)  return("Job is running")
  if (!FLAGjobRecover) return("Job is not yet started")
}


#' Fit model on cluster
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
#' @family Cluster
#' @importFrom conveniencefunctions dMod_files cf_as.parframe dMod_saveMstrust check_clusterTimeStamp
#' @importFrom dMod distributed_computing
#'
#' @examples
pd_cluster_mstrust <- function(pd, .outputFolder, n_startsPerNode = 16*3, n_nodes = 10, 
                               identifier = "mstrust", FLAGforcePurge = FALSE, opt.parameter_startpoints = "sample") {
  # .. General job handling -----
  jobnm <- paste0("mstrust_", identifier, "_", gsub("-","_",gsub("(S\\d+(-\\d+)?).*", "\\1", basename(.outputFolder))))
  
  fileJobDone    <- conveniencefunctions::dMod_files(.outputFolder, identifier)[["mstrust"]]
  fileJobPurged  <- file.path(dirname(fileJobDone), paste0(".", jobnm, "jobPurged"))
  fileJobRecover <- file.path(paste0(jobnm, "_folder"), paste0(jobnm, ".R"))
  
  FLAGjobDone    <- file.exists(fileJobDone)
  FLAGjobPurged  <- file.exists(fileJobPurged)
  FLAGjobRecover <- file.exists(fileJobRecover) | FLAGjobDone | FLAGjobPurged
  
  cat(clusterStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  if (!FLAGjobPurged) conveniencefunctions::check_clusterTimeStamp()
  
  # Assign Global variables: Important, in future, this might be a source of bugs, if other cluster-functions are written
  assign("n_startsPerNode",n_startsPerNode,.GlobalEnv)
  assign("opt.parameter_startpoints",opt.parameter_startpoints,.GlobalEnv)
  
  # Start mstrust job
  file.copy(file.path(pd$filenameParts$.currentFolder, pd$filenameParts$.compiledFolder, "/"), ".", recursive = TRUE)
  job <- dMod::distributed_computing(
    {
      loadDLL(pd$obj_data);
      
      seed <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1
      FLAGincludeCurrent <- seed == 1
      
      if (identical(opt.parameter_startpoints, "sample")){
        center <- pepy_sample_parameter_startpoints(pd$pe, n_starts = n_startsPerNode, seed = seed, 
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
              rinit = 0.1, rmax = 10, cores = 16,
              iterlim = 500, 
              optmethod = "trust", 
              output = TRUE, cautiousMode = TRUE,
              stats = FALSE, 
              parlower = parlower, parupper = parupper)
    },
    jobname = jobnm, 
    partition = "single", cores = 16, nodes = 1, walltime = "12:00:00",
    ssh_passwd = Sys.getenv("hurensohn"), machine = "cluster", 
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
      Sys.sleep(1) # avoid being blocked
      job$get()
      fitlist  <- if (exists("cluster_result")) do.call(c, cluster_result) else {NULL
        #cf_dMod_rescueFits()
        #   fitlist <- list.files(file.path(paste0(jobnm, "_folder"), "results","fit"), "\\.Rda$", recursive = TRUE, full.names = TRUE)
        #   fitlist <- lapply(fitlist, function(x) try(local(load(x))))
      }
      fits <- fitlist
      fits <- fits[vapply(fits, is.list, TRUE)]
      class(fits) <- "parlist"
      fits <- conveniencefunctions::cf_as.parframe(fits)
      conveniencefunctions::dMod_saveMstrust(fit = fits, path = .outputFolder, 
                                             identifier = identifier, FLAGoverwrite = TRUE)
      savedFits <- readRDS(fileJobDone)
      if (nrow(savedFits) != (n_startsPerNode * n_nodes)){
        return("Job done, but some fits had errors. You can check out the results by running `readPd` which will load the fit into pd$result$fits. Re-run this function once more to purge the job.")
      } else cat("Job done and all went fine.\n")
    }
  }
  
  FLAGjobDone    <- file.exists(fileJobDone)
  if (FLAGjobDone & !FLAGjobPurged) {
    if (readline("Purge job. Are you sure? Type yes: ") == "yes"){
      job$purge(purge_local = TRUE)
      writeLines("jobPurged", fileJobPurged)
      return("Job purged\n")
    }
  }
}


#' Run profiles on cluster
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
#' @family Cluster
#' @importFrom dMod profile_pars_per_node
#' @importFrom conveniencefunctions check_clusterTimeStamp
#'
#' @examples
pd_cluster_profile <- function(pd, .outputFolder, FLAGforcePurge = FALSE, FLAGfixParsOnBoundary = TRUE, 
                               profpars = names(pd$pars)) {
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
  
  cat(clusterStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  if (!FLAGjobPurged) conveniencefunctions::check_clusterTimeStamp()
  
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
    ssh_passwd = Sys.getenv("hurensohn"), machine = "cluster", 
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
      cat("run readPd() again to load the results into the pd")
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




#' Fit model on cluster
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
#' @family Cluster
#' @importFrom conveniencefunctions dMod_files cf_as.parframe dMod_saveMstrust cf_parf_metaNames0 check_clusterTimeStamp
#' @importFrom dMod distributed_computing profile_pars_per_node
#'
#' @examples
pd_cluster_L1 <- function(pd, .outputFolder, n_nodes = 6, lambdas = 10^(seq(log10(0.0001), log10(100), length.out = n_nodes*16-1)), 
                          identifier = "L1", FLAGforcePurge = FALSE) {
  # .. General job handling -----
  jobnm <- paste0("mstrust_", identifier, "_", gsub("(S\\d+).*", "\\1", basename(.outputFolder)))
  
  fileJobDone    <- conveniencefunctions::dMod_files(.outputFolder, identifier)$L1
  fileJobPurged  <- file.path(dirname(fileJobDone), paste0(".", jobnm, "jobPurged"))
  fileJobRecover <- file.path(paste0(jobnm, "_folder"), paste0(jobnm, ".R"))
  
  FLAGjobDone    <- file.exists(fileJobDone)
  FLAGjobPurged  <- file.exists(fileJobPurged)
  FLAGjobRecover <- file.exists(fileJobRecover) | FLAGjobDone | FLAGjobPurged
  if (!FLAGjobPurged) conveniencefunctions::check_clusterTimeStamp()
  
  cat(clusterStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  
  # Assign Global variables: Important, in future, this might be a source of bugs, if other cluster-functions are written
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
    ssh_passwd = Sys.getenv("hurensohn"), machine = "cluster", 
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
      fits <- cf_parframe(fits, metanames = conveniencefunctions::cf_parf_metaNames0$l1)
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


#' Fit model on cluster
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
#' @family Cluster
#' @importFrom conveniencefunctions dMod_files cf_as.parframe dMod_saveMstrust check_clusterTimeStamp
#' @importFrom dMod distributed_computing
#'
#' @examples
pd_cluster_L1_fitUnbiasedEachMstrust <- function(pd, .outputFolder, n_startsPerNode = 16*3, 
                                                 identifier = "L1UB", FLAGforcePurge = FALSE) {
  
  stop("implement saving the node id (see S311 for the problem)")
  # .. General job handling -----
  jobnm <- paste0("L1UB_", identifier, "_", gsub("(S\\d+).*", "\\1", basename(.outputFolder)))
  
  fileJobDone    <- conveniencefunctions::dMod_files(.outputFolder, paste0(identifier,1))[["mstrust"]]
  fileJobPurged  <- file.path(dirname(fileJobDone), paste0(".", jobnm, "jobPurged"))
  fileJobRecover <- file.path(paste0(jobnm, "_folder"), paste0(jobnm, ".R"))
  
  FLAGjobDone    <- file.exists(fileJobDone)
  FLAGjobPurged  <- file.exists(fileJobPurged)
  FLAGjobRecover <- file.exists(fileJobRecover) | FLAGjobDone | FLAGjobPurged
  if (!FLAGjobPurged) conveniencefunctions::check_clusterTimeStamp()
  
  cat(clusterStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  
  # Hacking var_list: want to have different nodes
  n_nodes <- nrow(L1_getModelCandidates(pd$result$L1))
  var_list <- list(node = 1:n_nodes)
  
  # Assign Global variables: Important, in future, this might be a source of bugs, if other cluster-functions are written
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
      
      pd$pars <- fit_par
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
      try(conveniencefunctions::cf_as.parframe(fit))
    },
    jobname = jobnm, 
    partition = "single", cores = 16, nodes = 1, walltime = "12:00:00",
    ssh_passwd = Sys.getenv("hurensohn"), machine = "cluster", 
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
      parfs <- cluster_result
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


#' Fit model on cluster
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
#' @family Cluster
#' @importFrom conveniencefunctions dMod_files cf_as.parframe dMod_saveMstrust check_clusterTimeStamp
#' @importFrom dMod distributed_computing
#'
#' @examples
pd_cluster_L1_fitUnbiasedEachOnce <- function(pd, .outputFolder, n_startsPerNode = 16*3, 
                                              identifier = "L1UBSingle", FLAGforcePurge = FALSE) {
  
  # .. General job handling -----
  jobnm <- paste0("L1UB_", identifier, "_", gsub("(S\\d+).*", "\\1", basename(.outputFolder)))
  
  fileJobDone    <- conveniencefunctions::dMod_files(.outputFolder, identifier)[["mstrust"]]
  fileJobPurged  <- file.path(dirname(fileJobDone), paste0(".", jobnm, "jobPurged"))
  fileJobRecover <- file.path(paste0(jobnm, "_folder"), paste0(jobnm, ".R"))
  
  FLAGjobDone    <- file.exists(fileJobDone)
  FLAGjobPurged  <- file.exists(fileJobPurged)
  FLAGjobRecover <- file.exists(fileJobRecover) | FLAGjobDone | FLAGjobPurged
  if (!FLAGjobPurged) conveniencefunctions::check_clusterTimeStamp()
  
  cat(clusterStatusMessage(FLAGjobDone, FLAGjobPurged, FLAGjobRecover), "\n")
  
  # Hacking var_list: want to have different nodes
  n_nodes <- nrow(L1_getModelCandidates(pd$result$L1))
  var_list <- dMod::profile_pars_per_node(1:n_nodes, 16)
  
  # Assign Global variables: Important, in future, this might be a source of bugs, if other cluster-functions are written
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
    ssh_passwd = Sys.getenv("hurensohn"), machine = "cluster", 
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
      fitlist  <- if (exists("cluster_result")) do.call(c, cluster_result) else {NULL
        #cf_dMod_rescueFits()
        #   fitlist <- list.files(file.path(paste0(jobnm, "_folder"), "results","fit"), "\\.Rda$", recursive = TRUE, full.names = TRUE)
        #   fitlist <- lapply(fitlist, function(x) try(local(load(x))))
      }
      fits <- fitlist
      fits <- fits[vapply(fits, is.list, TRUE)]
      class(fits) <- "parlist"
      fits <- conveniencefunctions::cf_as.parframe(fits)
      conveniencefunctions::dMod_saveMstrust(fit = fits, path = .outputFolder, 
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










# -------------------------------------------------------------------------#
# L1 ----
# -------------------------------------------------------------------------#

#' Prepare data for some L1 plots
#'
#' @param pd with pd$result$L1
#'
#' @return data.table(parameterId,parameterType,isL1Parameter,lambda,value,iterations,converged,muValue,estValue,estDeviance)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family L1
#' @importFrom data.table as.data.table data.table
#' 
#' @examples
pd_L1_preparePlottingData <- function(pd) {
  d <- pd$result$L1
  d <- data.table::as.data.table(as.data.frame(d))
  d <- melt(d, measure.vars = names(pd$pars), variable.name = "parameterId", variable.factor = FALSE, value.name = "estValue")
  d <- petab_getParameterType(pd$pe)[d, on = "parameterId"]
  d[,`:=`(isL1Parameter = parameterId %in% pd$L1$parametersL1)]
  d <- data.table::data.table(muValue = pd$pars, parameterId = names(pd$pars))[d, on = "parameterId"]
  d[,`:=`(estDeviance = muValue - estValue)]
  d[,`:=`(estDeviance = round(estDeviance, 4))]
  d
}

#' Title
#'
#' @param pd pd with pd$result$L1
#' @param tol Tolerance 
#'
#' @return data.table(lambdaL1 = lambda, Ntotal = Total parameters which were L1'd, NFree = Remaining free parameters, df = degrees of freedom, chis = The quantile value which corresponds to significance level alpha = 0.5)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family L1
pd_L1_getObjBounds <- function(pd,tol = 1e-4) {
  d <- pd_L1_preparePlottingData(pd)
  d <- d[isL1Parameter == TRUE]
  d <- d[,list(Ntotal = .N, NFree = sum(abs(estDeviance) > tol)), by = "lambdaL1"]
  # Best way would be to include a fit with lambda = 0 and compare against this fit
  # This way one could also exclude "non-informative" parameters, such as k_233,k_234,k_244 in JS1
  d[,`:=`(df = max(NFree) - NFree)]
  d[,`:=`(chis = qchisq(0.95, df))]
  d
}


#' Title
#'
#' @param pd with pd$result$L1
#' @param ... Arguments going to cf_outputFigure
#'
#' @return ggplot
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family L1
#' @importFrom conveniencefunctions cfggplot cf_outputFigure
pd_L1_plotDevianceVsLambda <- function(pd, ...) {
  d <- pd_L1_preparePlottingData(pd)
  pl <- conveniencefunctions::cfggplot(d, aes(lambdaL1, estDeviance, color = parameterId)) + 
    facet_wrap(~parameterType + isL1Parameter, scales = "free") + 
    geom_line() + 
    scale_color_viridis_d() + 
    scale_x_log10() + 
    theme(legend.position = "bottom") + 
    geom_blank()
  conveniencefunctions::cf_outputFigure(pl, ...)
}

#' Title
#'
#' @param pd 
#' @param ... 
#'
#' @return ggplot
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family L1
#' @importFrom conveniencefunctions cfggplot cf_outputFigure
#'
#' @examples
pd_L1_plotIsConvergedVsLambda <- function(pd, ...) {
  d <- pd_L1_preparePlottingData(pd)
  pl <- conveniencefunctions::cfggplot(unique(d[,list(isconverged = as.numeric(converged), lambdaL1)]), aes(lambdaL1, isconverged)) + 
    geom_point() + 
    scale_x_log10()
  conveniencefunctions::cf_outputFigure(pl, ...)
}

#' Title
#'
#' @param pd 
#' @param ... 
#'
#' @return ggplot
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family L1
#' @importFrom conveniencefunctions cfggplot cf_outputFigure
#'
#' @examples
pd_L1_plotValueVsLambda <- function(pd, ...) {
  d <- pd_L1_preparePlottingData(pd)
  dbounds <- pd_L1_getObjBounds(pd)
  dbounds <- unique(d[,list(value, lambdaL1)])[dbounds, on = "lambdaL1"]
  dbounds[,`:=`(valuePlot = min(value) + chis)]
  pl <- conveniencefunctions::cfggplot() + 
    geom_point(aes(lambdaL1, value), data = unique(d[,list(value, lambdaL1)])) + 
    geom_line(aes(lambdaL1, valuePlot), data = dbounds) + 
    scale_x_log10()
  conveniencefunctions::cf_outputFigure(pl, ...)
}


#' Title
#'
#' @param pd 
#' @param ... 
#'
#' @return ggplot
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family L1
#' @importFrom parallel mclapply
#' @importFrom dMod as.parvec
#'
#' @examples
pd_L1_getUnbiasedValues <- function(pd, ...) {
  ncores <- 1
  values <- parallel::mclapply(X = seq_len(nrow(pd$result$L1)), mc.cores = ncores, FUN = function(i) {
    # add L1-fixed parameters to "fixed"
    
    # Fit unbiased
    # Calculate unbiased obj value
    pd$obj(dMod::as.parvec(pd$result$L1,i), fixed = pd$fixed, deriv = FALSE)
  })
  pd$result$L1$valueUnbiased <- values
  pd
}


#' Title
#'
#' @param L1Scan L1-parframe Result from pd_cluster_L1
#'
#' @return Matrix indicating the model structure of L1 models
#' @export
#'
#' @examples
L1_getModelCandidates <- function(L1Scan) {
  fixed_L1 <- as.matrix(L1Scan)
  fixed_L1 <- fixed_L1[,grep("^L1_", colnames(fixed_L1)),drop=FALSE]
  fixed_L1 <- fixed_L1 == 0
  fixed_L1 <- unique(fixed_L1)
  fixed_L1
}




# -------------------------------------------------------------------------#
# Plotting ----
# -------------------------------------------------------------------------#

#' First version of plotting a pd
#'
#' Predicts with "pars" and "times"
#'
#' @param pd A pd
#' @param ... Going to [dMod::plotCombined()]
#' @param page,nrow,ncol Going to [ggforce::facet_wrap_paginate()]
#' @param filename,width,height,scale,units Going to [conveniencefunctions::cf_outputFigure()]
#'
#' @return ggplot
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Plotting
#' @importFrom dMod plotCombined
#' @importFrom ggforce facet_wrap_paginate
#' @importFrom conveniencefunctions theme_cf scale_color_cf cf_outputFigure
#'
#' @examples
#' pdx <- petab_exampleRead("01", "pd")
#' pd_plot(pdx)
pd_plot <- function(pd, ..., page = 1, nrow = 3, ncol = 4, filename = NULL, width = 29.7, height = 21, scale = 1, units = "cm"){
  pred <- pd$prd(pd$times, pd$pars)
  pl <- dMod::plotCombined(pred, pd$dModAtoms$data, ...) +
    ggforce::facet_wrap_paginate(~name, nrow = nrow, ncol = ncol, scales = "free", page = page) +
    conveniencefunctions::theme_cf() +
    conveniencefunctions::scale_color_cf()
  conveniencefunctions::cf_outputFigure(pl, filename, width = width, height = height, scale = scale, units = units)
}




#' Title
#'
#' @param pd
#' @param opt.base
#' @param opt.mstrust
#' @param FLAGsubsetPredictionToData
#' @param page
#' @param filename
#' @param width
#' @param height
#' @param scale
#' @param units
#' @param i "data.table i" to subset data and predictions
#' @param FLAGsummarizeProfilePredictions
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Plotting
#' @importFrom conveniencefunctions cf_predict cfggplot scale_color_cf cf_outputFigure
#' @importFrom ggforce facet_wrap_paginate
#'
#' @examples
pd_predictAndPlot <- function(pd, i,
                              opt.base = pd_parf_opt.base(),
                              opt.mstrust = pd_parf_opt.mstrust(),
                              opt.profile = pd_parf_opt.profile(FALSE),
                              NFLAGsubsetType = c(none = 0, strict = 1, keepInternal = 2)[2],
                              FLAGsummarizeProfilePredictions = TRUE,
                              nrow = 3, ncol = 4, scales = "free", page = 1,
                              n.breaks = 5,
                              filename = NULL, width = 29.7, height = 21, scale = 1, units = "cm") {
  # Catch i (see petab_mutateDCO for more ideas)
  mi <- missing(i)
  si <- substitute(i)
  
  # .. Collect parameters and predict -----
  parf <- pd_parf_collect(pd, opt.base = opt.base, opt.mstrust = opt.mstrust, opt.profile = opt.profile)
  if (nrow(parf) > 5) cat("Predicting for more than 5 parameter sets. Are you sure?")
  predictions <- conveniencefunctions::cf_predict(prd = pd$prd, times = pd$times, pars = parf)
  
  # .. Prepare data and prediction -----
  pplot <- copy(predictions)
  setnames(pplot,
           c("condition"  , "name"        , "value"),
           c("conditionId", "observableId", "measurement"))
  pplot[,`:=`(fitrank = as.factor(fitrank))]
  
  # Apply observableScale to data
  dplot <- petab_joinDCO(pd$pe)
  dplot[,`:=`(measurement = eval(parse(text = paste0(observableTransformation, "(", measurement, ")")))), by = 1:nrow(dplot)]
  
  # set Order: Put pure predictions to end of plot
  pplot[,`:=`(hasData = observableId %in% unique(dplot$observableId))]
  pplot <- pplot[order(-hasData)]
  pplot[,`:=`(observableId=factor(observableId, unique(observableId)))]
  dplot[,`:=`(observableId=factor(observableId, levels(pplot$observableId)))]
  
  # subset conditions and observables of predictions: mutually exclusive with previous...
  pplot <- subsetPredictionToData(pplot, dplot, NFLAGsubsetType = NFLAGsubsetType)
  
  # [ ] Idea for plotting predictions along a profile: Summarize by min and max?
  
  # Handle i
  if (!mi) {dplot <- dplot[eval(si)]; pplot <- pplot[eval(si)]}
  
  # .. Plot -----
  pl <- conveniencefunctions::cfggplot() +
    ggforce::facet_wrap_paginate(~observableId, nrow = nrow, ncol = ncol, scales = scales, page = page) +
    geom_line(aes(time, measurement, color = conditionId, linetype = parameterSetId), data = pplot) +
    geom_point(aes(time, measurement, color = conditionId), data = dplot) +
    scale_y_continuous(n.breaks = n.breaks) +
    conveniencefunctions::scale_color_cf()
  
  conveniencefunctions::cf_outputFigure(pl, filename = filename, width = width, height = height, scale = scale, units = units)
}



#' Title
#'
#' @param pd 
#' @param pe to draw other data
#' @param i,j data.table arguments working on dco, but note that j will be evaulated in a separate step *after* i
#' @param parf parframe to simulate with. if supplied, opt.base,opt.mstrust,opt.profile are meaningless
#' @param opt.base,opt.mstrust,opt.profile Options to get parameters from parframes for prediction purposes
#' @param NFLAGsubsetType subset *species*, *observableId*, *conditionId* and "time".
#'   * 0 none                 : no subsetting
#'   * 1 strict               : Only show predictions for conditions where there is data
#'   * 2 keepInternal         : For observableIds match conditionIds to data, for internal species, no subsetting
#'   * 3 strict_cutTimes      : as 1, but cut time to data for observableIds
#'   * 4 keepInternal_cutTimes: as 2, but cut time to data for observableIds
#' @param FLAGsummarizeProfilePredictions summarize predictions based on profile likelihoods by ribbon
#' @param FLAGmeanLine draw line for mean(data)
#' @param aeslist list of aestheticss
#' @param ggCallback list(ggplot2 calls), e.g. list(labs(title = "bla"), scale_y_log10())
#' @param filename,FLAGfuture,width,height,scale,units Agurments going to [conveniencefunctions::cf_outputFigure()]
#' @param FLAGreturnPlotData return list(data.tables) which go into plotting instead of ggplot
#'
#' @return ggplot
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Plotting
#' @importFrom conveniencefunctions cf_predict cfggplot cfgg_getAllAesthetics scale_color_cf cf_outputFigure
#' @importFrom data.table setnames
#' @importFrom ggforce n_pages
#'
#' @examples
#' # pd <- petab_exampleRead("02", "pd")
#' pe = pd$pe
#' opt.base = pd_parf_opt.base(FALSE,)
#' opt.mstrust = pd_parf_opt.mstrust(FALSE, fitrankRange = 1:2)
#' opt.profile <- pd_parf_opt.profile(T, parameters = c("geneSerpine1_act2", "geneSerpine1_inh1", "geneSerpine1_inh2", "geneSerpine1_inh3", "pRec_degind", "S_phos"))
#' NFLAGsubsetType = 0
#' FLAGsummarizeProfilePredictions = TRUE
#' FLAGmeanLine = FALSE
#' aeslist = petab_plotHelpers_aeslist()
#' ggCallback = list(facet_wrap_paginate(~observableId, nrow = 4, ncol = 4, scales = "free"), scale_y_continuous(n.breaks = 5))
#' opt.gg = list(ribbonAlpha = 0.2)
#' filename = NULL
#' FLAGfuture = TRUE
#' width = 29.7
#' height = 21
#' scale = 1
#' units = "cm"
pd_predictAndPlot2 <- function(pd, pe = pd$pe,
                               i,j,
                               opt.base = pd_parf_opt.base(),
                               opt.mstrust = pd_parf_opt.mstrust(),
                               opt.profile = pd_parf_opt.profile(FALSE),
                               opt.L1 = pd_parf_opt.L1(FALSE),
                               parf = NULL,
                               NFLAGsubsetType = c(none = 0, strict = 1, keepInternal = 2, strict_cutTimes = 3,keepInternal_cutTimes = 3)["strict_cutTimes"],
                               FLAGsummarizeProfilePredictions = TRUE,
                               FLAGmeanLine = FALSE,
                               aeslist = petab_plotHelpers_aeslist(),
                               ggCallback = list(facet_wrap_paginate(~observableId, nrow = 4, ncol = 4, scales = "free"),
                                                 scale_y_continuous(n.breaks = 5)),
                               opt.sim = list(Ntimes_gt5ParSetIds = 100, predtimes = NULL),
                               opt.gg = list(ribbonAlpha = 0.2), # would be nice to put this into opt.profile or maybe opt.gg?
                               filename = NULL, FLAGfuture = TRUE,
                               width = 29.7, height = 21, scale = 1, units = "cm",
                               FLAGreturnPlotData = FALSE
) {
  
  
  # .. Catch i (see petab_mutateDCO for more ideas) -----
  mi <- missing(i)
  si <- substitute(i)
  mj <- missing(j)
  sj <- substitute(j)
  
  # .. Data -----
  # observableTransformation
  dplot <- petab_joinDCO(pe)
  dplot[,`:=`(measurement = eval(parse(text = paste0(observableTransformation, "(", measurement, ")")))), by = 1:nrow(dplot)]
  dplot[,`:=`(observableId=factor(observableId, petab_plotHelpers_variableOrder(pd)))]
  
  # .. Prediction -----
  if (is.null(parf)) parf <- pd_parf_collect(pd, opt.base = opt.base, opt.mstrust = opt.mstrust, opt.profile = opt.profile, opt.L1 = opt.L1)
  if (nrow(parf) > 5) {
    pd$times <- pd_predtimes(pd, N = opt.sim$Ntimes_gt5ParSetIds)
    if (!opt.profile$include) cat("Predicting for more than 5 parameter sets. Are you sure?")
  }
  if (!is.null(opt.sim$predtimes)) pd$times <- opt.sim$predtimes
  simconds <- if (NFLAGsubsetType == 0) pd$pe$experimentalCondition else if (!mi && !NFLAGsubsetType%in%c(2,4)) dplot[eval(si)] else dplot
  simconds <- unique(simconds[,conditionId])
  pplot <- conveniencefunctions::cf_predict(prd = pd$prd, times = pd$times, pars = parf, fixed = pd$fixed, conditions = simconds)
  data.table::setnames(pplot, c("condition", "name", "value"), c("conditionId", "observableId", "measurement"))
  pplot[,`:=`(observableId=factor(observableId,  petab_plotHelpers_variableOrder(pd)))]
  pplot <- subsetPredictionToData(pplot, dplot, NFLAGsubsetType = NFLAGsubsetType)
  
  # .. Error model / prediction ribbon -----
  pplotRibbon <- NULL
  if (FLAGsummarizeProfilePredictions && opt.profile$include) {
    # pd_plotHelpers_summarizeProfilePredictions - in here not as separate function because of multiple return elements
    pplotRibbon <- pplot[profileDirection %in% c("left", "right")]
    pplotRibbon <- pplotRibbon[,list(
      measurementmin = min(measurement),
      measurementmax = max(measurement),
      parameterSetIdmin = parameterSetId[which.min(measurement)],
      parameterSetIdmax = parameterSetId[which.max(measurement)]
    ), by = c("conditionId", "observableId", "time")] # Last one is a bit hacky. Better use different aeslist functions for profiles and mstrusts
    pplot <- pplot[profileDirection == "optimum"]
  }
  
  # .. Make experimentalCondition Columns available -----
  # [ ] Idea: instead of merging experimentalCondition, one could also merge a boiled down version of dplot.
  #   Then one could have a j argument which first works on dplot before merging.
  #   However, this is not necessary because this is already possible because pe is supplied separately which allows to do exactly this thing.
  #   [ ] Document this use as example
  pplot <- pe$experimentalCondition[pplot, on = c("conditionId")]
  if (!is.null(pplotRibbon)) pplotRibbon <- pe$experimentalCondition[pplotRibbon, on = c("conditionId")]
  
  # .. HACK: Add parameterSetId to dplot do grouping can be performed with respect to it -----
  dplot[,`:=`(parameterSetId = "DATA")]
  
  # .. Handle i -----
  if (!mi) {
    dplot <- dplot[eval(si)]
    pplot <- pplot[eval(si)]
    if (!is.null(pplotRibbon)) pplotRibbon <- pplotRibbon[eval(si)]
  }
  if (!mj) {
    dplot[,eval(sj)]
    pplot[,eval(sj)]
    if (!is.null(pplotRibbon)) pplotRibbon[,eval(sj)]
  }
  
  # HACK: Return data, don't plot. Is this nice? Think about it.
  # Probably best to make a dedicated plotting function which takes this list as input
  if (FLAGreturnPlotData) {
    return(list(dplot = dplot, pplot = pplot, pplotRibbon = pplotRibbon))
  }
  
  # .. Plot -----
  pl <- conveniencefunctions::cfggplot()
  if (FLAGmeanLine) { # Add first so the lines don't mask the points
    dmean <- petab_plotHelpers_meanMeasurementsValue(dplot, aeslist)
    aesmeanlist <- list(linetype = ~conditionId, group = ~conditionId)
    aesmeanlist <- c(aeslist, aesmeanlist[setdiff(names(aesmeanlist), names(aeslist))])
    aesmeanlist <- aesmeanlist[setdiff(names(aesmeanlist), "size")]
    pl <- pl + geom_line(do.call(aes_q, aesmeanlist), data = dmean, size = 0.1) # make size available as parameter
  }
  if (nrow(dplot)) pl <- pl + geom_point(do.call(aes_q, aeslist[intersect(names(aeslist), conveniencefunctions::cfgg_getAllAesthetics()[["geom_point"]])]), data = dplot)
  if (nrow(pplot)) pl <- pl + geom_line( do.call(aes_q, aeslist[intersect(names(aeslist), conveniencefunctions::cfgg_getAllAesthetics()[["geom_line"]])]) , data = pplot)
  if (!is.null(pplotRibbon)) {
    aesl <- aeslist[intersect(names(aeslist), conveniencefunctions::cfgg_getAllAesthetics()[["geom_ribbon"]])]
    aesl <- aesl[setdiff(names(aesl), c("linetype", "lty", "y", "color", "colour"))]
    pl <- pl + geom_ribbon(do.call(aes_q, aesl), data = pplotRibbon, alpha = opt.gg$ribbonAlpha)}
  pl <- pl + conveniencefunctions::scale_color_cf(aesthetics = c("color", "fill"))
  for (plx in ggCallback) pl <- pl + plx
  
  # .. Print paginate message so user doesnt forget about additional pages -----
  message("Plot has ", ggforce::n_pages(pl), " pages\n")
  
  # Output
  conveniencefunctions::cf_outputFigure(pl = pl, filename = filename, width = width, height = height, scale = scale, units = units, FLAGFuture = FLAGfuture)
}




#' Plot multiple parameter vectors in parf as barplots next to each other
#'
#' @param pd
#' @param parf
#'
#' @return
#' @export
#'
#' @examples
pd_plot_compareParameters <- function(pd, parf,
                                      nrow = 1, ncol = 3, scales = "free", page = 1,
                                      ...
) {
  parameters <- attr(parf, "parameters")
  p <- data.table(parf)
  p <- melt(p, id.vars = "parameterSetId", measure.vars = parameters, variable.name = "parameterId", variable.factor = FALSE, value.name = "parameterValueEstScale")
  
  pt <- petab_getParameterType(pd$pe)
  p <- pt[p, on = "parameterId"]
  
  pl <- conveniencefunctions::cfggplot(p, aes(parameterId, parameterValueEstScale)) +
    conveniencefunctions::facet_wrap_paginate(~parameterType, nrow = nrow, ncol = ncol, scales = scales, page = page) +
    geom_col(aes(fill = parameterSetId), position = position_dodge2()) +
    conveniencefunctions::scale_color_cf(aesthetics = c("fill", "color")) +
    theme(axis.text.x = element_text(angle = 90))
  
  cf_outputFigure(pl, ...)
}


#' Compare two parameter vectors by plotting their difference as bar plot
#'
#' @param v1,v2 Parameters on estScale
#' @param pe Petab to infer parameterType and estimated parameters
#' @param filename for plotting
#'
#' @return ggplot
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @familyPlotting
#' @importFrom lemon facet_rep_grid
plotParameterComparison <- function(v1,v2, pe, filename, FLAGexcludeObsPars = TRUE, FLAGestimatedOnly = TRUE) {
  # match and calculate difference
  d1 <-  data.table(parameterId = names(v1), estValue = unclass(v1))
  d2 <-  data.table(parameterId = names(v2), estValue = unclass(v2))
  dx <- merge(d1,d2, by = "parameterId", suffixes = c("_1", "_2"), all = TRUE)
  dx[,`:=`(difference = estValue_2-estValue_1)]
  
  if (FLAGexcludeObsPars){
    pt <- petab_getParameterType(pe)
    dx <- merge(dx,pt, by =  "parameterId", all.x = TRUE, all.y = FALSE)
    dx <- dx[!parameterType %in% c("noiseParameters", "observableParameters")]
  }
  if (FLAGestimatedOnly) dx <- dx[parameterId %in% petab_getParametersToEstimate(pe)]
  
  pl <- cfggplot(dx, aes(parameterId, difference), FLAGbold = FALSE) + 
    lemon::facet_rep_grid(parameterType~., scales = "free", space = "free", repeat.tick.labels = TRUE) + 
    geom_col() + 
    coord_flip() +
    geom_blank()
  
  height <- nrow(dx) * 0.42 + 2
  cf_outputFigure(pl = pl, filename = filename, width = 15.5, height = height, scale = 1, units = "cm")
}





#' Plot resulting parameters from mstrust
#' 
#' This plot is inspired by the Plot from the Tutorial on modeling by Villaverde
#'
#' @param pd pd with field pd$result$mstrust
#' @param stepMax maximum step of waterfall
#' @param filename for plotting
#' @param i subset the data.table on their names c("parameterId", "parameterType", "fitrank", "step", "stepsize", 
#'                                                "index", "value", "converged", "iterations", "estValue")
#' @param ggCallback for plotting
#' @param ... output to [conveniencefunctions::cf_outputFigure()]
#'
#' @return ggplot
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family plotting
#' @importFrom data.table melt
#' @importFrom conveniencefunctions cfggplot cf_outputFigure
#'
#' @examples
pd_plotParsParallelLines <- function(pd, stepMax = 3, filename = NULL, i, ggCallback = NULL, ...) {
  
  si <- substitute(i)
  mi <- missing(i)
  
  parf <- pd$result$mstrust
  
  parameters <- attr(parf, "parameters")
  
  p <- as.data.table(as.data.frame(parf))
  p <- data.table::melt(p, measure.vars = parameters, variable.name = "parameterId", variable.factor = FALSE, value.name = "estValue", verbose = F)
  pt <- petab_getParameterType(pd$pe)
  p <- pt[p, on="parameterId"]
  p <- p[order(parameterType)]
  p[,`:=`(parameterId = factor(parameterId, unique(parameterId)))]
  p[,`:=`(parameterType = factor(parameterType, c("noiseParameters", "observableParameters", "other", "L1")))]
  p <- p[step <= stepMax]
  
  if (!mi) p <- p[eval(si)]
  
  pl <- conveniencefunctions::cfggplot(p, aes(parameterId, estValue, group = fitrank, color = parameterType)) + 
    facet_grid(step~parameterType, scales = "free_x", space = "free_x") + 
    geom_point(alpha = 0.1) + 
    geom_line( alpha = 0.1) + 
    guides(color = FALSE) +
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major.y = element_line(color="grey90"))
  
  for (plx in ggCallback) pl <- pl + plx
  
  conveniencefunctions::cf_outputFigure(pl, filename, 
                                        width = length(parameters) * 0.5 + 4, height = 21, units = "cm", ...)
  
  pl
  
}



#' Plot resulting parameters from mstrust
#' 
#' This plot is inspired by the Plot from the Tutorial on modeling by Villaverde
#'
#' @param pd pd with field pd$result$mstrust
#' @param stepMax maximum step of waterfall
#' @param filename for plotting
#' @param i subset the data.table on their names c("parameterId", "parameterType", "fitrank", "step", "stepsize", 
#'                                                "index", "value", "converged", "iterations", "estValue")
#' @param ggCallback for plotting
#' @param ... output to [conveniencefunctions::cf_outputFigure()]
#'
#' @return ggplot
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family plotting
#' @importFrom data.table melt
#' @importFrom conveniencefunctions cfggplot cf_outputFigure
#'
#' @examples
pd_plotParsParallelLines2 <- function(pd, stepMax = 3, filename = NULL, i, ggCallback = NULL, ...) {
  
  si <- substitute(i)
  mi <- missing(i)
  
  parf <- pd$result$mstrust
  
  parameters <- attr(parf, "parameters")
  
  p <- as.data.table(as.data.frame(parf))
  p <- data.table::melt(p, measure.vars = parameters, variable.name = "parameterId", variable.factor = FALSE, value.name = "estValue", verbose = F)
  pt <- petab_getParameterType(pd$pe)
  p <- pt[p, on="parameterId"]
  p <- p[order(parameterType)]
  p[,`:=`(parameterId = factor(parameterId, unique(parameterId)))]
  p[,`:=`(parameterType = factor(parameterType, c("noiseParameters", "observableParameters", "other", "L1")))]
  p <- p[step <= stepMax]
  
  if (!mi) p <- p[eval(si)]
  
  p[,`:=`(step = paste0(step,": ", stepsize))]
  p[,`:=`(step = factor(step))]
  p <- p[order(step, decreasing = TRUE)]
  
  
  
  pl <- conveniencefunctions::cfggplot(p, aes(parameterId, estValue, group = fitrank, color = step)) + 
    facet_grid(~parameterType, scales = "free_x", space = "free_x") + 
    geom_point(aes(alpha = step)) + 
    scale_color_cf() + 
    scale_alpha_manual(values = c(0.15,rep(0.08,20))) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major.y = element_line(color="grey95"),
          panel.grid.major.x = element_line(color="grey95")
          ) + 
    geom_blank()
  # Hack to draw lines in order: best step on top
  for (sx in sort(unique(p$step),decreasing = TRUE)) pl <- pl + geom_line(aes(alpha = step), data = p[step == sx])
  for (plx in ggCallback) pl <- pl + plx
  
  conveniencefunctions::cf_outputFigure(pl, filename, 
                                        width = length(parameters) * 0.5 + 5, height = 16, units = "cm", ...)
  
  pl
  
}





#' Title
#'
#' @param x
#' @param y
#' @param color
#' @param ...
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family plotting
#'
#' @examples
petab_plotHelpers_aeslist <- function(x = ~time,
                                      y = ~measurement,
                                      color = ~conditionId,
                                      fill = ~conditionId,
                                      linetype = ~parameterSetId,
                                      ymin = ~measurementmin,
                                      ymax = ~measurementmax,
                                      ...) {
  list(x = x, y = y, color = color, fill = fill, linetype = linetype, ymin = ymin, ymax = ymax, ...)
}


#' Title
#'
#' @param dplot
#' @param aeslist
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family plotting
#' @importFrom cOde getSymbols
#'
#' @examples
petab_plotHelpers_meanMeasurementsValue <- function(dplot, aeslist) {
  byvars <- lapply(aeslist, function(x) cOde::getSymbols(as.character(x)))
  byvars <- do.call(c, byvars)
  byvars <- setdiff(byvars, "measurement")
  byvars <- c(byvars, "observableId")
  byvars <- intersect(byvars, names(dplot))
  dmean <- copy(dplot)
  dmean <- dmean[,list(measurement = mean(measurement)), by = byvars]
}


#' Title
#'
#' @param pd
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family plotting
#'
#' @examples
petab_plotHelpers_variableOrder <- function(pd) {
  
  if (!is.null(pd$pe$meta$variableOrder)) return(variableOrder)
  
  states      <- pd$dModAtoms$symbolicEquations$reactions$states
  observables <- names(pd$dModAtoms$symbolicEquations$observables)
  
  c(observables, states)
}

#' Title
#'
#' @param pd 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family plotting
#' @importFrom dMod getParameters
#' @importFrom data.table data.table rbindlist
#'
#' @examples
petab_plotHelpers_parameterOrder <- function(pd) {
  pt <- petab_getParameterType(pd$pe)
  pt[,`:=`(parameterType = factor(parameterType, levels = c("other", "observableParameters", "noiseParameters")))]
  parameterOrder <- dMod::getParameters(pd$dModAtoms$symbolicEquations$reactions)
  parameterOrder <- data.table::data.table(parameterId = parameterOrder)
  parameterOrder <- parameterOrder[parameterId %in% pt$parameterId]
  pt <- data.table::rbindlist(list(pt[parameterOrder,on = .NATURAL][,`:=`(parameterOrder = 1)], pt[!parameterOrder,on = .NATURAL][,`:=`(parameterOrder = 2)]))
  pt <- pt[order(parameterOrder, parameterType)]
  pt[,`:=`(parameterOrder = 1:.N)]
  # subset to estimated
  pt <- pt[pd$pe$parameters[estimate == 1,list(parameterId)], on = .NATURAL][order(parameterOrder)]
  pt$parameterId
}


#' Title
#'
#' @param pd 
#' @param filename 
#' @param base_size 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family plotting
#' @importFrom dMod plotProfile
#' @importFrom conveniencefunctions cfggplot scale_color_cf theme_cf
#'
#' @examples
pd_plotProfile <- function(pd, filename = NULL, ggCallback = NULL) {
  dplot <- attr(dMod::plotProfile(pd$result$profiles), "data")
  dplot <- as.data.table(dplot)
  
  dplot <- dplot[mode %in% c("data", "prior", "total")]
  dplot[,`:=`(mode = factor(mode, levels = c("data", "prior", "total")))]
  
  nameOrder = petab_plotHelpers_parameterOrder(pd)
  nameOrder <- nameOrder[nameOrder %in% dplot$name]
  dplot[,`:=`(name = factor(name, levels = nameOrder))]
  
  dplot2 <- copy(dplot)
  dplot2 <- dplot2[mode == "data"]
  dplot2[,`:=`(isEndpoint = par %in% range(par)), by = "name"]
  dplot2 <- dplot2[isEndpoint == TRUE]
  dplot2[delta >= 3.84,`:=`(delta = NA)]
  
  dplot3 <- dplot[is.zero == TRUE & mode == "data"]
  
  pl <- conveniencefunctions::cfggplot(dplot, aes(par, delta)) + 
    facet_wrap_paginate(~name, nrow = 3, ncol = 4, scales = "free_x") + 
    geom_hline(aes(yintercept = yinter), data = data.frame(yinter = c(0, 3.84)), lty = 2, size = 0.2, color = "grey80") +
    geom_point(data = dplot2, size = 3, color = cfcolors[2], alpha = 0.7) + 
    geom_line(aes(color = mode, size = mode, linetype = mode)) + 
    geom_point(data = dplot3, size = 2, color = cfcolors[1]) + 
    scale_size_manual(values = c(1,0.5, 0.5)) + 
    scale_y_continuous(breaks=c(0, 1, 2.7, 3.84), labels = c("0", "68% / 1   ", "90% / 2.71", "95% / 3.84"), limits = c(-1, 5)) +
    conveniencefunctions::scale_color_cf() + 
    conveniencefunctions::theme_cf(base_size = 9) + 
    geom_blank()
  for (plx in ggCallback) pl <- pl + plx
  
  cf_outputFigure(pl, filename = filename, width = 29.7, height = 21, scale = 1, units = "cm")
}




# -------------------------------------------------------------------------#
# Simulate model ----
# -------------------------------------------------------------------------#

#' Subset predictions to observables and conditions observed in data
#'
#' @param pplot
#' @param dplot
#' @param NFLAGsubsetType
#'
#' @return
#' @export
#'
#' @examples
subsetPredictionToData <- function(pplot, dplot, NFLAGsubsetType = c(none = 0, strict = 1, keepInternal = 2, strict_cutTimes = 3,
                                                                     keepInternal_cutTimes = 3
)[2]) {
  pplot_out <- pplot
  if (NFLAGsubsetType == 0) return(pplot_out)
  
  # strict: Only combinations of observables and conditions which are present in data
  dplot_lookup <- copy(dplot)
  dplot_lookup <- dplot_lookup[,list(observableId, conditionId)]
  dplot_lookup <- unique(dplot_lookup)
  pplot_out <- pplot[dplot_lookup, on = c("observableId", "conditionId")]
  
  if (NFLAGsubsetType %in% c(3,4)) {
    dplot_lookup <- copy(dplot)
    dplot_lookup <- dplot_lookup[,list(maxtime = max(time)), by = c("observableId", "conditionId")]
    dplot_lookup <- unique(dplot_lookup)
    pplot_out <- pplot[dplot_lookup, on = c("observableId", "conditionId")]
    pplot_out <- pplot_out[time <= maxtime * 1.1]
    pplot_out[,`:=`(maxtime = NULL)]
  }
  
  if (NFLAGsubsetType %in% c(2,4)) {
    # keepInternal: Additionally keep ALL conditions of ALL internal states (with no data)
    pplot_keepInt <- pplot[!observableId %in% dplot$observableId]
    pplot_out <- rbindlist(list(pplot_out, pplot_keepInt))
  }
  
  pplot_out
}

# -------------------------------------------------------------------------#
# Test model ----
# -------------------------------------------------------------------------#


#' Some tests
#'
#' @param pd
#' @param page
#' @param cn
#'
#' @return
#' @export
#'
#' @importFrom dMod getDerivs plotPrediction
#' @importFrom ggforce facet_wrap_paginate
#'
#' @examples
pd_tests <- function(pd, page = 1, cn = 1, whichTests = c("plot" = 1, "objData" = 2, "derivs" = 3, "Rprof obj data" = 4)[2]) {
  
  cat("Number of implemented tests: 4")
  
  # Test prediction without derivs
  if (1 %in% whichTests){
    prediction <- pd$prd(objtimes(pd$pe$measurementData$time, 200), pd$pars)
    pl <- dMod::plotPrediction(prediction, name %in% pd$pe$observables$observableId) +
      ggforce::facet_wrap_paginate(~name, nrow = 4, ncol = 4, scales = "free", page = page)
    cat("\n===================================================", "\n")
    cat("01-plot: Plotting prediction page ", page, " / ", n_pages(pl), "\n")
    cat("===================================================", "\n")
    print(pl)
  }
  # Test obj
  if (2 %in% whichTests){
    cat("\n===================================================\n")
    cat("02-objData: Objective function\n")
    cat("===================================================\n")
    objval <- pd$objfns$obj_data(pd$pars, FLAGverbose = TRUE)
    print(objval)
  }
  # Test x for one condition
  if (3 %in% whichTests){
    pars <- pd$p(pd$pars)
    pars <- pars[[cn]]
    pred <- pd$dModAtoms$fns$x(objtimes(pd$pe$measurementData$time), pars)
    # Look at derivs
    derivs <- dMod::getDerivs(pred)
    cat("\n===================================================\n")
    cat("03-derivs: Derivs of x[", cn, "]\n")
    cat("===================================================\n")
    print(derivs[[1]][1:10, 1:10])
  }
  
  if (4 %in% whichTests){
    cat("\n===================================================\n")
    cat("04-summaryRprof of obj")
    cat("===================================================\n")
    rp <- tempfile()
    Rprof(rp)
    pd$obj(pd$pars)
    Rprof(NULL)
    rp <- summaryRprof(rp)
    print(rp$by.total)
    
  }
  
  
  
}


#' Mainly to have the code accessible
#'
#' @param pd
#' @param ID
#'
#' @return
#' @export
#'
#' @examples
pd_debug_p0 <- function(pd, ID = 1) {
  pars <- pd$pars
  fixed <- pd$fixed
  est.grid <- pd$dModAtoms$gridlist$est.grid
  fix.grid <- pd$dModAtoms$gridlist$fix.grid
  dummy <- make_pars(pars, fixed, est.grid, fix.grid, ID)
  pars_ <- dummy$pars
  fixed_ <- dummy$fixed
  
  cat("------------ Compare parameters ------------------------------------","\n",
      "x = pd$dModAtoms$fns$p0)", "\n", "y = c(names(pars_), names(fixed_))","\n",
      "setdiff(x,y) are MISSING in parameters\n"
  )
  conveniencefunctions::compare(getParameters(pd$dModAtoms$fns$p0), c(names(pars_), names(fixed_)))
  
  cat("------------ Parameter values ------------------------------------","\n")
  cat("------------ pars_ -------------","\n")
  print(pars_)
  cat("------------ fixed_ -------------","\n")
  print(fixed_)
  
  pd$dModAtoms$fns$p0(pars_, fixed = fixed_)
}







