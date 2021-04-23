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
#' @importFrom dMod loadDLL
#' @importFrom conveniencefunctions dMod_readProfiles dMod_readMstrust
#' @examples
readPd <- function(filename) {
  # 0 wd wrangling
  wd <- getwd()
  on.exit({setwd(wd)})
  # 1 Read RDS
  pd <- readRDS(filename)
  # 2 Load DLLs
  setwd(dirname(filename))
  dMod::loadDLL(pd$obj_data)
  setwd(wd)

  # 3 Read potential results if not yet in pd
  # If in pd already, some filenameParts variable might have been derived from them,
  #   so it would be dangerous to reload them if they were postprocessed
  path <- dirname(dirname(filename))
  if (is.null(pd$results$fits) && dir.exists(file.path(path, "Results", "mstrust")))
    pd$results$fits <- conveniencefunctions::dMod_readMstrust(path)
  if (is.null(pd$results$profile) && dir.exists(file.path(path, "Results", "profile")))
    pd$results$profile <- conveniencefunctions::dMod_readProfiles(path)

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
                                 paste0(".obsparsfitted")))
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
  if (!is.null(objtimes))    dMod::controls(pd$obj_data, "times") <- objtimes

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
  p <- dMod::P_indiv(pd$dModAtoms$fns$p0, pd$dModAtoms$gridlist$est.grid, pd$dModAtoms$gridlist$fix.grid)

  # Rebuild high-level prediction function
  prd0 <- Reduce("*", pd$dModAtoms$fns)
  prd <- PRD_indiv(prd0, pd$dModAtoms$gridlist$est.grid, pd$dModAtoms$gridlist$fix.grid)

  # Rebuild obj_data
  tobj <- if (!is.null(pd$obj_data)) dMod::controls(pd$obj_data, name = "times") else dMod::objtimes(pd$pe$measurementData$time, Nobjtimes = Nobjtimes)
  obj_data <- normL2_indiv(pd$dModAtoms$data, prd0,
                           pd$dModAtoms$e,
                           est.grid = pd$dModAtoms$gridlist$est.grid,
                           fix.grid = pd$dModAtoms$gridlist$fix.grid,
                           times = tobj)

  # Update p, prd and obj_data
  pd$p <- p
  pd$prd <- prd
  pd$obj_data <- obj_data

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
#'
#' @examples
pd_objtimes <- function(pd, N = 100) {
  datatimes <- unique(sort(pd$pe$measurementData$time))
  eventtimes <- NULL
  dMod::objtimes(datatimes,eventtimes,N)
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
#'
#' @examples
pd_updateEstPars <- function(pd, parsEst, FLAGupdatePE = TRUE, FLAGsavePd = FALSE) {
  pd$pars[names(parsEst)] <- parsEst
  if (FLAGupdatePE) {
    cat("pd$pe pars have been *set*")
    petab_setPars_estScale(pd$pe, parsEst)}
  if (FLAGsavePd) writePd(pd)
  pd
}



# -------------------------------------------------------------------------#
# Fitting functions ----
# -------------------------------------------------------------------------#

#' Run a fit only for observation parameters
#'
#' @param pd
#' @param NFLAGsavePd as in pd_importIndiv: 0: don't save, 1: save, but redo everything, 3: if it was done before and pd didn't change, just return the pd
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
pd_fitObsPars <- function(pd, NFLAGsavePd = 3) {
  logfile <- pd_files(pd$filenameParts)$obsParsFitted
  if (file.exists(logfile) &&
      !inputFileChanged(pd_files(pd$filenameParts)$rdsfile, logfile) &&
      NFLAGsavePd == 3)
    return(pd)

  obspars <- petab_getParameterType(pd$pe)
  obspars <- obspars[parameterType %in% c("observableParameters", "noiseParameters"), parameterId]

  fit_par <- pd$pars[obspars]
  fit_fix <- pd$pars[setdiff(names(pd$pars), obspars)]

  parlower <- petab_getParameterBoundaries(pd$pe, "lower")[obspars]
  parupper <- petab_getParameterBoundaries(pd$pe, "upper")[obspars]

  fit <- trust(pd$obj_data, fit_par, 1,10, iterlim = 1000, fixed = fit_fix, parlower = parlower, parupper = parupper)

  pd <- pd_updateEstPars(pd, parsEst = fit$argument, FLAGupdatePE = TRUE, FLAGsavePd = NFLAGsavePd > 0)
  if (NFLAGsavePd > 0) writeLines("obsPars fitted", logfile)
  pd
}


# idea: fit_hierarchical
# could be implemented as objective function which automatically fits observable parameters

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
#'
#' @examples
#' pdx <- petab_exampleRead("01", "pd")
#' pd_plot(pdx)
pd_plot <- function(pd, ..., page = 1, nrow = 3, ncol = 4, filename = NULL, width = 29.7, height = 21, scale = 1, units = "cm"){
  pred <- pd$prd(pd$times, pd$pars)
  pl <- plotCombined(pred, pd$dModAtoms$data, ...) +
    facet_wrap_paginate(~name, nrow = nrow, ncol = ncol, scales = "free", page = page) +
    theme_cf() +
    scale_color_cf()
  if (!is.null(filename)) cf_outputFigure(pl, filename, width = width, height = height, scale = scale, units = units)
  pl
}


# -------------------------------------------------------------------------#
# Simulate model ----
# -------------------------------------------------------------------------#


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
pd_tests <- function(pd, page = 1, cn = 1) {

  # Test prediction without derivs
  prediction <- pd$prd(objtimes(pd$pe$measurementData$time, 200), pd$pars)
  pl <- dMod::plotPrediction(prediction, name %in% pd$pe$observables$observableId) +
    ggforce::facet_wrap_paginate(~name, nrow = 4, ncol = 4, scales = "free", page = page)
  cat("\n===================================================", "\n")
  cat("Plotting prediction page ", page, " / ", n_pages(pl), "\n")
  cat("===================================================", "\n")
  print(pl)

  # Test obj
  objval <- pd$obj_data(pd$pars)
  cat("\n===================================================\n")
  cat("Objective function\n")
  cat("===================================================\n")
  print(objval)

  # Test x for one condition
  pars <- pd$p(pd$pars)
  pars <- pars[[cn]]
  pred <- pd$dModAtoms$fns$x(objtimes(pd$pe$measurementData$time), pars)
  # Look at derivs
  derivs <- dMod::getDerivs(pred)
  cat("\n===================================================\n")
  cat("Derivs of x[", cn, "]\n")
  cat("===================================================\n")
  print(derivs[[1]][1:10, 1:10])
}


