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
    pd$result$fits <- conveniencefunctions::dMod_readMstrust(path)
  if (is.null(pd$results$profile) && dir.exists(file.path(path, "Results", "profile")))
    pd$result$profile <- conveniencefunctions::dMod_readProfiles(path)

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
  p <- dMod::P_indiv((pd$dModAtoms$fns$p1 * pd$dModAtoms$fns$p0), pd$dModAtoms$gridlist$est.grid, pd$dModAtoms$gridlist$fix.grid)

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
  pd$p        <- p
  pd$prd      <- prd
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

#' Title
#'
#' @param pd
#' @param fitrankRange
#'
#' @return [dMod::parframe()]
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
pd_parf_collectMstrust <- function(pd, fitrankRange = 1:15, tol = 1) {
  parf0 <- pd$result$fits[fitrankRange]
  fitidxs = cf_parf_getStepRepresentatives(parf0, tol = tol)
  pars <- parf0[fitidxs]
  pars <- data.table::as.data.table(pars)
  pars[,`:=`(parameterSetId = paste0("step", step, ",", "rank", fitrank))]
  pars <- dMod::parframe(pars, parameters = names(pd$pars), metanames = setdiff(names(pars), names(pd$pars)))
  pars
}


#' Title
#'
#' @param conveniencefunctions
#' @param pars2parframe
#' @param pd
#'
#' @return
#' @export
#'
#' @examples
pd_parf_collectPars <- function(pd, parameterSetId = "Base") {
  conveniencefunctions::pars2parframe(pd$pars,parameterSetId = parameterSetId, pd$obj_data)
}


#' Title
#'
#' @param pd
#' @param rows
#' @param parameters
#'
#' @return
#' @export
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
    parsEnd[,`:=`(parameterSetId = "profile_endpoints")]
    setcolorder(parsEnd, c(names(pars), "parameterSetId"))
    }
  if ("optimum" %in% rows){
    parsOpt <- pars[which(profileDirection == "optimum")[1]]
    parsOpt[,`:=`(parameterSetId = "optimum")]
    setcolorder(parsOpt, c(names(pars), "parameterSetId"))
    }

  pars <- rbindlist(list(parsEnd, parsOpt))
  # [ ] Idea: return whichPar as well. Then, when summarizing profiles in pd_predictAndPlot2 by min and max,
  #   one could have a look in which profile the minimum/maximum values occured, i.e. which profiles contribute most to prediction uncertainty
  #   Don't know if this is informative, but for now I'd copy the prediction code which would lead to a lot of code duplication which I don't want
  #   And I shouldn't be distracted from my main goals now
  pars <- pars[,.SD, .SDcols = c("parameterSetId", "profileDirection", "value", parameternames)]
  pars <- dMod::parframe(pars, parameters = parameternames, metanames = setdiff(names(pars), parameternames))
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
pd_parf_collect <- function(pd,
                            opt.base = pd_parf_opt.base(include = TRUE, parameterSetId = "Base"),
                            opt.mstrust = pd_parf_opt.mstrust(include = TRUE, fitrankRange = 1:20, tol = 1),
                            opt.profile = pd_parf_opt.profile(include = FALSE, rows = "profile_endpoints", parameters = NULL)) {

  parf_base <- parf_fit <- parf_profile <- NULL

  if (opt.base$include) {
    args <- c(list(pd = pd), opt.base[setdiff(names(opt.base), "include")])
    parf_base <- do.call(pd_parf_collectPars, args)}

  if (opt.mstrust$include && !is.null(pd$result$fits)) {
    args <- c(list(pd = pd), opt.mstrust[setdiff(names(opt.mstrust), "include")])
    parf_fit <- do.call(pd_parf_collectMstrust, args)}

  if (opt.profile$include) {
    args <- c(list(pd = pd), opt.profile[setdiff(names(opt.profile), "include")])
    parf_profile <- do.call(pd_parf_collectProfile, args)}

  parf <- cf_parf_rbindlist(list(parf_base, parf_fit, parf_profile))
  parf <- parf[order(parf$value)]

  parf
}

#' @export
pd_parf_opt.base <- function(include = TRUE, parameterSetId = "Base") {
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
  if (!identical(names(parlower), names(pd$pars))) stop("parameter names and boundaries do not match")
  parlower <- pd$pars - parlower
  parlower <- parlower[parlower <= tol]

  parupper <- petab_getParameterBoundaries(pd$pe, "upper")
  parupper <- -(pd$pars - parupper)
  parupper <- parupper[parupper <= tol]

  fixedOnBoundary <- pd$pars[c(names(parlower), names(parupper))]
  fixedOnBoundary
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


# fit_hierarchical -----
# could be implemented as objective function which automatically fits observable parameters
pd_normHierarchical <- function(pd){
pd <- petab_exampleRead("01", "pd")
# normIndiv_hierarchical <- function (pd) {
  force(pd)

  # Get Objects from pd
  est.grid <- data.table(pd$dModAtoms$gridlist$est.grid)
  fix.grid <- data.table(pd$dModAtoms$gridlist$fix.grid)
  setkeyv(est.grid, c("ID", "condition"))
  setkeyv(fix.grid, c("ID", "condition"))
  xp0      <- (pd$dModAtoms$fns$x * pd$dModAtoms$fns$p0)
  gxp0     <- (pd$dModAtoms$fns$g * pd$dModAtoms$fns$x * pd$dModAtoms$fns$p0)
  g0       <- pd$dModAtoms$fns$g
  data     <- pd$dModAtoms$data
  errmodel <- pd$dModAtoms$e

  # Parameter types
  pepa <- pd$pe$parameters
  pepa <- pepa[petab_getParameterType(pd$pe), on = "parameterId"]
  pepa <- pepa[estimate == 1]
  parametersObsErr <- pepa[parameterType != "other", parameterId]
  parametersOther  <- pepa[parameterType != "other", parameterId]
  parsObsErr <- pd$pars[parametersObsErr]

  # Salat from dMod
  timesD          <- pd$times
  x.conditions    <- est.grid$condition
  data.conditions <- names(data)
  e.conditions    <- names(attr(errmodel, "mappings"))
  controls        <- list(times = timesD, attr.name = attr.name, conditions = intersect(x.conditions, data.conditions))

  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = controls$conditions,
                   simcores = 1,
                   NFLAGbrowser = 0,
                   FLAGverbose = FALSE,
                   FLAGNaNInfwarnings = FALSE) {

    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pars <- arglist[[1]]

    # Logic
    # * pars_optObs <- c(parsObsErr); fixed <- c(pars, fixed)
    # * predictionXP0 <- xp0(times, pars_optObs)
    # * normObsParsErr <- function(parsObsParsErr, predictionXP0) norm(g(predictionXP0, pars_optObs), data)
    # * fit <- trust(normObsParsErr, parsObsErr)
    # * parsObsErr <<- fit$argument
    # * pars_optOuter <- c(pars, parsObsErr); fixed <- fixed
    # * norm(g*x*p0, pars_optOuter), data)


    cn <- (setNames(nm = conditions))[[1]]
    objlists <- lapply(setNames(nm = conditions), function(cn) {

      if (FLAGbrowser) browser()

      ID <- est.grid[condition == cn, ID]
      if (FLAGverbose) cat(ID, cn, "\n", sep = " ---- ")
      dummy <- make_pars(pars, fixed, est.grid, fix.grid, ID)
      pars_ <- dummy$pars
      fixed_ <- dummy$fixed

      if (!length(pars_)) return(init_empty_objlist(pars, deriv = deriv, FLAGchisquare = TRUE)) # No pars_ can happen if one fits only condition specific parameters and in this condition there are none

      prediction <- try(prd0(times = controls$times, pars = pars_, fixed = fixed_, deriv = deriv))

      if (inherits(prediction, "try-error"))
        stop("Prediction failed in \n>>>condition = ", cn, "\n>>>ID = ", ID, "\n\nTry iterating p(pars), (x*p)(pars), ... to find the problem.")

      prediction <- prediction[[1]]
      prediction <- check_and_sanitize_prediction(prediction, data, cn, FLAGNaNInfwarnings)

      err <- NULL
      if (any(is.na(data[[cn]]$sigma))) {
        err <- errmodel(out = prediction, pars = getParameters(prediction), conditions = cn, deriv=deriv)
        mywrss <- nll(res(data[[cn]], prediction, err[[1]]), deriv = deriv, pars = pars)
      } else {
        mywrss <- nll(res(data[[cn]], prediction), deriv = deriv, pars = pars)
      }

      if (deriv) mywrss <- renameDerivParsInObjlist(mywrss, dummy$parnames)

      mywrss
    })

    # Sum all objlists
    out <- Reduce("+", objlists)

    # Consider fixed: return only derivs wrt pouter
    out$gradient <- out$gradient[names(pars)]
    out$hessian <- out$hessian[names(pars), names(pars)]

    # Populate attributes
    attr(out, controls$attr.name) <- out$value
    ll_conditions <- data.frame(
      logl = vapply(setNames(objlists, conditions), function(.x) .x$value, 1),
      chi2 = vapply(setNames(objlists, conditions), function(.x) attr(.x, "chisquare"), 1))
    ll_sum <- data.frame(logl = sum(ll_conditions$logl),
                         chi2 = sum(ll_conditions$chi2))
    attributes(out) <- c(attributes(out), list(ll_cond_df = ll_conditions))
    attributes(out) <- c(attributes(out), list(ll_sum_df = ll_sum))
    # attr(out, "AIC") <- out$value + length(pars) * 2
    # attr(out, "BIC") <- out$value + length(pars) * log(nrow(as.data.frame(data)))
    return(out)
  }

  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- data.conditions
  attr(myfn, "parameters") <- attr(prd0, "parameters")
  attr(myfn, "modelname") <- modelname(prd0, errmodel)
  return(myfn)
}


# -------------------------------------------------------------------------#
# Profiles ----
# -------------------------------------------------------------------------#





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
  cf_outputFigure(pl, filename, width = width, height = height, scale = scale, units = units)
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

  cf_outputFigure(pl, filename = filename, width = width, height = height, scale = scale, units = units)
}



#' Title
#'
#' @param pd
#' @param pe
#' @param i
#' @param opt.base
#' @param opt.mstrust
#' @param opt.profile
#' @param NFLAGsubsetType
#' @param FLAGsummarizeProfilePredictions
#' @param FLAGmeanLine
#' @param aeslist
#' @param ggCallback
#' @param filename
#' @param FLAGfuture
#' @param width
#' @param height
#' @param scale
#' @param units
#'
#' @return
#' @export
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
                               i,
                               opt.base = pd_parf_opt.base(),
                               opt.mstrust = pd_parf_opt.mstrust(),
                               opt.profile = pd_parf_opt.profile(FALSE),
                               NFLAGsubsetType = c(none = 0, strict = 1, keepInternal = 2)[2],
                               FLAGsummarizeProfilePredictions = TRUE,
                               FLAGmeanLine = FALSE,
                               aeslist = petab_plotHelpers_aeslist(),
                               ggCallback = list(facet_wrap_paginate(~observableId, nrow = 4, ncol = 4, scales = "free"),
                                                 scale_y_continuous(n.breaks = 5)),
                               opt.gg = list(ribbonAlpha = 0.2), # would be nice to put this into opt.profile or maybe opt.gg?
                               filename = NULL, FLAGfuture = TRUE,
                               width = 29.7, height = 21, scale = 1, units = "cm",
                               FLAGreturnPlotData = FALSE
) {


  # .. Catch i (see petab_mutateDCO for more ideas) -----
  mi <- missing(i)
  si <- substitute(i)

  # .. Data -----
  # observableTransformation
  dplot <- petab_joinDCO(pe)
  dplot[,`:=`(measurement = eval(parse(text = paste0(observableTransformation, "(", measurement, ")")))), by = 1:nrow(dplot)]
  dplot[,`:=`(observableId=factor(observableId, petab_plotHelpers_variableOrder(pd)))]


  # .. Prediction -----
  parf <- pd_parf_collect(pd, opt.base = opt.base, opt.mstrust = opt.mstrust, opt.profile = opt.profile)
  if (nrow(parf) > 5) {
    pd$times <- pd_predtimes(pd, N = 60)
    if (!opt.profile$include) cat("Predicting for more than 5 parameter sets. Are you sure?")
  }
  pplot <- conveniencefunctions::cf_predict(prd = pd$prd, times = pd$times, pars = parf, fixed = pd$fixed)
  setnames(pplot, c("condition"  , "name"        , "value"), c("conditionId", "observableId", "measurement"))
  pplot[,`:=`(observableId=factor(observableId,  petab_plotHelpers_variableOrder(pd)))]
  pplot <- subsetPredictionToData(pplot, dplot, NFLAGsubsetType = NFLAGsubsetType)

  # .. Error model / prediction ribbon -----
  pplotRibbon <- NULL
  if (FLAGsummarizeProfilePredictions && opt.profile$include) {
    # pd_plotHelpers_summarizeProfilePredictions - in here not as separate function because of multiple return elements
    pplotRibbon <- pplot[profileDirection %in% c("left", "right")]
    pplotRibbon <- pplotRibbon[,list(
      measurementmin = min(measurement),
      measurementmax = max(measurement)
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
  pl <- pl + geom_point(do.call(aes_q, aeslist[intersect(names(aeslist), conveniencefunctions::cfgg_getAllAesthetics()[["geom_point"]])]), data = dplot)
  pl <- pl + geom_line( do.call(aes_q, aeslist[intersect(names(aeslist), conveniencefunctions::cfgg_getAllAesthetics()[["geom_line"]])]) , data = pplot)
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



#' Title
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
                                      filename = NULL, width = 29.7, height = 21, scale = 1, units = "cm"
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

  cf_outputFigure(pl, filename = filename, width = width, height = height, scale = scale, units = units)

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
#'
#' @examples
petab_plotHelpers_variableOrder <- function(pd) {

  if (!is.null(pd$pe$meta$variableOrder)) return(variableOrder)

  states      <- pd$dModAtoms$symbolicEquations$reactions$states
  observables <- names(pd$dModAtoms$symbolicEquations$observables)

  c(observables, states)
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
subsetPredictionToData <- function(pplot, dplot, NFLAGsubsetType = c(none = 0, strict = 1, keepInternal = 2)[2]) {
  if (NFLAGsubsetType == 0) return(pplot)

  # strict: Only combinations of observables and conditions which are present in data
  dplot_lookup <- copy(dplot)
  dplot_lookup <- dplot_lookup[,list(observableId, conditionId)]
  dplot_lookup <- unique(dplot_lookup)
  pplot_out <- pplot[dplot_lookup, on = c("observableId", "conditionId")]

  # keepInternal: Additionally keep ALL conditions of ALL internal states (with no data)
  if (NFLAGsubsetType == 2) {
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
  objval <- pd$obj_data(pd$pars)
  cat("\n===================================================\n")
  cat("02-objData: Objective function\n")
  cat("===================================================\n")
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
  cat("04-summaryRprof of obj_data")
  cat("===================================================\n")
  rp <- tempfile()
  Rprof(rp)
  pd$obj_data(pd$pars)
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
  pd$dModAtoms$fns$p0(pars_, fixed = fixed_)
}



