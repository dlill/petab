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
                            opt.profile = list(include = FALSE, rows = "profile_endpoints", parameters = NULL)) {

  # [ ] refactor: each opt.* should get its own collector function: collect_opt.base(include, parameterSetId), collect_opt.fit(...)
  # [ ] refactor: each collect step should get individual function
  parf_base <- parf_fit <- NULL

  if (opt.base$include) {
    args <- c(list(pd = pd), opt.base[setdiff(names(opt.base), "include")])
    parf_base <- do.call(pd_parf_collectPars, args)}

  if (opt.mstrust$include && !is.null(pd$result$fits)) {
    args <- c(list(pd = pd), opt.mstrust[setdiff(names(opt.mstrust), "include")])
    parf_fit <- do.call(pd_parf_collectMstrust, args)}

  if (opt.profile$include) cat("collect profiles is not implemented yet")

  parf <- cf_parf_rbindlist(list(parf_base, parf_fit))
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


# .. fit_hierarchical -----
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
pd_predictAndPlot <- function(pd, i, opt.base = pd_parf_opt.base(), opt.mstrust = pd_parf_opt.mstrust(),
                              FLAGsubsetPredictionToData = TRUE,
                              FLAGsummarizeProfilePredictions = TRUE,
                              nrow = 3, ncol = 4, scales = "free", page = 1,
                              n.breaks = 5,
                              filename = NULL, width = 29.7, height = 21, scale = 1, units = "cm") {
  # Catch i (see petab_mutateDCO for more ideas)
  mi <- missing(i)
  si <- substitute(i)

  # .. Collect parameters and predict -----
  parf <- pd_parf_collect(pd, opt.base = opt.base, opt.mstrust = opt.mstrust)
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
  if (FLAGsubsetPredictionToData) pplot <- subsetPredictionToData(pplot, dplot)

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



# -------------------------------------------------------------------------#
# Simulate model ----
# -------------------------------------------------------------------------#

#' Subset predictions to observables and conditions observed in data
#'
#' @param pplot
#' @param dplot
#'
#' @return
#' @export
#'
#' @examples
subsetPredictionToData <- function(pplot, dplot) {
  dplot_lookup <- copy(dplot)
  dplot_lookup <- dplot_lookup[,list(observableId, conditionId)]
  dplot_lookup <- unique(dplot_lookup)
  pplot <- pplot[dplot_lookup, on = c("observableId", "conditionId")]
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
pd_tests <- function(pd, page = 1, cn = 1, whichTests = c("01-plot", "02-objData", "03-derivs")[1:2]) {

  cat("The following tests are implemented: \n", 'c("01-plot", "02-objData", "03-derivs")')

  # Test prediction without derivs
  if ("01-plot" %in% whichTests){
    prediction <- pd$prd(objtimes(pd$pe$measurementData$time, 200), pd$pars)
    pl <- dMod::plotPrediction(prediction, name %in% pd$pe$observables$observableId) +
      ggforce::facet_wrap_paginate(~name, nrow = 4, ncol = 4, scales = "free", page = page)
    cat("\n===================================================", "\n")
    cat("01-plot: Plotting prediction page ", page, " / ", n_pages(pl), "\n")
    cat("===================================================", "\n")
    print(pl)
  }
  # Test obj
  if ("02-objData" %in% whichTests){
  objval <- pd$obj_data(pd$pars)
  cat("\n===================================================\n")
  cat("02-objData: Objective function\n")
  cat("===================================================\n")
  print(objval)
  }
  # Test x for one condition
  if ("03-derivs" %in% whichTests){
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
}


