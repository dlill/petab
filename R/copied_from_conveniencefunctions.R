# -------------------------------------------------------------------------#
# Predict ----
# -------------------------------------------------------------------------#

#' predict.prdfn with options bowser and verbose
#'
#' @param times
#' @param pars parframe
#' @param ...
#' @param FLAGverbose
#' @param FLAGbrowser
#' @param prd 
#' @param keep_names 
#' @param FLAGverbose2 
#'
#' @return data.table
#' @export
#'
#' @importFrom data.table rbindlist data.table 
#' @importFrom purrr imap
#' @importFrom dMod as.parvec
#'
#' @examples
cf_predict <- function (prd, times, pars, keep_names = NULL, ncores = 4, FLAGverbose = FALSE, FLAGverbose2 = FALSE, FLAGbrowser = FALSE,deriv = FALSE, ...) {
  if (FLAGverbose2) cat("Simulating", "\n")
  
  i <- 1
  out <- parallel::mclapply(X = 1:nrow(pars), mc.cores = ncores, FUN = function(i) {
    if (FLAGverbose) cat("Parameter set", i, "\n")
    if (FLAGbrowser) browser()
    mypar <- dMod::as.parvec(pars[i],1)
    prediction <- try(prd(times, mypar, deriv = deriv, ...))
    if (inherits(prediction, "try-error")) {
      warning("parameter set ", i, " failed\n")
      return(NULL)
    }
    prediction <- purrr::imap(prediction, function(.x,.y){
      .x <- data.table(.x)
      if (!is.null(keep_names))
        .x[, (setdiff(names(.x), c(keep_names, "time"))) := NULL]
      .x[, `:=`(condition = .y, parframe_rowid = i)]
      .x
    })
    melt(rbindlist(prediction), variable.name = "name", value.name = "value", id.vars = c("time", "condition", "parframe_rowid"))
  })
  if (FLAGverbose2) cat("postprocessing", "\n")
  out <- rbindlist(out[!is.null(out)])
  
  
  # Make this available as a FLAG
  pars <- cf_parf_getMeta(pars)
  if (!is.null(pars)){
    pars <- data.table(pars)
    pars[, `:=`(parframe_rowid = 1:.N)]
    out <- merge(pars, out, by = "parframe_rowid")
    out$parframe_rowid <- NULL
  }
  out
}



# ---------------------------------------------------------- #
# Parframe-class ----
# ---------------------------------------------------------- #

#' Add step column and fitrank
#'
#' @param myparframe a parframe
#' @param tol integer for steps.
#'
#' @return the parframe with columns fitrank and step
#' @export
#' @importFrom purrr map_dbl
add_stepcolumn <- function(myparframe, tol = 1) {
  steps <- dMod:::stepDetect(myparframe$value, tol)
  bla <- 1:nrow(myparframe)
  stepcol <- cumsum(bla%in%steps)
  
  fitrank <- 1:length(stepcol)
  stepsize <- purrr::map_dbl(stepcol, ~sum(stepcol == .x)) 
  mydf <- as.data.frame(myparframe)
  mydf <- mydf[!names(mydf)%in%c("fitrank", "step", "stepsize")]
  mydf <- cbind(fitrank = fitrank, step  = stepcol, stepsize = stepsize, mydf)
  
  return(parframe(mydf, parameters = attr(myparframe, "parameters")))
}

#' get Parameter names of a parframe
#' @param x parframe
#' @export
cf_parf_parNames <- function(x) {
  attr(x, "parameters")
}

#' Title
#'
#' @param pars 
#'
#' @return data.frame of meta columns
#' @export
cf_parf_getMeta <- function(pars){
  pars <- as.data.frame(pars)[cf_parf_metaNames(pars)]
  if (length(pars) == 0)
    return(NULL)
  pars <- cbind(pars, parframe_rowid = 1:nrow(pars))
  if ("value" %in% names(pars))
    pars <- dplyr::rename(pars, objvalue = value)
  pars
}

#' Title
#'
#' @param pars 
#'
#' @return
#' @export
cf_parf_metaNames <- function(pars){
  setdiff(names(pars), attr(pars, "parameters"))
}


#' @export
cf_parf_metaNames0 <- list(
  mstrust = c("index", "value", "converged", "iterations", "fitrank", "step", "stepsize"),
  other = c("AIC", "BIC",  "valueData", "valueObj"),
  profile = c("constraint", "stepsize", "gamma", "whichPar", "value"),
  l1 = c("value", "converged", "iterations", "lambda")
)



#' Select parameter columns of parframe
#'
#' @param parf parframe
#' @param parameters keine ahnung mehr?
#'
#' @export
cf_parf_getPars <- function(parf) {
  as.data.frame(parf)[attr(pars, "parameters")] 
}



#' Better as_parframe
#' 
#' Adds AIC and BIC automatically, adds stepcolumn automatically
#'
#' @param x 
#' @param sort.by 
#' @param ... 
#'
#' @return
#' @export
#' @md
#' @importFrom data.table rbindlist
cf_as.parframe <- function (x, sort.by = "value", ...) {
  m_stat <- dMod:::stat.parlist(x)
  m_metanames <- c("index", "value", "converged", "iterations")
  m_idx <- which("error" != m_stat)
  m_parframe <- data.frame(index = m_idx, 
                           value = vapply(x[m_idx], function(.x) .x$value, 1), 
                           converged = vapply(x[m_idx], function(.x) .x$converged, TRUE), 
                           iterations = vapply(x[m_idx], function(.x) as.integer(.x$iterations), 1L))
  
  if (!is.null(attr(x[[m_idx[[1]]]], "BIC"))){
    m_parframe <- cbind(m_parframe, 
                        AIC = vapply(x[m_idx], function(.x) attr(.x, "AIC"), 1),
                        BIC = vapply(x[m_idx], function(.x) attr(.x, "BIC"), 1))
    m_metanames <- c(m_metanames, c("AIC", "BIC"))
  }
  
  parameters <- lapply(x[m_idx], function(x) as.data.table(as.list(x$argument)))
  parameters <- data.table::rbindlist(parameters, use.names = TRUE)
  m_parframe <- cbind(m_parframe, parameters)
  
  m_parframe <- m_parframe[order(m_parframe[sort.by]), ]
  
  cf_parframe(m_parframe, parameters = names(x[[m_idx[1]]]$argument), 
              metanames = m_metanames)
}

#' Improved version of parframe
#' 
#' Fits coerces to data.frame and guesses metanames. Adds stepcolumn automatically
#' 
#' @param x 
#' @param parameters 
#' @param metanames 
#' @param obj.attributes 
#' @param tol 
#'
#' @return
#' @export
#'
#' @examples
cf_parframe <- function(x = NULL, parameters = NULL, metanames = NULL, 
                        obj.attributes = NULL, tol = 1) {
  x <- as.data.frame(x) 
  if (is.null(metanames))
    metanames <- intersect(names(x), Reduce(union, cf_parf_metaNames0))
  if (is.null(parameters))
    parameters <- setdiff(names(x),metanames)
  x <- dMod::parframe(x,parameters, metanames, obj.attributes)
  if ("converged" %in% metanames & !"fitrank" %in% metanames)
    x <- add_stepcolumn(x, tol)
  x
}


#' Title
#'
#' @param parf 
#' @param tol 
#'
#' @return
#' @export
#'
#' @examples
cf_parf_getStepRepresentatives <- function(parf, tol = 1) {
  which(as.logical(c(1, diff(parf$step) != 0)))
}


#' Turn a pars vector into a single-row parframe
#'
#' @param obj Objective function like normL2
#' @param pars setNames(outervalues, parnames)
#' @param parameterSetId Identifier for this parameterset
#'
#' @return parframe
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
#' obj <- function(pars) {list(value = sum(pars^2))}
#' pars = c(a = 1, b = 2)
#' pars2parframe(pars, "base", obj)
pars2parframe <- function(pars, parameterSetId = "Base", obj = NULL) {
  value <- if (!is.null(obj)) obj(pars)$value else NA
  parf0 <- data.frame(parameterSetId = parameterSetId, value = value, as.data.frame(as.list(pars)), index = 0, converged = FALSE)
  cf_parframe(parf0, metanames = c("parameterSetId", "value", "index", "converged"))
}


#' Rbind some parframes
#'
#' @param parflist list(parframes)
#'
#' @return parframe
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom data.table rbindlist
#'
#' @examples
#' parflist <- list(cf_parframe(data.frame(value = 1, a = 1, b = 2), metanames = "value"),
#'                  cf_parframe(data.frame(value = 1, a = 2, parameterSetId = "par2"), metanames = c("value", "parameterSetId")))
#' parflist <- list(cf_parframe(data.frame(value = 1, a = 1), metanames = "value"),
#'                  cf_parframe(data.frame(value = 1, a = 2, parameterSetId = "par2"), metanames = c("value", "parameterSetId")))
#' cf_parf_rbindlist(parflist)
cf_parf_rbindlist <- function(parflist) {
  
  metanames <- lapply(parflist, function(x) attr(x, "metanames"))
  metanames <- do.call(c, metanames)
  metanames <- unique(metanames)
  
  parnames <- lapply(parflist, function(x) attr(x, "parameters")) # could implement check to see that parameter names are the same
  parnames <- do.call(c, parnames)
  parnames <- unique(parnames)
  
  mixedCol <- intersect(metanames, parnames)
  if (length(mixedCol)) stop("The following columns are parameters and metanames: ", paste0(mixedCol, collapse = ", "))
  
  parf <- data.table::rbindlist(parflist, use.names = TRUE, fill = TRUE)
  parf <- cf_parframe(parf, parameters = parnames, metanames = metanames)
  parf
}


# -------------------------------------------------------------------------#
# Files ----
# -------------------------------------------------------------------------#
#' Consistent dMod filenames
#'
#' @param path 
#' @param identifier 
#'
#' @return
#' @export
#'
#' @examples
dMod_files <- function(path, identifier = "") {
  list(
    mstrust   = file.path(path, "Results", "mstrust", paste0("mstrustList-",identifier,".rds")),
    profile   = file.path(path, "Results", "profile", paste0("profiles-",identifier,".rds")),
    L1        = file.path(path, "Results", "L1", paste0("L1-",identifier,".rds")),
    PPL   = file.path(path, "Results", "PPL", paste0("PPL-",identifier,".rds")),
    petabdMod = file.path(path, paste0("pd",".rds"))
  )
}



#' Consistently save profiles
#'
#' @param profiles 
#' @param path 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
dMod_saveProfiles <- function(profs, path) {
  profsplit <- split(profs, profs$whichPar)
  lapply(profsplit, function(prof) {
    whichNm <- as.character(unique(prof$whichPar))
    filename <- dMod_files(path, whichNm)$profile
    dir.create(dirname(filename), FALSE, TRUE)
    saveRDS(prof, filename)
    prof
  })
  invisible()
}

#' Title
#'
#' @param path 
#'
#' @return Profiles
#' @export
#'
#' @examples
dMod_readProfiles <- function(path = .outputFolder) {
  profPath <- file.path(path, "Results", "profile")
  files <- list.files(profPath, full.names =  TRUE)
  profs <- lapply(files, readRDS)
  
  hasError <- vapply(profs, function(x) inherits(x, "tre-error"), FALSE)
  if (any(hasError)) 
    cat("The following profiles have errors: \n",
        "* ", paste0(basename(files[hasError]), collapse = "\n* "))
  
  profs <- profs[!hasError]
  profs <- do.call(rbind, profs)
  profs
}


#' Title
#'
#' @param fit parlist or parframe, will be coerced to parframe with [cf_as.parframe()]
#' @param path,identifier result will be written to [dMod_files()](path, identifier)$mstrust
#' @param FLAGoverwrite TRUE or FALSE
#'
#' @return fit as parframe, invisibly
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
dMod_saveMstrust <- function(fit, path, identifier = "1", FLAGoverwrite = FALSE) {
  if (!is.parframe(fit)) fit <- cf_as.parframe(fit)
  
  filename <- dMod_files(path, identifier)$mstrust
  if (!FLAGoverwrite && file.exists(filename)) 
    stop("FLAGoverwrite is FALSE and file.exists: ", filename)
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  saveRDS(fit, filename)
}


#' Title
#'
#' @param path 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
dMod_readMstrust <- function(path, identifier = NULL) {
  filename <- dMod_files(path, identifier)$mstrust
  fits <- readRDS(filename)
  fits
}


#' Title
#'
#' @param fit 
#' @param path 
#' @param identifier 
#' @param FLAGoverwrite 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
dMod_saveL1 <- function(L1, path, identifier = "1", FLAGoverwrite = FALSE) {
  filename <- dMod_files(path, identifier)$L1
  if (!FLAGoverwrite && file.exists(filename)) 
    stop("FLAGoverwrite is FALSE and file.exists: ", filename)
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  saveRDS(L1, filename)
}


#' Title
#'
#' @param path 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
dMod_readL1 <- function(path) {
  filename <- dMod_files(path)$L1
  fits <- list.files(dirname(filename), "rds$", full.names = TRUE)
  fits <- lapply(fits, readRDS)
  fits <- do.call(rbind, fits)
  fits
}

# -------------------------------------------------------------------------#
# Profiles ----
# -------------------------------------------------------------------------#


#' getProfiles with optimum below the fit$argument
#' 
#' Sometimes, a profile finds a better optimum
#' Get only profiles of parameters which show this behaviour
#' 
#' @param profiles parframe
#'
#' @return data.table
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
cf_profile_getProfilesOptimumBelowFit <- function(profiles){
  dp <- as.data.table(profiles)
  value0 <- dp[constraint == 0, unique(value)]
  parsBelow <- dp[value < value0, unique(whichPar)]
  
  profiles <- profiles[profiles$whichPar %in% parsBelow]
  profiles
}


#' Title
#'
#' @param profiles dMod::profile
#' @param value_name Choose from numerical attributes of obj, e.g. "value", "data", "prior" ...
#'
#' @return data.table(whichPar, profileDirection, isOpen)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family profiles
#'
#' @examples
cf_profile_getOpenProfiles <- function(profiles, value_column = "data") {
  dp <- as.data.table(profiles)
  dp[,`:=`(profileDirection = if(constraint <= 0) "left" else "right"), by = 1:nrow(dp)]
  dp[,`:=`(VALUECOLUMN = eval(parse(text = value_column)))]
  valueOptimum <- dp[constraint == 0, unique(VALUECOLUMN)]
  dp[,`:=`(valueOptimum = valueOptimum)]
  dp[,`:=`(valueDifference = VALUECOLUMN - valueOptimum)]
  
  dpOpen <- dp[,list(isOpen = all(valueDifference < 3.84)), by = c("whichPar", "profileDirection")]
  dpOpen <- dpOpen[isOpen == TRUE]
  dpOpen
}



#' Prepare parameters affected by a profile
#' 
#' Its not helpful to plot all paths in plotPaths, 
#' but only of those parameters which are affected
#'
#' @param profiles Profiles from dMod::profile
#' @param tol tolerance to measure if a parameter is affected
#' @param FLAGnormalizeYParameters Let them all pass through the same point: subtract the value at constraint = 0
#'
#' @return data.table(value, constraint, PARAMETER1, PARAMETER2, PARVALUE1, PARVALUE2)]
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
cf_profile_prepareAffectedPaths <- function(profiles, tol = 1e-1, FLAGnormalizeYParameters = TRUE,SDcols = c("value", "constraint", "PARAMETER1", "PARAMETER2", "PARVALUE1", "PARVALUE2")) {
  parnames <- cf_parf_parNames(profiles)
  dp <- data.table(profiles)
  dp <- split(dp, dp$whichPar)
  d <- (dp)[[1]]
  dp <- lapply(dp, function(d) {
    wP <- unique(d$whichPar)
    d[,`:=`(PARAMETER1 = whichPar)]
    d[,`:=`(PARVALUE1  = .SD[[1]]), .SDcols = wP]
    d <- melt(d, measure.vars = setdiff(parnames, wP), 
              variable.name = "PARAMETER2", variable.factor = FALSE, 
              value.name = "PARVALUE2")
    d <- d[,.SD, .SDcols = SDcols ]
    d <- d[,`:=`(AFFECTED = diff(range(PARVALUE2)) > tol), by = c("PARAMETER1", "PARAMETER2")]
    d <- d[AFFECTED == TRUE]
    d[,`:=`(AFFECTED = NULL)]
    d
  })
  dp <- rbindlist(dp, use.names = TRUE)
  if (FLAGnormalizeYParameters) dp[,`:=`(PARVALUE2 = PARVALUE2-PARVALUE2[which.min(abs(constraint))]), by = c("PARAMETER1", "PARAMETER2")]
  dp
}

#' Title
#'
#' @param profiles 
#'
#' @return
#' @export
#'
#' @examples
cf_profile_getStartPar <- function(profiles) {
  parnames <- cf_parf_parNames(profiles)
  profopt <- profiles[profiles$constraint == 0, parnames]
  profopt <- as.data.frame(profopt)
  profopt <- unique(profopt)
  profopt <- unlist(profopt)
  profopt <- data.table(PARAMETEROPT = names(profopt), PARVALUEOPT = profopt)
  profopt
}

#' Title
#'
#' @param profiles 
#'
#' @return
#' @export
#'
#' @examples
cf_profile_getOptimum <- function(profiles) {
  parnames <- cf_parf_parNames(profiles)
  profopt <- profiles[which.min(profiles$value), parnames]
  profopt <- as.data.frame(profopt)
  profopt <- unique(profopt)
  if (nrow(profopt) > 1) {
    warning("Multiple parameter sets with optimum value. Taking the first one")
    profopt <- profopt[1,]
  }
  profopt <- unlist(profopt)
  profopt <- data.table(PARAMETEROPT = names(profopt), PARVALUEOPT = profopt)
  profopt
}

#' Title
#'
#' @param profiles profiles to plotPaths
#' @param page pagination page
#' @param tol 
#'
#' @return paginated ggplot
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' 
cf_profile_plotPathsAffected <- function(profiles, tol = 1e-1, 
                                         FLAGnormalizeYParameters = TRUE,
                                         nrow = 4,ncol = 5,
                                         ggCallback = list(),
                                         ...
) {
  
  dp <- cf_profile_prepareAffectedPaths(profiles, tol = tol, FLAGnormalizeYParameters = FLAGnormalizeYParameters)
  dfit <- dp[cf_profile_getStartPar(profiles), on = c("PARAMETER1" = "PARAMETEROPT", PARVALUE1 = "PARVALUEOPT")]
  dfit <- dfit[!is.na(PARVALUE2)]
  dopt <- dp[cf_profile_getOptimum(profiles), on = c("PARAMETER1" = "PARAMETEROPT", PARVALUE1 = "PARVALUEOPT")]
  dopt <- dopt[!is.na(PARVALUE2)]
  
  pl <- cfggplot(dp, aes(PARVALUE1, PARVALUE2, group = PARAMETER2, color = PARAMETER2)) + 
    facet_wrap_paginate(~PARAMETER1, nrow = nrow, ncol = ncol, scales = "free", page = 1)+
    geom_line() + 
    geom_point(data = dfit, size = 2) + 
    geom_point(data = dopt, shape = 4, size = 2) +
    scale_color_cf() + 
    guides(color = guide_legend())
  for (plx in ggCallback) pl <- pl + plx
  
  message("Plot has ", n_pages(pl), " pages.\n")
  
  cf_outputFigure(pl = pl, ...)
}

# -------------------------------------------------------------------------#
# Plotting ----
# -------------------------------------------------------------------------#

#' Title
#'
#' @param fp,type see code of [getPaginateInfo()]
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @examples
#' pl <- ggplot(diamonds) +
#'   geom_point(aes(carat, price), alpha = 0.1) +
#'   facet_grid_paginate(color ~ cut:clarity, ncol = 3, nrow = 3, page = 4)
#' fp <- pl$facet$params
#' type <- "facet_grid_paginate"
#' getFacets(fp, type)
getFacets <- function(fp, type) {
  facets <- NULL
  if (type == "facet_grid_paginate"){
    nmr <- ifelse(length(names(fp$rows)), paste0(names(fp$rows), collapse = " + "), ".")
    nmc <- ifelse(length(names(fp$cols)), paste0(names(fp$cols), collapse = " + "), ".")
    facets <- as.formula(paste0(nmr, " ~ ", nmc))
  } else {
    facets <- as.formula(paste0(". ~ ", paste0(names(fp$facets), collapse = " + ")))
  }
  facets
}


#' get the pagination info from a plot
#'
#' @param pl ggplot
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
#' pl <- ggplot(diamonds) +
#'   geom_point(aes(carat, price), alpha = 0.1) +
#'   facet_grid_paginate(color ~ cut:clarity, ncol = 3, nrow = 3, page = 4)
#' getPaginateInfo(pl)
#' pl <- ggplot(diamonds) +
#'   geom_point(aes(carat, price), alpha = 0.1) +
#'   facet_grid_paginate(color ~ ., nrow = 2, ncol = 2, page = 1)
#' getPaginateInfo(pl)
#' pl <- ggplot(diamonds) +
#'   geom_point(aes(carat, price), alpha = 0.1) +
#'   facet_grid_paginate( ~ color, nrow = 2, ncol = 2, page = 1)
#' getPaginateInfo(pl)
#' pl <- ggplot(diamonds) +
#'   geom_point(aes(carat, price), alpha = 0.1) +
#'   facet_wrap_paginate( ~ color, nrow = 2, ncol = 2, page = 1)
#' getPaginateInfo(pl)
getPaginateInfo <- function(pl) {
  
  type <- switch(class(pl$facet)[1], FacetGridPaginate = "facet_grid_paginate", 
                 FacetWrapPaginate = "facet_wrap_paginate", NA)
  if (is.na(type)) return(NA)
  
  facet_paginate <- utils::getFromNamespace(type, "ggforce")
  
  fp <- pl$facet$params
  # recover arguments "facets", "scales", "space"
  facets <- getFacets(fp = fp, type = type)
  scales <- switch(as.character(as.numeric(fp$free$x + 2*fp$free$y)), "0" = "fixed", "1" = "free_x", "2" = "free_y", "3" = "free")
  space <- NULL
  if ("space_free" %in% names(fp))
    space  <- switch(as.character(as.numeric(fp$space_free$x + 2*fp$space_free$y)), "0" = "fixed", "1" = "free_x", "2" = "free_y", "3" = "free")
  
  # Assemble final arglist
  paginateInfo <- c(list(facets = facets, scales = scales, space = space), fp)
  paginateInfo <- paginateInfo[intersect(names(paginateInfo), names(formals(facet_paginate)))]
  paginateInfo <- paginateInfo[setdiff(names(paginateInfo), "page")]
  # "shrink" is not supported for grid
  # "shrink" and switch is not supported for wrap
  list(facet_paginate = facet_paginate, paginateInfo = paginateInfo)
}

#' Get list of paginated plots
#' 
#' use ggforce
#' 
#' @param pl (not yet facetted) ggplot
#' @param paginateInfo output from [cf_paginateInfo]
#'
#' @return list of ggplots
#' 
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' 
#' @importFrom ggforce facet_wrap_paginate facet_grid_paginate
#'
#' @examples
#' # adapted from ggforce examples
#' pl <- ggplot(diamonds) +
#'   geom_point(aes(carat, price), alpha = 0.1) +
#'   facet_grid_paginate(color ~ cut:clarity, ncol = 3, nrow = 3, page = 4)
#' plotlist <- cf_applyPaginate(pl)
#' plotlist[[1]]
#' plotlist[[5]]
#' n_pages(pl)
#' length(plotlist)
cf_applyPaginate <- function(pl) {
  
  pi <- getPaginateInfo(pl)
  
  if (!is.list(pi)) {
    if ("gg" %in% class(pl) || "ggplot" %in% class(pl)) {
      return(list(pl))
    } else if (is.list(pl))
      return(pl)
  }
  facet_paginate <- pi$facet_paginate
  paginateInfo <- pi$paginateInfo
  n <- ggforce::n_pages(pl)
  lapply(1:n, function(i) {
    pl + {do.call(facet_paginate, c(paginateInfo,list(page = i)))}
  })
}


#' Output figures
#' 
#' Handle paginate automatically, save asynchronously with future
#'
#' @param pl plot like ggplot, potentially created with facet_grid_paginate or list of ggplots
#' @param filename,width,height,scale,units,dpi,limitsize,device,... see [ggplot2::ggsave()]
#' @param FLAGFuture Export asynchronously with the future-package
#'
#' @return nothing
#' @export
#' 
#' @importFrom future plan "%<-%"
#' @importFrom grDevices dev.off dev.cur dev.set
#' @importFrom utils capture.output
#' 
#' @details # importFrom ggplot2 plot_dev parse_dpi plot_dim
#' 
#' @example 
cf_outputFigure <- function(pl, filename = NULL, 
                            width = 16, height = 10, scale = 1, 
                            units = c("cm", "mm", "in")[1], 
                            dpi = 300, limitsize = TRUE, 
                            device = NULL,
                            heightrel = NULL,
                            FLAGFuture = TRUE,
                            FLAGoverwrite = TRUE,
                            ...) {
  
  # Handle overwrite etc
  if (is.null(filename)) return(pl)
  if (!FLAGoverwrite & file.exists(filename)) {
    cat("FLAGoverwrite = FALSE. Plot is not written to disk\n")
    return(invisible(pl))
  }
  
  # Handle heightrel
  if (!is.null(heightrel)) height <- width * heightrel
  
  # Handle paginate: Wraps plot in list of length n_pages. 
  # For unpaginated plots, length(pl)=1
  pl <- cf_applyPaginate(pl) 
  if (length(pl)>1) {
    if(!grepl("pdf$", filename) && !grepl("%03d.png$", filename)) {
      filename <- gsub(".png", "%03d.png", filename)
    }
  } 
  
  # handle pdf
  # Was buggy for multipage pdf
  if (is.null(device) && tools::file_ext(filename) == "pdf") {
    if (length(pl) == 1) {
      cat("1 page, using cairo for export\n")
      device <- grDevices::cairo_pdf
    } else {
      cat("Multiple pages, using normal pdf device")
    }
    
  }
  # device wrestling (from ggsave)
  dpi <- ggplot2:::parse_dpi(dpi)
  dev <- ggplot2:::plot_dev(device, filename, dpi = dpi)
  dim <- ggplot2:::plot_dim(c(width, height), scale = scale, units = units, 
                            limitsize = limitsize)
  
  doPlot <- function() {
    old_dev <- grDevices::dev.cur()
    dev(filename = filename, width = dim[1], height = dim[2], ...)
    on.exit(utils::capture.output({
      grDevices::dev.off()
      if (old_dev > 1) grDevices::dev.set(old_dev)
    }))
    for (p in pl) print(p)
    "done"
  }
  
  if (FLAGFuture) {
    if (!"multisession" %in% class(future::plan())) {
      future::plan("multisession")
      cat("Future plan is now 'multisession'\n")
    }
    # message("Temporarily affecting Sys.getenv('OMP_NUM_THREADS')'\n")
    Sys.setenv(OMP_NUM_THREADS = 2)
    future::`%<-%`(.dummy, {doPlot()})
    Sys.setenv(OMP_NUM_THREADS = 1)
  } else {
    doPlot()
  }
  
  invisible(pl)
}


#' Grab a base plot into grid object
#'
#' @return gtree, can be exported with [ggplot2::ggsave()]
#' @export
#' @importFrom grid grid.grab
#' @importFrom gridGraphics grid.echo
cfgrab <- function() {
  gridGraphics::grid.echo()
  grid::grid.grab()
}


#' Nice ggplot theme
#' 
#' Taken from dMod
#' 
#' @param base_size,base_family see ?theme_bw
#' @param FLAGbold Make text bold faced
#'
#' @export
theme_cf <- function(base_size = 11, base_family = "", FLAGbold = TRUE) {
  colors <- list(
    medium = c(gray = '#737373', red = '#F15A60', green = '#7AC36A', blue = '#5A9BD4', orange = '#FAA75B', purple = '#9E67AB', maroon = '#CE7058', magenta = '#D77FB4'),
    dark = c(black = '#010202', red = '#EE2E2F', green = '#008C48', blue = '#185AA9', orange = '#F47D23', purple = '#662C91', maroon = '#A21D21', magenta = '#B43894'),
    light = c(gray = '#CCCCCC', red = '#F2AFAD', green = '#D9E4AA', blue = '#B8D2EC', orange = '#F3D1B0', purple = '#D5B2D4', maroon = '#DDB9A9', magenta = '#EBC0DA')
  )
  gray <- colors$medium["gray"]
  black <- colors$dark["black"]
  
  out <- theme_bw(base_size = base_size, base_family = base_family) +
    theme(line = element_line(colour = "black"),
          rect = element_rect(fill = "white", colour = NA),
          text = element_text(colour = "black"),
          axis.text = element_text(size = rel(1.0), colour = "black"),
          axis.text.x = element_text(margin=unit(c(3, 3, 0, 3), "mm")),
          axis.text.y = element_text(margin=unit(c(3, 3, 3, 0), "mm")),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(-2, "mm"),
          legend.key = element_rect(colour = NA),
          panel.border = element_rect(colour = "black"),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(size = rel(1.0))
    ) 
  if (FLAGbold) out <- out + theme(text = element_text(face = "bold"))
  out
}


#' cfggplot
#'
#' @param data,mapping see ?ggplot
#' @param FLAGbold Make text bold faced
#'
#' @return ggplot
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
cfggplot <- function(data = NULL, mapping = aes(), FLAGbold = TRUE) {
  ggplot(data,mapping) + 
    theme_cf(FLAGbold = FLAGbold)
}

#' scalecolorcf
#' Copied from dMod
#' @param ... see scale_color_manual
#'
#' @return added to plots
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
scale_color_cf <- function(...) {
  scale_color_manual(..., values = cfcolors)
}

#' Colors Copied from dMod
#' @export
cfcolorsFULL <- list(
  medium = c(black = '#010202', gray = '#737373', red = '#F15A60', green = '#7AC36A', blue = '#5A9BD4', orange = '#FAA75B', purple = '#9E67AB', maroon = '#CE7058', magenta = '#D77FB4'),
  dark   = c(black = '#010202', gray = '#333333', red = '#EE2E2F', green = '#008C48', blue = '#185AA9', orange = '#F47D23', purple = '#662C91', maroon = '#A21D21', magenta = '#B43894'),
  light  = c(black = '#010202', gray = '#CCCCCC', red = '#F2AFAD', green = '#D9E4AA', blue = '#B8D2EC', orange = '#F3D1B0', purple = '#D5B2D4', maroon = '#DDB9A9', magenta = '#EBC0DA')
)

#' Colors
#' copied and adapted from dMod
#' @examples 
#' ggplot(data.frame(x = 1:42,color = factor(cfcolors[1:42], unique(cfcolors[1:42])))) + 
#'   geom_tile(aes(x=x, y = 1, fill = color)) + 
#'   scale_color_manual(values= setNames(nm = cfcolors),
#'   aesthetics = c("color", "fill"))
#' @export
cfcolors <- c("#000000", "#C5000B", "#0084D1", "#579D1C", "#FF950E", 
              "#4B1F6F", "#CC79A7","#006400", "#F0E442", "#8B4513",
              "salmon", "slateblue1", "chocolate3", "firebrick", 
              "cyan3", "chartreuse4", "gold", "ivory4", "seagreen3", "dodgerblue",
              RColorBrewer::brewer.pal(9, "Reds"  )[c(3,5,7,9)],
              RColorBrewer::brewer.pal(9, "Blues" )[c(3,5,7,9)],
              RColorBrewer::brewer.pal(9, "Greens")[c(3,5,7,9)],
              RColorBrewer::brewer.pal(8, "Accent")[c(3,5,7,9)],
              rep("gray", 100))



# -------------------------------------------------------------------------#
# Output table ----
# -------------------------------------------------------------------------#



#' Output table as markdown
#'
#' @param dt data.table
#' @param split_by Charactor of colnames to split the table by inserting separator lines
#' @param filename 
#'
#' @return object of class knitr::kable
#' 
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @importFrom stringr str_pad
#'
#' @examples
#' tableName <- "bla"
#' dt <- data.table(iris)[1:5]
#' cfoutput_MdTable(dt)
#' cfoutput_MdTable(data.table(iris), split_by = "Species")
#' split_by <- NULL
#' filename <- NULL
#' format <- c("markdown", "pandoc")
#' caption <- NULL
#' na.strings <- "NA"
#' FLAGsummaryRow <- TRUE
#' NFLAGtribble <- 2
cfoutput_MdTable <- function(dt, split_by = NULL, filename = NULL, format = c("markdown", "pandoc"), 
                             caption = NULL, na.strings = ifelse(NFLAGtribble == 2, "NA", "-"), FLAGsummaryRow = TRUE, 
                             NFLAGtribble = 0,
                             ...) {
  tableName <- as.character(substitute(dt))
  
  options(knitr.kable.NA = na.strings)
  
  # kt <- knitr::kable(dt,format = format[1], caption = caption)
  kt <- knitr::kable(dt,format = format[1], caption = caption, ...)
  
  seprow <- gsub(":","-",kt[2 + 2*(!is.null(caption))])
  widths <- nchar(strsplit(seprow, "|", fixed = TRUE)[[1]][-1])
  
  if (NFLAGtribble) {
    types <- vapply(dt, class, "double")
    for (fact in names(types)[types=="factor"]) dt[[fact]] <- as.character(dt[[fact]])
    hasDoubleQuote <- vapply(names(types)[types%in%c("character","factor")], function(nm) any(grepl('"', dt[[nm]])), FUN.VALUE = TRUE)
    hasSingleQuote <- vapply(names(types)[types%in%c("character","factor")], function(nm) any(grepl("'", dt[[nm]])), FUN.VALUE = TRUE)
    if (any(hasDoubleQuote & hasSingleQuote)) 
      warning("This table has double quotes (\") AND single quotes ('). Using (\"), but the output must be checked manually.")
    quoteSymbol <- ifelse(any(hasDoubleQuote), "'", '"')
    for (fact in names(types)[types%in%c("character","factor")]) dt[[fact]] <- paste0(quoteSymbol, dt[[fact]], quoteSymbol)
    kt <- knitr::kable(dt,format = "markdown")
    kt <- kt[-(1:2)]
    FLAGpastedt <- Sys.info()["nodename"] != "IQdesktop"
    row0 <- ifelse(FLAGpastedt, "data.table(","")
    row0 <- paste0(row0, "tibble::tribble(")
    row1 <- paste0(stringr::str_pad(paste0("~", names(dt), ","), width = widths, side = "left"), collapse = "")
    kt <- substr(kt, 2, nchar(kt))
    kt <- gsub("\\|", ",", kt)
    rowN <- kt[length(kt)]
    rowN <- substr(rowN, 1, nchar(rowN)-1)
    rowN <- paste0(rowN, ")")
    rowN <- ifelse(FLAGpastedt,paste0(rowN, ")"),paste0(rowN, "")) # not nice but too lazy to do properly
    kt <- kt[-length(kt)]
    kt <- c(row0, row1, kt, rowN, "")
    
    if (length(tableName) == 1) 
      kt[1] <- paste0(tableName, " <- ", kt[1])
    
    if (NFLAGtribble == 1) cat(kt, sep = "\n")
    if (NFLAGtribble == 2) {
      e <- rstudioapi::getSourceEditorContext()
      rstudioapi::documentSave(e$id)
      rstudioapi::insertText(e$selection[[1]]$range$start, paste0(kt, collapse = "\n"))
    }
    return(invisible(kt))
  }
  
  if (!is.null(split_by)){
    dt <- copy(as.data.table(dt))
    dn <- dt[,list(nlines = .N), by = split_by]
    dn[,`:=`(rowid = cumsum(nlines))]
    dn[,`:=`(rowid = rowid + 2 + 2*(!is.null(caption)))]
    atr <- attributes(kt)
    for (i in rev(dn$rowid)[-1]) {
      kt <- c(kt[1:i], seprow, kt[(i+1):length(kt)])
    }
    attributes(kt) <- atr
  }
  
  if (FLAGsummaryRow) {
    summaryrow <- vapply(dt, function(x) {
      if (is.numeric(x))   return(paste0("Sum=", sum(x)))
      if (is.character(x)) return(paste0("Lvl=", length(unique(x))))
    }, FUN.VALUE = "N=nunique")
    summaryrow <- vapply(seq_along(summaryrow), function(i) sprintf(paste0("%",widths[i],"s"), summaryrow[i]), "bla")
    summaryrow <- paste0("|", paste0(summaryrow, collapse = "|"), "|")
    kt <- c(kt, seprow, summaryrow)
  }
  
  if(!is.null(filename))
    writeLines(kt, filename)
  
  cat(kt, sep = "\n")
  invisible(kt)
}



# -------------------------------------------------------------------------#
# Cluster time stamps ----
# -------------------------------------------------------------------------#
#' Write the time stamp when you logged into the cluster
#'
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Cluster login helpers
write_clusterTimeStamp <- function() {
  clusterTimeStampFile <- file.path("~", ".clusterTimeStamp.rds")
  if (file.exists(clusterTimeStampFile)) {
    rstudioapi::insertText(c(1,1), text = "check_clusterTimeStamp()", "#console")
    stop("An old clusterTimeStamp file exists. Please run check_clusterTimeStamp()")
  }
  saveRDS(Sys.time(), file = clusterTimeStampFile)
}

#' Check for the time stamp
#'
#' If no recent login-stamp is available or if stamp is invalidated (> 1 hour){
#' if stamp is invalidated: {remove stamp}
#' stop, force creation of new stamp}
#'
#'
#' @param FLAGforcePurge Force the removal of a time stamp
#'
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Cluster login helpers
check_clusterTimeStamp <- function(FLAGforcePurge = FALSE) {
  clusterTimeStampFile <- file.path("~", ".clusterTimeStamp.rds")
  
  if (file.exists(clusterTimeStampFile)) {
    clusterTimeStamp <- readRDS(clusterTimeStampFile)
    dt <- difftime(Sys.time(), clusterTimeStamp, units = "mins")
    cat("Last login: ", round(dt,2), " minutes ago\n")
    
    if (dt >= 59 || FLAGforcePurge) {
      cat("Removing old time stamp\n")
      unlink(clusterTimeStampFile)
    }
  }
  
  if (!file.exists(clusterTimeStampFile)) {
    rstudioapi::insertText(c(1,1), text = "write_clusterTimeStamp();Q", "#console")
    stop("Please login again manually via console. Then, in R call write_clusterTimeStamp()")
  }
  "All good, login is still active :)"
}



