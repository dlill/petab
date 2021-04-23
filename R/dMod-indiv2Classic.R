#' Title
#'
#' @param eg est.grid
#' @param fg fix.grid
#'
#' @return condition.grid(condition, parsMixedCols...)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
indiv2Classic_merge_grids <- function(eg,fg) {
  cg <- copy(eg)
  cg <- cbind(cg, fg[,.SD,.SDcols = setdiff(names(fg), names(eg))])
  mixed <- intersect(names(fg), names(eg))
  mixed <- setdiff(mixed, c("condition", "ID"))
  for (m in mixed) cg[[m]][is.na(cg[[m]])] <- fg[[m]][is.na(cg[[m]])]
  cg
}

#' Title
#'
#' @param gridlist gridlist
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
indiv2Classic_gridlist2cond.grid <- function(gridlist) {
  cg <- indiv2Classic_merge_grids(gridlist$est.grid,gridlist$fix.grid)
  cg[,`:=`(ID = NULL)]
  cg <- as.data.frame(cg)
  rownames(cg) <- cg$condition
  cg
}

#' Title
#'
#' @param trafo base trafo
#' @param cg output of [indiv2Classic_gridlist2cond.grid()]
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom dMod branch insert
#'
#' @examples
indiv2Classic_trafo <- function(trafo, cg) {
  trafoL <- dMod::branch(trafo, cg)
  assign(".pars_to_insert", setdiff(names(cg), "condition"), .GlobalEnv) # hacky
  trafoL <- dMod::insert(trafoL, "name ~ value", value = unlist(mget(.pars_to_insert)), name = .pars_to_insert)
  trafoL
}

#' Compile classic trafo
#'
#' @param trafoL output of [indiv2Classic_trafo]
#' @param .compiledFolder Folder for compiled objects
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom dMod P
#'
#' @examples
indiv2Classic_compileTrafo <- function(trafoL, .compiledFolder = file.path("CompiledObjects")) {
  w <- getwd()
  on.exit({setwd(w)})
  setwd(.compiledFolder)
  p <- P(trafoL, compile = TRUE, modelname = "pClassic")
  setwd(w)
  p
}

#' Title
#'
#' @param pd
#'
#' @return a `pd` object but with classic fns instead of indiv fns
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom dMod normL2
#'
#' @examples
indiv2Classic <- function(pd,
                          NFLAGcompile = c(Recompile = 0, LoadPrevious = 1)[3],
                          Nobjtimes = 50) {

  pd$filenameParts$type                        <- "classic"
  rdsfile <- pd_files(pd$filenameParts)$rdsfile

  if (NFLAGcompile > 0 && file.exists(rdsfile))
    return(readPd(rdsfile))

  cg <- indiv2Classic_gridlist2cond.grid(pd$dModAtoms$gridlist)

  # Calculate objects
  trafoL   <- indiv2Classic_trafo(trafo = pd$dModAtoms$symbolicEquations$trafo, cg = cg)
  p        <- indiv2Classic_compileTrafo(trafoL = trafoL, .compiledFolder = .compiledFolder)
  prd      <- (pd$dModAtoms$fns$g * pd$dModAtoms$fns$x * p)
  obj_data <- dMod::normL2(pd$dModAtoms$data, prd, pd$dModAtoms$e, dMod::objtimes(pd$pe$measurementData$time, Nobjtimes = Nobjtimes))

  # Re-populate pd
  pd$dModAtoms$symbolicEquations$trafo <- trafoL
  pd$dModAtoms$fns$p0                  <- p
  pd$prd                               <- prd
  pd$obj_data                          <- obj_data

  # Save and return
  saveRDS(pd, rdsfile)

  pd
}



#' Get the faster version of classic and indiv
#'
#' @param pd pd_indiv
#' @param NFLAGcompile recompile the parameter trafo of pdc?
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
indiv2classic_getFasterVersion <- function(pd, NFLAGcompile = 1) {

  # build classic with comparable objtimes
  objtimes <- dMod::controls(pd$obj_data, name = "times")
  pdc <- indiv2Classic(pd, NFLAGcompile = NFLAGcompile)
  dMod::controls(pdc$obj_data, name = "times") <- objtimes

  # compare
  t1 <- system.time(pd$obj_data(pd$pars))
  t2 <- system.time(pdc$obj_data(pd$pars))
  cat("indiv: ", round(t1,3), "\n")
  cat("classic: ", round(t2,3), "\n")

  # return faster version
  if (t1 < t2) return(pd)
  pdc
}


