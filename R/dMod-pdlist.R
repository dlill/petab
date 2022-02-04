
#' Slightly modified default aesthetics for multimodel plotting
#'
#' @param x 
#' @param y 
#' @param color 
#' @param fill 
#' @param linetype 
#' @param ymin 
#' @param ymax 
#' @param shape 
#' @param ... 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family pdlist
#'
#' @examples
pdlist_petab_plotHelpers_aeslist <- function(x = ~time,
                                             y = ~measurement,
                                             color = ~modelId,
                                             fill = ~modelId,
                                             linetype = ~parameterSetId,
                                             ymin = ~measurementmin,
                                             ymax = ~measurementmax,
                                             shape = ~modelId,
                                             ...) {
  list(x = x, y = y, color = color, fill = fill, linetype = linetype, ymin = ymin, ymax = ymax, shape = shape, ...)
}

#' Plot with list of pds
#'
#' @param pdlist named list of pds
#' @inheritParams pd_predictAndPlot2
#' @param ... Arguments going to pd_predictAndPlot2
#'
#' @return ggplot or plotData = list(dplot,pplot,pplotRibbon)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family pdlist
#' @family Plotting
#' @importFrom data.table rbindlist
#' @importFrom conveniencefunctions cfggplot cfgg_getAllAesthetics scale_color_cf
#' 
#' @examples
pdlist_predictAndPlot2 <- function(pdlist, 
                                   FLAGreturnPlotData = FALSE,
                                   nrow = 4, ncol = 5,
                                   aeslist = pdlist_petab_plotHelpers_aeslist(),
                                   ggCallback = list(),
                                   opt.gg = list(ribbonAlpha = 0.2), # would be nice to put this into opt.profile or maybe opt.gg?
                                   ...
) {
  
  plotData <- lapply(pdlist, function(pd) pd_predictAndPlot2(pd, FLAGreturnPlotData = TRUE, ...))
  dplot <- data.table::rbindlist(lapply(plotData, function(x) x$dplot), use.names = TRUE, fill = TRUE, idcol = "modelId")
  pplot <- data.table::rbindlist(lapply(plotData, function(x) x$pplot), use.names = TRUE, fill = TRUE, idcol = "modelId")
  pplotRibbon <- data.table::rbindlist(lapply(plotData, function(x) x$pplotRibbon), use.names = TRUE, fill = TRUE, idcol = "modelId")
  pplotRibbon <- if(!length(pplotRibbon)) NULL else pplotRibbon
  if (FLAGreturnPlotData) return(list(dplot = dplot, pplot = pplot, pplotRibbon = pplotRibbon))
  
  if (length(unique(c(dplot$conditionId, pplot$conditionId)))>1) cat("More than one condition, you might want to map it to an aesthetic via aeslist = pdlist_petab_plotHelpers_aeslist()")
  
  pl <- conveniencefunctions::cfggplot()
  if (nrow(dplot)) pl <- pl + geom_point(do.call(aes_q, aeslist[intersect(names(aeslist), conveniencefunctions::cfgg_getAllAesthetics()[["geom_point"]])]), data = dplot)
  if (nrow(pplot)) pl <- pl + geom_line( do.call(aes_q, aeslist[intersect(names(aeslist), conveniencefunctions::cfgg_getAllAesthetics()[["geom_line"]])]) , data = pplot)
  if (!is.null(pplotRibbon)) {
    aesl <- aeslist[intersect(names(aeslist), conveniencefunctions::cfgg_getAllAesthetics()[["geom_ribbon"]])]
    aesl <- aesl[setdiff(names(aesl), c("linetype", "lty", "y", "color", "colour"))]
    pl <- pl + geom_ribbon(do.call(aes_q, aesl), data = pplotRibbon, alpha = opt.gg$ribbonAlpha)}
  pl <- pl + conveniencefunctions::scale_color_cf(aesthetics = c("color", "fill")) + 
    facet_wrap_paginate(~observableId, nrow = nrow, ncol = ncol, scales = "free") + 
    scale_y_continuous(n.breaks = 5)
  for (plx in ggCallback) pl <- pl + plx
  pl
}
