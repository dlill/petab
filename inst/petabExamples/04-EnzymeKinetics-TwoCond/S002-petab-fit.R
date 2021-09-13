devtools::load_all("~/Promotion/Promotion/Projects/petab")
setwd("~")
petab_python_setup()
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# -------------------------------------------------------------------------#
# Fit ----
# -------------------------------------------------------------------------#
filename <- "petab"

pd <- importPEtabSBML_indiv("petab", NFLAGcompile = 2, .compiledFolder = "Compiled", SFLAGbrowser = "0")
# pd_fit(pd)
pd <- readPd(pd_files(pd$filenameParts)$rdsfile)


# fits <- mstrust(pd$obj, pepy_sample_parameter_startpoints(pd$pe, 200), studyname = "fit")
# fits <- cf_as.parframe(fits)
# dMod_saveMstrust(fits, path = ".", "mstrust", FLAGoverwrite = TRUE)

pd <- readPd(pd_files(pd$filenameParts)$rdsfile)
pl <- pd_plotParsParallelLines(pd)
cf_outputFigure(pl, filename = "wup.pdf")

# -------------------------------------------------------------------------#
# Collect report ----
# -------------------------------------------------------------------------#
pd_petabSelect_reportYaml(pd)

# Exit ----
