# library(petab)
devtools::load_all("~/Promotion/Promotion/Projects/petab")
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# -------------------------------------------------------------------------#
# Create enzyme kinetics model and data ----
# -------------------------------------------------------------------------#
filename <- "petab"
NFLAGcompile = 1
.compiledFolder = "Compiled"
SFLAGbrowser = "1"

testCases = FALSE
path2TestCases = "PEtabTests/"

# ..  -----

pd <- importPEtabSBML_indiv(filename, NFLAGcompile = 1, .compiledFolder = "Compiled", SFLAGbrowser = "1")

# Test model
pred <- pd$prd(seq(0,100), pd$pars)
pd$obj_data(pd$pars)

# Test fitting
myfit <- trust(pd$obj_data, pd$pars,1,10,iterlim = 1000)
plotCombined(pd$prd(seq(0,100), pd$pars), pd$dModAtoms$data)
plotCombined(pd$prd(seq(0,100), myfit$argument), pd$dModAtoms$data)

# # Mstrust
# center <- pepy_sample_parameter_startpoints(pd$pe, n_starts = 8)
# fits <- mstrust(pd$obj_data, center, "fit", fits = 5, iterlim = 20, cores = 8)
# fits <- as.parframe(fits)
# plotValues(fits)
# dMod_saveMstrust(fits, ".", FLAGoverwrite = TRUE)
# unlink("fit", T)


# Exit ----
