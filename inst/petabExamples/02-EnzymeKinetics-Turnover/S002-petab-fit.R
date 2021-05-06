library(petab)
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# -------------------------------------------------------------------------#
# Create enzyme kinetics model and data ----
# -------------------------------------------------------------------------#
pd <- importPEtabSBML_indiv("petab", NFLAGcompile = 3, .compiledFolder = "Compiled", SFLAGbrowser = "0")

# Test model
pred <- pd$prd(seq(0,100), pd$pars)
val <- pd$obj_data(pd$pars)
print(val, 20,20)

pd_predictAndPlot(pd)

pd <- pd_fitObsPars(pd, NFLAGsavePd = 3)

# Test fitting
myfit <- trust(pd$obj_data, pd$pars,1,10,iterlim = 1000)
plotCombined(pd$prd(seq(0,100), pd$pars), pd$dModAtoms$data)
plotCombined(pd$prd(seq(0,100), myfit$argument), pd$dModAtoms$data)

# Mstrust
center <- pepy_sample_parameter_startpoints(pd$pe, n_starts = 8)
fits <- mstrust(pd$obj_data, center, "fit", fits = 5, iterlim = 50, cores = 8)
fits <- conveniencefunctions::cf_as.parframe(fits)
plotValues(fits)
dMod_saveMstrust(fits, ".", FLAGoverwrite = TRUE)
unlink("fit", T)
# -------------------------------------------------------------------------#
# Indiv ----
# -------------------------------------------------------------------------#

# Exit ----
