library(petab)
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# -------------------------------------------------------------------------#
# Create enzyme kinetics model and data ----
# -------------------------------------------------------------------------#
pd <- importPEtabSBML_indiv("petab", NFLAGcompile = 3, .compiledFolder = "Compiled", SFLAGbrowser = "0")

# .. Test model -----
pred <- pd$prd(seq(0,100), pd$pars)
val <- pd$obj_data(pd$pars)
print(val, 20,20)

pd_predictAndPlot(pd)

pd <- pd_fitObsPars(pd, NFLAGsavePd = 3)

# Test fitting
# myfit <- trust(pd$obj_data, pd$pars,1,10,iterlim = 1000)
# plotCombined(pd$prd(seq(0,100), pd$pars), pd$dModAtoms$data)
# plotCombined(pd$prd(seq(0,100), myfit$argument), pd$dModAtoms$data)
pd_predictAndPlot2(pd)

# Mstrust
center <- pepy_sample_parameter_startpoints(pd$pe, n_starts = 8)
fits <- mstrust(pd$obj_data, center, paste0("fit", 1), 
                fits = 4, iterlim = 1000, cores = 4, output = TRUE,
                cautiousMode = TRUE)
fits <- conveniencefunctions::cf_as.parframe(fits)
plotValues(fits)
dMod_saveMstrust(fits, ".", FLAGoverwrite = TRUE)
unlink("fit", T)

# .. profiles -----
profs <- dMod::profile(pd$obj_data,as.parvec(pd$result$fit), names(pd$pars)[1:2], cautiousMode = TRUE)


# Exit ----
