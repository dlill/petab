library(petab)
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# -------------------------------------------------------------------------#
# Create enzyme kinetics model and data ----
# -------------------------------------------------------------------------#
pd <- importPEtabSBML("petab", path2model = "./")
# ..  -----
conveniencefunctions::compare(getParameters(pd$dModAtoms$fns$p0), names(pd$pars))
pd$dModAtoms$symbolicEquations
pd$dModAtoms$fns$p0(c(pd$pars, structure(rep(0,4), .Names = c("obsE", "obsES", "obsP", "obsS"))))

# ..  -----
# Test model
pred <- pd$prd(seq(0,100), c(pd$pars, structure(rep(0,4), .Names = c("obsE", "obsES", "obsP", "obsS"))))
val <- pd$obj_data(c(pd$pars, structure(rep(0,4), .Names = c("obsE", "obsES", "obsP", "obsS"))))
print(val, 20,20)
# ..  -----
pred <- pd$prd(seq(0,100), pd$pars)
pd$obj_data(pd$pars)

# Test fitting
myfit <- trust(pd$obj_data, pd$pars,1,10,iterlim = 1000)
plotCombined(pd$prd(seq(0,100), pd$pars), pd$dModAtoms$data)
plotCombined(pd$prd(seq(0,100), myfit$argument), pd$dModAtoms$data)

# Mstrust
center <- pepy_sample_parameter_startpoints(pd$pe, n_starts = 8)
fits <- mstrust(pd$obj_data, center, "fit", fits = 5, iterlim = 20, cores = 8)
fits <- as.parframe(fits)
plotValues(fits)
dMod_saveMstrust(fits, ".", FLAGoverwrite = TRUE)
unlink("fit", T)
# -------------------------------------------------------------------------#
# Indiv ----
# -------------------------------------------------------------------------#
# pd <- importPEtabSBML_indiv("petab", NFLAGcompile = 3, .compiledFolder = "Compiled")

# Exit ----
