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

pd <- importPEtabSBML_indiv(filename, NFLAGcompile = 0, .compiledFolder = "Compiled", SFLAGbrowser = "0")

# Test model
pd_predictAndPlot2(pd, i = time > 10)
pred <- pd$prd(seq(0,100), pd$pars)
pd$obj(pd$pars)

# Test fitting
myfit <- trust(pd$obj, pd$pars,1,10,iterlim = 1000)
plotCombined(pd$prd(seq(0,100), pd$pars), pd$dModAtoms$data)
plotCombined(pd$prd(seq(0,100), myfit$argument), pd$dModAtoms$data)

# # Mstrust
# center <- pepy_sample_parameter_startpoints(pd$pe, n_starts = 8)
# fits <- mstrust(pd$obj_data, center, "fit", fits = 5, iterlim = 20, cores = 8)
# fits <- as.parframe(fits)
# plotValues(fits)
# dMod_saveMstrust(fits, ".", FLAGoverwrite = TRUE)
# unlink("fit", T)

parf <- expand.grid(A = myfit$argument[["A"]],
                    B = myfit$argument[["B"]],
                    C = myfit$argument[["C"]],
                    k_AB = myfit$argument[["k_AB"]] + seq(-1,1,0.1),
                    k_BC = myfit$argument[["k_BC"]] + seq(-1,1,0.1),
                    sigma_obsC = myfit$argument[["sigma_obsC"]])
rownames(parf) <- NULL
parf <- parframe(parf, parameters = names(parf))                    

objvals <- lapply(seq_len(nrow(parf)), function(i) pd$obj(as.parvec(parf,i))$value)

p <- data.table(parf, objval = do.call(c, objvals))

cfggplot(p, aes(k_AB, k_BC, color = log10(objval - min(objval)+0.01), fill = log10(objval - min(objval)+0.01))) + geom_tile()


# Exit ----
