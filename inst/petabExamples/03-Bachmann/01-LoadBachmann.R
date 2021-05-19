# -------------------------------------------------------------------------#
# 0 Header ----
# -------------------------------------------------------------------------#
#
# 01-LoadBachmann.R
#
# [PURPOSE]
# 
# Load the Bachmann model
#
#
# [AUTHOR]
# Daniel Lill
#
# [Date]
# Mon May 17 23:08:36 2021
#
# [Template Info]
# .. Template name TGFB_PETABEXPORT -----
# .. Template name TGFB_MSTRUST -----
# .. Template name TGFB_PROFILE -----
# .. Template version 0.0.1 -----
#
rm(list = ls(all.names = TRUE))
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# library(petab)
devtools::load_all("~/Promotion/Promotion/Projects/petab")

.outputFolder <- paste0("Output/")
.compiledFolder <- file.path(.outputFolder, "Compiled")
for(folder in c(.outputFolder,.compiledFolder)) 
if(!dir.exists(folder)) dir.create(folder)

# -------------------------------------------------------------------------#
# 1 Import ----
# -------------------------------------------------------------------------#
pe_file <- "Bachmann_MSB2011"
pe <- readPetab(pe_file)


# Apparently, the scale is not found
pd <- importPEtabSBML_indiv(pe_file, .compiledFolder = .compiledFolder, NFLAGcompile = 0, SFLAGbrowser = "2")


# -------------------------------------------------------------------------#
# Debug ----
# -------------------------------------------------------------------------#

pars <- pd$p(pd$pars)
pars <- pars[[1]]
pars <- unclass_parvec(pars)
f <- as.eqnvec(pd$dModAtoms$symbolicEquations$reactions)
fnum <- lapply(f, function(fx) {
  with(as.list(pars), eval(parse(text = fx)))
})
fnum <- do.call(c, fnum)
na_states <- fnum[is.na(fnum)]

parameters_naState <- f[names(na_states)[1]] %>% getSymbols()
pars[parameters_naState]

# The problem is CISRNAEqc which is zero
pe$parameters[grep("CISR", parameterId)]
# It has log10 scale but is not estimated. This might be the root of the problem, as petab_lint already warned us :)
petab_lint(pe)
# Indeed, the value in the fix.grid is transformed, but the trafo is not applied in the trafo
pd$dModAtoms$gridlist$fix.grid$CISRNAEqc
pd$dModAtoms$symbolicEquations$trafo["CISRNAEqc"]


# -------------------------------------------------------------------------#
# 2 Predict and plot ----
# -------------------------------------------------------------------------#
pd_predictAndPlot2(pd)

pd$prd(pd$times, pd$pars)





# Exit ----
future::plan("sequential")
