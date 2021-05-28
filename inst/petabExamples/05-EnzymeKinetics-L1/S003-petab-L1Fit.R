# library(petab)
devtools::load_all("~/Promotion/Promotion/Projects/petab")
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# -------------------------------------------------------------------------#
# Arguments ----
# -------------------------------------------------------------------------#
filename <- "petab"

NFLAGcompile = 0
.compiledFolder = "Compiled"
SFLAGbrowser = "0"

testCases = FALSE
path2TestCases = "PEtabTests/"

# -------------------------------------------------------------------------#
# OLD ----
# -------------------------------------------------------------------------#
pd <- importPEtabSBML_indiv("petab", NFLAGcompile = 0, .compiledFolder = "Compiled", SFLAGbrowser = "0")
pd$pe

pd$prd(pd$times, pd$pars)


# Exit ----
