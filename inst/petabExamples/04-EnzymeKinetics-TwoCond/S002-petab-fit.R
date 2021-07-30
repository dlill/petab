library(petab)
setwd("~")
petab_python_setup()
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# -------------------------------------------------------------------------#
# Fit ----
# -------------------------------------------------------------------------#
filename <- "petab"

pd <- importPEtabSBML_indiv("petab", NFLAGcompile = 2, .compiledFolder = "Compiled", SFLAGbrowser = "0")
pd_fit(pd)
pd <- readPd(pd_files(pd$filenameParts)$rdsfile)

# -------------------------------------------------------------------------#
# Collect report ----
# -------------------------------------------------------------------------#
pd_petabSelect_reportYaml(pd)

# Exit ----
