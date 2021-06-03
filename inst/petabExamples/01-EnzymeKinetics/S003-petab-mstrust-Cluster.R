# library(petab)
devtools::load_all("~/Promotion/Promotion/Projects/petab")
setwd(tempdir())
petab_python_setup()
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

.outputFolder <- "."
# -------------------------------------------------------------------------#
# Load petab, create L1 problem, write L1 problem ----
# -------------------------------------------------------------------------#

SFLAGbrowser = "1"

testCases = FALSE
path2TestCases = "PEtabTests/"

# debugonce(readPd)
pd <- importPEtabSBML_indiv("petab", NFLAGcompile = 3, .compiledFolder = "Compiled", SFLAGbrowser = "0")

pd_cluster_mstrust(pd, .outputFolder = ".", 16, 2)
pd_cluster_mstrust(pd, .outputFolder = ".", 16, 2, identifier = "redoit")

pd <- readPd(pd_files(pd$filenameParts)$rdsfile)

# Exit ---- 
