# library(petab)
devtools::load_all("~/Promotion/Promotion/Projects/petab")
setwd(tempdir())
petab_python_setup()
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# -------------------------------------------------------------------------#
# Load petab, create L1 problem, write L1 problem ----
# -------------------------------------------------------------------------#
pe <- readPetab("basemodel")
priorStrength = c("kcat"= 0.1,"E" = 0.2)
pe <- petab_createqNLMEProblem(pe = pe, priorStrength = priorStrength, Nsub = 3, seed = 1)
writePetab(pe, "petab")

# subjectId is now present in metaInformation
petab_plotData(pe, ggCallback = list(facet_wrap(~subjectId + observableId)))

# Exit ---- 
