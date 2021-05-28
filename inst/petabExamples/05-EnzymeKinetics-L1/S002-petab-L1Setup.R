# library(petab)
devtools::load_all("~/Promotion/Promotion/Projects/petab")
setwd(tempdir())
petab_python_setup()
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# -------------------------------------------------------------------------#
# Load petab, create L1 problem, write L1 problem ----
# -------------------------------------------------------------------------#
pe <- readPetab("basemodel")
debugonce(pe_L1_createL1Problem)
pe <- pe_L1_createL1Problem(pe, parameterId_base = c("kcat","E"), conditionSpecL1_reference = "C1", j_conditionSpecL1 = conditionId)
writePetab(pe, "petab")

# Exit ---- 
