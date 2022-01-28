library(petab)
library(dMod)
library(dplyr)
library(stringr)
library(conveniencefunctions)

try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# pd <- importPEtabSBML(modelname = "petab", path2model = "")

pd <- importPEtabSBML_indiv(filename = "petab/petab.yaml")

# original value was -89.64576
pd$obj(pd$pars)
# new value is 85.36924 
# [ ] get rid of difference!

fitModelPEtabSBML(model_name = "petab")


# Exit ----
