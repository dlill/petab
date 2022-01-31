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
# new value is 84.80751
# [ ] get rid of difference!

fitModelPEtabSBML(model_name = "petab")


filename = "petab/petab.yaml"
testCases = FALSE
path2TestCases = "PEtabTests/"
.compiledFolder = file.path("CompiledObjects")
NFLAGcompile = c(Auto = 3, Recompile = 0, RebuildGrids = 1, LoadPrevious = 2)[1]
SFLAGbrowser = c("0None", "1Beginning", "2BuildGrids", "3Compilation", "4CollectList",
                 "5Scales", "6ParameterFormulaInjection")[1]


# Exit ----
