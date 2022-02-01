library(petab)
library(dMod)
library(dplyr)
library(stringr)
library(conveniencefunctions)

try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# pd <- importPEtabSBML(modelname = "petab", path2model = "")

pd <- importPEtabSBML_indiv(filename = "petab/petab.yaml")

# original value was -89.63
pd$obj(pd$pars)




.unbug <- F
if(.unbug){
  filename = "petab/petab.yaml"
  testCases = FALSE
  path2TestCases = "PEtabTests/"
  .compiledFolder = file.path("CompiledObjects")
  NFLAGcompile = c(Auto = 3, Recompile = 0, RebuildGrids = 1, LoadPrevious = 2)[1]
  SFLAGbrowser = c("0None", "1Beginning", "2BuildGrids", "3Compilation", "4CollectList",
                   "5Scales", "6ParameterFormulaInjection")[1]
  
  
  
  prd0 <- Reduce("*", list(pd$dModAtoms$fns$g, pd$dModAtoms$fns$x, pd$dModAtoms$fns$p0))
  obj_data <- normL2_indiv(pd$dModAtoms$data, prd0,
                           pd$dModAtoms$e,
                           est.grid = pd$dModAtoms$gridlist$est.grid,
                           fix.grid = pd$dModAtoms$gridlist$fix.grid,
                           times = mytimes)
  obj_prior <- constraintL2(bestfit, sigma = 12)
  
  # Rebuild obj
  obj <- Reduce("+", list(obj_data, obj_prior))
  obj(bestfit)$value
}


pd_predictAndPlot2(pd)

pd_fit(pd)
pd <- readPd(pd_files(pd$filenameParts)$rdsfile)
pd_predictAndPlot2(pd)

fitModelPEtabSBML(model_name = "petab")





# Exit ----
