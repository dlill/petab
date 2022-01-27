library(petab)
library(dMod)
library(dplyr)
library(stringr)
library(conveniencefunctions)

try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

pd <- importPEtabSBML(modelname = "petab", path2model = "")

pd <- importPEtabSBML_indiv(filename = "petab/petab.yaml")

fitModelPEtabSBML()


# Exit ----
