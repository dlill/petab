library(petab)
library(dMod)
library(dplyr)
library(stringr)
library(conveniencefunctions)

try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# -------------------------------------------------------------------------#
# Export Petab ----
# -------------------------------------------------------------------------#

# bestfit <- pouter
pe <- petab_dModmodel2PE(ODEmodel=reactions,
                         obsFun=observables,
                         errormodel=errors,
                         data=mydata,  
                         bestfit=bestfit,
                         trafo=trafo,
                         estGrid=est.grid,
                         fixedGrid=fixed.grid,
                         eventList=eventlist
                         )
filename <- "petab"
writePetab(pe, filename)
unlink(list.files(".", "\\.o$|\\.so$|\\.c$"))

              
ODEmodel=reactions
obsFun=observables
errormodel=errors
data=pred
bestfit=pouter
trafo=trafo
estGrid=est.grid
fixedGrid=fixed.grid
eventList=eventlist

# sbml_exportEquationList
equationList = args$equationList
events = args$events
parameterFormulaList = args$parameterFormulaList
modelname = "Model"
speciesInfo     = getSpeciesInfo(equationList)
parInfo         = getParInfo(equationList, eventList = events)
compartmentInfo = getCompartmentInfo(equationList)




# Exit ----
