library(petab)
library(dMod)
library(dplyr)
library(stringr)
library(conveniencefunctions)

try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# -------------------------------------------------------------------------#
# Export Petab ----
# -------------------------------------------------------------------------#
# .. Create petab tables -----
# combine est.grid und fixed.grid
pe_ex <- petab_experimentalCondition(conditionId = c("C1", "C2", "C3"), conditionName = c("C1", "C2", "C3"),
                                     k1 = c(-1000, "k1_C2", "k1_C3"),
                                     k2 = c(-1000, "k2_C2", "k2_C3"))

pe_ob <- petab_observables(observableId = c("pERK_obs"),
                           observableName = c("pERK_obs"),
                           observableFormula = c("scale_pERK*pERK+offset_pERK"),
                           observableTransformation = "log10",
                           noiseFormula = c("0.1"),
                           noiseDistribution = c("normal"))
 
pe_me <- petab_measurementData(observableId = pred$name,
                               simulationConditionId = pred$condition,
                               measurement = pred$value,
                               time = pred$time,
                               observableParameters = NA_character_,
                               noiseParameters = pred$sigma,
                               datasetId = "data1",
                               replicateId = 1,
                               preequilibrationConditionId = NA_character_,
                               datapointId = 1:nrow(pred))
# .. error model -----
error_model <- errors
for(i in 1:length(error_model)){
  myobs <- names(error_model)[i]
  pe_me[observableId == myobs,`:=`(noiseParameters =  error_model[myobs])]
}

pe_ob[,`:=`(noiseFormula = paste0("noiseParameter1_", observableId))]


# .. observation parameters -----
observables <- observables
for(i in 1:length(observables)){
  myobs <- names(observables)[i]
  pe_me[observableId == myobs,`:=`(observableParameters =  observables[myobs])]
}

pe_me[observableId == "obsE",`:=`(observableParameters = "offset_E")]
pe_mo <- petab_model(el,events = NULL,parInfo = parInfo, speciesInfo = speciesInfo)


# .. Create petab -----
pe <- petab(model = pe_mo,
            experimentalCondition = pe_ex,
            measurementData = pe_me,
            observables = pe_ob)
pe$parameters <- petab_create_parameter_df(pe)

pe$parameters$objectivePriorType <- NA_character_

filename <- "petab"
writePetab(pe, filename)
unlink(list.files(".", "\\.o$|\\.so$|\\.c$"))

# Exit ----
