library(petab)
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# -------------------------------------------------------------------------#
# Create enzyme kinetics model and data ----
# -------------------------------------------------------------------------#
# .. Eqnlist and objects -----
modelname <- "petab"

el <- NULL
el <- addReaction(el, from = "A", to = "B", rate = "k_AB * A", description = "A to B")
el <- addReaction(el, from = "B", to = "C", rate = "k_BC * B", description = "B to C")
el <- eqnlist_addDefaultCompartment(el, "cytoplasm") # Need compartment information for SBML

parInfo <- data.table(tibble::tribble(
  ~parName, ~parValue, ~parUnit,
  "k_AB"   ,     1  ,"per_second" , 
  "k_BC"  ,     0.1,"per_second" ))

speciesInfo <- data.table(tibble::tribble(
  ~speciesName, ~compName, ~initialAmount,
  "A"         ,"cytoplasm" ,             1,          # Amount, not concentration
  "B"         ,"cytoplasm" ,             1e-12,
  "C"         ,"cytoplasm" ,             1e-12))

# compartmentInfo is left as the default getCompartmentInfo(el)
# unitInfo is left as the default getUnitInfo(): If you need other units, you need to add them

# .. Simulate Data -----
compiled <- odemodel(f = el,modelname = modelname)
x <- Xs(compiled, condition = "C1", optionsOde = list(atol = 1e-12,rtol = 1e-12))
pars <- c(setNames(parInfo$parValue, parInfo$parName),
          setNames(speciesInfo$initialAmount, speciesInfo$speciesName))

pred <- x(seq(0,50, 1), pars)
pred <- data.table(as.data.frame(pred))
pred <- rbind(pred,pred,pred)
pred <- pred[time %in% seq(25,50,5)]
pred <- pred[name == "C"]
set.seed(1)
pred[,`:=`(sigma = 0.1)]
pred[,`:=`(value = exp(log(value) + rnorm(length(value), sd = sigma)))]
pred[,`:=`(name = paste0("obs", name))]

# -------------------------------------------------------------------------#
# Export Petab ----
# -------------------------------------------------------------------------#
# .. Create petab tables -----
pe_ex <- petab_experimentalCondition("C1", "C1")
pe_ob <- petab_observables(observableId = c("obsC"),
                           observableName = c("obsC"),
                           observableFormula = c("C"),
                           observableTransformation = "log",
                           noiseFormula = c("0.1"),
                           noiseDistribution = c("normal"))
pe_me <- petab_measurementData(observableId = pred$name,
                               simulationConditionId = "C1",
                               measurement = pred$value,
                               time = pred$time,
                               observableParameters = NA_character_,
                               noiseParameters = pred$sigma,
                               datasetId = "data1",
                               replicateId = rep(1:3, each = nrow(pred)/3),
                               preequilibrationConditionId = NA_character_,
                               datapointId = 1:nrow(pred))
# .. error model -----
pe_ob[,`:=`(noiseFormula = paste0("noiseParameter1_", observableId))]
pe_me[,`:=`(noiseParameters = paste0("sigma_", observableId))]

# ..  -----
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
