library(petab)
setwd(tempdir())
petab_python_setup()
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# -------------------------------------------------------------------------#
# Create enzyme kinetics model and data ----
# -------------------------------------------------------------------------#
# .. Eqnlist and objects -----
modelname <- "petab"

el <- NULL
el <- addReaction(el, from = "E + S", to = "ES", rate = "(kon)*E*S",
                  description = "production of complex")
el <- addReaction(el, from = "ES", to = "E + S", rate = "koff*ES",
                  description = "decay of complex")
el <- addReaction(el, from = "ES", to = "E + P", rate = "kcat*ES",
                  description = "production of product")
el <- eqnlist_addDefaultCompartment(el, "cytoplasm") # Need compartment information for SBML

parInfo <- data.table(tibble::tribble(
  ~parName, ~parValue, ~parUnit,
  "kon"   ,     1  ,"litre_per_mole_per_second" ,    # Because of compartment, all dMod-fluxes are multiplied with cytoplasm volume
  "koff"  ,     0.1,"per_second" ,
  "kcat"  ,     0.1,"per_second" ))

speciesInfo <- data.table(tibble::tribble(
  ~speciesName, ~compName, ~initialAmount,
  "E"         ,"cytoplasm" ,             1,          # Amount, not concentration
  "S"         ,"cytoplasm" ,             100,
  "ES"        ,"cytoplasm" ,             1e-12,
  "P"         ,"cytoplasm" ,             1e-12))

# compartmentInfo is left as the default getCompartmentInfo(el)
# unitInfo is left as the default getUnitInfo(): If you need other units, you need to add them

# .. Simulate Data -----
compiled <- odemodel(f = el,modelname = modelname)
x <- Xs(compiled, condition = "C1")
pars <- c(setNames(parInfo$parValue, parInfo$parName),
          setNames(speciesInfo$initialAmount, speciesInfo$speciesName))
# .... Condition 1 ------
parsC1 <- c(kon = 1, koff = 0.1, kcat = 0.1, E = 1, S = 100, ES = 1e-12, P = 1e-12)

pred <- x(seq(0,100, 10), parsC1)
pred <- data.table(as.data.frame(pred))
pred <- pred[time > 0]
pred[,`:=`(sigma = 0.1)]
pred <- rbind(pred,pred,pred)
pred[,`:=`(value = exp(log(value) + rnorm(length(value), sd = sigma)))]
pred[,`:=`(name = paste0("obs", name))]
predC1 <- pred

# .... Condition 2 ------
parsC1 <- c(kon = 1, koff = 0.1, kcat = 0.1, E = 1, S = 100, ES = 1e-12, P = 1e-12)
parsC2["kcat"] <- 0.5
parsC2["kon"]  <- 0.5
parsC2["koff"] <- 0.2
parsC2["E"]    <- 1.3

pred <- x(seq(0,100, 10), parsC2)
pred <- data.table(as.data.frame(pred))
pred <- pred[time > 0]
pred[,`:=`(sigma = 0.1)]
pred <- rbind(pred,pred,pred)
pred[,`:=`(value = exp(log(value) + rnorm(length(value), sd = sigma)))]
pred[,`:=`(name = paste0("obs", name))]
pred[,`:=`(condition = "C2")]
predC2 <- pred

# .. Combine -----
pred <- rbindlist(list(predC1,predC2))

# -------------------------------------------------------------------------#
# Export Petab ----
# -------------------------------------------------------------------------#
# .. Create petab tables -----
pe_ex <- petab_experimentalCondition(conditionId = c("C1", "C2"), conditionName = c("C1", "C2"))
pe_ob <- petab_observables(observableId = c("obsE","obsS","obsES","obsP"),
                           observableName = c("obsE","obsS","obsES","obsP"),
                           observableFormula = c("E","S","ES","P"),
                           observableTransformation = "log",
                           noiseFormula = c("0.1"),
                           noiseDistribution = c("normal"))
pe_me <- petab_measurementData(observableId = pred$name,
                               simulationConditionId = pred$condition,
                               measurement = pred$value,
                               time = pred$time,
                               observableParameters = NA_character_,
                               noiseParameters = pred$sigma,
                               datasetId = "data1",
                               replicateId = rep(1:3, each = nrow(pred)/3),
                               preequilibrationConditionId = NA_character_)

# One error parameter
pe_ob[,`:=`(noiseFormula = paste0("noiseParameter1_", observableId))]
pe_me[,`:=`(noiseParameters = paste0("sigma_", observableId))]
# Two error parameters
pe_ob[observableId == "obsE",`:=`(noiseFormula = paste0("noiseParameter1_", observableId, "*", observableId, " + ", "noiseParameter2_", observableId))]
pe_me[observableId == "obsE",`:=`(noiseParameters = paste0("sigmaRel_", observableId, ";","sigmaAbs_", observableId))]

# One observable parameter
pe_me[observableId == "obsE",`:=`(observableParameters = "offset_E")]

# Condition specific observable parameters and error parameters
pe_me[observableId == "obsS",`:=`(observableParameters = "offset_S")]
pe_me[observableId == "obsS" & simulationConditionId == "C2",`:=`(noiseParameters = paste0("sigma_", observableId, "_C2"),
                                                                  observableParameters = "offset_S_C2")]


pe_mo <- petab_model(el,events = NULL,parInfo = parInfo, speciesInfo = speciesInfo)

# .. Create petab -----
pe <- petab(model = pe_mo,
            experimentalCondition = pe_ex,
            measurementData = pe_me,
            observables = pe_ob)

pe$parameters <- petab_create_parameter_df(pe)
pe$parameters[grep("kon", parameterId),`:=`(estimate = 0)]

filename <- "basemodel"
writePetab(pe, filename)
unlink(list.files(".", "\\.o$|\\.so$|\\.c$"))
# Exit ----
