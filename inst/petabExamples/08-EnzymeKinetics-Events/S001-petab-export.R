devtools::load_all("~/Promotion/Promotion/Projects/petab")
library(petab)
# try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

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


speciesInfo <- data.table(tibble::tribble(
  ~speciesName, ~compName, ~initialAmount,
  "E"         ,"cytoplasm" ,             1,          # Amount, not concentration
  "S"         ,"cytoplasm" ,             100,
  "ES"        ,"cytoplasm" ,             1e-12,
  "P"         ,"cytoplasm" ,             1e-12))

# compartmentInfo is left as the default getCompartmentInfo(el)
# unitInfo is left as the default getUnitInfo(): If you need other units, you need to add them

eventList <- eventlist(var = "ES", time = 20, value = 0, root = NA, method = "replace")
eventList <- addEvent(eventList, var = "ES", time = 50, value = "addES", root = NA, method = "replace")

parInfo <- data.table(tibble::tribble(
  ~parName, ~parValue, ~parUnit,
  "kon"   ,     1  ,"litre_per_mole_per_second" ,    # Because of compartment, all dMod-fluxes are multiplied with cytoplasm volume
  "koff"  ,     0.1,"per_second" ,
  "kcat"  ,     0.1,"per_second" ,
  "addES", 30, "mole_per_litre"))



# .. Simulate Data -----
compiled <- odemodel(f = el,modelname = modelname, events = eventList)
x <- Xs(compiled, condition = "C1")
pars <- c(setNames(parInfo$parValue, parInfo$parName),
          setNames(speciesInfo$initialAmount, speciesInfo$speciesName))

pred <- as.data.table(x(seq(0,100), pars)[[1]])
pred <- x(seq(0,100, 10), pars)
pred <- data.table(as.data.frame(pred))
pred <- pred[time > 0]
pred[,`:=`(sigma = 0.1)]
pred <- rbind(pred,pred,pred)
pred[,`:=`(value = exp(log(value) + rnorm(length(value), sd = sigma)))]
pred[,`:=`(name = paste0("obs", name))]
cfggplot(pred, aes(time, value)) + geom_point() + facet_wrap(~name, scales = "free")
# -------------------------------------------------------------------------#
# Export Petab ----
# -------------------------------------------------------------------------#
# .. Create petab tables -----
pe_ex <- petab_experimentalCondition("C1", "C1")
pe_ob <- petab_observables(observableId = c("obsE","obsS","obsES","obsP"),
                           observableName = c("obsE","obsS","obsES","obsP"),
                           observableFormula = c("E","S","ES","P"),
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
                               preequilibrationConditionId = NA_character_)
# .. error model -----
pe_ob[,`:=`(noiseFormula = paste0("noiseParameter1_", observableId))]
pe_me[,`:=`(noiseParameters = paste0("sigma_", observableId))]

# ..  -----
pe_me[observableId == "obsE",`:=`(observableParameters = "offset_E")]
pe_mo <- petab_model(el,events = eventList,parInfo = parInfo, speciesInfo = speciesInfo)


# .. Create petab -----
pe <- petab(model = pe_mo,
            experimentalCondition = pe_ex,
            measurementData = pe_me,
            observables = pe_ob)
pe$parameters <- petab_create_parameter_df(pe)

filename <- "petab"
# debugonce(sbml_exportEquationList)
writePetab(pe, filename)
unlink(list.files(".", "\\.o$|\\.so$|\\.c$"))


# Exit ----
