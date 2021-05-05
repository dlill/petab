library(petab)
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

.outputFolder <- paste0("Outputs")
for(folder in c(.outputFolder))
if(!dir.exists(folder)) dir.create(folder)
# -------------------------------------------------------------------------#
# Create enzyme kinetics model and data ----
# -------------------------------------------------------------------------#

# >>>> Model description <<<<<<<<<<< ----
# * Substrate has turnover (production and degradation)
# * Enzyme is added at time = 0 by an event (preequilibration with E=0)
# * Apart from that, normal enzyme kinetics

# .. Eqnlist and objects -----
modelname <- "petab"

el <- NULL
el <- addReaction(el, from = "", to = "S", rate = "kproS",
                  description = "Production of substrate")
el <- addReaction(el, from = "S", to = "", rate = "kdegS*S",
                  description = "Degradation of substrate")
el <- addReaction(el, from = "E + S", to = "ES", rate = "(kon)*E*S",
                  description = "production of complex")
el <- addReaction(el, from = "ES", to = "E + S", rate = "koff*ES",
                  description = "decay of complex")
el <- addReaction(el, from = "ES", to = "E + P", rate = "kcat*ES",
                  description = "production of product")
el <- eqnlist_addDefaultCompartment(el, "cytoplasm") # Need compartment information for SBML

events <- eventlist()
events <- addEvent(events, var = "E", time = 0, value = "Eadd", method = "replace")

# .. SteadyStates  -----

# # Trafos
# # est -> beforeSS -> SS -> x -> g
#
# # condition noE
# # noE
# trafoSS <- c(S = "kproS / kdegS", E = 0, C = 0, P = 0)
# # E
# trafoSS <- c(S = 0,...)

# .. Infos -----

parInfo <- data.table(tibble::tribble(
  ~parName, ~parValue, ~parUnit,
  "kon"   , 1        ,  "litre_per_mole_per_second" ,    # Because of compartment, all dMod-fluxes are multiplied with cytoplasm volume
  "koff"  , 0.1      ,  "per_second" ,
  "kcat"  , 1        ,  "per_second" ,
  "kproS" , 2        ,  "mole_per_litre_per_second" ,
  "kdegS" , 0.1      ,  "per_second",
  "Eadd"  , 1        ,  "mole_per_litre"
  ))

speciesInfo <- data.table(tibble::tribble(
  ~speciesName, ~compName, ~initialAmount,
  "E"         ,"cytoplasm" ,             0,          # Amount, not concentration
  "S"         ,"cytoplasm" ,             100,
  "ES"        ,"cytoplasm" ,             0,
  "P"         ,"cytoplasm" ,             0))

# compartmentInfo is left as the default getCompartmentInfo(el)
# unitInfo is left as the default getUnitInfo(): If you need other units, you need to add them

# .. Compile and plot model -----
compiled <- odemodel(f = el,modelname = modelname, events = events)
x        <- Xs(compiled, condition = "C1")
pars     <- c(setNames(parInfo$parValue, parInfo$parName),
              setNames(speciesInfo$initialAmount, speciesInfo$speciesName))

pl <- plot(x(seq(-10,30, 10), pars))
ggsave(file.path(.outputFolder, "01-BasicSimulation.png"), pl, width = 15.5, height = 10, scale = 1, units = "cm")

# .. Simulate Data -----
pred <- x(unique(c(seq(0,10), seq(0,100, 10))), pars)
pred <- data.table(as.data.frame(pred))
pred <- pred[time > 0]
pred[,`:=`(sigma = 0.1)]
pred <- rbind(pred,pred,pred)
pred[,`:=`(value = exp(log(value) + rnorm(length(value), sd = sigma)))]
pred[,`:=`(name = paste0("obs", name))]

pl <- ggplot(pred, aes(time, value, color = name)) + geom_point() + theme_bw() + scale_y_log10()
ggsave(file.path(.outputFolder, "02-Data.png"), pl, width = 15.5, height = 10, scale = 1, units = "cm")

# -------------------------------------------------------------------------#
# Export Petab ----
# -------------------------------------------------------------------------#
# .. Create petab tables -----
pe_ex <- petab_experimentalCondition(conditionId = c("NoEnzyme", "Enzyme"), conditionName = c("NoEnzyme", "Enzyme"), E = c(0,"Eadd"))
pe_ob <- petab_observables(observableId = c("obsE","obsS","obsES","obsP"),
                           observableName = c("obsE","obsS","obsES","obsP"),
                           observableFormula = c("E","S","ES","P"),
                           observableTransformation = "lin", # For simplicity
                           noiseFormula = c("0.1"),
                           noiseDistribution = c("normal"))
pe_ob[,`:=`(noiseFormula = paste0("noiseParameter1_", observableId, "*", observableId))]
pe_me <- petab_measurementData(observableId = pred$name,
                               simulationConditionId = "Enzyme",
                               measurement = pred$value,
                               time = pred$time,
                               observableParameters = NA_character_,
                               noiseParameters = pred$sigma,
                               datasetId = "data1",
                               replicateId = rep(1:3, each = nrow(pred)/3),
                               preequilibrationConditionId = "NoEnzyme")
pe_me[,`:=`(noiseParameters = paste0("sigma_", observableId))]
pe_me[observableId == "obsE",`:=`(observableParameters = "offset_E")]


pe_mo <- petab_model(el,events = events,parInfo = parInfo, speciesInfo = speciesInfo)


# .. Create petab -----
pe <- petab(model = pe_mo,
            experimentalCondition = pe_ex,
            measurementData = pe_me,
            observables = pe_ob)
pe$parameters <- petab_create_parameter_df(pe)

filename <- "petab"
writePetab(pe, filename)
unlink(list.files(".", "\\.o$|\\.so$|\\.c$"))
# Exit ----