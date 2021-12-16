try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

devtools::load_all("~/Promotion/Promotion/Projects/petab")
# # .. simpleSBML -----
equationList <- NULL
equationList <- addReaction(equationList, from = "E + S", to = "ES", rate = "(kon)*E*S",
                  description = "production of complex")
equationList <- addReaction(equationList, from = "ES", to = "E + S", rate = "koff*ES",
                  description = "decay of complex")
equationList <- addReaction(equationList, from = "ES", to = "E + P", rate = "kcat*ES",
                  description = "production of product")
equationList <- addReaction(equationList, "P", "", "P * vmax / (km + P)", "decay of p")
equationList <- eqnlist_addDefaultCompartment(equationList, "c1")


eventList <- eventlist(var = "E", time = 0, value = 1, root = NA, method = "replace")
eventList <- addEvent(eventList, var = "E", time = 1, value = 0, root = NA, method = "replace")

filename <- file.path(tempdir(), "model.xml")

parInfo <- data.table(tibble::tribble(
  ~parName, ~parValue, ~parUnit,
  "kon"   ,     0.10,"litre_per_mole_per_second" ,
  "koff"  ,     0.55,"per_second" ,
  "kcat"  ,     1.00,"per_second" ))

modelname = "Model"
# unitInfo = getUnitInfo()
compartmentInfo = getCompartmentInfo(equationList)
speciesInfo = getSpeciesInfo(equationList)
parInfo = getParInfo(equationList)
reactionInfo <- getReactionInfo(equationList,parInfo = parInfo)



pys <- petab_python_setup(FLAGreturnSimpleSBML = TRUE)


# Collect arguments

# setnames(unitInfo, old = c("speciesName", "compName", "initialAmount"))
setnames(compartmentInfo, old = c("compName", "compSize"), c("comp_id", "vol"))
setnames(speciesInfo, old = c("speciesName", "compName", "initialAmount"), new = c("species_id", "comp", "amt"))
setnames(parInfo, old = c("param_id", "val", "units"))
reactionInfo <- reactionInfo[,list(reactants, products = product, expression = equation, rxn_id = reactionName)]


# unitInfoList        <- purrr::transpose(unitInfo)
compartmentInfoList <- purrr::transpose(compartmentInfo)
speciesInfoList     <- purrr::transpose(speciesInfo)
parInfoList         <- purrr::transpose(parInfo)
reactionInfoList    <- purrr::transpose(reactionInfo)

# Start SBML document
model <- pys$SbmlModel()

# Populate with content
x <- (speciesInfoList)[[1]]
for (x in speciesInfoList)     do.call(model$addSpecies,       c(x))
for (x in parInfoList)         do.call(model$addParameter,   c(x))

# >>>> Fail: This parses the product ES as E and S. Conclusion: SimpleSBML is unusable for us <<<<<<<<<<< ----
model$addReaction("E", c('ES'), 'kon*ES')

