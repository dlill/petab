# -------------------------------------------------------------------------#
# dMod helper functions ----
# -------------------------------------------------------------------------#


#' Title
#'
#' @param el 
#' @param compName 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
#' # see getReactionInfo
eqnlist_addDefaultCompartment <- function(equationList, compName) {
  as.eqnlist(as.data.frame(equationList), 
             volumes = setNames(rep(compName, length(equationList$states)), equationList$states))
}


#' Extract parInfo data.table from equationlist
#'
#' @param el equationList
#'
#' @return data.table(parName,parValue,parUnit)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom data.table data.table
#' @family SBML export
#'
#' @examples
#' library(dMod)
#' example(eqnlist)
#' getParInfo(f)
getParInfo <- function(equationList) {
  parName <- setdiff(getParameters(equationList), c(equationList$states, equationList$volumes))
  parInfo <- data.table::data.table(parName = parName)
  parInfo[,`:=`(parValue = seq(0.1,1,length.out = .N))]
  parInfo[,`:=`(parUnit = "per_second")]
  parInfo
}



#' Extract speciesInfo data.table from equationlist
#'
#' @param el equationList
#'
#' @return data.table(speciesName, compName, initialAmount)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
#' library(dMod)
#' example(eqnlist)
#' getSpeciesInfo(f)
getSpeciesInfo <- function(equationList){
  data.table(speciesName = equationList$states,
             compName    = equationList$volumes,
             initialAmount = 0)
}

#' get compartmentInfo
#'
#' @param el equationlist
#'
#' @return data.table(compName, compSize)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
#' library(dMod)
#' example(eqnlist)
#' getCompartmentInfo(f)
getCompartmentInfo <- function(equationList) {
  data.table(compName = unique(equationList$volumes),
                   compSize = 1)
}


#' Title
#'
#' @param el 
#' @param p 
#'
#' @return list of reactions
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#' @importFrom cOde getSymbols
#' 
#' @examples 
#' 
#' # simple
#' el <- NULL
#' el <- addReaction(el, from = "E + S", to = "ES", rate = "(kon)*E*S",
#'                   description = "production of complex")
#' el <- addReaction(el, from = "ES", to = "E + S", rate = "koff*ES",
#'                   description = "decay of complex")
#' el <- addReaction(el, from = "ES", to = "E + P", rate = "kcat*ES",
#'                   description = "production of product")
#' el <- eqnlist_addDefaultCompartment(el, "cytoplasm") # Need compartment information for SBML
#' equationList <- el
#' parInfo = getParInfo(equationList)
#' 
#' # complex stoichiometries, modifier in production of complex
#' el <- NULL
#' el <- addReaction(el, from = "2*E + S", to = "ES", rate = "(kon)*E*S * 1/(1+kinh*P)",
#'                   description = "production of complex")
#' el <- addReaction(el, from = "ES", to = "2*E + S", rate = "koff*ES",
#'                   description = "decay of complex")
#' el <- addReaction(el, from = "ES", to = "E + P", rate = "kcat*ES",
#'                   description = "production of product")
#' el <- eqnlist_addDefaultCompartment(el, "cytoplasm") # Need compartment information for SBML
#' equationList <- el
#' parInfo = getParInfo(equationList)
#' 
#' getReactionInfo(equationList,parInfo)
getReactionInfo <- function(equationList, parInfo = getParInfo(equationList)) {
  re <- getReactions(equationList)
  
  s <- equationList$smatrix
  rp <- lapply(1:nrow(s), function(i) {
    si <- s[i,,drop = TRUE]
    si <- si[!is.na(si)]
    reactants <- abs(si[si<0])
    products  <- abs(si[si>0])
    list(reactants = list(name = names(reactants), stoichiometry = reactants), 
         products = list(name  = names(products), stoichiometry = products))
  })
  rp <- purrr::transpose(rp)
  
  re <- data.table(re)
  re[,`:=`(reactants = rp$reactants)] # reactants is a list
  re[,`:=`(product  = rp$products)]
  
  fl <- getFluxes(equationList, type = "amount")
  fl <- rbindlist(lapply(fl, function(f) data.table(Description = names(f), Flux = f)), idcol = "Species")
  fl <- fl[,list(Description, Flux)]
  fl[,`:=`(Flux = gsub("^[+-]?\\d+\\*", "", Flux))]
  fl <- unique(fl)
  if (nrow(fl) != nrow(re)) stop("fluxes and reactionlist don't match")
  re <- fl[re, on = c("Description")]
  
  re[,`:=`(reactionName = gsub(" ", "", Description))]
  re[,`:=`(parName  = lapply(Rate,    function(x)  setdiff(cOde::getSymbols(x), equationList$states)))]
  re[,`:=`(parValue = lapply(parName, function(pn) parInfo[parName %in% pn, parValue]))]
  re[,`:=`(parUnit  = lapply(parName, function(pn) parInfo[parName %in% pn, parUnit]))]
  re[,`:=`(reversible = 0)]
  re[,`:=`(equation = copy(Flux))]
  re[,`:=`(modifiers = lapply(seq_along(Flux), function(i) {
    modifiers <- intersect(equationList$states, cOde::getSymbols(Flux[i]))
    modifiers <- setdiff(modifiers, reactants[[i]]$name)
    modifiers <- setdiff(modifiers, product[[i]]$name) # questionable if products are not modifiers?
    modifiers
  }))]
  
  
  re <- re[,list(reactants,product,reactionName,parName,parValue,parUnit,reversible,equation,modifiers)] 
  re
}

#' Assemble unitInfo
#' 
#' Some default units, some others
#' 
#' @param unitName character vector
#' @param unitKind list of character vectors
#' @param exponent list of numeric vectors
#'
#' @return data.table(unitName, unitKind, exponent)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
#' getUnitInfo()
getUnitInfo <- function(unitName = NULL, unitKind = NULL, exponent = NULL) {
  # Default list of units
  unitInfo <- data.table(tibble::tribble(
    ~unitName, ~unitKind, ~exponent,
    "per_second"                      ,c("UNIT_KIND_SECOND")                                           ,c(-1)     ,
    "per_mole"                        ,c("UNIT_KIND_MOLE")                                             ,c(-1)     ,
    "per_mole_per_second"             ,c("UNIT_KIND_MOLE"  , "UNIT_KIND_SECOND")                       ,c(-1,-1) ,
    "litre_per_mole_per_second"       ,c("UNIT_KIND_LITRE",  "UNIT_KIND_MOLE"  , "UNIT_KIND_SECOND")   ,c(1,-1,-1),
    "litre2_per_mole2_per_second"     ,c("UNIT_KIND_LITRE",  "UNIT_KIND_MOLE"  , "UNIT_KIND_SECOND")   ,c(2,-2,-1),
    "mole_per_litre_per_second"       ,c("UNIT_KIND_MOLE",  "UNIT_KIND_LITRE"  , "UNIT_KIND_SECOND")   ,c(1,-1,-1),
    "mole_per_litre"                  ,c("UNIT_KIND_MOLE",  "UNIT_KIND_LITRE")                         ,c(1,-1),
    "mole_per_second"                  ,c("UNIT_KIND_MOLE",  "UNIT_KIND_SECOND")                       ,c(1,-1),
    "identity"                        ,c("UNIT_KIND_DIMENSIONLESS")                                    ,c(1)
  ))
  
  # Add custom
  unitInfo <- rbindlist(list(
    unitInfo,
    data.table(unitName = unitName, unitKind = unitKind, exponent = exponent)
  ))
  unitInfo
}

#' Title
#'
#' @param eventList [dMod::eventList()]
#'
#' @return data.table()
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family 
#' @importFrom data.table data.table
#'
#' @examples
getEventInfo <- function(eventList) {
  if (any(!is.na(eventList$root))) stop("root in events not yet supported")
  if (any(eventList$method != "replace")) stop("only replace events are supported")
  eventInfo <- data.table::data.table(eventList)
  eventInfo <- eventInfo[, list(eventSpecies = var, eventTrigger = paste0("time >= ", time), eventFormula = value)]
  eventInfo
}




# -------------------------------------------------------------------------#
# libSBML helper functions ----
# -------------------------------------------------------------------------#

#' Initialize the sbml document
#'
#' @param modelname String
#' @param sbmlDoc return value from a call like `SBMLDocument(level = 2, version = 4)`
#'
#' @return a libsbml model representation
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @family SBML export
#'
#' @examples
sbml_initialize <- function(modelname = "EnzymaticReaction", sbmlDoc) {
  model = SBMLDocument_createModel(sbmlDoc)
  Model_setId(model, modelname)
  model
}



#' Add unit
#'
#' @param model the sbml model
#' @param unitName string, used as identifier
#' @param unitKind vector of involved units
#' @param exponent vector of exponents associated to units
#'
#' @return called for side-effect, but `model` is returned invisibly
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
sbml_addOneUnit <- function(model, unitName = "litre_per_mole_per_second", 
                            unitKind =c("UNIT_KIND_LITRE", "UNIT_KIND_MOLE", "UNIT_KIND_SECOND"),
                            exponent = c(1,-1,-1)) {
  unitdef = Model_createUnitDefinition(model)
  UnitDefinition_setId(unitdef, unitName)
  for (i in seq_along(unitKind)){
    unit = UnitDefinition_createUnit(unitdef)
    Unit_setKind(unit, unitKind[i])
    Unit_setExponent(unit,exponent[i])
  }
  invisible(model)
}



#' Add compartment
#'
#' @param model  the sbml model
#' @param compName string, e.g. "cytoplasm"
#' @param compsize 
#'
#' @return called for side-effect
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
sbml_addOneCompartment <- function(model, compName, compSize) {
  # [ ] Compartment units?
  comp = Model_createCompartment(model)
  Compartment_setId(comp,compName)
  Compartment_setSize(comp,compSize)
}


#' Title
#'
#' @param model 
#' @param speciesName 
#' @param compName 
#' @param speciesInitialAmount 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
sbml_addOneSpecies <- function(model, speciesName, compName, initialAmount) {
  sp = Model_createSpecies(model)  
  Species_setCompartment(sp,compName)
  Species_setId(sp,speciesName)
  Species_setName(sp,speciesName)
  Species_setInitialAmount(sp,initialAmount)
}

#' Title
#'
#' @param model 
#' @param speciesName 
#' @param compName 
#' @param speciesInitialAmount 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
sbml_addOneParameter <- function(model, parName, parValue, parUnit) {
  parm = Model_createParameter(model)  
  Parameter_setId(parm, parName)
  Parameter_setValue(parm, parValue)
  Parameter_setUnits(parm, parUnit)
}

#' Title
#'
#' @param reaction 
#' @param speciesNames 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
sbml_reactionAddReactants <- function(reaction, reactants) {
  reactants <- purrr::transpose(reactants)
  for (x in reactants){
    spr = Reaction_createReactant(reaction)
    SimpleSpeciesReference_setSpecies(spr, x$name)
    SpeciesReference_setStoichiometry(spr, x$stoichiometry)
  }
}

#' Title
#'
#' @param reaction 
#' @param products list(name = c("a", "b"), stoichiometry = c(1,2))
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
sbml_reactionAddProducts <- function(reaction, products) {
  products <- purrr::transpose(products)
  for (x in products){
    spr = Reaction_createProduct(reaction)
    SimpleSpeciesReference_setSpecies(spr, x$name)
    SpeciesReference_setStoichiometry(spr, x$stoichiometry)
  }
}
#' Title
#'
#' @param reaction 
#' @param products list(name = c("a", "b"), stoichiometry = c(1,2))
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
sbml_reactionAddModifiers <- function(reaction, modifiers) {
  for (x in modifiers){
    spr = Reaction_createModifier(reaction)
    SimpleSpeciesReference_setSpecies(spr, x)
  }
}

#' Title
#'
#' @param kl 
#' @param parName 
#' @param parValue 
#' @param unitName 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
sbml_kineticLawAddParameters <- function(kl, parName, parValue, unitName) {
  para = KineticLaw_createParameter(kl)
  for (i in seq_along(parName)){
    Parameter_setId(para, parName[i])
    Parameter_setValue( para, parValue[i])
    Parameter_setUnits( para, unitName[i])}
}

#' Title
#'
#' @param reaction 
#' @param equation 
#' @param parName 
#' @param parValue 
#' @param parUnit 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
sbml_reactionAddKineticLaw <- function(reaction, equation, parName, parValue, parUnit) {
  kl = Reaction_createKineticLaw(reaction)
  astMath <- parseFormula(equation)
  KineticLaw_setMath( kl, astMath)  
  # sbml_kineticLawAddParameters(kl, parName, parValue, parUnit) # solved by global parameters?
}

#' Title
#'
#' @param model 
#' @param reactionName 
#' @param reversible 
#' @param reactants 
#' @param products 
#' @param equation 
#' @param parName 
#' @param parValue 
#' @param parUnit 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
sbml_addOneReaction <- function(model, reactionName,
                                reversible = 1, reactants, products,
                                equation, parName, parValue, parUnit,
                                modifiers
) {
  reaction = Model_createReaction(model)
  Reaction_setId(reaction,reactionName)
  if (reversible==0) Reaction_setReversible(reaction, reversible)
  sbml_reactionAddReactants(reaction, reactants)
  sbml_reactionAddProducts(reaction, products)
  sbml_reactionAddModifiers(reaction, modifiers)
  sbml_reactionAddKineticLaw(reaction, equation, parName, parValue, parUnit)
}


#' Title
#' 
#' Rather generic already with eventTrigger and eventFormula,
#' but I don't know if import can handle anything else but "replace"-events at the moment
#'
#' @param model 
#' @param eventSpecies 
#' @param eventTrigger 
#' @param eventFormula 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML export
#'
#' @examples
sbml_addOneEvent <- function(model, eventSpecies, eventTrigger, eventFormula) {
  # Create event
  ev <- Event(3,2) # sbml level,version, I guess ??
  
  # Trigger
  Event_createTrigger(ev)
  tr <- Trigger(3,2)
  Trigger_setInitialValue(tr, FALSE)
  astMath <- parseL3Formula(eventTrigger)
  Trigger_setMath(tr, astMath)
  Event_setTrigger(ev, tr)
  
  # Delay # Possibility to add in future
  # delay <- Delay(3,2)
  # Delay_setMath(delay, parseL3Formula("0"))
  # Event_setDelay(ev,delay)
  
  # Assignment
  ea <- EventAssignment(3,2)
  EventAssignment_setVariable(ea,eventSpecies)
  EventAssignment_setMath(ea, parseL3Formula(eventFormula))
  Event_addEventAssignment(ev,ea)
  
  Model_addEvent(model, ev)
}


#' From libSBML example
#' 
#' @author Frank Bergmann
#'
#' @param sbmlDoc 
#'
#' @return
#' @export
#' @md
#' @family SBML export
#'
#' @examples
sbml_validateSBML <- function(sbmlDoc)
{
  noProblems             = 1;
  numCheckFailures       = 0;
  numConsistencyErrors   = 0;
  numConsistencyWarnings = 0;
  numValidationErrors    = 0;
  numValidationWarnings  = 0;
  
  # LibSBML 3.3 is lenient when generating models from scratch using the
  # API for creating objects.  Once the whole model is done and before it
  # gets written out, it's important to check that the whole model is in
  # fact complete, consistent and valid.
  
  numCheckFailures = SBMLDocument_checkInternalConsistency( sbmlDoc);
  if ( numCheckFailures > 0 ) {
    noProblems = 0;
    for (i in 0:(numCheckFailures-1)) {
      sbmlErr = SBMLDocument_getError( sbmlDoc, i);
      if ( XMLError_isFatal(sbmlErr) || XMLError_isError(sbmlErr) ) {
        numConsistencyErrors = numConsistencyErrors + 1;
      } else {
        numConsistencyWarnings = 1 + numConsistencyWarnings;
      }      
    } 
    
    SBMLDocument_printErrors(sbmlDoc);
  }
  
  # If the internal checks fail, it makes little sense to attempt
  # further validation, because the model may be too compromised to
  # be properly interpreted.
  
  if (numConsistencyErrors > 0) {
    cat("Further validation aborted.\n"); 
  } else {
    numCheckFailures = SBMLDocument_checkConsistency( sbmlDoc );
    if ( numCheckFailures > 0 ) {
      noProblems = 0;
      for (i in 0:(numCheckFailures-1)) {
        sbmlErr = SBMLDocument_getError( sbmlDoc, i);
        if ( XMLError_isFatal( sbmlErr) || XMLError_isError( sbmlErr) ) {
          numValidationErrors = 1+ numValidationErrors;
        } else {
          numValidationWarnings = 1+ numValidationWarnings;
        }      
      } 
      SBMLDocument_printErrors(sbmlDoc);
    }
  }
  
  if (noProblems) {
    return (1);
  } else {
    if (numConsistencyErrors > 0) {
      cat("ERROR: encountered ",numConsistencyErrors," consistency error(s) in model '",Model_getId(SBMLDocument_getModel( sbmlDoc)),"'.\n")
    }
    if (numConsistencyWarnings > 0) {
      cat( "Notice: encountered ",numConsistencyWarnings," consistency warning(s) in model '",Model_getId(SBMLDocument_getModel( sbmlDoc)),"'.\n")
    }
    
    if (numValidationErrors > 0) {
      cat("ERROR: encountered ",numValidationErrors," validation error(s) in model '",Model_getId(SBMLDocument_getModel( sbmlDoc)),"'.\n")
    }
    if (numValidationWarnings > 0) {
      cat( "Notice: encountered ",numValidationWarnings," validation warning(s) in model '",Model_getId(SBMLDocument_getModel( sbmlDoc)),"'.\n")
    }
    
    return (numConsistencyErrors == 0 && numValidationErrors == 0);
  }
}


# -------------------------------------------------------------------------#
# sbml export function ----
# -------------------------------------------------------------------------#


#' Write sbml model to filename
#'
#' @param modelname 
#' @param filename 
#' @param equationList 
#' @param unitInfo 
#' @param speciesInfo 
#' @param parInfo 
#' @param compartmentInfo 
#'
#' @return NULL, called for side effect 
#' @export
#' @md
#' @family SBML export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @importFrom tibble tribble
#' @importFrom purrr transpose
#'
#' @examples
#' equationList <- NULL
#' equationList <- addReaction(equationList, from = "E + S", to = "ES", rate = "(kon)*E*S",
#'                             description = "production of complex")
#' equationList <- addReaction(equationList, from = "ES", to = "E + S", rate = "koff*ES",
#'                             description = "decay of complex")
#' equationList <- addReaction(equationList, from = "ES", to = "E + P", rate = "kcat*ES",
#'                             description = "production of product")
#' equationList <- eqnlist_addDefaultCompartment(equationList, "cytoplasm")
#' 
#' filename <- modelname <- file.path(tempdir(), "model.xml")
#' 
#' parInfo <- data.table(tibble::tribble(
#'   ~parName, ~parValue, ~parUnit,
#'   "kon"   ,     0.10,"litre_per_mole_per_second" ,
#'   "koff"  ,     0.55,"per_second" ,
#'   "kcat"  ,     1.00,"per_second" ))
#' unitInfo <- getUnitInfo()
#' compartmentInfo <- getCompartmentInfo(equationList)
#' speciesInfo <- getSpeciesInfo(equationList)
#' 
#' eventList <- eventlist(var = "E", time = 0, value = 1, root = NA, method = "replace")
#' eventList <- addEvent(eventList, var = "E", time = 1, value = 0, root = NA, method = "replace")
#' eventInfo <- getEventInfo(eventList)
#' 
#' sbml_exportEquationList(equationList, filename, parInfo = parInfo)
#' sbml_exportEquationList(equationList, filename, parInfo = parInfo, eventList = eventList)
sbml_exportEquationList <- function(equationList,
                                    filename,
                                    modelname = "Model",
                                    unitInfo        = getUnitInfo(),
                                    speciesInfo     = getSpeciesInfo(equationList),
                                    parInfo         = getParInfo(equationList),
                                    compartmentInfo = getCompartmentInfo(equationList),
                                    eventList = NULL) {
  
  # Load libSBML
  library(libSBML)
  
  # Collect arguments
  reactionInfo <- getReactionInfo(equationList,parInfo = parInfo)
  if (!is.null(eventList)) eventInfo <- getEventInfo(eventList)
  unitInfoList        <- purrr::transpose(unitInfo)
  speciesInfoList     <- purrr::transpose(speciesInfo)
  parInfoList         <- purrr::transpose(parInfo)
  compartmentInfoList <- purrr::transpose(compartmentInfo)
  reactionInfoList    <- purrr::transpose(reactionInfo)
  if (!is.null(eventList)) eventInfoList    <- purrr::transpose(eventInfo)
  
  # Start SBML document
  sbmlDoc = SBMLDocument(level = 2, version = 4) 
  model <- sbml_initialize(modelname, sbmlDoc)
  
  # Populate with content
  for (x in unitInfoList)        do.call(sbml_addOneUnit,        c(list(model = model),x))
  for (x in compartmentInfoList) do.call(sbml_addOneCompartment, c(list(model = model),x))
  for (x in speciesInfoList)     do.call(sbml_addOneSpecies,     c(list(model = model),x))
  for (x in parInfoList)         do.call(sbml_addOneParameter,   c(list(model = model),x))
  for (x in reactionInfoList)    do.call(sbml_addOneReaction,    c(list(model = model),x))
  if (!is.null(eventList)) for (x in eventInfoList)       do.call(sbml_addOneEvent,       c(list(model = model),x))
  
  # Promote parameters to global parameters - this is an sbml thingy
  props = ConversionProperties();
  ConversionProperties_addOption(props, "promoteLocalParameters", TRUE, "Promotes all Local Parameters to Global ones");
  
  # Validate and write to file
  sbml_validateSBML(sbmlDoc) 
  writeSBML(sbmlDoc, filename = filename)
  invisible(NULL)
}
