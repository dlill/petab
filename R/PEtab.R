# -------------------------------------------------------------------------#
# Convenience functions ----
# -------------------------------------------------------------------------#

#' Original petab function didn't work
#'
#' @param pe petab without a parameter_df
#' @param observableParameterScale "lin", "log" or "log10"
#'
#' @return [petab_parameters()] data.table
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family 
#' @importFrom cOde getSymbols
#' @importFrom data.table rbindlist
#'
#' @examples
#' pe <- petab_exampleRead("01", "pe")
#' pe$parameters <- NULL
#' petab_create_parameter_df(pe)
petab_create_parameter_df <- function(pe, observableParameterScale = "log10") {
  
  model                 <- pe$model
  measurementData       <- pe$measurementData
  experimentalCondition <- pe$experimentalCondition
  
  # Species
  speciesInfo <- model$speciesInfo
  par_sp <- petab_parameters(parameterId =   speciesInfo$speciesName,
                             parameterName = speciesInfo$speciesName,
                             nominalValue =  speciesInfo$initialAmount,
                             estimate = as.numeric(speciesInfo$initialAmount > 0),
                             parameterScale = ifelse(speciesInfo$initialAmount > 0, "log10", "lin")
  )
  # Kinetic parameters in ODEs
  parInfo <- model$parInfo
  par_pa <- petab_parameters(parameterId =   parInfo$parName,
                             parameterName = parInfo$parName,
                             nominalValue =  parInfo$parValue)
  
  # HACK: Should this be here or further down? I put it here so the dynamic model is not interrupted
  if (length(cOde::getSymbols(pe$meta$events$value))){
    par_ev <- petab_parameters(parameterId   = cOde::getSymbols(pe$meta$events$value),
                               parameterName = cOde::getSymbols(pe$meta$events$value),
                               nominalValue  = 1,
                               estimate      = 0,
                               parameterScale = "lin") # up to debate
    par_pa <- data.table::rbindlist(list(par_pa, par_ev))
  }
  
  # Observable parameters
  par_ob <- NULL
  if (length(cOde::getSymbols(measurementData$observableParameters)))
    par_ob <- petab_parameters(parameterId =  cOde::getSymbols(measurementData$observableParameters),
                               parameterName = cOde::getSymbols(measurementData$observableParameters),
                               parameterScale = observableParameterScale)
  
  # MeasurementErrors
  par_meErr <- NULL
  if (length(cOde::getSymbols(measurementData$noiseParameters)))
    par_meErr <- petab_parameters(parameterId =   cOde::getSymbols(measurementData$noiseParameters),
                                  parameterName = cOde::getSymbols(measurementData$noiseParameters),
                                  nominalValue = 0.1)
  
  # Get all base-parameters
  par <- data.table::rbindlist(list(par_sp, par_pa, par_ob, par_meErr))
  
  
  # Parameters from experimentalConditions
  # More complicated, need also to exclude colnames
  #   from previously collected parameters
  par_ec <- NULL
  parnamesOuter <- petab_getParametersExperimentalCondition(experimentalCondition)
  if (length(parnamesOuter)){
    # Collect ec_pars
    par_ec <- petab_parameters(parameterId = parnamesOuter,
                               parameterName = parnamesOuter)
    
    # Adjust scales according to inner parameters
    parnamesInner <- setdiff(colnames(experimentalCondition), c("conditionId", "conditionName"))
    for (pxinner in parnamesInner) {
      pxouter <- cOde::getSymbols(experimentalCondition[[pxinner]])
      parameterScales <- par[parameterId == pxinner,parameterScale]
      if (!length(parameterScales)) next # In this case it's probably a L1-parameter which is dealt with later
      par_ec[parameterId %in% pxouter,`:=`(parameterScale = parameterScales)]
    }
    
    # Remove base parameters from par
    par <- par[!parameterId %in% parnamesInner]
    
    # Append par_ec
    par <- data.table::rbindlist(list(par, par_ec))
  }
  
  pfi <- pe$meta$parameterFormulaInjection
  if (!is.null(pfi)) {
    
    # parameterFormulaInjection by L1
    # Ensure parameterScale and nominalValue is set correctly for L1 parameters
    pfi_L1 <- pfi[trafoType == "L1"]
    if (nrow(pfi_L1)) {
      for (px in pfi_L1$parameterId) {
        # Set scale of L1 parameters
        pxscale <- par[parameterId == px, parameterScale]
        if (!length(pxscale)) pxscale <- par[parameterId == paste0("L1Ref_", px), parameterScale] # Might be that a parameter was already renamed
        if (!length(pxscale)) warning(px, " has no associated scale")
        par[grepl(paste0("L1_", px), parameterId),`:=`(parameterScale = pxscale)] # why grepl? fragile...
        par[grepl(paste0("L1_", px), parameterId) & parameterScale != "lin",`:=`(nominalValue = 1)]
        par[grepl(paste0("L1_", px), parameterId) & parameterScale == "lin",`:=`(nominalValue = 0)]
        par[grepl(paste0("L1_", px), parameterId),`:=`(objectivePriorType = "parameterScaleNormal")] # [ ] Should be parameterScaleLaplace
        par[grepl(paste0("L1_", px), parameterId),`:=`(objectivePriorParameters = "0;10")] 
        # Update names of inner parameters to L1Ref_par, so their names are set correctly in parameters
        par[parameterId==px,`:=`(parameterId = paste0("L1Ref_", parameterId))]
      }}
    
    parnamesPFI <- setdiff(cOde::getSymbols(pfi$parameterFormula), c(par$parameterId, names(pe$experimentalCondition)))
    if (length(parnamesPFI)){
      warning("The following parameters appear in a parameterFormulaInjection but are otherwise unset: ", paste0(parnamesPFI, collapse = ","))
      par_pfi <- petab_parameters(parameterId = parnamesPFI)
      par <- rbindlist(list(par, par_pfi))
    }
    
    # Remove inner parameters which are set by injection from the petab_parameters
    par <- par[!parameterId %in% pfi$parameterId]
    
  }
  
  par
}


#' dput for petab variables
#'
#' @param pe [petab()]
#' @param variable (single) character denoting a column, e.g. "observableId"
#'
#' @return dput output
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
petab_dput <- function(pe, variable) {
  dx <- copy(pe$measurementData)
  dx <- pe$experimentalCondition[dx, on = c("conditionId" = "simulationConditionId")]
  dx <- pe$observables[dx, on = c("observableId")]
  dput(unique(dx[[variable]]))
}


#' Get list of column names in petab
#'
#' @param pe NULL: Default names, [petab()] Actual names present in petab
#'
#' @return list of column names
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
petab_columns <- function(pe = NULL) {
  
  if(!is.null(pe)) {
    pc <- lapply(pe[c("experimentalCondition","measurementData","observables","parameters")], names)
    
    # If a "pattern" is given, add those hypothetical column names as well
    metaInformation <- petab_metaInformation_getPatterns(pe$meta$metaInformation)
    metaInformation <- lapply(metaInformation, function(mi) strsplit(mi, "_")[[1]])
    metaInformation <- do.call(c, metaInformation)
    metaInformation <- unname(metaInformation)
    pc <- c(pc, list(metaInformation = metaInformation))
    pc <- c(pc, list(dco = names(petab_joinDCO(pe))))
    
    return(pc)
  }
  
  # Default information
  m <- c("observableId","simulationConditionId","measurement","time","observableParameters","noiseParameters","datasetId","replicateId","preequilibrationConditionId", "datapointId", "lloq")
  o <- c("observableId","observableName","observableFormula","observableTransformation","noiseFormula","noiseDistribution")
  e <- c("conditionId","conditionName")
  p <- c("parameterId","parameterName","parameterScale","lowerBound","upperBound","nominalValue","estimate","initializationPriorType","initializationPriorParameters","objectivePriorType","objectivePriorParameters")
  metaInformation <- NULL 
  
  list(experimentalCondition = e,
       measurementData = m,
       observables = o,
       parameters = p,
       metaInformation = metaInformation)
}



# -------------------------------------------------------------------------#
# DCO ----
# -------------------------------------------------------------------------#


#' Create one big table containing measurementData, observables and experimentalCondition
#'
#' conditionId is joined with simulationConditionId 
#'
#' @param petab [petab]
#'
#' @return data.table
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @family dco
#' 
#' @examples 
#' # Without metaInformation
#' pe <- petab_exampleRead("01")
#' dco <- petab_joinDCO(pe)
#' bla <- petab_unjoinDCO(dco,pe)
#' 
#' # With metaInformation including a "pattern" entry
#' pe <- petab_exampleRead("04")
#' pe <- petab_mutateDCO(pe, i = conditionId == "C1", j = `:=`(conditionId = "cell1_dose1"))
#' pe <- petab_mutateDCO(pe, i = conditionId == "C2", j = `:=`(conditionId = "cell1_dose2"))
#' pe$meta$metaInformation <- list(experimentalCondition = list(conditionId = list(pattern = "celltype_dose")))
#' dco <- petab_joinDCO(pe)
#' dco
petab_joinDCO <- function(pe, FLAGincludeMetaInformation = TRUE) {
  
  if (length(pe$measurementData$preequilibrationConditionId) &&
      any(!is.na(pe$preequilibrationConditionId)))
    warning("conditionId is joined with simulationConditionId")
  
  if (is.null(pe$experimentalCondition)) stop("no experimentalCondition table to join on. You can create a dummy table with \npetab_experimentalCondition(unique(pe$measurementData$simulationConditionId))")
  if (is.null(pe$observables)) stop("no observables table to join on. You can create a dummy table with \npetab_observables(unique(pe$measurementData$observableId))")
  
  dco <- copy(pe$measurementData)
  dco <- pe$experimentalCondition[dco, on = c("conditionId" = "simulationConditionId")]
  dco <- pe$observables[dco, on = c("observableId")]
  
  if (FLAGincludeMetaInformation & !is.null(pe$meta$metaInformation)) {
    # Add units
    # units <- pe$meta$metaInformation$units
    # names(units) <- paste0("unit_", names(units))
    # dco <- data.table(dco, as.data.table(units))
    
    # Expand petab_columns
    dco <- dco_expandPatterns(dco, pe)
  }
  dco
}

#' Unjoin DataConditonObs
#'
#' @param DCO A dco: a DataConditionObs data.table, from petab_joinDCO
#' @param pe [petab()] object
#'
#' @return pe with measurementData,experimentalCondition,Observables replaced by the columns created from DCO
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @family dco
#' 
#' @examples 
#' # see ?petab_joinDCO
petab_unjoinDCO <- function(DCO, pe = NULL) {
  # Get standard column names
  pc <- petab_columns(pe = pe)
  
  # Extract the tables
  # Excols: A bit overcautious with the names, when there is the potential to
  #   pass a petab to pe it should be able to get the names
  #   from that, but let's keep it for now.
  excols <- unique(c(intersect(names(DCO), pc$experimentalCondition),
                     setdiff(names(DCO), c(pc$measurementData, pc$observables, pc$metaInformation))))
  experimentalCondition = do.call(petab_experimentalCondition, DCO[,..excols])
  experimentalCondition <- unique(experimentalCondition)
  
  mecols <- intersect(names(DCO), c("conditionId", pc$measurementData))
  measurementData <- DCO[,..mecols]
  setnames(measurementData, "conditionId", "simulationConditionId") # need to name back
  measurementData = do.call(petab_measurementData, measurementData)
  
  obcols <- intersect(names(DCO), pc$observables)
  observables = do.call(petab_observables, DCO[,..obcols])
  observables <- unique(observables)
  
  # Do some sanity checks
  lostConditions <- setdiff(unique(pe$experimentalCondition$conditionId),
                            unique(experimentalCondition$conditionId))
  if (length(lostConditions))
    warning("Conditions got lost in the transformation: ", paste0(lostConditions, collapse = ", "))
  
  # Modify and return petab
  pe$experimentalCondition <- experimentalCondition
  pe$measurementData       <- measurementData
  pe$observables           <- observables
  pe
}


#' Consistently mutate the dco
#'
#' If j is supplied, the full dco is always returned and no subsetting is done!
#'
#' @param pe [petab()]
#' @param i expression to subset the dco, can be missing
#' @param j Must be a call of the form `:=`(hsdhfgjs). If missing, pe is subsetted by i
#'
#' @return pe with modified DCO
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @family dco
#'
#' @examples
#' # 1. supply i and j, mutate in rows defined by i. Important: Need to assign!
#' pe <- petab_exampleRead("02")
#' pe <- petab_mutateDCO(pe, replicateId == 2 & observableId == "obsE", `:=`(observableId = "obsE_special", observableFormula = "log(E)/log(7)"))
#' pe
#' # 2. supply i only. Subset to rows
#' pe <- petab_exampleRead("02")
#' pe <- petab_mutateDCO(pe, replicateId == 2 & observableId == "obsE")
#' pe
#' # 3. supply j only. Mutate all rows
#' pe <- petab_exampleRead("02")
#' pe <- petab_mutateDCO(pe, j = `:=`(observableId = "obsE_special", observableFormula = "log(E)/log(7)", observableName = "obsE_special", noiseFormula = "1"))
#' pe
petab_mutateDCO <- function(pe, i, j) {
  
  mi <- missing(i)
  mj <- missing(j)
  if ( mi &&  mj) return(pe)
  
  # Catch the arguments
  si <- substitute(i)
  sj <- substitute(j)
  
  # Probably unnecessary check that j is supplied properly
  if (!mj){
    pj <- getParseData(parse(text=deparse(sj)))
    is_set <- "`:=`" %in% pj[,"text"] # is the command a "set_*" command in the data.table sense?
    if (!is_set) stop("j should be called with `:=`")
  }
  
  # Perform the data transformation on the DCO
  dco <- petab_joinDCO(pe)
  if (!mi &&  mj) {dco <- dco[eval(si)];   cat("subsetted dco\n")}       # filter
  if ( mi && !mj) {dco[,eval(sj)];         cat("mutated dco\n")}         # mutate
  if (!mi && !mj) {dco[eval(si),eval(sj)]; cat("mutated dco in rows\n")} # mutate_in_row
  
  # Return modified petab
  pe_out <- petab_unjoinDCO(dco, pe)
  pe_out
}






#' Expand columns which have a specified pattern
#'
#' @param dco dco without meta columns
#' @param pe pe. If meta$emtaInformation$petab_columns specification present, this information will be used to generate new columns
#' 
#' Example:
#' pe$meta$emtaInformation$petab_columns = list(replicateId = "replicate_gel_experiment", conditionId = "cellline_stimulus")
#' will split 
#' a) the replicateId column into three columns called "replicate", "gel" and "experiment" which are added to the dco
#' b) the conditionId column into cellline and stimulus
#' 
#' @return dco
#' @export
#' @importFrom data.table tstrsplit
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family dco
#' @family metaInformation
dco_expandPatterns <- function(dco, pe) {
  
  patterns <- petab_metaInformation_getPatterns(pe$meta$metaInformation)
  
  for (nm in names(patterns)) {
    varnames <- strsplit(patterns[[nm]], "_")[[1]]
    keep <- which(!varnames %in% names(dco)) # don't overwrite existing columns
    if (length(keep)) 
      dco[,(varnames[keep]):=(eval(parse(text = paste0('tstrsplit(',nm,', "_", keep = keep)'))))]
  }
  dco
}








# -------------------------------------------------------------------------#
# Initializers of core objects ----
# -------------------------------------------------------------------------#

#' Constructor for Conditions
#'
#' https://petab.readthedocs.io/en/stable/documentation_data_format.html#condition-table
#'
#' @param conditionId character Condition ID
#' @param conditionName character Condition ID for plotting
#' @param ... parameterOrSpeciesOrCompartmentId1. Numeric (value) or string (different parameterId)
#'
#' @return data.table
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
petab_experimentalCondition <- function(
  conditionId,
  conditionName = NA,
  ...) {
  pe_ex_Pars <- list(...)
  d <- data.table(conditionId =   as.character(conditionId),
                  conditionName = as.character(conditionName),
                  as.data.table(pe_ex_Pars))
  
  # Sorting tables with base::order() because of some specific use case I had once.
  d[base::order(conditionId)]
}


#' Constructor for Measurements
#'
#' @param observableId character
#' @param simulationConditionId character
#' @param measurement numeric
#' @param time numeric
#' @param observableParameters numeric string or NA
#' @param datasetId character
#' @param replicateId character
#' @param preequilibrationConditionId character
#' @param noiseParameters numeric, string or NA: Measurement noise or parameter name
#' @param datapointId: Deviating from original petab, a unique identifier of each data point
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @family measurementData
#' 
petab_measurementData <- function(
  observableId,
  simulationConditionId,
  measurement,
  time,
  observableParameters        = list(NA, "1;1", "scale_obsi;offset_obsi")[[1]],
  noiseParameters             = list(1, "error_ADD_obsi;error_REL_obsi")[[1]],
  datasetId                   = NA,
  replicateId                 = NA,
  preequilibrationConditionId = NA,
  datapointId                 = NA,
  lloq                        = -Inf
) {
  
  d <- data.table(
    observableId                = as.character(observableId),
    preequilibrationConditionId = as.character(preequilibrationConditionId),
    simulationConditionId       = as.character(simulationConditionId),
    measurement                 = as.double(measurement),
    time                        = as.double(time),
    observableParameters        = as.character(observableParameters),
    noiseParameters             = as.character(noiseParameters),
    datasetId                   = as.character(datasetId),
    replicateId                 = as.character(replicateId),
    datapointId                 = as.character(datapointId),
    lloq                        = as.double(lloq)
  )
  
  if (any(is.na(d$measurement))) {
    cat("Removing missing measurements: \n")
    print(d[is.na(measurement), c("observableId", "simulationConditionId", "datapointId")])
    d <- d[!is.na(measurement)]
  }
  
  # Handle BLOQ
  bloq <- d$measurement < lloq
  if (any(bloq)) {
    cat("Replacing BLOQ measuremnts by LLOQ")
    d$measurement[bloq] <- d$lloq[bloq]
  }
  
  
  
  d[base::order(simulationConditionId, observableId, time)]
}

#' Constructor for Observables
#'
#' @param observableId character
#' @param observableName character
#' @param observableFormula character
#' @param observableTransformation "lin", "log", "log10"
#' @param noiseFormula character
#' @param noiseDistribution "normal", "laplace"
#'
#' @return data.table
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
petab_observables <- function(
  observableId,
  observableName           = NA,
  observableFormula        = "observableParameter1_${observableId} * state1",
  observableTransformation = c("lin", "log", "log10")[[1]],
  noiseFormula             = c(1, "noiseParameter${n}_${observableId} + noiseParameter${n}_${observableId}*${observableId}")[[1]], # aka errormodel
  noiseDistribution        = c("normal", "laplace")[[1]]) {
  d <- data.table(
    observableId             = as.character(observableId),
    observableName           = as.character(observableName),
    observableFormula        = as.character(observableFormula),
    observableTransformation = as.character(observableTransformation),
    noiseFormula             = as.character(noiseFormula),
    noiseDistribution        = as.character(noiseDistribution)
  )
  d[base::order(observableId)]
}

#' Constructor for Parameters
#'
#' @param parameterId character
#' @param parameterName character
#' @param parameterScale "lin", "log", "log10"
#' @param lowerBound,upperBound,nominalValue numeric
#' @param estimate 0 or 1
#' @param initializationPriorType see argument suggestions
#' @param initializationPriorParameters see argument suggestions
#' @param objectivePriorType see argument suggestions
#' @param objectivePriorParameters see argument suggestions
#'
#' @return data.table
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
petab_parameters <- function(
  parameterId,
  parameterName                 = NA,
  parameterScale                = c("log10", "log", "lin")[[1]],
  lowerBound                    = 0.00001, # given on linear scale
  upperBound                    = 1000,    # given on linear scale
  nominalValue                  = 1,       # given on linear scale
  estimate                      = c(1,0)[[1]],
  initializationPriorType       = c("parameterScaleUniform","uniform","normal","laplace","logNormal","logLaplace","parameterScaleNormal","parameterScaleLaplace")[[1]],
  initializationPriorParameters = "-1;1",
  objectivePriorType            = c("parameterScaleNormal","parameterScaleUniform","uniform","normal","laplace","logNormal","logLaplace","parameterScaleLaplace")[[1]],
  objectivePriorParameters      = "0;2") {
  
  d <- data.table(
    parameterId                   = as.character(parameterId),
    parameterName                 = as.character(parameterName),
    parameterScale                = as.character(parameterScale),
    lowerBound                    = as.double(lowerBound),
    upperBound                    = as.double(upperBound),
    nominalValue                  = as.double(nominalValue),
    estimate                      = as.double(estimate),
    initializationPriorType       = as.character(initializationPriorType),
    initializationPriorParameters = as.character(initializationPriorParameters),
    objectivePriorType            = as.character(objectivePriorType),
    objectivePriorParameters      = as.character(objectivePriorParameters)
  )
  d[base::order(parameterId)]
}

#' PEtab structural model without sbml
#'
#' @param equationList eqnlist
#' @param events eventlist
#' @param ... not used, but could be used in the future for imitating assignment rules etc
#' @param parInfo [getParInfo()]
#' @param speciesInfo [getSpeciesInfo()]
#'
#' @return list
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
petab_model <- function(equationList, events = NA,
                        parInfo = getParInfo(equationList, events),
                        speciesInfo = getSpeciesInfo(equationList),
                        ...) {
  if (anyDuplicated(equationList$description)) stop("equationList description entries cannot contain duplicates")
  list(equationList = equationList, events = events,
       parInfo = parInfo, speciesInfo = speciesInfo,
       ...)
}

# -------------------------------------------------------------------------#
# Initializer of meta objects ----
# -------------------------------------------------------------------------#


#' Additional Information
#'
#' @param parameterFormulaInjection data.table(parameterId, parameterFormula)
#' @param variableOrder vector of stateIds and observableIds in the order which you want for plotting
#' @param ...
#' @param metaInformation ordered list potentially containing descriptions of columns, e.g. list(experimentalCondition = list(conditionId = list(patterns = "cellline_dose")), measurementData = list(time = list(unit = "hours"), measurement = list(unit = "value", lloq = 1))) see example of petab_joinDCO
#'
#' @return named list of supplied arguments
#' @export
petab_meta <- function(parameterFormulaInjection = NULL, variableOrder = NULL, metaInformation = NULL,...) {
  list(parameterFormulaInjection = parameterFormulaInjection, variableOrder = variableOrder, metaInformation = metaInformation, ...)
}


#' Create a parameterFormulaInjection
#'
#' See example-02
#'
#' Inject a trafo between pEst and x
#'
#' This is a somewhat elegant, somewhat hacky implementation of "complicated" parameter
#' trafos, e.g. custom steady state trafos.
#' It would be nicer to have this stuff directly in the SBML file, but I don't have the time for that now.
#' Furthermore, some functionality of the import will be reusable when this is moved to SBML as assignmentRule or sth like that
#'
#' @param parameterId e.g. "kprodS"
#' @param parameterFormula "kdegS * S"
#'
#' @return data.table(parameterId, parameterFormula)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
petab_parameterFormulaInjection = function(parameterId, parameterFormula, trafoType = "SteadyState") {
  data.table(
    parameterId      = as.character(parameterId),
    parameterFormula = as.character(parameterFormula),
    trafoType        = as.character(trafoType)
  )
}


# -------------------------------------------------------------------------#
# PEtab representation ----
# -------------------------------------------------------------------------#

#' Collector function for petab files
#'
#' @param model [dMod::eqnlist()]
#' @param condition see [petab_condition()]
#' @param measurements see [petab_measurements()]
#' @param observables see [petab_observables()]
#' @param parameters see [petab_parameters()]
#'
#' @return list of the input arguments
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
#' # [ ] Todo
#' # filename <- system.file("examples/petab/enyzmeKinetics", package = "conveniencefunctions")
#' # readPetab(filename)
petab <- function(
  model = NULL,
  experimentalCondition = NULL,
  measurementData = NULL,
  observables = NULL,
  parameters = NULL,
  meta = NULL,
  ...
) {
  
  # Do type coercion and initialize list
  if(!is.null(model))                 model                 = do.call(petab_model, model)
  if(!is.null(experimentalCondition)) experimentalCondition = do.call(petab_experimentalCondition, experimentalCondition)
  if(!is.null(measurementData))       measurementData       = do.call(petab_measurementData, measurementData)
  if(!is.null(observables))           observables           = do.call(petab_observables, observables)
  if(!is.null(parameters))            parameters            = do.call(petab_parameters, parameters)
  if(!is.null(meta))                  meta                  = do.call(petab_meta, meta)
  petab <- list(
    model                 = model,
    experimentalCondition = experimentalCondition,
    measurementData       = measurementData,
    observables           = observables,
    parameters            = parameters,
    meta                  = meta,
    ...
  )
  
  petab
}

#' Consistent file names
#'
#' @param filename path, potentially ending in .petab for backwards compatibility
#'
#' @return list(modelname, path)
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @importFrom tools file_path_sans_ext
#'
#' @examples
#' # same:
#' petab_modelname_path("Models/Example")
#' petab_modelname_path("Models/Example/Example.yaml")
petab_modelname_path <- function(filename) {
  modelname <- basename(tools::file_path_sans_ext(filename))
  path <- if (grepl("yaml",filename)) dirname(filename) else filename
  list(modelname = modelname, path = path)
}


#' List petab files
#'
#' @param filename "path/to/modelname.petab". Will generate filenames like
#'        "path/to/modelname/model_modelname.xml"
#' @param FLAGreturnList return list or vector?
#'
#' @return list or character vector of file paths
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
#' filename <- "Models/petabFolder/mypetab.yaml"
#' petab_files("Models/Example")
petab_files <- function(filename, FLAGoverwrite = FALSE) {
  
  FLAGisYaml   <- grepl(".yaml$", filename)  
  FLAGFromYaml <- FLAGisYaml  & file.exists(filename) & !FLAGoverwrite
  
  modelname <- petab_modelname_path(filename)$modelname
  path      <- petab_modelname_path(filename)$path
  
  out <- NULL
  
  if (FLAGFromYaml) {
    out <- petab_files_fromYaml(filename)
  } else {
    out <- petab_files_default(modelname, path)
  }
  out
}


#' Get default filenames
#'
#' @param modelname 
#' @param path 
#'
#' @return list of file.paths
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab_files
#' @importFrom stats setNames
#'
#' @examples
#' petab_files_default("modelname", "path")
petab_files_default <- function(modelname, path) {
  out <- list(
    yaml                       = paste0(modelname, ".yaml"),
    experimentalCondition      = paste0("experimentalCondition_"     , modelname, ".tsv"),
    measurementData            = paste0("measurementData_"           , modelname, ".tsv"),
    modelXML                   = paste0("model_"                     , modelname, ".xml"),
    # [ ] not very elegant. Remove rds when sbml is stable
    model                      = paste0("model_"                     , modelname, ".rds"),
    observables                = paste0("observables_"               , modelname, ".tsv"),
    parameters                 = paste0("parameters_"                , modelname, ".tsv"),
    simulatedData              = paste0("simulatedData_"             , modelname, ".tsv"),
    visualizationSpecification = paste0("visualizationSpecification_", modelname, ".tsv"),
    meta                       = paste0("meta_"                      , modelname, ".rds"),
    metaInformation            = paste0("metaInformation_"           , modelname, ".yaml"),
    reportYaml                 = paste0(modelname, "_report"          , ".yaml")
  )
  out <- lapply(out, function(x) file.path(path, x))
  out
}

#' Read list of filenames from a petabYaml
#'
#' @param filename file.path to petab.yaml
#'
#' @return list of filenames
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab files
#' @importFrom yaml read_yaml
#' @importFrom tools file_path_sans_ext
#'
#' @examples
#' filename <- file.path(petab_examplePath("01"), "petab.yaml")
#' petab_files_fromYaml(filename)
petab_files_fromYaml <- function(filename) {
  path <- dirname(filename)
  yaml_content <- yaml::read_yaml(filename)
  if (length(yaml_content$problems) > 1) 
    stop("More than one 'problems' section is not yet supported. If you want to implement it, consider that several petab_combine* functions exist already.")
  out <- list(
    yaml                       = basename(filename),
    experimentalCondition      = yaml_content$problems[[1]]$condition_files[[1]],
    measurementData            = yaml_content$problems[[1]]$measurement_files[[1]],
    modelXML                   = yaml_content$problems[[1]]$sbml_files[[1]],
    # [ ] not very elegant. Remove rds when sbml is stable
    model                      = yaml_content$problems[[1]]$model_files[[1]],
    observables                = yaml_content$problems[[1]]$observable_files[[1]],
    parameters                 = yaml_content$parameter_file,
    # simulatedData              = paste0("_simulatedData"             , ".tsv"),
    visualizationSpecification = yaml_content$problems[[1]]$visualization_files[[1]],
    meta                      = yaml_content$meta_file[[1]],
    metaInformation           = yaml_content$metaInformation_file[[1]],
    reportYaml                = paste0(tools::file_path_sans_ext(basename(filename)), "_report.yaml")
  )
  out <- lapply(out, function(x) file.path(path, x))
  out
}



#' Read PEtab files
#'
#' @param modelname
#' @param path
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom yaml read_yaml
#' @importFrom data.table fread
#'
#' @examples
readPetab <- function(filename) {
  
  files <- petab_files(filename = filename)
  # might break if more than one "problem" is in the yaml
  files <- unlist(files, recursive = F)
  files <- files[file.exists(files)]
  # tables
  files_tsv <- grep("tsv", files, value = TRUE)
  files_tsv <- lapply(files_tsv, data.table::fread)
  # model
  # files_model <- grep("xml", files, value = TRUE) # Do nothing, read rds
  # model+meta
  files_model <- grep("rds", files, value = TRUE)
  files_model <- lapply(files_model, readRDS)
  
  if ("metaInformation" %in% names(files)) files_model$meta$metaInformation <- yaml::read_yaml(files["metaInformation"])
  
  pe <- do.call(petab, c(files_model, files_tsv))
  petab_lint(pe)
  pe
}

#' Title
#'
#' @param petab
#' @param filename
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @importFrom data.table fwrite
#' @importFrom yaml write_yaml
#'
#' @examples
writePetab <- function(pe, filename = "petab/model") {
  
  # run linter
  petab_lint(pe)
  
  # Create folder
  dir.create(petab_modelname_path(filename)$path, FALSE, TRUE)
  
  # Get filenames of objects to write except for yaml
  files0 <- petab_files(filename = filename, FLAGoverwrite = TRUE)
  modelname <- gsub(".petab$","", basename(filename))
  
  # Hack: Duplicate model, export as RDS as well. 
  # This makes it possible to "read" a petab instead of "import"ing it. 
  # To resolve this hack, one would need to parse the sbml model into its original dMod objects. 
  # But this is currently nested too deeply into the import and not refactored nicely.
  if ("model" %in% names(pe)) pe$modelXML <- pe$model
  
  # Select files to write: Ugly, this removes metaInformation
  files <- files0[names(pe)]
  files <- files[vapply(pe, function(x) !is.null(x), TRUE)]
  
  # Write tables
  files_tsv <- grep("tsv", files, value = TRUE)
  if (length(files_tsv))
    lapply(names(files_tsv), function(nm) {
      data.table::fwrite(pe[[nm]], files[[nm]], sep = "\t")})
  
  # Write metaInformation.yaml
  file_metaInformation <- NULL
  if (!is.null(pe$meta$metaInformation)) {
    # Ugly to split up meta into meta and metaInformation.yaml
    file_metaInformation <- petab_files(filename = filename)[["metaInformation"]]
    yaml::write_yaml(pe$meta$metaInformation, file_metaInformation) 
    pe$meta$metaInformation <- NULL
  }
  
  # Write rds files: model+meta
  files_rds <- grep("rds", files, value = TRUE)
  if (length(files_rds))
    lapply(names(files_rds), function(nm) {
      saveRDS(pe[[nm]], files[[nm]])})
  
  # Write model xml
  files_xml <- grep("xml", files, value = TRUE)
  if (length(files_xml)) {
    args <- c(pe$model, list(filename = files_xml,
                             modelname = modelname))
    # args <- args[setdiff(names(args), "events")] # [ ] Todo: Events
    do.call(sbml_exportEquationList, args)
  }
  
  
  # Write petabYaml
  files <- c(files, metaInformation = file_metaInformation) # metaInformation due to hacky split up of meta
  files <- lapply(files, basename)                          # Only use basenames
  
  petabYaml(parameter_file = files$parameters,
            condition_files = files$experimentalCondition,
            measurement_files = files$measurementData,
            observable_files = files$observables,
            sbml_files = files$modelXML,
            visualization_files = files$visualizationSpecification,
            meta_file = files$meta,
            metaInformation_file = files$metaInformation,
            filename = files0$yaml)
  
  invisible(pe)
}


# -------------------------------------------------------------------------#
# petab yaml ----
# -------------------------------------------------------------------------#

#' Gather all filenames in the list structure of a petab.yaml
#'
#' @param condition_files,measurement_files,observable_files,sbml_files,visualization_files,meta_file,metaInformation_file,parameter_file file.path of length 1 or NULL
#' @param filename file.path to write yaml to
#' 
#' 
#' @return list of these objects, with sublist "problems"
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab_files
#' @importFrom yaml write_yaml
#'
#' @examples
#' petabYaml(parameter_file = "wup", 
#' condition_files = NULL, 
#' measurement_files = "bla", 
#' observable_files = NULL, 
#' sbml_files = NULL, 
#' visualization_files = NULL,
#' meta_file = NULL, 
#' metaInformation_file = "bam", filename = "bla.yaml")
petabYaml <- function(parameter_file = NULL, 
                      condition_files = NULL, measurement_files = NULL, observable_files = NULL, sbml_files = NULL, visualization_files = NULL,
                      meta_file = NULL, metaInformation_file = NULL,
                      filename = NULL
) {
  
  if (length(condition_files) > 1 || length(measurement_files) > 1 || length(observable_files) > 1 || length(sbml_files) > 1 || length(visualization_files) > 1)
    stop("Only one problem allowed")
  
  yaml_content <- list(format_version = 1, # Update if necessary
                       parameter_file = parameter_file, 
                       problems = list(list( # One problem
                         condition_files = condition_files, 
                         measurement_files = measurement_files, 
                         observable_files = observable_files, 
                         sbml_files = sbml_files, 
                         visualization_files = visualization_files)), 
                       meta_file = meta_file, 
                       metaInformation_file = metaInformation_file)
  
  if (!is.null(filename)) {
    yaml::write_yaml(yaml_content, filename) 
    return(invisible(yaml_content))
  }
  yaml_content
}


# -------------------------------------------------------------------------#
# Combining two petabs ----
# -------------------------------------------------------------------------#

#' Title
#'
#' @param ec1
#' @param ec2
#'
#' @return
#' @export
#' 
#' @family Combine petabs
#' @examples
petab_combine_experimentalCondition <- function(ec1, ec2, NFLAGconflict = c("stop" = 0, "use_pe1" = 1, "use_pe2" = 2)[2]) {
  s12 <- setdiff(names(ec1), names(ec2))
  if (length(s12)) stop("The following names are in experimentalCondition 1, but not in 2:\n",
                        paste0(s12, collapse = ", "), "\n",
                        "Please fix manually before combining.\n")
  s21 <- setdiff(names(ec1), names(ec2))
  if (length(s21)) stop("The following names are in experimentalCondition 1, but not in 2:\n",
                        paste0(s21, collapse = ", "), "\n",
                        "Please fix manually before combining.\n")
  
  i12 <- intersect(ec1$conditionId,ec2$conditionId)
  if (length(i12)) {
    cat("Overlapping conditionIds =================\n-----------------pe1----------------------\n")
    print(ec1[conditionId %in% i12])
    cat("-----------------pe2----------------------\n")
    print(ec2[conditionId %in% i12])
    if (NFLAGconflict == 0) {
      stop("The following conditionId is in both petabs: ", paste0(i12, collapse = ","))
    } else if (NFLAGconflict == 1) {
      ec2 <- ec2[!conditionId %in% i12]
    } else if (NFLAGconflict == 2) {
      ec1 <- ec1[!conditionId %in% i12]
    }
  }
  
  rbindlist(list(ec1,ec2), use.names = TRUE)
}

#' Title
#'
#' @param md1
#' @param md2
#'
#' @return
#' @export
#' @family Combine petabs
#'
#' @examples
petab_combine_measurementData <- function(md1, md2) {
  rbindlist(list(md1,md2), use.names = TRUE)
}

#' Title
#'
#' @param o1
#' @param o2
#'
#' @return
#' @export
#' @family Combine petabs
#' @importFrom data.table rbindlist
#'
#' @examples
petab_combine_observables <- function(o1,o2, NFLAGconflict = c("stop" = 0, "use_pe1" = 1, "use_pe2" = 2)[2]) {
  i12 <- intersect(o1$observableId,o2$observableId)
  if (length(i12)) {
    cat("Overlapping observables =================\n-----------------pe1----------------------\n")
    print(o1[observableId %in% i12])
    cat("-----------------pe2----------------------\n")
    print(o2[observableId %in% i12])
    if (NFLAGconflict == 0) {
      stop("The following observableId is in both petabs: ", paste0(i12, collapse = ","))
    } else if (NFLAGconflict == 1) {
      o2 <- o2[!observableId %in% i12]
    } else if (NFLAGconflict == 2) {
      o1 <- o1[!observableId %in% i12]
    }
  }
  data.table::rbindlist(list(o1,o2), use.names = TRUE)
}



#' Title
#'
#' @param p1
#' @param p2
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Combine petabs
#'
#' @importFrom data.table rbindlist
#' @examples
#' 
petab_combine_parameters <- function(p1,p2) {
  
  if (is.null(p1)) return(p2)
  if (is.null(p2)) return(p1)
  
  i12 <- intersect(p2$parameterId,p1$parameterId)
  if (length(i12)) {
    message("The following parameterId are in both petabs. Using parameters from pe1. \n",
            paste0(i12, collapse = ","))
  }
  p2 <- p2[!parameterId %in% p1$parameterId]
  data.table::rbindlist(list(p1,p2), use.names = TRUE)
}

#' Combine two petabs
#' 
#' * Row-bind the relevant tables
#' * If conflicts occur, inform user about it (e.g. different column names in experimentalCondition)
#'
#' @param pe1,pe2 petabs
#' @param NFLAGconflict one of c("stop" = 0, "use_pe1" = 1, "use_pe2" = 2)
#'
#' @return pe
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Combine petabs
#'
#' @examples
petab_combine <- function(pe1,pe2, NFLAGconflict = c("stop" = 0, "use_pe1" = 1, "use_pe2" = 2)[2]) {
  if (NFLAGconflict>0) cat("Preferring pe", NFLAGconflict, "\n")
  
  petab(
    model                 = pe1$model,
    experimentalCondition = petab_combine_experimentalCondition(pe1$experimentalCondition, pe2$experimentalCondition),
    measurementData       = petab_combine_measurementData(      pe1$measurementData      , pe2$measurementData),
    observables           = petab_combine_observables(          pe1$observables          , pe2$observables, NFLAGconflict = NFLAGconflict),
    parameters            = petab_combine_parameters(           pe1$parameters           , pe2$parameters),
    meta                  = list(pe1,pe2)[[NFLAGconflict]]$meta
  )
  
}




# -------------------------------------------------------------------------#
# Interface to useful PEtab.py functions ----
# -------------------------------------------------------------------------#

#' Own little linter until petab.py lint is stable
#'
#' @param pe petab
#'
#' @return list of errors
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @importFrom stringr str_count
#'
#' @examples
#' petab_lint(petab_exampleRead("01", "pe"))
#' # Todo: Implement an example where the petab 
petab_lint <- function(pe) {
  
  # Basic linting, collect error messages here
  
  # [ ] Implement access to petab.lint
  
  # Collect all errors in a list
  errlist <- list()
  
  # experimentalCondition
  dupes <- which(duplicated(pe$experimentalCondition$conditionId))
  if(length(dupes)) {
    warning("These rows are duplicates in conditionId :", paste0(head(dupes,10), collapse = ","), "...")
    errlist <- c(errlist, list(conditionIdDupes = dupes))}
  
  nm_fixed <- setdiff(names(pe$experimentalCondition), c("conditionId", "conditionName"))
  if (length(nm_fixed)) {
    nm_fixed <- vapply(setNames(nm = nm_fixed), function(nm) {suppressWarnings(any(!is.na(as.numeric(pe$experimentalCondition[[nm]]))))}, FUN.VALUE = TRUE)
    if (any(nm_fixed)) warning("The numeric values in experimentalCondition are assumed on estScale. The following parameters are affected: ", paste0(names(nm_fixed)[nm_fixed], collapse = ", "))
  }
  
  # measurementData
  dupes <- which(duplicated(pe$measurementData))
  if(length(dupes)) {
    warning("These rows are duplicates in measurementData: ", paste0(dupes, collapse = ","))
    errlist <- c(errlist, list(measurementDataDupes = dupes))}
  
  if (any(is.na(pe$measurementData$time))){
    cat("Missing times here: \n")
    print(pe$measurementData[is.na(time), unique(.SD), .SDcols = c("observableId", "simulationConditionId")])
    stop("Petab contains missing times")
  }
  
  # check for existing observable parameters
  if (!is.null(pe$observables) & !is.null(pe$measurementData)) {
    dco <- pe$observables[pe$measurementData, on = .NATURAL]
    dco <- unique(dco[,list(observableId, observableFormula, observableParameters)])
    dco[,`:=`(nobspars = stringr::str_count(observableFormula,   "observableParameter"))]
    dco[,`:=`(nobsparsAvalailable = if(is.na(observableParameters)|observableParameters == "") 0 else 1+stringr::str_count(observableParameters, ";")), by = 1:nrow(dco)]
    hasMissingObsPar <- dco[,which(nobsparsAvalailable < nobspars)]
    if(length(hasMissingObsPar)) stop("The following observableIds have no observableParameters specified: ", paste0(dco[hasMissingObsPar, observableId], collapse = ", "))
    # check for existing noise parameters
    dco <- pe$observables[pe$measurementData, on = .NATURAL]
    unique(dco[,list(observableId, noiseFormula, noiseParameters)])
    dco <- unique(dco[,list(observableId, noiseFormula, noiseParameters)])
    dco[,`:=`(nobspars = stringr::str_count(noiseFormula,   "noiseParameter"))]
    dco[,`:=`(nobsparsAvalailable = if(is.na(noiseParameters)|noiseParameters == "") 0 else 1+stringr::str_count(noiseParameters, ";")), by = 1:nrow(dco)]
    hasMissingNoisePar <- dco[,which(nobsparsAvalailable < nobspars)]
    if(length(hasMissingNoisePar)) stop("The following observableIds have no noiseParameters specified: ", paste0(dco[hasMissingNoisePar, observableId], collapse = ", "))
  }
  
  
  # observables
  dupes <- which(duplicated(pe$observables$observableID))
  if(length(dupes)) {
    warning("These rows are duplicates in observableId: ", paste0(head(dupes,10), collapse = ","), "...")
    errlist <- c(errlist, list(observableIdDupes = dupes))}
  
  # parameters
  if (!is.null(pe$parameters)) {
    pars_NA <- pe$parameters[is.na(nominalValue), parameterId]
    if (length(pars_NA)) warning("These parameters have no correct nominal value: ", paste0(pars_NA, collapse = ","))
    
    parsNotEstimatedNotLin <- pe$parameters[estimate == 0 & parameterScale != "lin"]
    if (nrow(parsNotEstimatedNotLin)) {
      logscale_but_zero <- parsNotEstimatedNotLin$nominalValue == 0
      if (any(logscale_but_zero)) stop("Fixed parameters on log-scale, but their nominal value is zero: ",
                                       paste0(parsNotEstimatedNotLin$parameterId[logscale_but_zero], collapse = ","))
    }
    parsNotLinNominal0 <- pe$parameters[nominalValue == 0 & parameterScale != "lin"]
    if (nrow(parsNotLinNominal0)) stop("Parameters on log-scale, but their nominal value is zero: ",
                                       paste0(parsNotLinNominal0$parameterId, collapse = ","))
  }
  
  # meta
  if (!is.null(pe$meta$parameterFormulaInjection)) {
    names_overwritten <- pe$meta$parameterFormulaInjection$parameterId %in% names(pe$experimentalCondition)
    if (any(names_overwritten)) stop("Parameters given in est.grid are overwritten by parameterFormulaInjection: ", 
                                     paste0(pe$meta$parameterFormulaInjection$parameterId[names_overwritten], ", "))}
  
  errlist
}


#' Sample parameter start points
#'
#' Uses Petab.py, is based on pe$parameters
#'
#' @param pe
#' @param n_starts number of starting points
#' @param seed
#' @param FLAGincludeCurrent Include values of pe$parameters in returned start parameters
#'
#' @return parframe
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family pepy
#' 
#' @importFrom dMod parframe
#'
#'
#' @examples
#' n_starts = 100L
#' seed = 1L
#' FLAGincludeCurrent = TRUE
pepy_sample_parameter_startpoints <- function(pe, n_starts = 100L, seed = 1L, FLAGincludeCurrent = TRUE) {
  n_starts <- as.integer(n_starts-as.numeric(FLAGincludeCurrent))
  seed     <- as.integer(seed)
  
  pepy <- petab_python_setup()
  
  pars <- pepy$sample_parameter_startpoints(
    parameter_df = pe$parameters,
    n_starts = n_starts,
    seed = seed)
  pars <- `colnames<-`(pars, pe$parameters$parameterId[pe$parameters$estimate==1])
  if (FLAGincludeCurrent) {
    pars <- pars[,names(pd$pars)]
    pars <- rbind(pd$pars, pars)
  }
  pars <- dMod::parframe(pars)
  pars
}

# -------------------------------------------------------------------------#
# Python setup ----
# -------------------------------------------------------------------------#

#' Reinstall the python setup
#'
#' @return nothing, installs the virtualenv again
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family pepy
#' @importFrom reticulate virtualenv_create use_virtualenv virtualenv_install install_python use_python
#' 
petab_python_installPackages <- function(FLAGcleanInstall = FALSE, FLAGforcePip = FALSE, FLAGforcePythonVersion = FALSE) {
  message("If initialization fails, try running with FLAGcleanInstall or FLAGforcePythonVersion. Might help")
  
  if (FLAGcleanInstall){
    # Hacky version for Linux only
    message("Please restart RStudio. If that doesn't help call this function again")
    if (readline("enter 'yes' to continue") != "yes") return("Nothing was done")
    unlink("~/.virtualenvs/petab", T)
    unlink("~/.local/share/r-reticulate/", T)
  }
  
  if (FLAGforcePythonVersion){
    pyversion <- "3.9.7"
    pyver <- reticulate::install_python(pyversion)
    reticulate::use_python(pyver, TRUE)
  }
  
  reticulate::use_virtualenv("petab")
  reticulate::virtualenv_install("petab", "petab", ignore_installed = TRUE)
  reticulate::virtualenv_install("petab", "petab-select", ignore_installed = TRUE)
  "installed petab in virtual environment"
}

#' Setup the connection to python petab
#'
#' * if petab virtualenv not present:
#'   * sets up a virtualenv called "petab"
#'   * pip installs petab
#' * uses "petab" virtual env
#' * imports and returns petab
#'
#' use as pe <- petab_python_setup()
#'  
#' If this function fails, try one of these: 
#'   1. Restart RStudio from your terminal
#'   2. petab::petab_python_installPackages(TRUE)
#'   3. Recreate the *.Rproj file of your current project
#'  
#' @param FLAGreturnpetabSelect 
#'
#' @return python module, see [reticulate::import()]
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @family pepy
#' 
#' @importFrom reticulate virtualenv_list virtualenv_install use_virtualenv import
#'
#' @examples
#' pepy <- petab_python_setup()
#' # pepy$lint(pe)
#' 
#' peps <- petab_python_setup(FLAGreturnpetabSelect = TRUE)
petab_python_setup <- function(FLAGreturnpetabSelect = FALSE) {
  
  # Necessary due to reticulate/rstudio interaction 
  mywd <- getwd()
  on.exit(setwd(mywd))
  setwd("~")
  
  if (!"petab" %in% reticulate::virtualenv_list()){
    petab_python_installPackages(FLAGcleanInstall = FALSE)
  }
  
  message("Using petab virtualenv\n")
  reticulate::use_virtualenv("petab")
  if (FLAGreturnpetabSelect) return(reticulate::import("petab_select"))
  reticulate::import("petab")
}



# -------------------------------------------------------------------------#
# PEtab import for dMod ----
# -------------------------------------------------------------------------#
#' Parse the parameters columns in measurementData
#'
#'
#' @param measurementData petab$measurementData. Uses only columns observableId,simulationConditionId and "column"
#' @param column character(1L). one of c("observableParameters", "noiseParameters")
#'
#' @return data.table(condition, Innerpars... = Outerpars...) of local parameters.
#'   Can be added to a gridlist with [dMod::indiv_addLocalParsToGridList()]
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
#' measurementData <- data.table(
#' observableId = c("obsE"),
#' simulationConditionId = c("C1", "C2"),
#' observableParameters = c("scale;offset", "1;offset"),
#' noiseParameters = "sigmaAbs"
#' )
#' experimentalCondition <- data.table(conditionId = c("C1", "C2", "C3"))
#' pe <- list(measurementData = measurementData, experimentalCondition = experimentalCondition)
#' petab_getMeasurementParsMapping(pe, "noiseParameters")
#' petab_getMeasurementParsMapping(pe, "observableParameters")
#' # suggested use: indiv_addLocalParsToGridlist ...
petab_getMeasurementParsMapping <- function(pe, column = c("observableParameters", "noiseParameters")[1]) {
  
  measurementData <- copy(pe$measurementData)
  
  # Select column and name
  parameters <- measurementData[[column]]
  parameterString <- gsub("Parameters", "Parameter", column)
  
  if (all(is.na(parameters))) return(data.table(conditionId = unique(measurementData$simulationConditionId)))
  
  # Pipeline of death
  mp <- strsplit(parameters, ";")
  mp <- lapply(mp, function(x) {if(length(x)) return(as.data.table(as.list(x))) else data.table(NA)})
  mp <- rbindlist(mp, fill = TRUE)
  setnames(mp, paste0(parameterString, 1:length(mp), "_"))
  mp <- data.table(observableId = measurementData$observableId, conditionId = measurementData$simulationConditionId, mp)
  mp <- melt(mp, id.vars = c("observableId", "conditionId"), variable.name = "INNERPARAMETER", variable.factor = FALSE, value.name = "OUTERPARAMETER")
  mp <- unique(mp)
  mp <- mp[!is.na(OUTERPARAMETER)]
  mp[,`:=`(INNERPARAMETER = paste0(INNERPARAMETER, observableId))]
  mp <- mp[,list(conditionId, INNERPARAMETER, OUTERPARAMETER)]
  dupes <- duplicated(mp[,list(conditionId,INNERPARAMETER)])
  if (any(dupes)) {
    print(mp[dupes, list(conditionId,INNERPARAMETER)])
    stop(column, "contain non-unique values in some conditionIds (see print output above).\n
         Error models with multiple erros per observable are *not* supported as of today (2022-02-01)")}
  mp <- dcast(mp, conditionId ~ INNERPARAMETER, value.var = "OUTERPARAMETER")
  
  # Check that all conditionIds are specified, and if not, create empty row
  experimentalCondition <- copy(pe$experimentalCondition)
  if (length(setdiff(experimentalCondition$conditionId, mp$conditionId))) {
    mp <- mp[data.table(conditionId = experimentalCondition$conditionId), on = "conditionId"]
  }
  
  # Replace NAs with dummy value 1
  colsWithNa <- vapply(mp, function(x) any(is.na(x)), FALSE)
  colsWithNa <- names(colsWithNa)[colsWithNa]
  if (length(colsWithNa)) {
    # warning("The following parameters are not specified in all conditionIds. ",
    #         "In the unspecified conditionIds, they are set to 1 (should not make a difference):\n",
    #         paste0(colsWithNa, collapse = ","))
    mp[,(colsWithNa):=lapply(.SD, function(x) replace(x, is.na(x),1)), .SDcols = colsWithNa]
  }
  
  mp
}



#' Create the prior function specified in parameters_df
#'
#' @param pe petab
#' @param FLAGuseNominalCenter Should be deprecated
#'
#' @return [dMod::constraintL2()]
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family asdf
#' @importFrom dMod constraintL2
#' 
#' @examples
petab_createObjPrior <- function(pe, FLAGuseNominalValueAsCenter = FALSE) {
  # [ ] Todo: Move this function to dMod-*.R file
  
  p <- copy(pe$parameters)
  p <- p[estimate == 1&!is.na(objectivePriorType)]
  
  notImplemented <- which(p$objectivePriorType != "parameterScaleNormal")
  if (length(notImplemented))
    warning("objectivePriorType not implemented: ", paste0(unique(p$objectivePriorType[notImplemented]),collapse = ","), "\n",
            "The respective parameters are not going to be priored: ", paste0(unique(p$parameterId[notImplemented]),collapse = ","))
  p <- p[objectivePriorType == "parameterScaleNormal"]
  
  
  p[,`:=`(mu = as.numeric(gsub(";.*", "", objectivePriorParameters)))]
  if (FLAGuseNominalValueAsCenter) {
    # Apply transformation before using as center
    cat("using pe$parameters$nominalValue as center for prior\n")
    p[,`:=`(estValue =
              eval(parse(
                text = paste0(parameterScale, "(", nominalValue, ")")))),
      by = 1:nrow(p)]
    p[,`:=`(mu = estValue)]
  }
  p[,`:=`(sd = as.numeric(gsub(".*;", "", objectivePriorParameters)))]
  
  
  dMod::constraintL2(mu = setNames(p$mu, p$parameterId),
                     sigma = setNames(p$sd, p$parameterId))
  
}


#' Title
#'
#' @param pe petab
#' @param whichBoundary upper or lower
#'
#' @return c(name = value)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
petab_getParameterBoundaries <- function(pe, whichBoundary = c("upper", "lower")[1]) {
  p <- copy(pe$parameters)
  p <- p[estimate == 1]
  lin <- function(x) x
  p[,`:=`(boundary =
            eval(parse(
              text = paste0(parameterScale, "(", whichBoundary, "Bound", ")")))),
    by = 1:nrow(p)]
  setNames(p$boundary, p$parameterId)
}





petab_getMeasurementParsScales <- function(measurementData,parameters) {
  
  # >>>>>>> obsolete function!
  # >>>>>>> functionality covered currently by
  # >>>>>>> dMod::updateParscalesToBaseTrafo <<<<<<<<
  
  parameters <- pe$parameters
  measurementData <- pe$measurementData
  
  # measurementPars
  mp <- lapply(c("observableParameters", "noiseParameters"), function(column) {
    petab_getMeasurementParsMapping(measurementData, column)})
  mp <- merge(mp[[1]],mp[[2]])
  mp <- melt(mp, id.vars = "condition",
             variable.name = "INNERPARAMETER",
             variable.factor = FALSE, value.name = "parameterId")
  
  # parScales
  ps <- parameters[mp, on = c("parameterId") ]
  ps <- ps[,list(INNERPARAMETER, parameterScale)]
  ps <- unique(ps)
  
  # check that trafo is the same in all conditions
  dupeScale <- duplicated(ps$INNERPARAMETER)
  if (any(dupeScale))
    stop("Model parameter has two different scales mapped to it in different conditions: ",
         ps$INNERPARAMETER[dupeScale])
  
  list(parScales = setNames(ps$parameterScale, ps$INNERPARAMETER),
       coveredPars = unique(mp$parameterId))
}


# -------------------------------------------------------------------------#
# Plotting ----
# -------------------------------------------------------------------------#

#' Plot Petab data
#'
#' Automatic multipage export in the background with the future package and ggforce::facet_wrap_paginate
#'
#' @param petab A [petab] object
#' @param aeslist list of formulas used to call [ggplot2::aes_q()].
#'  Default aesthetics are list(x = ~time, y = ~measurement, color = ~conditionId).
#'  Can be overriden by aeslist
#' @param ggCallback additional stuff to add to the ggplot, such as a call to [ggplot2::labs()] or scales
#' @param FLAGmeanLine draw line connecting the means of groups defined by c("observableId", "time", "conditionId", names(aeslist)=
#' @param FLAGfuture export asynchronously with the future package
#' @param ... Arguments to [conveniencefunctions::cf_outputFigure()]
#'
#' @return ggplot
#'
#' @family plotting
#'
#' @importFrom ggforce n_pages facet_wrap_paginate
#' @importFrom conveniencefunctions cf_outputFigure
#' @importFrom cOde getSymbols
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' 
#' @examples 
#' 
#' pe <- petab_exampleRead(exampleName = "03", "pe")
#' 
#' # Plot data
#' petab_plotData(pe)
#' 
#' # Plot dose response: subset and remap aes(x)
#' petab_plotData(pe, i = conditionName == "dose response 30 min", aeslist = list(x = ~Epo))
#' 
#' # Change layout, scale and theme via ggCallback
#' petab_plotData(pe, i = conditionName == "dose response 30 min", aeslist = list(x = ~Epo),
#'                ggCallback = list(facet_wrap(~observableId, scales = "free"), 
#'                                  scale_x_log10(breaks =unique(pe$experimentalCondition[conditionName=="dose response 30 min", Epo])),
#'                                      theme_bw()))
#' 
#' # Add datapointId, plot data points with labels
#' pe <- petab_exampleRead(exampleName = "01", "pe")
#' pe <- petab_mutateDCO(pe,j=`:=`(datapointId = paste0(datasetId, "_", 1:.N)))
#' petab_plotData(pe, ggCallback = list(facet_wrap(~observableId, scales = "free"),
#'                                      ggrepel::geom_text_repel(aes(x= time, y = measurement, label = datapointId), color = "grey", size = 2)))
#'
#' # Save plot to multipage pdf
#' td <- tempdir()
#' petab_plotData(pe,
#'                ggCallback = list(facet_wrap_paginate(~observableId, nrow = 1, ncol = 2, scales = "free"),
#'                                  ggrepel::geom_text_repel(aes(x= time, y = measurement, label = datapointId), color = "grey", size = 1)),
#'                filename = file.path(td, "plot.pdf"),
#'                width = 15.5, heightrel = 8/16, scale = 1, units = "cm")
#' system(paste0("nautilus ", td), wait = FALSE)
#'
#' # Display second page of multipage plot interactively
#' petab_plotData(pe,
#'                ggCallback = list(facet_wrap_paginate(~observableId, nrow = 1, ncol = 2, scales = "free", page = 2),
#'                                  ggrepel::geom_text_repel(aes(x= time, y = measurement, label = datapointId), color = "grey", size = 1)))
#'
petab_plotData <- function(petab,
                           i,j,
                           aeslist = NULL,
                           FLAGUseObservableTransformation = TRUE,
                           FLAGmeanLine = TRUE,
                           ggCallback = list(),
                           ...
) {
  
  # create plotting data.table
  dplot <- petab_joinDCO(petab)
  # apply log transformation to data when applicable
  if (FLAGUseObservableTransformation) {
    dplot[,`:=`(measurement = eval(parse(text = paste0(observableTransformation, "(measurement)")))#,
                # observableId = paste0(observableTransformation, " ", observableId) # HACK for Jamboree, where I didn't want the scale in the plot: Make this accessible as option
    ),
    by = 1:nrow(dplot)]
  }
  
  # Handle i and j
  mi <- missing(i)
  mj <- missing(j)
  si <- substitute(i)
  sj <- substitute(j)
  if (!mi) {dplot <- dplot[eval(si)]}
  if (!mj) {dplot[,eval(sj)]}
  
  # mean values to draw lines
  if (FLAGmeanLine){
    byvars <- lapply(aeslist, function(x) cOde::getSymbols(as.character(x)))
    byvars <- do.call(c, byvars)
    byvars <- unique(c("observableId", "time", "conditionId", byvars))
    dmean <- copy(dplot)
    dmean <- dmean[,list(measurement = mean(measurement)), by = byvars]
  }
  
  # handle aesthetics
  aes0 <- list(x = ~time, y = ~measurement, color = ~conditionId)
  aeslist <- c(aeslist, aes0[setdiff(names(aes0), names(aeslist))])
  
  # Create plot
  pl <- cfggplot(dplot)
  pl <- pl + ggforce::facet_wrap_paginate(~observableId, nrow = 4, ncol = 4, scales = "free", page = 1)
  if (FLAGmeanLine) { # Add first so te lines don't mask the points
    aesmean0 <- list(linetype = ~conditionId, group = as.formula(paste0("~ interaction(", paste0(setdiff(byvars, c("time")), collapse = ", "), ")")))
    aesmeanlist <- c(aeslist, aesmean0[setdiff(names(aesmean0), names(aeslist))])
    pl <- pl + geom_line(do.call(aes_q, aesmeanlist), data = dmean)
  }
  pl <- pl + geom_point(do.call(aes_q, aeslist), data = dplot)
  for (plx in ggCallback) pl <- pl + plx
  
  # Print paginate message so user doesnt forget about additional pages
  npages <- ggforce::n_pages(pl)
  if (length(npages) && npages > 1) message("Plot has ", npages, " pages. View them via passing a ggCallback with facet_wrap_paginate(..., page=N) (see examples) \n")
  
  # output
  cf_outputFigure(pl = pl, ...)
}



# #' Mark data points in a plot
# #'
# #' * Opens shiny app 
# #' * Writes clicked points to a file
# #' 
# #' @param plot A plot from petab_plotData. Must have aeslist(customdata~datapointId) mapped
# #' @param fileCSV filename for csv file
# #'
# #' @return Called for side-effect
# #' @export
# #' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
# #' @md
# #' @family plotting
# #'
# #' @examples
# #' # Plot data
# #' pe <- petab_exampleRead(exampleName = "01", "pe")
# #' pe$measurementData[,`:=`(datapointId = sprintf(paste0("%s_%0",nchar(as.character(.N)), "i"),  datasetId, 1:.N))]
# #' plot <- petab_plotData(pe, aeslist = list(customdata=~datapointId))
# #' 
# #' # Mark outliers
# #' fileCSV <- "~/wup.txt"
# #' petab_markDataPointsShiny_step1(plot, fileCSV)
# #' outliers <- petab_markDataPointsShiny_step2(fileCSV)
# #' 
# #' # Flag them in the original plot
# #' pe$measurementData[,`:=`(outlier = FALSE)]
# #' pe$measurementData[datapointId %in% outliers$datapointId, `:=`(outlier = TRUE)]
# #' petab_plotData(pe, aeslist = list(customdata=~datapointId, shape = ~outlier, size = ~outlier))
# petab_markDataPointsShiny_step1 <- function(plot, fileCSV) {
#   
#   library(shiny)
#   library(plotly)
#   
#   if (file.exists(fileCSV)) stop("file exists, please remove file")
#   
#   ui <- fluidPage(
#     plotlyOutput("plot", height = "960px"),
#     verbatimTextOutput("click")
#   )
#   
#   server <- function(input, output, session) {
#     
#     
#     output$plot <- renderPlotly({
#       p <- ggplotly(plot)
#       p %>% 
#         layout(dragmode = "select") %>%
#         event_register("plotly_selecting")
#     })
#     
#     output$click <- renderPrint({
#       d <- event_data("plotly_click")
#       write.table(d, file = fileCSV, append=TRUE, sep=",", row.names = FALSE)
#       if (is.null(d)) "Click events appear here (double-click to clear)" else d
#     })
#   }
#   shinyApp(ui, server)
# }
# 
# 
# #' Post process marked data points
# #' @rdname petab_markDataPointsShiny_step1
# #' @export
# #' @importFrom data.table fread fwrite
# petab_markDataPointsShiny_step2 <- function(fileCSV, filename = NULL, NFLAGtribble = 0) {
#   wup <- readLines(fileCSV)
#   wup <- wup[grepl("\\w",wup)]
#   wup <- wup[grep(",{4}", gsub("[^,]","",wup))] # only keep entries with 5 fields (i.e. entries which have customdata)
#   wup <- wup[!(duplicated(wup) & grepl("pointNumber|curveNumber|customdata", wup))]
#   selectedPoints <- data.table::fread(text = wup)
#   selectedPoints <- selectedPoints[,list(datapointId = customdata)]
#   selectedPoints <- unique(selectedPoints)
#   if (!is.null(filename)) data.table::fwrite(selectedPoints, file = filename)
#   if (NFLAGtribble) conveniencefunctions::cfoutput_MdTable(selectedPoints, NFLAGtribble = NFLAGtribble)
#   selectedPoints
# }


# -------------------------------------------------------------------------#
# Overview tables ----
# -------------------------------------------------------------------------#

#' Generate overview table which observables are in which condition
#'
#' @param pe [petab()] object
#' @param Ntruncate truncate pasted observables at this many characters
#' @param FLAGincludedatasetId summarize per conditionId, datasetId and replicateId
#' @param ... arguments going to [conveniencefunctions::cfoutput_MdTable()]
#'
#' @return prints table to console or writes it to disk
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Overview Tables
#' @importFrom conveniencefunctions cfoutput_MdTable
petab_overviewObsPerCond <- function(pe, Ntruncate = 1000, FLAGincludedatasetId = TRUE, ...) {
  dco <- petab_joinDCO(pe)
  if ("conditionName" %in% names(dco)) dco[,`:=`(conditionName = conditionId)]
  bycols <- if (FLAGincludedatasetId) c("conditionId", "datasetId", "replicateId") else c("conditionId")
  dco <- dco[,list(observableId = paste0(sort(unique(observableId)), collapse = ",")), by = bycols]
  dco <- dco[,`:=`(observableId = substr(observableId, 1, Ntruncate))]
  conveniencefunctions::cfoutput_MdTable(dco, ...)
}


#' Print available names in the DCO
#'
#' @param pe 
#'
#' @return prints list of available columns
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Overview Tables
#'
#' @examples
petab_overviewDCONames <- function(pe) {
  columns <- petab_columns(pe)
  columns$parameters <- NULL
  
  metaNames <- names(petab_joinDCO(pe))
  metaNames <- setdiff(metaNames, do.call(c, columns))
  columns$metaInformation <- metaNames
  
  print(columns)
}


#' Print some overview stuff
#'
#' @param pe pe
#' @param FLAGcolumns,FLAGmetaInformation,FLAGobsPerCond include this and that info?
#' @param ... arguments going to [petab_overviewObsPerCond()]
#'
#' @return Prints some stuff to console
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Overview Tables
#' @importFrom yaml as.yaml
#' @importFrom conveniencefunctions cfoutput_MdTable
#'
#' @examples
#' petab_overview(petab_exampleRead("01", "pe"))
petab_overview <- function(pe, FLAGcolumns = TRUE, FLAGmetaInformation = TRUE, FLAGobsPerCond = TRUE, ...) {
  if (FLAGcolumns){
    cat("Available columns\n==================\n")
    pc <- petab_columns(pe)
    pc <- lapply(pc, function(x) paste0(x, collapse = ", "))
    conveniencefunctions::cfoutput_MdTable(data.table(Table = names(pc), Columns = unlist(pc, F,F)), FLAGsummaryRow = FALSE)
  }
  if (FLAGmetaInformation) {
    cat("\n")
    cat("MetaInformation\n=================\n")
    cat(yaml::as.yaml(pe$meta$metaInformation))
  }
  if (FLAGobsPerCond) {
    cat("\nObservableId per condition\n======================\n")
    petab_overviewObsPerCond(pe, ...)
  }
  # To be continued
}



# -------------------------------------------------------------------------#
# Very small helpers ----
# -------------------------------------------------------------------------#

# Might delete some of them again

#' Get a mapping of parameter types
#'
#' @param pe [petab()] (with parameters defined)
#'
#' @return data.table(parameterId, parameterType)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom cOde getSymbols
#'
#' @examples
petab_getParameterType <- function(pe) {
  
  # Parameter types
  observableParameters <- cOde::getSymbols(pe$measurementData$observableParameters)
  noiseParameters      <- cOde::getSymbols(pe$measurementData$noiseParameters)
  
  parameters <- data.table(parameterId = pe$parameters$parameterId)
  parameters[,`:=`(parameterType = "other")]
  parameters[parameterId %in% observableParameters,`:=`(parameterType = "observableParameters")]
  parameters[parameterId %in% noiseParameters     ,`:=`(parameterType = "noiseParameters")]
  parameters[grep("^L1_",parameterId)             ,`:=`(parameterType = "L1")]
  
  parameters
}



#' Title
#'
#' @param pe
#'
#' @return character vector
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
petab_getParametersToEstimate <- function(pe) {
  pe$parameters[estimate==1,parameterId]
}



#' Get parameter names from pe$experimentalCondition
#'
#' @param experimentalCondition
#'
#' @return character vector
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' 
#' @importFrom data.table copy
#'
#' @examples
petab_getParametersExperimentalCondition <- function(experimentalCondition) {
  ec <- data.table::copy(experimentalCondition)
  ec[,`:=`(conditionId = NULL, conditionName = NULL)]
  ec <- lapply(ec, getSymbols)
  ec <- do.call(c,ec)
  ec <- unique(ec)
  ec
}


#' Hash a petab
#' 
#' The difficulty is that data.tables have attributes specific to R sessions
#' 
#' @param pe petab
#'
#' @return a digest
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab helpers
#' @importFrom digest digest
#'
#' @examples
#' # Should evaluate to TRUE, else the petab has changed
#' petab_hash(petab_exampleRead("01","pe")) == "45281114c1d5f92ca00f9f2171379ba2"
petab_hash <- function(pe) {
  undatatable <- function(x) {
    if(is.data.table(x)) data.frame(x) else x}
  pe$model <- lapply(pe$model, undatatable)
  pe$meta <- lapply(pe$meta, undatatable)
  pe <- lapply(pe, undatatable)
  digest::digest(pe)
}




# -------------------------------------------------------------------------#
# Pars ----
# -------------------------------------------------------------------------#

#' Get parameters on outer scale
#'
#' @param pe
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Parameter wrangling
#'
#' @importFrom data.table copy
#' @examples
#' pe <- petab_exampleRead("01")
#' parsEst <- petab_getPars_estScale(pe)
petab_getPars_estScale <- function(pe) {
  p <- data.table::copy(pe$parameters)
  p[,`:=`(estValue = eval(parse(text = paste0(parameterScale, "(", nominalValue, ")")))),by = 1:nrow(p)]
  setNames(p$estValue, p$parameterId)
}

#' Get parameters on nominal scale
#'
#' @param pe
#'
#' @return parameter vector
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Parameter wrangling
#'
#' @examples
#' pe <- petab_exampleRead("01")
#' parsLin <- petab_getPars_linScale(pe)
petab_getPars_linScale <- function(pe) {
  p <- copy(pe$parameters)
  setNames(p$nominalValue, p$parameterId)
}



#' Title
#'
#' @param pe
#' @param parsLin
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Parameter wrangling
#'
#' @examples
#' pe <- petab_exampleRead("01")
#' parsEst <- petab_getPars_estScale(pe)
#' parsLin <- petab_getPars_linScale(pe)
#' parsLin <- parsLin[sample(1:length(parsLin), 5)] + rnorm(5)
petab_setPars_linScale <- function(pe, parsLin) {
  # Careful: This is technically a set* function in data.table parlance
  parsLin <- parsLin[order(names(parsLin))]
  pe$parameters[parameterId %in% names(parsLin),`:=`(nominalValue = parsLin)]
  invisible(pe)
}

#' Title
#'
#' @param pe
#' @param parsEst
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Parameter wrangling
#'
#' @examples
#' pe <- petab_exampleRead("01")
#' parsEst <- petab_getPars_estScale(pe)
#' parsLin <- petab_getPars_linScale(pe)
#' parsLin <- parsLin[sample(1:length(parsLin), 5)] + rnorm(5)
petab_setPars_estScale <- function(pe, parsEst) {
  petab_setPars_linScale(pe, petab_transformPars_est2Lin(pe, parsEst))
}

#' Title
#'
#' @param pe
#' @param parsLin
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Parameter wrangling
#'
#' @examples
#' pe <- petab_exampleRead("01")
#' parsEst <- petab_getPars_estScale(pe)
#' parsLin <- petab_getPars_linScale(pe)
#' parsLin <- parsLin[sample(1:length(parsLin), 5)] + rnorm(5)
petab_transformPars_lin2Est <- function(pe, parsLin) {
  parsLinMerge <- parsLin[order(names(parsLin))]
  p <- copy(pe$parameters)
  p <- p[parameterId %in% names(parsLinMerge)]
  p[,`:=`(nominalValue = parsLinMerge)]
  p[,`:=`(estValue = eval(parse(text = paste0(parameterScale, "(", nominalValue, ")")))),by = 1:nrow(p)]
  parsEst <- setNames(p$estValue, p$parameterId)
  parsEst <- parsEst[names(parsLin)]
  parsEst
}

#' Title
#'
#' @param pe
#' @param parsEst
#'
#' @return vector
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Parameter wrangling
#'
#' @examples
#' pe <- petab_exampleRead("01")
#' parsEst <- petab_getPars_estScale(pe)
#' parsLin <- petab_getPars_linScale(pe)
#' parsLin <- parsLin[sample(1:length(parsLin), 5)] + rnorm(5)
petab_transformPars_est2Lin <- function(pe, parsEst) {
  parsEstMerge <- parsEst[order(names(parsEst))]
  p <- copy(pe$parameters)
  p <- p[parameterId %in% names(parsEstMerge)]
  p[,`:=`(estValue = parsEstMerge)]
  p[,`:=`(nominalValue = eval(parse(text = paste0(inverse_scale(parameterScale), "(", estValue, ")")))),by = 1:nrow(p)]
  parsLin <- setNames(p$nominalValue, p$parameterId)
  parsLin <- parsLin[names(parsEst)]
  parsLin
}

#' Merge parameter values from one petab into another petab
#'
#' @param pe1 parameters from pe1
#' @param pe2 parameters from pe2
#' @param mergeCols columns to merge into pe_pa1
#'
#' @return updated pe1
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Parameter wrangling
#'
#' @examples
petab_mergeParameters <- function(pe1,pe2, mergeCols = c("nominalValue", "parameterScale", "estimate")) {
  # Get values
  pe_pa1 <- pe1$parameters
  pe_pa2 <- pe2$parameters
  # Merge
  pe_pa <- petab_parameters_mergeParameters(pe_pa1, pe_pa2, mergeCols)
  # Update pe1
  pe1$parameters <- pe_pa
  pe1
}

#' Title
#'
#' @param pe_pa1 parameters_df to merge values into
#' @param pe_pa2 parameters_df to merge values from
#' @param mergeCols Columns you want to update in pe_pa1
#'
#' @return updated pe_pa1
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @family Parameter wrangling
#' @importFrom data.table merge.data.table
#'
#' @examples
petab_parameters_mergeParameters <- function(pe_pa1, pe_pa2, mergeCols = setdiff(intersect(names(pe_pa1), names(pe_pa2)), "parameterId")) {
  message("Merging columns: ", paste0(mergeCols, collapse = ","), "\n")
  
  mergeColsPA2 <- paste0(mergeCols, "PA2")
  pe_pa2 <- pe_pa2[,c("parameterId", ..mergeCols)]
  # Merge, update, clean
  pe_pa <- data.table::merge.data.table(pe_pa1, pe_pa2, by = "parameterId", all.x = TRUE, all.y = FALSE, suffixes = c("", "PA2"))
  pe_pa[!is.na(nominalValuePA2),(mergeCols) := lapply(.SD, function(x) x), .SDcols = mergeColsPA2]
  pe_pa[,(mergeColsPA2) := NULL]
  pe_pa
}


# -------------------------------------------------------------------------#
# MeasurementData ----
# -------------------------------------------------------------------------#

#' Duplicate control data points
#'
#' Mainly for plotting, e.g. such that mean-lines start at zero
#'
#' @param pe
#' @param i logical expression operating on the DC (not DCO) indicating which data points are to be duplicated
#' @param conditionMapping data.table(conditionId, conditionIdNew) where conditionId is the control condition and conditionIdNew is the condition where to the data point shall be copied
#' @param FLAGDuplicatesFirst put duplicated data on top or bottom of measurementData? Relevant for plotting, which points are getting drawn on top of each other
#'
#' @return pe with updated
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @family plotData
#' @family measurementData
#'
#' @importFrom data.table rbindlist setnames
#'
#' @examples
petab_duplicateControls <- function(pe, i, conditionMapping, FLAGDuplicatesFirst = TRUE) {
  si <- substitute(i)
  # join_DC (couldn't I just have used the DCO?)
  DC <- pe$experimentalCondition[pe$measurementData, on = c("conditionId" = "simulationConditionId")]
  # copy controls
  dc_ctrl <- DC[eval(si)]
  # rename
  dc_ctrl <- lapply(1:nrow(conditionMapping), function(idx) {
    dc_ctrlx <- dc_ctrl[conditionMapping[idx], on =c("conditionId")]
    dc_ctrlx[,`:=`(conditionId = NULL)]
    data.table::setnames(dc_ctrlx, "conditionIdNew", "simulationConditionId")
    dc_ctrlx
  })
  dc_ctrl <- data.table::rbindlist(dc_ctrl)
  pe_measurementData_ctrl <- do.call(petab_measurementData,dc_ctrl[,.SD, .SDcols = names(pe$measurementData)])
  if (FLAGDuplicatesFirst){
    # HACK: Make nicer. I needed this in one plot because controls were light-colored crosses and stimulated
    # were big dark points and the crosses within the points looked ugly (because of the order, duplicates were plotted first, then the true data got overlaid)
    #
    pe$measurementData <- data.table::rbindlist(list(pe_measurementData_ctrl,pe$measurementData))
  } else {
    pe$measurementData <- data.table::rbindlist(list(pe$measurementData,pe_measurementData_ctrl))
  }
  pe
}


#' Title
#'
#' @param pe 
#'
#' @return pe
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family measurementData
#'
#' @examples
petab_applyObservableTransformation <- function(pe) {
  dco <- petab_joinDCO(pe)
  dco[,`:=`(measurement = eval(parse(text = paste0(observableTransformation, "(measurement)"))))]
  petab_unjoinDCO(dco, pe)
}

#' Title
#'
#' @param pe 
#'
#' @return pe
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family measurementData
#' 
#' @examples
petab_applyInverseObservableTransformation <- function(pe) {
  dco <- petab_joinDCO(pe)
  dco[observableTransformation == "log10",`:=`(measurement = eval(parse(text = paste0("10^(measurement)"))))]
  dco[observableTransformation == "log",`:=`(measurement = eval(parse(text = paste0("exp(measurement)"))))]
  petab_unjoinDCO(dco, pe)
}

#' Replace error parameters in measurementData by sensible estimates
#' 
#' 1. Calculates sd per time, observableId and conditionId. 
#' 2. Then the mean of these sd's is taken per observableId
#' 
#' All calculations are done on the est-scale of each observable
#' 
#' Assumes constant error model (per observable)
#' 
#' @param pe 
#'
#' @return pe
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family measurementData
#'
#' @examples
petab_fixErrorModel <- function(pe) {
  pe <- petab_applyObservableTransformation(pe)
  pe_me <- pe$measurementData
  idx <- !grepl(";",pe_me$noiseParameters)&suppressWarnings(is.na(as.numeric(pe_me$noiseParameters)))
  pe_me[idx,`:=`(noiseParameters = sd(measurement)), by = c("time", "observableId", "simulationConditionId")]
  pe_me[idx,`:=`(noiseParameters = as.character(mean(as.numeric(noiseParameters), na.rm = TRUE))), by = c("observableId")]
  pe$measurementData <- pe_me
  petab_applyInverseObservableTransformation(pe)
}

# -------------------------------------------------------------------------#
# metaInformation ----
# -------------------------------------------------------------------------#

#' Get "patterns" from metaInformation
#'
#' @param pe pe with pe$meta$metaInformation
#'
#' @return list(column_name = "pattern")
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family metaInformation
#'
#' @examples
#' # In two tables
#' mi <- yaml::read_yaml(text = "measurementData:\n  measurement:\n    unit: absolute_conc\n    lloq: ~\n  time:\n    unit: hours\n  replicateId:\n    pattern: experiment_replicate_technicalReplicate\nexperimentalCondition:\n  conditionId:\n    pattern: cellline_TGFb_Statin\n  TGFb:\n    unit: ng/ul\n  Statin:\n    unit: ug/ml\n")
#' petab_metaInformation_getPatterns(mi)
#' 
#' # In one table
#' mi <- yaml::read_yaml(text = "measurementData:\n  measurement:\n    unit: absolute_conc\n    lloq: ~\n  time:\n    unit: hours\n  experimentalCondition:\n  conditionId:\n    pattern: cellline_TGFb_Statin\n  TGFb:\n    unit: ng/ul\n  Statin:\n    unit: ug/ml\n")
#' petab_metaInformation_getPatterns(mi)
#' 
#' # No pattern
#' mi <- yaml::read_yaml(text = "measurementData:\n  measurement:\n    unit: absolute_conc")
#' petab_metaInformation_getPatterns(mi)
petab_metaInformation_getPatterns <- function(mi) {
  mi <- mi[intersect(names(mi), c("experimentalCondition", "measurementData", "observables"))]
  patterns <- lapply(unname(mi), function(x) {
    patterns_x <- lapply(x, function(y) y$pattern)
    patterns_x <- do.call(c, patterns_x)
  }) 
  patterns <- do.call(c, patterns)
  patterns <- as.list(patterns)
  patterns
}



# -------------------------------------------------------------------------#
# Scales ----
# -------------------------------------------------------------------------#

#' Title
#'
#' @param parameterScale vector with entries "lin", "log", "log10"
#'
#' @return vector with entries "lin", "exp", "10^"
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family scales
#' @examples
inverse_scale <- function(parameterScale) {
  parameterScale[parameterScale == "log"] <- "exp"
  parameterScale[parameterScale == "log10"] <- "10^"
  parameterScale
}


#' Identity
#'
#' @param x
#'
#' @return x
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family scales
#' @examples
lin <- function(x) {
  x
}



# -------------------------------------------------------------------------#
# L1 ----
# -------------------------------------------------------------------------#

#' Create L1 Trafo, depending on parameterScale
#' 
#' Trafo injection
#' 
#' * For logpars: par * L1_par
#' * For linpars: par + L1_par
#' 
#' @param pe 
#' @param parameterId_base vector of parameterId which should be L1'd.
#'
#' @return data.table(parameterId, parameterFormula, trafoType = "L1")
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab L1
#'
#' @examples
#' pe <- petab_exampleRead("05")
#' L1_getParameterFormulaInjection(pe, c("kcat", "E"))
L1_getParameterFormulaInjection <- function(pe, parameterId_base) {
  p <- pe$parameters[parameterId %in% parameterId_base]
  p[,`:=`(parameterIdL1 = paste0("L1_", parameterId))]
  p[,`:=`(parameterIdBase = paste0("L1Ref_", parameterId))]
  p[,`:=`(parameterFormula = paste0(parameterIdBase, ifelse(grepl("lin", parameterScale), " + ", " * "), parameterIdL1))]
  p[,list(parameterId, parameterFormula, trafoType = "L1")]
}


#' Title
#'
#' @param l1 pe$meta$L1
#'
#' @return
#' @export
#'
#' @examples
L1_augmentL1Info <- function(l1) {
  # [ ] Should be standard in pe$meta$L1
  conditionSpecL1_L1Conds <- setdiff(unique(l1$L1Spec$L1Spec), l1$conditionSpecL1_reference)
  parameters_L1 <- lapply(conditionSpecL1_L1Conds, function(cn) paste0("L1_",l1$parameterId_base, "_", cn))
  parameters_L1 <- do.call(c, parameters_L1)
  parameters_L1reference <- paste0("L1Ref_",l1$parameterId_base)
  
  list(conditionSpecL1_L1Conds = conditionSpecL1_L1Conds,parameters_L1 = parameters_L1,parameters_L1reference = parameters_L1reference)
}



#' Title
#'
#' @inheritParams L1_getParameterFormulaInjection
#' @return petab with updated parameterFormulaInjection, see [L1_getParameterFormulaInjection()]
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab L1
#' @importFrom data.table rbindlist
#'
#' @examples
#' pe <- petab_exampleRead("05")
#' pe_L1_updateParameterFormulaInjection(pe, c("kcat", "E"))
pe_L1_updateParameterFormulaInjection <- function(pe, parameterId_base) {
  pfi <- L1_getParameterFormulaInjection(pe, parameterId_base)
  pe$meta$parameterFormulaInjection <- data.table::rbindlist(list(pe$meta$parameterFormulaInjection, pfi), use.names = TRUE, fill = TRUE)
  pe
}


#' Add parameter columns to experimentalCondition
#' 
#'                 
#' @param pe petab with non-Null pe$experimentalCondition$L1Spec column which serves a) subsetting and b) creating parameter names
#' @inheritParams pe_L1_createL1Problem
#' 
#' @return petab with updated experimentalCondition: L1Spec is removed, L1 parameter columns are added
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab L1
#'
#' @examples
#' pe <- petab_exampleRead("04")
#' pe$experimentalCondition$L1Spec <- pe$experimentalCondition$conditionId
#' pe_L1_addL1ParsToExperimentalCondition(pe, parameterId_base = c("kcat", "E"), conditionSpecL1_reference = "C1")
pe_L1_addL1ParsToExperimentalCondition <- function(pe, parameterId_base, conditionSpecL1_reference) {
  # Checks
  if (is.null(pe$meta$L1$L1Spec)) stop("Please add pe$meta$L1Spec")  
  
  # Rename parameters in experimentalCondition
  pe_ex_colnames <- parameterId_base[parameterId_base %in% names(pe$experimentalCondition)]
  pe_ex_parameters <- parameterId_base[parameterId_base %in% petab_getParametersExperimentalCondition(pe$experimentalCondition)]
  if (length(union(pe_ex_colnames,pe_ex_parameters))) warning("These L1'd parameters are already condition specific: ", paste0(union(pe_ex_colnames,pe_ex_parameters), collapse = ","))
  if (length(pe_ex_colnames)) setnames(pe$experimentalCondition, pe_ex_colnames, paste0("L1Ref_", pe_ex_colnames))
  if (length(pe_ex_parameters)) pe$experimentalCondition[,(names(pe$experimentalCondition)):=lapply(.SD, function(x) replaceSymbols(pe_ex_parameters,
                                                                                                                                    paste0("L1Ref_", pe_ex_parameters),
                                                                                                                                    x))]
  # Add parameters to experimentalCondition
  pe$experimentalCondition <- pe$meta$L1$L1Spec[pe$experimentalCondition, on = "conditionId"]
  pe$experimentalCondition[,(paste0("L1_", parameterId_base)) := lapply(parameterId_base, function(px) paste0("L1_", px, "_", L1Spec))]
  pe$experimentalCondition[L1Spec %in% conditionSpecL1_reference,(paste0("L1_", parameterId_base)) := 0]
  pe$experimentalCondition[,`:=`(L1Spec = NULL)]
  pe
}


#' From a petab, create a petab which encodes an L1 problem
#'
#' @param pe petab
#' @param parameterId_base parameterIds to be L1'd
#' @param conditionSpecL1_reference Vector of entries referring to L1Spec which are deemed the base condition
#' @param j_conditionSpecL1 data.table-j argument to create L1Spec within pe$experimentalCondition. E.g. conditionId (as symbol, will be substitute-eval'd)
#'
#' @return [petab()] with L1 parameters
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab L1
#' @importFrom data.table copy
#'
#' @examples
#' pe <- petab_exampleRead("04")
#' pe_L1_createL1Problem(pe, c("kcat", "E"), conditionSpecL1_reference = "C1", j_conditionSpecL1 = conditionId)
pe_L1_createL1Problem <- function(pe, parameterId_base, conditionSpecL1_reference, j_conditionSpecL1 = conditionId) {
  pe <- copy(pe)
  
  # 1 Create parameterFormulaInjection
  pe <- pe_L1_updateParameterFormulaInjection(pe, parameterId_base)
  
  # 2 Add L1 spec and add other information to pe$meta$L1
  pe$meta$L1$L1Spec <- data.table::copy(pe$experimentalCondition[,list(conditionId = conditionId)])
  sj <- substitute(j_conditionSpecL1)
  pe$meta$L1$L1Spec[,`:=`(L1Spec = eval(sj))]
  
  pe$meta$L1$parameterId_base <- parameterId_base
  pe$meta$L1$conditionSpecL1_reference <- conditionSpecL1_reference
  pe$meta$L1 <- c(pe$meta$L1, L1_augmentL1Info(pe$meta$L1))
  
  # 3 Add L1 parameters to experimentalCondition
  pe <- pe_L1_addL1ParsToExperimentalCondition(pe, parameterId_base, conditionSpecL1_reference)
  
  # 4 Re-create parameters_df including L1 parameters
  pepaL1 <- petab_create_parameter_df(pe)
  pepaOld <- pe$parameters
  pepaOld[parameterId %in% parameterId_base,`:=`(parameterId = paste0("L1Ref_", parameterId))]
  pe$parameters <- petab_parameters_mergeParameters(pepaL1, pepaOld)
  pe
}

# -------------------------------------------------------------------------#
# Todolist/Wishlist ----
# -------------------------------------------------------------------------#

# [ ] sbml - getSpeciesInfo: speciesUnit
# [ ] sbml - initialAssignments
# [ ] sbml - assignmentRules
# [ ] sbml - events

