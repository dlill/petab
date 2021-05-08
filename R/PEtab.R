# -------------------------------------------------------------------------#
# Convenience functions ----
# -------------------------------------------------------------------------#

#' Original petab function didn't work
#'
#' @param model
#' @param measurementData
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
petab_create_parameter_df <- function(pe, observableParameterScale = "log10") {

  model                 <- pe$model
  measurementData       <- pe$measurementData
  experimentalCondition <- pe$experimentalCondition

  # Species
  message("Do we need to prepend init_ before species?")
  speciesInfo <- model$speciesInfo
  par_sp <- petab_parameters(parameterId =   speciesInfo$speciesName,
                             parameterName = speciesInfo$speciesName,
                             nominalValue =  speciesInfo$initialAmount,
                             estimate = as.numeric(speciesInfo$initialAmount > 0))
  # Kinetic parameters in ODEs
  parInfo <- model$parInfo
  par_pa <- petab_parameters(parameterId =   parInfo$parName,
                             parameterName = parInfo$parName,
                             nominalValue =  parInfo$parValue)
  # Observable parameters
  par_ob <- NULL
  message("Scale for all observable parameters: ", observableParameterScale, "\n")
  if (length(getSymbols(measurementData$observableParameters)))
    par_ob <- petab_parameters(parameterId =  getSymbols(measurementData$observableParameters),
                               parameterName = getSymbols(measurementData$observableParameters),
                               parameterScale = observableParameterScale)

  # MeasurementErrors
  par_meErr <- NULL
  if (length(getSymbols(measurementData$noiseParameters)))
    par_meErr <- petab_parameters(parameterId =   getSymbols(measurementData$noiseParameters),
                                  parameterName = getSymbols(measurementData$noiseParameters),
                                  nominalValue = 0.1)

  # Get all base-parameters
  par <- rbindlist(list(par_sp, par_pa, par_ob, par_meErr))


  # Parameters from experimentalConditions
  # More complicated, need also to exclude colnames
  #   from previously collected parameters
  par_ec <- NULL
  parnamesOuter <- petab_getParametersExperimentalCondition(experimentalCondition)
  if (length(parnamesOuter)){
    # Collect ec_pars
    par_ec <- petab_parameters(parameterId = parnamesOuter,
                               parameterName = parnamesOuter)
    message("Please check parameter scale of parameter names occuring in experimentalCondition or implement the scale matching in this function: ",
            paste(parnamesOuter, collapse = ", "))
    # [ ] Note: If you implement the scale-matching, refactor this function. Put all par_ec handling into another function

    # Remove base parameters from par
    parnamesInner <- setdiff(colnames(experimentalCondition), c("conditionId", "conditionName"))
    par <- par[!parameterId %in% parnamesInner]

    # Append par_ec
    par <- rbindlist(list(par, par_ec))
  }


  # parameterFormulaInjection
  pfi <- pe$meta$parameterFormulaInjection
  if (!is.null(pfi))
    par[parameterId %in% pfi$parameterId,`:=`(estimate=0, parameterScale = "lin")] # Ensure all parameters set by injection are set to estimate

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


#' Title
#'
#' @param pe NULL: Default names, [petab()] Actual names in petab
#'
#' @return list of column names
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
petab_columns <- function(pe = NULL) {
  if(!is.null(pe))
    return(lapply(pe[c("experimentalCondition",
                       "measurementData",
                       "observables",
                       "parameters")], names))

  m <- c("observableId","simulationConditionId","measurement","time","observableParameters","noiseParameters","datasetId","replicateId","preequilibrationConditionId")
  o <- c("observableId","observableName","observableFormula","observableTransformation","noiseFormula","noiseDistribution")
  e <- c("conditionId","conditionName")
  p <- c("parameterId","parameterName","parameterScale","lowerBound","upperBound","nominalValue","estimate","initializationPriorType","initializationPriorParameters","objectivePriorType","objectivePriorParameters")

  list(experimentalCondition = e,
       measurementData = m,
       observables = o,
       parameters = p)
}

#' Create one big table containing measurementData, observables and experimentalCondition
#'
#' @param petab [petab]
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
petab_joinDCO <- function(pe) {
  if (length(pe$measurementData$preequilibrationConditionId) &&
      any(!is.na(pe$preequilibrationConditionId)))
    warning("DCO might not be able to handle preequilibrationConditionId")

  dx <- copy(pe$measurementData)
  dx <- pe$experimentalCondition[dx, on = c("conditionId" = "simulationConditionId")]
  dx <- pe$observables[dx, on = c("observableId")]
  dx
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
petab_unjoinDCO <- function(DCO, pe = NULL) {
  # Get standard column names
  pc <- petab_columns(pe = pe)

  # Extract the tables
  # Excols: A bit overcautious with the names, when there is the potential to
  #   pass a petab to pe it should be able to get the names
  #   from that, but let's keep it for now.
  excols <- unique(c(intersect(names(DCO), pc$experimentalCondition),
                     setdiff(names(DCO), c(pc$measurementData, pc$observables))))
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
#'
#' @examples
#' # todo!!
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
#' @return
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
petab_experimentalCondition <- function(
  conditionId,
  conditionName = NA,
  ...) {
  d <- data.table(conditionId =   as.character(conditionId),
                  conditionName = as.character(conditionName),
                  as.data.table(list(...)))
  d[base::order(conditionId)]
}


#' Constructor for Measurements
#'
#' @param observableId
#' @param simulationConditionId
#' @param measurement
#' @param time
#' @param observableParameters numeric string or NA
#' @param datasetId
#' @param replicateId
#' @param preequilibrationConditionId
#' @param noiseParameters numeric, string or NA: Measurement noise or parameter name
#'
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
petab_measurementData <- function(
  observableId,
  simulationConditionId,
  measurement,
  time,
  observableParameters        = list(NA, "1;1", "scale_obsi;offset_obsi")[[1]],
  noiseParameters             = list(1, "error_ADD_obsi;error_REL_obsi")[[1]],
  datasetId                   = NA,
  replicateId                 = NA,
  preequilibrationConditionId = NA
) {
  d <- data.table(
    observableId                = as.character(observableId),
    preequilibrationConditionId = as.character(preequilibrationConditionId),
    simulationConditionId       = as.character(simulationConditionId),
    measurement                 = measurement,
    time                        = time,
    observableParameters        = as.character(observableParameters),
    noiseParameters             = as.character(noiseParameters),
    datasetId                   = as.character(datasetId),
    replicateId                 = as.character(replicateId)
  )
  d[base::order(simulationConditionId, observableId, time)]
}

#' Constructor for Observables
#'
#' @param observableId
#' @param observableName
#' @param observableFormula
#' @param observableTransformation
#' @param noiseFormula
#' @param noiseDistribution
#'
#' @return
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
#' @param parameterId
#' @param parameterName
#' @param parameterScale
#' @param lowerBound
#' @param upperBound
#' @param nominalValue
#' @param estimate
#' @param initializationPriorType
#' @param initializationPriorParameters
#' @param objectivePriorType
#' @param objectivePriorParameters
#'
#' @return
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
    lowerBound                    = lowerBound,
    upperBound                    = upperBound,
    nominalValue                  = nominalValue,
    estimate                      = estimate,
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
#' @param parInfo
#' @param speciesInfo
#' @param parameterFormulaInjection
#'
#' @return list
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
petab_model <- function(equationList, events = NA,
                        parInfo = getParInfo(equationList),
                        speciesInfo = getSpeciesInfo(equationList),
                        ...) {
  list(equationList = equationList, events = events,
       parInfo = parInfo, speciesInfo = speciesInfo,
       ...)
}

# -------------------------------------------------------------------------#
# Initializer of meta objects ----
# -------------------------------------------------------------------------#


#' Title
#'
#' @param parameterFormulaInjection
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
petab_meta <- function(parameterFormulaInjection = NULL, variableOrder = NULL, ...) {
       list(parameterFormulaInjection = parameterFormulaInjection, variableOrder = variableOrder, ...)
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
petab_parameterFormulaInjection = function(parameterId, parameterFormula) {
  data.table(
    parameterId      = as.character(parameterId),
    parameterFormula = as.character(parameterFormula)
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
  cat("Sorting tables with base::order()\n")
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

  petab_lint(petab)

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
#' @importFrom tools file_ext
#'
#' @examples
#' # same:
#' petab_modelname_path("Models/Example")
#' petab_modelname_path("Models/Example.petab")
petab_modelname_path <- function(filename) {
  if (tools::file_ext(filename) == "petab")
    filename <- gsub(".petab$", "", filename)
  modelname <- basename(filename)
  path <- filename
  list(modelname = modelname, path = path)
}


#' List petab files
#'
#' @param FLAGTestCase generate TestCases filename
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
#' petab_files("Models/Example")
petab_files <- function(filename, FLAGTestCase = FALSE, FLAGreturnList = FALSE) {

  modelname <- petab_modelname_path(filename)$modelname
  path      <- petab_modelname_path(filename)$path

  # [ ] warning("model refers to rds instead of xml\n")
  out <- NULL
  if (FLAGTestCase) {
    out <- c(
      yaml                       = paste0(modelname, ".yaml"),
      experimentalCondition      = paste0("_experimentalCondition"     , ".tsv"),
      measurementData            = paste0("_measurementData"           , ".tsv"),
      modelXML                   = paste0("_model"                     , ".xml"),
      # [ ] not very elegant. Remove rds when sbml is stable
      model                      = paste0("_model"                     , ".rds"),
      observables                = paste0("_observables"               , ".tsv"),
      parameters                 = paste0("_parameters"                , ".tsv"),
      simulatedData              = paste0("_simulatedData"             , ".tsv"),
      visualizationSpecification = paste0("_visualizationSpecification", ".tsv"),
      meta                       = paste0("_meta"                      , ".rds"))
  } else {
    out <- c(
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
      meta                       = paste0("meta_"                      , modelname, ".rds"))
  }
  nm <- names(out)
  out <- setNames(file.path(path, out), nm)
  if (FLAGreturnList) out <- as.list(out)
  out
}


#' Read PEtab files
#'
#' @param modelname
#' @param path
#' @param FLAGTestCase
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @importFrom data.table fread
#'
#' @examples
readPetab <- function(filename, FLAGTestCase = FALSE) {

  files <- petab_files(filename = filename, FLAGTestCase = FLAGTestCase)
  files <- files[file.exists(files)]
  # tables
  files_tsv <- grep("tsv", files, value = TRUE)
  files_tsv <- lapply(files_tsv, data.table::fread)
  # model
  # files_model <- grep("xml", files, value = TRUE) # Do nothing, read rds
  # model+meta
  files_model <- grep("rds", files, value = TRUE)
  files_model <- lapply(files_model, readRDS)

  do.call(petab, c(files_model, files_tsv))
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
#'
#' @examples
writePetab <- function(petab, filename = "petab/model") {

  # Create folder, load petab
  dir.create(petab_modelname_path(filename)$path, FALSE, TRUE)
  pe <- petab_python_setup()

  # Get filenames
  files <- petab_files(filename = filename)
  modelname <- gsub(".petab$","", basename(filename))

  # Write yaml
  pe$create_problem_yaml(sbml_files        = basename(files["modelXML"]),
                         condition_files   = basename(files["experimentalCondition"]),
                         measurement_files = basename(files["measurementData"]),
                         parameter_file    = basename(files["parameters"]),
                         observable_files  = basename(files["observables"]),
                         yaml_file         = files["yaml"])

  # [ ] Hack: Remove once sbml export is stable
  if ("model" %in% names(petab)) petab$modelXML <- petab$model

  # Select files to write
  files <- files[names(petab)]
  files <- files[vapply(petab, function(x) !is.null(x), TRUE)]

  # Write tables
  files_tsv <- grep("tsv", files, value = TRUE)
  if (length(files_tsv))
    lapply(names(files_tsv), function(nm) {
      data.table::fwrite(petab[[nm]], files[[nm]], sep = "\t")})

  # Write model+meta rds
  files_model <- grep("rds", files, value = TRUE)
  if (length(files_model))
    lapply(names(files_model), function(nm) {
      saveRDS(petab[[nm]], files[[nm]])})

  # Write model xml
  files_model <- grep("xml", files, value = TRUE)
  if (length(files_model)) {
    args <- c(petab$model, list(filename = files_model,
                                modelname = modelname))
    args <- args[setdiff(names(args), "events")] # [ ] Todo: Events
    do.call(sbml_exportEquationList, args)
  }

  cat("Success?")
  invisible(petab)
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
  rbindlist(list(o1,o2), use.names = TRUE)
}



#' Title
#'
#' @param p1
#' @param p2
#'
#' @return
#' @export
#'
#' @examples
petab_combine_parameters <- function(p1,p2) {

  if (is.null(p1)) return(p2)
  if (is.null(p2)) return(p1)

  i12 <- intersect(p2$parameterId,p1$parameterId)
  if (length(i12)) {
    message("The following parameterId are in both petabs. Using parameters from pe1. \n",
            paste0(i12, collapse = ","))
  }
  p2 <- p2[!parameterId %in% p1$parameterId]
  rbindlist(list(p1,p2), use.names = TRUE)
}

#' Title
#'
#' @param pe1
#' @param pe2
#' @param NFLAGconflict one of c("stop" = 0, "use_pe1" = 1, "use_pe2" = 2)
#'
#' @return
#' @export
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
    meta                  = pe1$meta
  )

}




# -------------------------------------------------------------------------#
# Interface to useful PEtab.py functions ----
# -------------------------------------------------------------------------#

#' Own little linter until petab.py lint is stable
#'
#' @param petab
#'
#' @return list of errors
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
petab_lint <- function(pe) {

  # [ ] Implement access to petab.lint
  errlist <- list()

  # Some quick own checks
  dupes <- which(duplicated(pe$measurementData))
  if(length(dupes)) {
    warning("These rows are duplicates in measurementData: ", paste0(head(dupes,10), collapse = ","), "...")
    errlist <- c(errlist, list(measurementDataDupes = dupes))}

  dupes <- which(duplicated(pe$observables$observableID))
  if(length(dupes)) {
    warning("These rows are duplicates in observableId: ", paste0(head(dupes,10), collapse = ","), "...")
    errlist <- c(errlist, list(observableIdDupes = dupes))}

  dupes <- which(duplicated(pe$experimentalCondition$conditionId))
  if(length(dupes)) {
    warning("These rows are duplicates in conditionId :", paste0(head(dupes,10), collapse = ","), "...")
    errlist <- c(errlist, list(conditionIdDupes = dupes))}

  if (!is.null(pe$parameters)) {
    pars_NA <- pe$parameters[is.na(nominalValue), parameterId]
    if (length(pars_NA)) warning("These parameters have no correct nominal value: ", paste0(pars_NA, collapse = ","))}

  if (any(is.na(pe$measurementData$time)))
    stop("Petab contains missing times")
  if (any(is.na(pe$measurementData$measurement)))
    stop("Petab contains missing times")


  if (!is.null(pe$meta$parameterFormulaInjection)) {
    names_overwritten <- pe$meta$parameterFormulaInjection$parameterId %in% names(pe$experimentalCondition)
    if (any(names_overwritten)) stop("Parameters given in est.grid are overwritten by parameterFormulaInjection: ", paste0(names_overwritten, ", "))}

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
#' @importFrom dMod parframe
#'
#'
#' @examples
#' n_starts = 100L
#' seed = 1L
#' FLAGincludeCurrent = TRUE
pepy_sample_parameter_startpoints <- function(pe, n_starts = 100L, seed = 1L, FLAGincludeCurrent = TRUE) {
  n_starts <- as.integer(n_starts)
  seed     <- as.integer(seed)

  pepy <- petab_python_setup()

  pars <- pepy$sample_parameter_startpoints(
    parameter_df = pe$parameters,
    n_starts = n_starts,
    seed = seed)
  if (FLAGincludeCurrent) pars <- rbind(pd$pars, pars)
  pars <- `colnames<-`(pars, pe$parameters$parameterId[pe$parameters$estimate==1])
  pars <- dMod::parframe(pars)
  pars
}

# -------------------------------------------------------------------------#
# Python setup ----
# -------------------------------------------------------------------------#

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
#' @return python module, see [reticulate::import()]
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @importFrom reticulate virtualenv_list virtualenv_install use_virtualenv import
#'
#' @examples
#' pepy <- petab_python_setup()
#' # pepy$lint(pe)
petab_python_setup <- function() {
  if (!"petab" %in% reticulate::virtualenv_list()){
    reticulate::virtualenv_install("petab", "petab", ignore_installed = TRUE)
  }

  # Stupid RStudio "ich mach mein eigenes environment variables ding"
  # PATH <- Sys.getenv("PATH")
  # PATH <- paste0(file.path(Sys.getenv("HOME"), ".virtualenvs/petab/bin"), ":", PATH)
  # Sys.setenv(PATH = PATH)
  message("======================\nIf this function fails, restart RStudio from your terminal\n==============")

  message("Using petab virtualenv\n")
  reticulate::use_virtualenv("petab")
  reticulate::import("petab")
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
petab_python_reinstall <- function() {
  message("Please restart RStudio. If that doesn't help call this function again")
  if (readline("enter 'yes' to continue") != "yes") return("Nothing was done")
  # Hacky version for Linux only
  unlink("~/.virtualenvs/petab", T)
  unlink("~/.local/share/r-reticulate/", T)
  reticulate::use_virtualenv("petab")
  reticulate::virtualenv_install("petab", "petab", ignore_installed = TRUE)
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

  if (all(is.na(parameters))) return(data.table(condition = unique(measurementData$simulationConditionId)))

  # Pipeline of death
  mp <- strsplit(parameters, ";")
  mp <- lapply(mp, function(x) {if(length(x)) return(as.data.table(as.list(x))) else data.table(NA)})
  mp <- rbindlist(mp, fill = TRUE)
  setnames(mp, paste0(parameterString, 1:length(mp), "_"))
  mp <- data.table(observableId = measurementData$observableId, condition = measurementData$simulationConditionId, mp)
  mp <- melt(mp, id.vars = c("observableId", "condition"), variable.name = "INNERPARAMETER", variable.factor = FALSE, value.name = "OUTERPARAMETER")
  mp <- unique(mp)
  mp <- mp[!is.na(OUTERPARAMETER)]
  mp[,`:=`(INNERPARAMETER = paste0(INNERPARAMETER, observableId))]
  mp <- mp[,list(condition, INNERPARAMETER, OUTERPARAMETER)]
  mp <- dcast(mp, condition ~ INNERPARAMETER, value.var = "OUTERPARAMETER")

  # Check that all conditions are specified, and if not, create empty row
  experimentalCondition <- copy(pe$experimentalCondition)
  if (length(setdiff(experimentalCondition$conditionId, mp$condition))) {
    mp <- mp[data.table(condition = experimentalCondition$conditionId), on = "condition"]
  }

  # Replace NAs with dummy value 1
  colsWithNa <- vapply(mp, function(x) any(is.na(x)), FALSE)
  colsWithNa <- names(colsWithNa)[colsWithNa]
  if (length(colsWithNa)) {
    warning("The following parameters are not specified in all conditions. ",
            "In the unspecified conditions, they are set to 1 (should not make a difference):\n",
            paste0(colsWithNa, collapse = ","))
    mp[,(colsWithNa):=lapply(.SD, function(x) replace(x, is.na(x),1)), .SDcols = colsWithNa]
  }

  mp
}



#' Create the prior function specified in parameters_df
#'
#'
#'
#' @param pe
#' @param FLAGuseNominalCenter
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
petab_createObjPrior <- function(pe, FLAGuseNominalValueAsCenter) {
  # [ ] Todo: Move this function to dMod-*.R file

  p <- copy(pe$parameters)
  p <- p[estimate == 1]


  notImplemented <- which(p$objectivePriorType != "parameterScaleNormal")
  if (length(notImplemented))
    stop("objectivePriorType not implemented: ", paste0(unique(p$objectivePriorType[notImplemented]),collapse = ","))

  p[,`:=`(mu = as.numeric(gsub(";.*", "", objectivePriorParameters)))]
  if (FLAGuseNominalValueAsCenter) {
    # Apply transformation before using as center
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
#' @param pe
#' @param whichBoundary
#'
#' @return
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
#' @param filename,width,height,scale,units see [ggplot2::ggsave()]
#'
#' @return ggplot
#'
#' @importFrom ggforce n_pages
#' @importFrom conveniencefunctions cf_outputFigure
#' @importFrom cOde getSymbols
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
petab_plotData <- function(petab,
                           aeslist = NULL,
                           FLAGUseObservableTransformation = TRUE,
                           FLAGmeanLine = TRUE,
                           ggCallback = list(
                             facet_wrap_paginate(~observableId, nrow = 4, ncol = 4, scales = "free")),
                           filename = NULL,
                           FLAGfuture = TRUE,
                           width = 29.7, height = 21, scale = 1, units = "cm"
) {

  # create plotting data.table
  dplot <- petab_joinDCO(petab)
  # apply log transformation to data when applicable
  if (FLAGUseObservableTransformation) {
    dplot[,`:=`(measurement = eval(parse(text = paste0(observableTransformation, "(measurement)"))),
               observableId = paste0(observableTransformation, " ", observableId)),
          by = 1:nrow(dplot)]
  }

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
  pl <- cfggplot()
  if (FLAGmeanLine) { # Add first so te lines don't mask the points
    aesmean0 <- list(linetype = ~conditionId, group = as.formula(paste0("~ interaction(", paste0(setdiff(byvars, c("time")), collapse = ", "), ")")))
    aesmeanlist <- c(aeslist, aesmean0[setdiff(names(aesmean0), names(aeslist))])
    pl <- pl + geom_line(do.call(aes_q, aesmeanlist), data = dmean)
  }
  pl <- pl + geom_point(do.call(aes_q, aeslist), data = dplot)
  for (plx in ggCallback) pl <- pl + plx

  # Print paginate message so user doesnt forget about additional pages
  message("Plot has ", ggforce::n_pages(pl), " pages\n")

  # output
  cf_outputFigure(pl = pl, filename = filename, width = width, height = height, scale = scale, units = units, FLAGFuture = FLAGfuture)
}

# -------------------------------------------------------------------------#
# Overview tables ----
# -------------------------------------------------------------------------#

#' Generate overview table which observables are in which condition
#'
#' @param pe [petab()] object
#' @param Ntruncate truncate pasted observables at this many characters
#' @param ...
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' @importFrom conveniencefunctions cfoutput_MdTable
petab_overviewObsPerCond <- function(pe, Ntruncate = 1000, ...) {
  dx <- petab_joinDCO(pe)
  if ("conditionName" %in% names(dx)) dx[,`:=`(conditionName = conditionId)]
  dx <- dx[,list(observableId = paste0(sort(unique(observableId)), collapse = ",")),
           by = c("conditionId", "conditionName")]
  dx <- dx[,`:=`(observableId = substr(observableId, 1, Ntruncate))]
  cfoutput_MdTable(dx, ...)
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
#' @examples
petab_getParametersExperimentalCondition <- function(experimentalCondition) {
  ec <- copy(experimentalCondition)
  ec[,`:=`(conditionId = NULL, conditionName = NULL)]
  ec <- lapply(ec, getSymbols)
  ec <- do.call(c,ec)
  ec <- unique(ec)
  ec
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
#' @examples
#' pe <- petab_exampleRead("01")
#' parsEst <- petab_getPars_estScale(pe)
petab_getPars_estScale <- function(pe) {
  p <- copy(pe$parameters)
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
#' @return updated pe_pa1
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Parameter wrangling
#'
#' @examples
petab_mergeParameterValues <- function(pe1,pe2, mergeCols = c("nominalValue", "parameterScale")) {
  # Get values
  pe_pa1 <- pe1$parameters
  pe_pa2 <- pe2$parameters
  mergeColsPA2 <- paste0(mergeCols, "PA2")

  # Merge, update, clean
  pe_pa <- merge(pe_pa1, pe_pa2[,c("parameterId", ..mergeCols)], by = "parameterId", all.x = TRUE, all.y = FALSE, suffixes = c("", "PA2"))
  pe_pa[!is.na(nominalValuePA2),(mergeCols) := lapply(.SD, function(x) x), .SDcols = mergeColsPA2]
  pe_pa[,(mergeColsPA2) := NULL]

  # Update pe1
  pe1$parameters <- pe_pa
  pe1
}


# -------------------------------------------------------------------------#
# Scales ----
# -------------------------------------------------------------------------#

#' Title
#'
#' @param parameterScale vector with entries "lin", "log", "log10"
#'
#' @return
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
# Todolist/Wishlist ----
# -------------------------------------------------------------------------#

# [ ] petab_meta = list(meta_measurementData, meta_observables,
#                       meta_experimentalConditions, meta_parameters)
#   * Tables contain e.g. units, other annotations etc
#   * Could be used for programming on the DCO,
#     e.g. with meta_measurementData$cellline, it's easy to parameterize by cell line
# [ ] petab_simulate,
# [ ] petab_createObj_prior: Can be called either in petab_fit or in importPetab*
# [ ] petab_fit(FLAGrunbg,FLAGslurm)
# [ ] petab_updateParameters Incorporate fit values
# [ ] petab_addFitResults
#   * Elegant way to deal with waterfalls needed
# [ ] petab_profile(opt.runbg = NULL or runbg_list())
#   * petab_addProfileResults Elegant way to deal with profiles needed
# [ ] runbg_list()
#   * Machine: c("knecht1", "knecht2") or "cluster"
#   * nodes
#   * cores
#   * jobname
# [ ] sbml - getSpeciesInfo: speciesUnit
# [ ] sbml - initialAssignments
# [ ] sbml - assignmentRules
# [x] petab_import: append original petab

# * Sample from prior, ...
# * petablint ...


