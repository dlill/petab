#' Collect Parameters from PD
#'
#' @param pd a pd
#'
#' @return named list(fixedPar = value, estimatedPar = "estimate")
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#' @importFrom data.table copy
#'
#' @examples 
#' pd <- petab_exampleRead("01", "pd")
#' petab_select_collectParametersFromPD(pd)
petab_select_collectParametersFromPD <- function(pd) {
  pars <- data.table::copy(pd$pe$parameters)
  parameters <- setNames(pars$nominalValue, pars$parameterId)
  parameters <- as.list(parameters)
  idx_estimate <- pars[,estimate==1]
  parameters[idx_estimate] <- "estimate"
  parameters
}

#' Collect list of estimated parameters
#' 
#' Chooses based on availability
#' 1. If mstrust is available, take from mstrust
#' 2. else take from "base" (which might be a single fit) CAREFUL: base might be worse than base_obsParsFitted, as base is not necessarily fitted
#' 3. else take from "base_obsParsFitted"
#' 4. else take from pd$pars (which might be unfitted)
#' 
#' @param pd a pd
#'
#' @return list(parameter = value)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#' 
#' @examples 
#' pd <- petab_exampleRead("01", "pd")
#' petab_select_collectEstimatedParametersFromPD(pd)
petab_select_collectEstimatedParametersFromPD <- function(pd) {
  estimated_parameters <- if (!is.null(pd$result$mstrust)) as.parvec(pd$result$mstrust) else if (!is.null(pd$result$base)) as.parvec(pd$result$base) else if (!is.null(pd$result$base_obsParsFitted)) as.parvec(pd$result$base_obsParsFitted) else pd$pars
  as.list(estimated_parameters)
}


#' Title
#'
#' @param AIC,BIC  numeric(1L), the respective criteria
#'
#' @return list(AIC, BIC)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
petab_select_collectCriteria <- function(AIC, BIC) {
  list(AIC = AIC, BIC = BIC)
}

#' Title
#'
#' @param pd a pd
#'
#' @return list(AIC = numeric(1L), BIC = numeric(1L))
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#'
#' @examples
petab_select_collectCriteriaFromPD <- function(pd) {
  LL <- if (!is.null(pd$result$mstrust)) min(pd$result$mstrust$value) else if (!is.null(pd$result$base)) min(pd$result$base$value) else if (!is.null(pd$result$base_obsParsFitted)) min(pd$result$base_obsParsFitted$value) else pd$obj(pd$pars)
  AIC <- LL + sum(pd$pe$parameters$estimate) * 2
  BIC <- LL + sum(pd$pe$parameters$estimate) * log(nrow(pd$pe$measurementData))
  
  petab_select_collectCriteria(AIC = AIC, BIC = BIC)
}

#' Collect all information for the reportYaml
#'
#' @param model_id string
#' @param petab_yaml filename of yaml
#' @param sbml filename of sbml
#' @param parameters named list of numeric entries for fixed parameters and "estimate" for estimated parameters
#' @param estimated_parameters named list of numeric parameters
#' @param criteria list(LL, AIC, BIC)
#'
#' @return simply a list of the supplied objects
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
petab_select_createReportYaml <- function(model_id, petab_yaml, sbml, parameters, estimated_parameters, criteria) {
  list(model_id             = model_id,
       petab_yaml           = petab_yaml,
       sbml                 = sbml,
       parameters           = parameters,
       estimated_parameters = estimated_parameters,
       criteria             = criteria)
}

#' Title
#'
#' @param reportYaml reportYaml as list, see [petab_select_createReportYaml]
#' @param filename filename to write to
#'
#' @return NULL, called for side-effect
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#' @importFrom yaml write_yaml
petab_select_writeReportYaml <- function(reportYaml, filename = petab_files(reportYaml$model_id)["reportYaml"]) {
  yaml::write_yaml(reportYaml, file = filename)
  invisible()
}


#' Create the reportYaml from a PD
#'
#' @param pd a pd
#' @param FLAGwriteYaml TRUE: Writes reportYaml into the original petab folder
#'
#' @return reportYaml as list, see [petab_select_createReportYaml]
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#'
#' @examples
#' # create reportYaml
#' pd <- petab_exampleRead("01", "pd")
#' reportYaml <- petab_select_createReportYamlFromPD(pd, FLAGwriteYaml = FALSE)
#' 
#' # write it to disk and look at it
#' tf <- tempfile(fileext = ".yaml")
#' petab_select_writeReportYaml(reportYaml, tf)
#' file.edit(tf)
petab_select_createReportYamlFromPD <- function(pd, FLAGwriteYaml = TRUE) {
  
  # [ ] Paths are debatable. Current choice: 
  #   * In the reportYaml, the petab_yaml should be a relative file.path from the directory of reportYaml
  #   * Currently, the reportYaml is written to the same directory as the petab itself and the relative file.paths are implemented by not walking any directories.
  
  model_id             <- pd$filenameParts$modelname
  petab_yaml           <- basename(petab_files(pd$filenameParts$modelname)["yaml"])
  sbml                 <- basename(petab_files(pd$filenameParts$modelname)["modelXML"])
  parameters           <- petab_select_collectParametersFromPD(pd)
  estimated_parameters <- petab_select_collectEstimatedParametersFromPD(pd)
  criteria             <- petab_select_collectCriteriaFromPD(pd)
  
  reportYaml <- petab_select_createReportYaml(model_id, petab_yaml, sbml, parameters, estimated_parameters, criteria)
  
  if (FLAGwriteYaml) {
    filename <- petab_files(file.path(dirname(pd$filenameParts$.compiledFolder), pd$filenameParts$modelname))["reportYaml"]
    petab_select_writeReportYaml(reportYaml, filename = filename)
    cat("Report yaml written to ", filename, "\n")
    return(invisible(reportYaml)) # if writing to disk, don't print
  }
  reportYaml
}
