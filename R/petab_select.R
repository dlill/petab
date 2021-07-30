# -------------------------------------------------------------------------#
# Tool-agnostic reporting ----
# -------------------------------------------------------------------------#

#' Ensure consistent naming of Criteria
#'
#' @param AIC,BIC  numeric(1L), the respective criteria
#'
#' @return list(AIC, BIC)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
petabSelect_collectReportYamlCriteria <- function(AIC, BIC) {
  list(AIC = AIC, BIC = BIC)
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
petabSelect_reportYaml <- function(model_id, petab_yaml, sbml, parameters, estimated_parameters, criteria) {
  
  # Some checks. Add more to increase code safety
  if (!is.list(parameters)) stop("parameters should be a list (is ", class(parameters), ")")
  if (!is.list(estimated_parameters)) stop("parameters should be a list (is ", class(estimated_parameters), ")")
  
  list(model_id             = model_id,
       petab_yaml           = petab_yaml,
       sbml                 = sbml,
       criteria             = criteria,
       parameters           = parameters,
       estimated_parameters = estimated_parameters
       )
}

#' Title
#'
#' @param reportYaml reportYaml as list, see [petabSelect_reportYaml]
#' @param filename filename to write to
#'
#' @return NULL, called for side-effect
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#' @importFrom yaml write_yaml
petabSelect_writeReportYaml <- function(reportYaml, filename = petab_files(reportYaml$model_id)["reportYaml"]) {
  yaml::write_yaml(reportYaml, file = filename)
  invisible(reportYaml)
}

# -------------------------------------------------------------------------#
# petab-dMod based reporting ----
# -------------------------------------------------------------------------#

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
#' pd_petabSelect_collectParameters(pd)
pd_petabSelect_collectParameters <- function(pd) {
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
#' pd_petabSelect_collectEstimatedParameters(pd)
pd_petabSelect_collectEstimatedParameters <- function(pd) {
  estimated_parameters <- if (!is.null(pd$result$mstrust)) as.parvec(pd$result$mstrust) else if (!is.null(pd$result$base)) as.parvec(pd$result$base) else if (!is.null(pd$result$base_obsParsFitted)) as.parvec(pd$result$base_obsParsFitted) else pd$pars
  as.list(estimated_parameters)
}




#' Collect AIC and BIC from pd
#' 
#' Chooses LL based on availability
#' 1. If mstrust is available, take from mstrust
#' 2. else take from "base" (which might be a single fit) CAREFUL: base might be worse than base_obsParsFitted, as base is not necessarily fitted
#' 3. else take from "base_obsParsFitted"
#' 4. else take from pd$pars (which might be unfitted)
#' 
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
pd_petabSelect_collectReportYamlCriteria <- function(pd) {
  LL <- if (!is.null(pd$result$mstrust)) min(pd$result$mstrust$value) else if (!is.null(pd$result$base)) min(pd$result$base$value) else if (!is.null(pd$result$base_obsParsFitted)) min(pd$result$base_obsParsFitted$value) else pd$obj(pd$pars)
  AIC <- LL + sum(pd$pe$parameters$estimate) * 2
  BIC <- LL + sum(pd$pe$parameters$estimate) * log(nrow(pd$pe$measurementData))
  
  petabSelect_collectReportYamlCriteria(AIC = AIC, BIC = BIC)
}


#' Create the reportYaml from a PD
#'
#' @param pd a pd
#' @param FLAGwriteYaml TRUE: Writes reportYaml into the original petab folder
#'
#' @return reportYaml as list, see [petabSelect_reportYaml]
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#'
#' @examples
#' # create reportYaml
#' pd <- petab_exampleRead("01", "pd")
#' reportYaml <- pd_petabSelect_reportYaml(pd, FLAGwriteYaml = FALSE)
#' 
#' # write it to disk and look at it
#' tf <- tempfile(fileext = ".yaml")
#' petabSelect_writeReportYaml(reportYaml, tf)
#' file.edit(tf)
pd_petabSelect_reportYaml <- function(pd, FLAGwriteYaml = TRUE) {
  
  # [ ] Paths are debatable. Current choice: 
  #   * In the reportYaml, the petab_yaml should be a relative file.path from the directory of reportYaml
  #   * Currently, the reportYaml is written to the same directory as the petab itself and the relative file.paths are implemented by not walking any directories.
  
  model_id             <- pd$filenameParts$modelname
  petab_yaml           <- basename(petab_files(pd$filenameParts$modelname)["yaml"])
  sbml                 <- basename(petab_files(pd$filenameParts$modelname)["modelXML"])
  parameters           <- pd_petabSelect_collectParameters(pd)
  estimated_parameters <- pd_petabSelect_collectEstimatedParameters(pd)
  criteria             <- pd_petabSelect_collectReportYamlCriteria(pd)
  
  reportYaml <- petabSelect_reportYaml(model_id, petab_yaml, sbml, parameters, estimated_parameters, criteria)
  
  if (FLAGwriteYaml) {
    filename <- petab_files(file.path(dirname(pd$filenameParts$.compiledFolder), pd$filenameParts$modelname))["reportYaml"]
    petabSelect_writeReportYaml(reportYaml, filename = filename)
    cat("Report yaml written to ", filename, "\n")
    return(invisible(reportYaml)) # if writing to disk, don't print
  }
  reportYaml
}
