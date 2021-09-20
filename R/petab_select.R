# -------------------------------------------------------------------------#
# Tool-agnostic reporting ----
# -------------------------------------------------------------------------#

#' Ensure consistent naming of Criteria
#'
#' @param AIC,BIC,AICc  numeric(1L), the respective criteria
#'
#' @return list(AIC, BIC, AICc)
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
petabSelect_collectReportYamlCriteria <- function(AIC, BIC, AICc) {
  list(AIC = AIC, BIC = BIC, AICc = AICc)
}



#' Collect all information for the reportYaml
#'
#' @param model_id string
#' @param petab_yaml filename of yaml
#' @param sbml filename of sbml
#' @param parameters named list of numeric entries for fixed parameters and "estimate" for estimated parameters
#' @param estimated_parameters named list of numeric parameters
#' @param criteria list(AIC, BIC, AICc)
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


#' Read report yaml content
#'
#' @param filename file path
#'
#' @return the yaml contents as list
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#' @importFrom yaml read_yaml
#'
#' @examples
#' pd <- petab_exampleRead("04", "pd")
#' reportYaml <- pd_petabSelect_reportYaml(pd,F)
#' tf <- tempfile(fileext = ".yaml")
#' petabSelect_writeReportYaml(reportYaml, tf)
#' reportYaml_fromFile <- petabSelect_readReportYaml(tf)
#' 
petabSelect_readReportYaml <- function(filename) {
  yaml::read_yaml(filename)
}

#' Title
#'
#' @param reportYaml 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#' @importFrom data.table data.table
#'
#' @examples
as.data.table.reportYaml <- function(reportYaml) {
  data.table::data.table(model_id             = reportYaml$model_id,
                         nEstimatedPars       = length(reportYaml$estimated_parameters),
                         `-2LL`               = reportYaml$criteria$AIC - 2 * length(reportYaml$estimated_parameters),
                         AIC                  = reportYaml$criteria$AIC, 
                         AICc                 = reportYaml$criteria$AICc,
                         BIC                  = reportYaml$criteria$BIC 
  )
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
  # load petab_select as peps
  peps <- petab_python_setup(FLAGreturnpetabSelect = T)
  
  LL <- if (!is.null(pd$result$mstrust)) min(pd$result$mstrust$value) else if (!is.null(pd$result$base)) min(pd$result$base$value) else if (!is.null(pd$result$base_obsParsFitted)) min(pd$result$base_obsParsFitted$value) else pd$obj(pd$pars)
  nllh <- LL/2
  
  n_estimated <- sum(pd$pe$parameters$estimate)
  n_measurements <- nrow(pd$pe$measurementData)
  
  AIC <- peps$calculate_aic(n_estimated, nllh)
  BIC <- peps$calculate_bic(n_estimated, nllh, n_measurements, n_priors = 0)
  AICc <- peps$calculate_aicc(n_estimated, nllh, n_measurements, n_priors = 0)

  petabSelect_collectReportYamlCriteria(AIC = AIC, BIC = BIC, AICc = AICc)
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
  
  path                 <- dirname(dirname(pd$filenameParts$.compiledFolder))
  # >>>> This is not very intuitive and only specific to petabSelect problems as set up in petabSelect(). 
  # However, I'd like to use the reportYamls outside of petabSelect, as it's very useful to compare models also in a manual model reduction. 
  # If we want to differentiate "problem_id" and "model_id", we should have different fields in the petab_yaml 
  # And, why do we define problem_id but don't use it afterwards?
  # <<<<<<<<<<< ----
  # for me this wrote model_id as ".."
  # model_id             <- basename(dirname(path))
  # problem_id           <- pd$filenameParts$modelname
  
  model_id <- pd$filenameParts$modelname
  
  # for backwards compatibility
  petab_yaml           <- pd_guessPetabYaml(pd)
  sbml                 <- petab_files(petab_yaml)[["modelXML"]]
  parameters           <- pd_petabSelect_collectParameters(pd)
  estimated_parameters <- pd_petabSelect_collectEstimatedParameters(pd)
  criteria             <- pd_petabSelect_collectReportYamlCriteria(pd)
  
  reportYaml <- petabSelect_reportYaml(model_id, petab_yaml, sbml, parameters, estimated_parameters, criteria)
  
  if (FLAGwriteYaml) {
    filename <- petab_files(petab_yaml)[["reportYaml"]]
    petabSelect_writeReportYaml(reportYaml, filename = filename)
    cat("Report yaml written to ", filename, "\n")
    return(invisible(reportYaml)) # if writing to disk, don't print
  }
  reportYaml
}



# -------------------------------------------------------------------------#
# Comparing models ----
# -------------------------------------------------------------------------#

#' Compare models
#'
#' @param reportYamlFilenames 
#' @param sort_by column to sort models
#'
#' @return data.table
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#' @importFrom data.table rbindlist
#'
#' @examples
#' # TODO
petabSelect_compareModels <- function(reportYamlFilenames, sort_by = "BIC") {
  reportYamls <- lapply(reportYamlFilenames, petabSelect_readReportYaml)
  tab <- data.table::rbindlist(lapply(reportYamls, as.data.table.reportYaml))
  tab <- tab[eval(parse(text = paste0("order(", sort_by, ")")))]
  tab
}

#' List yaml files in subdirectories based on a pattern
#'
#' @param path 
#' @param pattern 
#'
#' @return
#' @export
#'
#' @examples
petabSelect_searchYamls <- function(path, pattern) {
  pattern <- paste0(pattern, ".*_report.yaml$")
  list.files(path, pattern, recursive = T, full.names = TRUE)
}


#' Likelihood ratio test
#'
#' @param comparisonTable output from [petabSelect_compareModels()] with two nested models
#' @param alpha Significance level of statistical test
#'
#' @return Prints a message, returns the printed table invisbly
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family petab select
#' @importFrom tibble tribble
#' @importFrom data.table data.table
#'
#' @examples
#' comparisonTable <- data.table(tibble::tribble(
#' ~model_id, ~nEstimatedPars   , ~`-2LL`, 
#' "Model_full"   ,            2, 5,
#' "Model_reduced",            1, 9))
petabSelect_LRT <- function(comparisonTable, alpha = 0.05) {
  comparisonTable <- comparisonTable[order(`-2LL`)]
  likelihoodRatio <- diff(comparisonTable[["-2LL"]])
  df <- -diff(comparisonTable[["nEstimatedPars"]])
  pValue <- 1-pchisq(likelihoodRatio, df)
  if (round(pValue,16) == 0) pValue <- "<1e-16"
  
  tab <- data.table::data.table(df = df, chisq = likelihoodRatio, `p(X>=chisq)` = pValue)
  cat("H0: Reduced model equally good as full model\n-------------------------------------------\n")
  print(tab)
  if (pValue <= alpha) {
    cat("-------------------------------------------\nH0 is rejected, reduced model is worse than full model (p =",alpha,")\n")
  } else cat("H0 is not rejected, reduction is plausible (p =",alpha,")\n")
  invisible(tab)
}

