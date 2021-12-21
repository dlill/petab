#' dMod Observables to PEtab converter
#' 
#' @description Function that takes an eqnvec object from dMod and converts it
#' to the desired PEtab format including renaming of scaling and offset-
#' parameters.
#' 
#' @param obs \link{eqnvec}
#' @param param_identifier vector of strings identifying parameters that should
#' be replaced in the PEtab format. The passed phrases must be appear in all 
#' parameters that should be replaced and in none that not. Per default it is
#' assumed that there are (only) \code{scale} and \code{offset} parameters
#' identifiable by these exact phrases
#' 
#' @return data.table formatted as required by PEtab
#' 
#' @author Severin Bang <severin.bang@physik.uni-freiburg.de>
#' @examples 
#' 
#' ## generate 'eqnvec' object with four parameters variably formatted 
#' observ <- eqnvec(test = "parafoo + meter_bar * cyberglarb(zorbhyper)")
#' 
#' ## define what the parameter identiyfyers are
#' identifiers <- c("foo", "bar", "glarb", "zorb")
#' 
#' ## execute function
#' out <- dmod_obs_2_petab_obs(obs = observ, param_identifier = identifiers)
#' 
#' @export


dmod_obs_2_petab_obs <- function (
  obs, param_identifier = c("scale", "offset")
) {
  if(!is.eqnvec(obs)) {
    stop("Parameter 'obs' is not of type 'eqnvec'.")
  }
  
  out <- data.table::data.table(
    observableId = attr(obs, "names"),
    observableFormula = as.character(obs)
  )
  vsub <- Vectorize(sub)

  for (i in seq_along(param_identifier)) {
    prefix <-  paste0("observableParameter", i ,"_")
    out[
      ,
      observableFormula := vsub(
        paste0("\\b\\w*",param_identifier[i],"\\w*"),
        paste0(prefix, observableId),
        observableFormula)
    ]
  }
  out <- data.table::copy(out)
  
  return(out)
}
