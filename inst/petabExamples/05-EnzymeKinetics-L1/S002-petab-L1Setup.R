# library(petab)
devtools::load_all("~/Promotion/Promotion/Projects/petab")
setwd(tempdir())
petab_python_setup()
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# -------------------------------------------------------------------------#
# Load petab ----
# -------------------------------------------------------------------------#
pe <- pe0 <- readPetab("basemodel")

# -------------------------------------------------------------------------#
# Define inputs ----
# -------------------------------------------------------------------------#

# * Parameter formula injection
#   * pbase_scale == log10,exp: pinner := pnominal * fc
#   * pbase_scale == lin:       pinner := pnominal + fc
#   * => function(pouter, pouterscale)
# * Set parameter values
#   * Base conditions fc_outerscale = 0
#   * Affected conditions (fc = parametername_conditionIdentifier)
#   * Need "augmented" condition grid to summarize e.g. multiple conditions for a cell line
#   * => function(conditionIds = list(baseConditions, L1Conditions, unaffectedConditions (treated like base conditions)),
#                 parameterAppendixColumn,
#                 pouter)
# * Add parameters to pe_pa
#   * => Augment petab_create_parameters_df: recognize fcL1_*-parameters as L1-parameters, look up the base parameter for scale, use priorScale Laplace, set sensible defaults


#' Create L1 Trafo, depending on parameterScale
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
  p[,`:=`(parameterFormula = paste0(parameterId, ifelse(grepl("lin", parameterScale), " + ", " * "), parameterIdL1))]
  p[,list(parameterId, parameterFormula, trafoType = "L1")]
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
  if (!"L1Spec" %in% names(pe$experimentalCondition)) stop("Please add column 'L1Spec' to pe$experimentalCondition before calling this function")
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
#'
#' @examples
#' pe <- petab_exampleRead("04")
#' pe_L1_createL1Problem(pe, c("kcat", "E"), conditionSpecL1_reference = "C1", j_conditionSpecL1 = conditionId)
pe_L1_createL1Problem <- function(pe, parameterId_base, conditionSpecL1_reference, j_conditionSpecL1 = conditionId) {
  # 1 Create parameterFormulaInjection
  pe <- pe_L1_updateParameterFormulaInjection(pe, parameterId_base)
  # 2 Add columns of L1 parameters to experimentalCondition
  # Don't know how 
  sj <- substitute(j_conditionSpecL1)
  pe$experimentalCondition[,`:=`(L1Spec = eval(eval(sj)))]
  pe <- pe_L1_addL1ParsToExperimentalCondition(pe, parameterId_base, conditionSpecL1_reference)
  # 3 Re-create parameters_df including L1 parameters
  pepaL1 <- petab_create_parameter_df(pe)
  pe$parameters <- petab_parameters_mergeParameters(pepaL1, pe$parameters)
  pe
}

parameterId_base <- c("kcat","E")
conditionSpecL1_reference <- "C1"
pe_L1_createL1Problem(pe, parameterId_base, conditionSpecL1_reference, j_conditionSpecL1 = conditionId)


# Exit ---- 
