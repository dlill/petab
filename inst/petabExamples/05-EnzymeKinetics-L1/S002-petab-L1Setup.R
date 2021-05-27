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


parameterId_base <- c("kcat","E")
L1_getParameterFormulaInjection <- function(pe, parameterId_base) {
  p <- pe$parameters[parameterId %in% parameterId_base]
  p[,`:=`(parameterIdL1 = paste0("L1_", parameterId))]
  p[,`:=`(parameterFormula = paste0(parameterId, ifelse(grepl("lin", parameterScale), " + ", " * "), parameterIdL1))]
  p[,list(parameterId, parameterFormula, trafoType = "L1")]
}

L1_pe_updateParameterFormulaInjection <- function(pe, parameterId_base) {
  pfi <- L1_getParameterFormulaInjection(pe, parameterId_base)
  pe$meta$parameterFormulaInjection <- data.table::rbindlist(list(pe$meta$parameterFormulaInjection, pfi), use.names = TRUE, fill = TRUE)
  pe
}


L1_addL1ParsToExperimentalCondition <- function(pe, parameterId_base, conditionSpecL1_reference) {
  pe$experimentalCondition[,(paste0("L1_", parameterId_base)) := lapply(parameterId_base, function(px) paste0("L1_", px, "_", L1Spec))]
  pe$experimentalCondition[L1Spec %in% conditionSpecL1_reference,(paste0("L1_", parameterId_base)) := 0]
  pe$experimentalCondition[,`:=`(L1Spec = NULL)]
  pe
}




L1_createL1Problem <- function(pe, parameterId_base, conditionSpecL1_reference, j_conditionSpecL1 = conditionId) {
  # 1 Create parameterFormulaInjection
  pe <- L1_pe_updateParameterFormulaInjection(pe, parameterId_base)
  # 2 Add columns of L1 parameters to experimentalCondition
  sj <- substitute(j_conditionSpecL1)
  pe$experimentalCondition[,`:=`(L1Spec = eval(eval(sj)))]
  pe <- L1_addL1ParsToExperimentalCondition(pe, parameterId_base, conditionSpecL1_reference)
  # 3 Re-create parameters_df including L1 parameters
  pepaL1 <- petab_create_parameter_df(pe)
  pe$parameters <- petab_parameters_mergeParameters(pepaL1, pe$parameters)
  pe
}

parameterId_base <- c("kcat","E")
conditionSpecL1_reference <- "C1"
L1_createL1Problem(pe, parameterId_base, conditionSpecL1_reference, j_conditionSpecL1 = conditionId)
# Exit ---- 
