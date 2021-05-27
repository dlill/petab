library(petab)
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






# Exit ---- 
