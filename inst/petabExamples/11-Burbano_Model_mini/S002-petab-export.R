# -------------------------------------------------------------------------#
# PEtab export ----
# -------------------------------------------------------------------------#
#
# S200-Export_PEtab.R
#
# [PURPOSE]
# 
# Export established dMod models as petab
#
#
# [AUTHOR]
# Svenja Kemmer
#
# [Date]
# Fri Feb 11 11:37:42 2022
#
# [Info]
# Run after having sourced the model to export
#
# rm(list = ls(all.names = TRUE))
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# .. Load libraries -----
library(petab)

.currentWD <- getwd()
# -------------------------------------------------------------------------#
# 1 Load bestfit to export ----
# -------------------------------------------------------------------------#

# .. 1 load bestfit -----
myframe <- readRDS(file = file.path(.fitFolder, paste0("myframe.rds")))
step <- 1
bestfit <- as.parvec(myframe, step)
fitvalue <- round(myframe[step]$value, digits = 2)
print(fitvalue)


# .. eventually rename pars -----

# names(bestfit) <- str_replace(names(bestfit), "degrad_", "degrade_")
# names(bestfit) <- str_replace(names(bestfit), "mutation_off", "mutation")
# names(bestfit) <- str_replace(names(bestfit), "drug_b_", "drug_")
# names(bestfit) <- str_replace(names(bestfit), "cRAF_on", "pcRAF")

# .. Collect PE -----

pe <- petab_dModmodel2PE(ODEmodel=reactions,
                         obsFun=observables,
                         # errormodel=errors,
                         data=mydata,  
                         bestfit=bestfit,
                         trafo=trafo,
                         estGrid=est.grid,
                         fixedGrid=fixed.grid,
                         # eventList=events,
                         priorSigma = 12,
                         priorCenter = -1
)

pe$parameters <- pe$parameters[!is.na(parameterId)]
# .. Write petab -----

filename <- "Model"
writePetab(pe, file.path(.resultsFolder, "05-PEtab",filename))

# .. Check obj value -----
obj(bestfit)$value


# -------------------------------------------------------------------------#
# 2 Debugging ----
# -------------------------------------------------------------------------#


ODEmodel=reactions
obsFun=observables
errormodel=NULL
data=mydata
estGrid=est.grid
fixedGrid=fixed.grid
eventList=NULL
lb = 6.14e-06
ub = 162754.8
priorSigma = 12
priorCenter = -1




# Exit ----
future::plan("sequential")


