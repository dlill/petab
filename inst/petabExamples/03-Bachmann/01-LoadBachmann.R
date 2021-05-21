# -------------------------------------------------------------------------#
# 0 Header ----
# -------------------------------------------------------------------------#
#
# 01-LoadBachmann.R
#
# [PURPOSE]
# 
# Load the Bachmann model
#
#
# [AUTHOR]
# Daniel Lill
#
# [Date]
# Mon May 17 23:08:36 2021
#
# [Template Info]
# .. Template name TGFB_PETABEXPORT -----
# .. Template name TGFB_MSTRUST -----
# .. Template name TGFB_PROFILE -----
# .. Template version 0.0.1 -----
#
rm(list = ls(all.names = TRUE))
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# library(petab)
devtools::load_all("~/Promotion/Promotion/Projects/petab")

.outputFolder <- paste0("Output/")
.compiledFolder <- file.path(.outputFolder, "Compiled")
for(folder in c(.outputFolder,.compiledFolder)) 
if(!dir.exists(folder)) dir.create(folder)

# -------------------------------------------------------------------------#
# 1 Import ----
# -------------------------------------------------------------------------#
pe_file <- "Bachmann_MSB2011"

pe <- readPetab(pe_file)
petab_overviewObsPerCond(pe)

pd <- importPEtabSBML_indiv(pe_file, .compiledFolder = .compiledFolder, NFLAGcompile = 3, SFLAGbrowser = "0")
# pd <- pd_updateEstPars(pd, pd$pars, FLAGupdatePE = FALSE, FLAGsavePd = TRUE)

# -------------------------------------------------------------------------#
# 2 Predict and plot ----
# -------------------------------------------------------------------------#
pd_predictAndPlot2(pd, opt.base = pd_parf_opt.base(T))
pd$obj_data(pd$pars)




# Exit ----
future::plan("sequential")
