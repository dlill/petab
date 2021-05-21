# library(petab)
devtools::load_all("~/Promotion/Promotion/Projects/petab")
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# -------------------------------------------------------------------------#
# Arguments ----
# -------------------------------------------------------------------------#
filename <- "petab"

# filename <- "/home/daniel/Promotion/Promotion/Software/PEtab/Benchmark-Models-PEtab/Benchmark-Models/Bachmann_MSB2011/"

NFLAGcompile = 0
.compiledFolder = "Compiled"
SFLAGbrowser = "0"

testCases = FALSE
path2TestCases = "PEtabTests/"

# -------------------------------------------------------------------------#
# import_pd_Indiv ----
# -------------------------------------------------------------------------#

# .. Set exit behaviour -----
mywd <- getwd()
on.exit({setwd(mywd)})

# .. Define path to SBML, PEtab and pd files -----
filename  <- path.expand(filename)
path      <- petab_modelname_path(filename)$path
modelname <- petab_modelname_path(filename)$modelname
files     <- petab_files(filename, FLAGTestCase = testCases, FLAGreturnList = TRUE)
filenameParts = list(modelname = modelname,.currentFolder = mywd,.compiledFolder = .compiledFolder,type = "indiv")
rdsfile   <- pd_files(filenameParts)$rdsfile

# .. Read previously imported files -----
dir.create(.compiledFolder, showWarnings = FALSE)
if (NFLAGcompile == 3)
  NFLAGcompile <- as.numeric(!inputFileChanged(files[[1]], rdsfile))*2 # Um die Ecke wegen suboptimaler NFLAGcompile Definition

if(NFLAGcompile > 0) {
  pd <- readPd(rdsfile)
  if (NFLAGcompile == 2) return(pd)
}

## load required packages
require(libSBML) # => Not very nice, better explicitly import the required functions

# .. Read PEtab tables -----
pe <- readPetab(filename)

# .. Model Definition - Equations -----
dummy           <- getReactionsSBML(files$modelXML, files$experimentalCondition)
myreactions      <- dummy$reactions
myreactions_orig <- dummy$reactions_orig
myevents         <- dummy$events
mypreeqEvents    <- dummy$preeqEvents
myobservables    <- getObservablesSBML(files$observables)

# [ ] Pre-Equi Events
if(!is.null(mypreeqEvents))
  stop("Pre-Equilibration is not yet implemented")
# [ ] Need example for preeqEvents


if (!is.null(myevents)){
  stop("events not yet implemented")
  # set remaining event initials to 0
  inits_events <- setdiff(unique(myevents$var), unique(mypreeqEvents$var))
  inits_events <- setNames(rep(0, length(inits_events)), inits_events)
}

# .. Get Data -----
mydataSBML <- getDataPEtabSBML(files$measurementData, files$observables)
mydata     <- mydataSBML$data
myerrors   <- mydataSBML$errors
myerr <- NULL

# .. Define constraints, initials, parameters and compartments -----
cg <- pe$experimentalCondition[,!"conditionName"]

# Get parameter information from all possible sources
# Estimated parameters. Possible sources
# * pe$parameters
parsEst <- pe$parameters[estimate == 1] # replaces myfit_values
parsEst[,`:=`(estValue = eval(parse(text = paste0(parameterScale, "(", nominalValue, ")")))),by = 1:nrow(parsEst)]

# Fixed parameters. Possible sources 
# * pe$parameters
# * SBML parameters
# * SBML compartments
# * SBML inits
# * Events
parsFix <- pe$parameters[estimate == 0] # replaces myconstraints part1
parsFix_SBMLPars <- sbmlImport_getParametersFromSBML(files$modelXML)   # replaces myconstraints part2
parsFix_SBMLPars <- parsFix_SBMLPars[!parameterId %in% c(names(cg), parsFix$parameterId, parsEst$parameterId)] 
parsFix_SBMLComp <- sbmlImport_getCompartmentsFromSBML(files$modelXML) # replaces myconstraints part3
parsFix_SBMLComp <- parsFix_SBMLComp[!parameterId %in% c(names(cg), parsFix$parameterId, parsEst$parameterId)] 
parsFix_SBMLInit <- sbmlImport_getInitialsFromSBML(files$modelXML)$inits_num # Should be handled as well
parsFix_SBMLInit <- parsFix_SBMLInit[!parameterId %in% c(names(cg), parsFix$parameterId, parsEst$parameterId)] 
parsFix <- rbindlist(list(parsFix, parsFix_SBMLPars, parsFix_SBMLComp, parsFix_SBMLInit))
parsFix[,`:=`(estValue = eval(parse(text = paste0(parameterScale, "(", nominalValue, ")")))),by = 1:nrow(parsFix)]

# .. Initialize gridlist and trafo -----
if (grepl(SFLAGbrowser, "2BuildGrids")) browser()
cg <- copy(pe$experimentalCondition)

# .. Add measurement parameters -----
obsParMapping <- petab_getMeasurementParsMapping(pe, column = "observableParameters")
cg <- merge(cg, obsParMapping, by = "conditionId", all.x = TRUE)

errParMapping <- petab_getMeasurementParsMapping(pe, column = "noiseParameters")
cg <- merge(cg, errParMapping, by = "conditionId", all.x = TRUE)

# obsErrPars were already added => remove
parametersObsErr <- petab_getParameterType(pe)
parametersObsErr <- parametersObsErr[parameterType %in% c("observableParameters", "noiseParameters"), parameterId]
parsEstNoObsErr <- parsEst[!parameterId %in% parametersObsErr]
parsFixNoObsErr <- parsFix[!parameterId %in% parametersObsErr]

# .. Add estPars and fixPars symbolically -----
# Determine parameters already specified in cg. Assume they are fully mapped to inner parameters => remove
parsEstNoObsErrNoCS <- parsEstNoObsErr[!parameterId %in% pdIndiv_getParameters_conditionGrid(cg)]
# Add remaining parsEst globally
parametersEstGlobal <- as.data.table(as.list(setNames(nm = parsEstNoObsErrNoCS$parameterId)))
cg <- data.table(cg, parametersEstGlobal)

# Add fixPars symbolically, they will later be replaced by actual values
# Determine parameters already specified in cg. Assume they are fully mapped to inner parameters => remove
parsFixNoObsErrNoCS <- parsFixNoObsErr[!parameterId %in% pdIndiv_getParameters_conditionGrid(cg)]
parametersFixGlobal <- as.data.table(as.list(setNames(nm = parsFixNoObsErrNoCS$parameterId)))
cg <- data.table(cg, parametersFixGlobal)

# .. Parameter scales -----
# Ensure that parameter scales are consistent for each CS <-> BASE mapping
scalesOuter <- c(setNames(parsEst$parameterScale, parsEst$parameterId), setNames(parsFix$parameterScale, parsFix$parameterId))
parametersBase <- setdiff(names(cg), c("conditionId", "conditionName", "ID"))
scalesBase <- vapply(setNames(nm = parametersBase), function(par) {
  px <- cg[[par]]
  px <- suppressWarnings(px[is.na(as.numeric(px))])
  
  if (!length(px)) 
    return("lin")
  
  sx <- scalesOuter[px]
  if (length(unique(sx)) > 1)
    stop("The following parameters refer to the same structural model parameter, but have different ",
         "scales in different conditions. This is not allowed, please fix manually. \n",
         "Base parameter: ", par , "\n",
         "Condition specific parameters: ", paste0(paste0(names(sx), " (", sx, ")"), collapse = ", "), "\n")
  unique(sx)
}, "lin")

# .. Replace symbolic fixPars by their values (on est scale) -----
cg <- as.matrix(cg)
for (par in parsFix$parameterId) cg[cg == par] <- parsFix[parameterId == par, estValue]
cg <- as.data.table(cg)

# .. Split into gridlist -----
gl <- petab::pdIndiv_initializeGridlist(cg)

# .. Build trafo -----
# Initialize
trafo <- setNames(nm = unique(c(getParameters(myreactions),
                                getSymbols(myobservables),
                                setdiff(getSymbols(myerrors), names(myobservables)),
                                getSymbols(as.character(myevents$value)))))
trafo <- trafo[trafo != "time"]

# Insert inits
trafo_inits <- sbmlImport_getInitialsFromSBML(files$modelXML)$inits_sym
trafo[trafo_inits$parameterId] <- trafo_inits$parameterFormula
trafo_assignmentRules <- NULL # [] Todo, if there is anything that needs to be done...

# Insert scales
trafo <- repar("x ~ 10**(x)", trafo = trafo, x = names(which(scalesBase=="log10")))
trafo <- repar("x ~ exp(x)" , trafo = trafo, x = names(which(scalesBase=="log")))


# .. ParameterFormulaInjection -----
if (grepl(SFLAGbrowser, "6ParameterFormulaInjection")) browser()
trafoInjected <- NULL
pfi <- pe$meta$parameterFormulaInjection
if (!is.null(pfi)) {
  # Probably over-cautios: This shouldn't happen. Can probably be removed
  check_pfiEstimated <- intersect(pfi$parameterId, names(gl$est.grid))
  if (length(check_pfiEstimated)) stop("These trafoInjected parameter are in est.grid: ", paste0(check_pfiEstimated, collapse = ", "))
  # Remove injected Parameter from est-trafo
  trafo <- trafo[setdiff(names(trafo), pfi$parameterId)]
  gl$fix.grid[,(pfi$parameterId) := NULL,]
  # trafoInjected
  trafoInjected <- setNames(pfi$parameterFormula, pfi$parameterId)
}







# -------------------------------------------------------------------------#
# OLD ----
# -------------------------------------------------------------------------#
# pd <- importPEtabSBML_indiv("petab", NFLAGcompile = 1, .compiledFolder = "Compiled", SFLAGbrowser = "1")
# 
# # Test model
# pred <- pd$prd(seq(0,100), pd$pars)
# pd$obj_data(pd$pars)
# 
# # Test fitting
# myfit <- trust(pd$obj_data, pd$pars,1,10,iterlim = 1000)
# plotCombined(pd$prd(seq(0,100), pd$pars), pd$dModAtoms$data)
# plotCombined(pd$prd(seq(0,100), myfit$argument), pd$dModAtoms$data)
# 
# # Mstrust
# center <- pepy_sample_parameter_startpoints(pd$pe, n_starts = 8)
# fits <- mstrust(pd$obj_data, center, "fit", fits = 5, iterlim = 20, cores = 8)
# fits <- as.parframe(fits)
# plotValues(fits)
# dMod_saveMstrust(fits, ".", FLAGoverwrite = TRUE)
# unlink("fit", T)


# Exit ----
