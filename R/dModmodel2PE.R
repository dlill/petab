#' Write PE from dMod objects
#'
#' @param ODEmodel eqnlist with the list of equations
#' @param obsFun 
#' @param errormodel 
#' @param data DT with obligatory columns name, time, value, sigma, condition
#' @param parameters 
#' @param trafo 
#' @param estGrid
#' @param fixedGrid
#' @param eventList
#'
#' @return pe: list with PEtab content that can be exported with writePetab()
#' @export
#' @author Svenja Kemmer
#' @md
#'
#' @importFrom data.table data.table
#' @family SBML export
petab_dModmodel2PE <- function(ODEmodel,
                               obsFun,
                               data,  
                               bestfit,
                               trafo,
                               estGrid,
                               fixedGrid,
                               errormodel = NULL,
                               eventList = NULL,
                               priorSigma = NULL,
                               priorCenter = NULL,
                               lb = 6.14e-06, 
                               ub = 162754.8){
  
  
  cat("Writing model ...\n")
  ti <- petab_getTrafoInfo(trafo)
  
  pfl <- petab_getParameterFormulaList(ti)
  pfil <- petab_getParameterFixedList(ti)
  
  # Create pe_mo
  if(is.null(ODEmodel$volumes)) ODEmodel <- eqnlist_addDefaultCompartment(ODEmodel, "cytoplasm")
  
  parInfo <- getParInfo(equationList = ODEmodel, 
                        eventList = eventList, 
                        parameterFormulaList = pfl,
                        parameterFixedList = pfil) 
  speciesInfo <- getSpeciesInfo(equationList = ODEmodel,
                                parameterFormulaList = pfl,
                                parameterFixedList = pfil)
  
  pe_mo <- petab_model(ODEmodel,
                       events = eventList, 
                       parInfo = parInfo, 
                       speciesInfo = speciesInfo,
                       parameterFormulaList = pfl)
  
  
  cat("Writing conditions ...\n")
  pe_ex <- getEXgrid(estGrid, fixedGrid, trafo)
  
  
  cat("Writing observables ...\n")
  obsDF <- as.data.frame(obsFun)
  obsDF$obsFun <- as.character(obsDF$obsFun)
  formula <- NULL
  noiseformula <- NULL
  obsscale <- NULL
  obsParMatch <- NULL
  for(i in 1:length(obsDF$obsFun)){
    obs <- obsDF$obsFun[i]
    obs_name <- rownames(obsDF)[i]
    if(str_detect(obs, "log10\\(")){
      formula <- c(formula, gsub(" ", "", substr(obs,7,length(strsplit(obs, "")[[1]])-1)))
      obsscale <- c(obsscale, "log10")
    } else if (str_detect(obs, "log\\(")){
      formula <- c(formula, gsub(" ", "", substr(obs,5,length(strsplit(obs, "")[[1]])-1)))
      obsscale <- c(obsscale, "log")
    } else {
      formula <- c(formula, gsub(" ", "", obs))
      obsscale <- c(obsscale, "lin")
    }
    
    # replace obspars by standard nomenclature
    obsPars <- setdiff(getSymbols(obs), ODEmodel$states)
    if(length(obsPars) != 0) {
      names(obsPars) <- paste0("observableParameter", 1:length(obsPars), "_", obs_name)
      obsStr <- paste(obsPars, collapse= ";")
      for(i in 1:length(obsPars)){
        op <- obsPars[i]
        op_name <- names(op)
        formula <- gsub(op, op_name, formula)
      }
    } else obsStr <- NA
    
    noise <- errormodel[obs_name]
    # replace noisepars by standard nomenclature
    noisePars <- setdiff(getSymbols(noise), ODEmodel$states)
    if(length(noisePars) != 0) {
      names(noisePars) <- paste0("noiseParameter", 1:length(noisePars), "_", obs_name)
      noiseStr <- paste(noisePars, collapse= ";")
      for(i in 1:length(noisePars)){
        np <- noisePars[i]
        np_name <- names(np)
        noise <- gsub(np, np_name, noise)
      }
    } else noiseStr <- NA
    noiseformula <- rbind(noiseformula, data.table(observableId = obs_name, noiseFormula = noise))
    
    # create match table for observable parameters
    obsParMatch <- rbind(obsParMatch, data.table(observableId = obs_name, observableParameters = obsStr, noiseParameters = noiseStr))
  }
  
  pe_ob <- petab_observables(observableId = rownames(obsDF), 
                             observableName = rownames(obsDF), 
                             observableFormula = formula, 
                             observableTransformation = obsscale)
  
  if(!is.null(errormodel)) pe_ob[noiseformula, ":=" (noiseFormula = i.noiseFormula, 
                                                     noiseDistribution = "normal"), on = .(observableId)]
  
  
  cat("Writing measurements ...\n")
  if(!is.data.frame(data)) data <- as.data.table(as.data.frame(data))
  pe_me <- petab_measurementData(observableId = data$name,
                                 simulationConditionId = data$condition,
                                 measurement = data$value,
                                 time = data$time,
                                 observableParameters = NA_character_,
                                 noiseParameters = data$sigma,
                                 datasetId = "data1", # is adjusted below
                                 replicateId = 1, #[] could be adjusted
                                 preequilibrationConditionId = NA_character_,
                                 datapointId = 1:nrow(data)
                                 
  )
  # transform data on lin scale
  obs2log <- pe_ob$observableId[which(pe_ob$observableTransformation=="log")]
  pe_me$measurement[which(pe_me$observableId%in%obs2log)] <- exp(pe_me$measurement[which(pe_me$observableId%in%obs2log)])
  obs2log10 <- pe_ob$observableId[which(pe_ob$observableTransformation=="log10")]
  pe_me$measurement[which(pe_me$observableId%in%obs2log10)] <- 10^(pe_me$measurement[which(pe_me$observableId%in%obs2log10)])
  
  # add observable and noise parameters
  pe_me[obsParMatch, ":=" (observableParameters = i.observableParameters), on = .(observableId)]
  if(!is.null(errormodel)) pe_me[obsParMatch, ":=" (noiseParameters = i.noiseParameters), on = .(observableId)]
  
  # replace condition specific observable and noise parameters in pe_ex
  obspars <- getSymbols(obsParMatch$observableParameters)
  selpars <- c("conditionId", obspars)
  for(c in unique(pe_me$simulationConditionId)) {
    replacements <- pe_ex[conditionId == c, ..selpars][,-1]
    for(r in 1:ncol(replacements)){
      replpar <- replacements[,..r]
      pe_me[simulationConditionId == c, observableParameters := gsub(names(replpar), replpar, observableParameters)]
    }
  }
  noisepars <- getSymbols(obsParMatch$noiseParameters)
  if(length(noisepars) > 0) {
    selpars <- c("conditionId", noisepars)
    for(c in unique(pe_me$simulationConditionId)) {
      replacements <- pe_ex[conditionId == c, ..selpars][,-1]
      for(r in 1:ncol(replacements)){
        replpar <- replacements[,..r]
        pe_me[simulationConditionId == c, noiseParameters := gsub(names(replpar), replpar, noiseParameters)]
      }
    }
  }
  
  
  # exclude obs and noise pars and all unchanged pars from pe_ex
  pe_ex <- pe_ex[, !..obspars]
  pe_ex <- pe_ex[, !..noisepars]
  pe_ex_orig <- copy(pe_ex)
  pfix <- NULL
  for (c in 1:ncol(pe_ex_orig)){
    mycol <- pe_ex_orig[,..c]
    if(all(names(mycol) == mycol[[1]])) {
      pe_ex[, names(mycol) := NULL] 
      
      # exclude parameters that are fixed to the same value in all conditions
    } 
    else if (nrow(unique(mycol))==1 & all(suppressWarnings(!is.na(as.numeric(mycol[[1]])))) ){
      fixed_in_all_condis <- data.table(parameterId = names(mycol), nominalValue = as.numeric(unique(mycol[[1]])))
      pfix <- rbind(pfix, fixed_in_all_condis)
      pe_ex[, names(mycol) := NULL]
    }
  }
  
  # exclude fixed pars already defined in the model from pfix
  pfix <- pfix[!parameterId %in% pfil$parameterId]
  if(nrow(pfix)==0) pfix <- NULL
  
  # adjust datasetId according to scale and offset
  count <- 1
  for (d in unique(pe_me$observableParameters)){
    pe_me[observableParameters == d, datasetId := paste0("dataset", count)]
    count <- count + 1
  }
  
  
  cat("Initialize PE ...\n")
  pe <- petab(model = pe_mo,
              experimentalCondition = pe_ex,
              measurementData = pe_me,
              observables = pe_ob)
  
  
  cat("Writing parameters ...\n")
  pe$parameters <- petab_create_parameter_df(pe, priorPars = paste0(priorCenter, ";", priorSigma))
  if(!is.null(priorSigma)) {
    pe$parameters$objectivePriorType <- "parameterScaleNormal"
  } else pe$parameters$objectivePriorType <- NA_character_
  
  # adjust bounds
  pe$parameters$lowerBound <- lb
  pe$parameters$upperBound <- ub
  
  # add bestfit
  pe$parameters$estimate <- 0
  if(attr(pfl, "generalScale")=="log") bestfit <- exp(bestfit)
  if(attr(pfl, "generalScale")=="log10") bestfit <- 10^(bestfit)
  bestfitDT <- data.table(parameterId = names(bestfit), 
                          parameterScale = attr(pfl, "generalScale"),
                          nominalValue = bestfit, 
                          estimate = 1)
  pe$parameters <- petab_parameters_mergeParameters(pe$parameters, bestfitDT)
  
  # adjust scale of pfil pars
  pe$parameters[parameterId %in% pfil$parameterId, parameterScale := "lin"]
  
  # assign parameters that are fixed to the same value in all conditions 
  if(!is.null(pfix)) {
    pe$parameters <- petab_parameters_mergeParameters(pe$parameters, pfix)
    pe$parameters[parameterId %in% pfix$parameterId, parameterScale := "lin"]
  }
  
  pe
}


#' Get condition specific parameters
#'
#' @param est.grid as output by dMod::getParGrids()[[1]]
#' @param fixed.grid as output by dMod::getParGrids()[[2]]
#' @param scale 
#'
#' @return merge of est.grid and fixed.grid
#' @export
#' @author Svenja Kemmer
#' @md
#' @family Parameter wrangling
#'
#' @importFrom data.table data.table 
getEXgrid <- function(est.grid, fixed.grid, trafo){
  
  est.grid <- as.data.table(est.grid)
  est.grid[est.grid == "dummy"] <- NA
  fixed.grid <- as.data.table(fixed.grid)
  fixed.grid[fixed.grid == "NA"] <- NA
  
  # bring values to lin scale when defined on log or log10 scale
  mycols <- names(fixed.grid)[3:ncol(fixed.grid)]
  pars_on_lin_scale <- names(grep("exp\\(|10\\^\\(", trafo[intersect(mycols, names(trafo))], value = T, invert = T))
  pars_on_log_scale <- setdiff(mycols, pars_on_lin_scale)
  
  ti <- petab_getTrafoInfo(trafo)
  if(setdiff(ti$trafoScale, "lin") == "log") fixed.grid[,(pars_on_log_scale) := lapply(.SD, function(x) exp(as.numeric(x))), .SDcols = pars_on_log_scale]
  if(setdiff(ti$trafoScale, "lin") == "log10") fixed.grid[,(pars_on_log_scale) := lapply(.SD, function(x) 10^(as.numeric(x))), .SDcols = pars_on_log_scale]
  
  shared_columns <- intersect(names(est.grid), names(fixed.grid))
  shared_pars <- setdiff(shared_columns, c("ID", "condition"))
  
  # merge shared columns
  mixed.grid <- NULL
  for(c in unique(est.grid$condition, fixed.grid$condition)){
    sub_fixed <-  fixed.grid[condition == c, fixed.grid[condition == c, !is.na(.SD)], with=FALSE]
    sub_est <-  est.grid[condition == c, est.grid[condition == c, !is.na(.SD)], with=FALSE]
    sub_merge <- merge(sub_fixed, sub_est, by = c("ID", "condition"))
    
    mixed.grid <- rbind(mixed.grid, sub_merge)
  }
  
  mixed.grid[,ID:= condition]
  setnames(mixed.grid, c("ID", "condition"), c("conditionId", "conditionName"))
  
  mixed.grid
}

#' Get parameter infos from dMod trafo
#'
#' @param trafo as output by dMod::define()
#'
#' @return DT with columns parameterId, parameterFormula and trafoScale
#' @export
#' @author Svenja Kemmer
#' @md
#' @family Parameter wrangling
#'
#' @importFrom data.table data.table
petab_getTrafoInfo <- function(trafo){
  
  trafoDF <- as.data.frame(trafo)
  trafoDF$trafo <- as.character(trafoDF$trafo)
  trafoDF$name <- NA
  trafoDF$scale <- NA
  parscale <- NULL
  for(i in 1:length(trafoDF$trafo)){
    par_value <- trafoDF$trafo[i]
    par_name <- rownames(trafoDF)[i]
    
    if (suppressWarnings(!is.na(as.numeric(par_value)))){
      par_scale <- "lin"
    } else if(str_detect(par_value, "exp\\(")){
      
      par_value_spl <- strsplit2(par_value, "exp\\(", type = "before")[[1]]
      par_value_spl <- unlist(strsplit2(par_value_spl, "\\)", type = "after"))
      par_value <- NULL
      for(el in par_value_spl){
        if (str_detect(el, "exp\\(")){
          el <- gsub("exp\\(", "", el)
          el <- gsub("\\)", "", el)
          par_value <- paste0(par_value, el)
          par_scale <- "log"
        } else par_value <- paste0(par_value, el)
      }
      
      
    } else if (str_detect(par_value, "10\\^\\(")){
      
      par_value_spl <- strsplit2(par_value, "10\\^\\(", type = "before")[[1]]
      par_value_spl <- unlist(strsplit2(par_value_spl, "\\)", type = "after"))
      par_value <- NULL
      for(el in par_value_spl){
        if (str_detect(el, "10\\^\\(")){
          el <- gsub("10\\^\\(", "", el)
          el <- gsub("\\)", "", el)
          par_value <- paste0(par_value, el)
          par_scale <- "log10"
        } else par_value <- paste0(par_value, el)
      }
      
    } 
    
    trafoDF$name[i] <- par_name
    trafoDF$trafo[i] <- par_value
    trafoDF$scale[i] <- par_scale
  }
  trafoDF <- as.data.table(trafoDF)
  trafoDF <- trafoDF[, list(parameterId = name, parameterFormula = trafo, trafoScale = scale)]
  trafoDF <- trafoDF[parameterId!=parameterFormula]
  
  trafoDF
}

#' Get parameter formulas from dMod trafo
#'
#' @param trafoInfo as output by petab_getTrafoInfo()
#'
#' @return DT with columns parameterId, parameterFormula and the attribute generalScale
#' @export
#' @author Svenja Kemmer
#' @md
#' @family Parameter wrangling
#'
#' @importFrom data.table data.table
petab_getParameterFormulaList <- function(trafoInfo){
  
  pfl <- trafoInfo[suppressWarnings(is.na(as.numeric(parameterFormula)))]
  trafoDF <- pfl[, list(parameterId, parameterFormula)]
  
  attr(trafoDF, "generalScale") <- setdiff(unique(pfl$trafoScale), "lin")[1]
  trafoDF
}


#' Get parameter formulas from dMod trafo
#'
#' @param trafoInfo as output by petab_getTrafoInfo()
#'
#' @return DT with columns parameterId, parameterFormula and the attribute generalScale
#' @export
#' @author Svenja Kemmer
#' @md
#' @family Parameter wrangling
#'
#' @importFrom data.table data.table
petab_getParameterFixedList <- function(trafoInfo){
  
  pfl <- trafoInfo[suppressWarnings(!is.na(as.numeric(parameterFormula)))]
  trafoDF <- pfl[, list(parameterId, parameterValue = parameterFormula)]
  trafoDF
}


#' strsplit with additional position argument
#'
#' @param x character to split
#' @param split pattern as split criterion
#' @param type position of split c("remove", "before", "after")
#' @param perl
#' @param ...
#'
#' @return character vector with split x
#' @author Svenja Kemmer
#' @md
#' @family Parameter wrangling
#'
strsplit2 <- function(x,
                      split,
                      type = "remove",
                      perl = FALSE,
                      ...) {
  if (type == "remove") {
    # use base::strsplit
    out <- base::strsplit(x = x, split = split, perl = perl, ...)
  } else if (type == "before") {
    # split before the delimiter and keep it
    out <- base::strsplit(x = x,
                          split = paste0("(?<=.)(?=", split, ")"),
                          perl = TRUE,
                          ...)
  } else if (type == "after") {
    # split after the delimiter and keep it
    out <- base::strsplit(x = x,
                          split = paste0("(?<=", split, ")"),
                          perl = TRUE,
                          ...)
  } else {
    # wrong type input
    stop("type must be remove, after or before!")
  }
  return(out)
}
