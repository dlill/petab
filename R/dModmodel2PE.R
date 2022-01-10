#' Write PE from dMod objects
#'
#' @param modelname 
#' @param ODEmodel 
#' @param obs_fun 
#' @param errormodel 
#' @param data DF with obligatory columns name, time, value, sigma, condition
#' @param parameters 
#' @param trafo 
#' @param est_grid
#' @param fixed_grid
#' @param parUnit
#'
#' @return pe list
#' @export
#' @author Svenja Kemmer
#' @md
#'
#' @importFrom parallel mclapply
petab_dModmodel2PE <- function(modelname,
                               ODEmodel,
                               obs_fun,
                               errormodel,
                               data,  
                               parameters,
                               trafo,
                               est_grid,
                               fixed_grid){
                    
  
  cat("Writing model ...\n")
  # Create parameterFormulaInjection from trafo
  
  # pfi <- petab_parameterFormulaInjection()
  # Create pe_mo
  parInfo <- getParInfo(equationList = reactions, 
                        eventList = eventlist, 
                        unit = "identity") # include parameterFormulaList
  speciesInfo <- getSpeciesInfo(equationList = reactions) # include parameterFormulaList
  
  pe_mo <- petab_model(reactions,
                       events = eventlist, 
                       parInfo = parInfo, 
                       speciesInfo = speciesInfo)
  
  
  cat("Writing conditions ...\n")
  pe_ex <- getEXgrid(est_grid, fixed_grid)
  
  
  cat("Writing observables ...\n")
  obsDF <- as.data.frame(obs_fun)
  obsDF$obs_fun <- as.character(obsDF$obs_fun)
  formula <- NULL
  obsscale <- NULL
  obsParMatch <- NULL
  for(i in 1:length(obsDF$obs_fun)){
    obs <- obsDF$obs_fun[i]
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
    
    # replace scale and offset by standard nomenclature
    obsPar1 <- grep("scale", getSymbols(obs), value = T)
    obsPar2 <- grep("offset", getSymbols(obs), value = T)
    if(str_detect(obs, "scale") & str_detect(obs, "offset")){
      formula <- gsub(obsPar1, paste0("observableParameter1_", obs_name), formula)
      formula <- gsub(obsPar2, paste0("observableParameter2_", obs_name), formula)
      obsStr <- paste0(obsPar1,";", obsPar2)
    } else if (str_detect(obs, "scale")){
      formula <- gsub(obsPar1, paste0("observableParameter1_", obs_name), formula)
      obsStr <- obsPar1
    } else if (str_detect(obs, "offset")){
      formula <- gsub(obsPar2, paste0("observableParameter1_", obs_name), formula)
      obsStr <- obsPar2
    }
    
    # create match table for observable parameters
    obsParMatch <- rbind(obsParMatch, data.table(observableId = obs_name, observableParameters = obsStr))
  }
  
  pe_ob <- petab_observables(observableId = rownames(obsDF), 
                             observableName = rownames(obsDF), 
                             observableFormula = formula, 
                             observableTransformation = obsscale)
                      
  #[] adjust for multiple noise parameters
  pe_ob[,`:=`(noiseFormula = paste0("noiseParameter1_", observableId),
              noiseDistribution = "normal")]
  
  noiseParMatch <- data.table(observableId=names(errormodel), noiseParameters=errormodel)
  
  
  cat("Writing measurements ...\n")
  pe_me <- petab_measurementData(observableId = data$name,
                                 simulationConditionId = data$condition,
                                 measurement = data$value,
                                 time = data$time,
                                 observableParameters = NA_character_,
                                 noiseParameters = NA_character_,
                                 datasetId = "data1", # is adjusted below
                                 replicateId = 1, #[] could be adjusted
                                 preequilibrationConditionId = NA_character_,
                                 datapointId = 1:nrow(data)
                      
  )
  # add observable parameters
  pe_me[obsParMatch, observableParameters := i.observableParameters, on = .(observableId)]
  # replace condition specific obspars
  obspars <- getSymbols(obsParMatch$observableParameters)
  selpars <- c("conditionId", obspars)
  for(c in unique(pe_me$simulationConditionId)) {
    replacements <- pe_ex[conditionId == c, ..selpars][,-1]
    for(r in 1:ncol(replacements)){
      replpar <- replacements[,..r]
      pe_me[simulationConditionId == c, observableParameters := gsub(names(replpar), replpar, observableParameters)]
    }
  }
  
  # add noise parameters (to be extended for several noise pars)
  if (!is.null(errormodel)){
    pe_me[noiseParMatch, noiseParameters := i.noiseParameters, on = .(observableId)]
  } else if (!is.na(data$sigma)){
    pe_me[, noiseParameters := data$sigma]
  } else print("Warning: No errors provided!")
  
  # exclude obs and noise pars and all unchanged pars from pe_ex
  pe_ex <- pe_ex[, !..obspars]
  noisepars <- getSymbols(pe_me$noiseParameters)
  pe_ex <- pe_ex[, !..noisepars]
  pe_ex_orig <- copy(pe_ex)
  for (c in 1:ncol(pe_ex_orig)){
    mycol <- pe_ex_orig[,..c]
    if(all(names(mycol) == mycol[[1]])) pe_ex[, names(mycol) := NULL] 
  }
  
  # adjust datasetId according to scale and offset
  count <- 1
  for (d in unique(pe_me$observableParameters)){
    pe_me[observableParameters == d, datasetId := paste0("dataset", count)]
    count <- count + 1
  }
  
  # if(!is.null(attr(parameters, "parscales"))){
  #   out <- data.frame(parameterId = names(parameters), parameterScale=attr(parameters, "parscales"))
  # } else out <- data.frame(parameterId = names(parameters), parameterScale="log")
  # if(!is.null(attr(parameters, "lowerBound"))){
  #   out <- merge(out, lowerBound = attr(parameters, "lowerBound"))  
  # } else out <- cbind(out, data.frame(parameterId = names(parameters), lowerBound="-12"))
  # if(!is.null(attr(parameters, "upperBound"))){
  #   out <- cbind(out, upperBound = attr(parameters, "upperBound"))  
  # } else out <- cbind(out, data.frame(parameterId = names(parameters), upperBound="12"))
  # out <- cbind(out, estimate=parameters)
  # write_tsv(out, path = paste0(exportwd,"/parameters_",modelname,".tsv"))
  
  
  cat("Initialize PE ...\n")
  pe <- petab(model = pe_mo,
              experimentalCondition = pe_ex,
              measurementData = pe_me,
              observables = pe_ob)
  
  
  cat("Writing parameters ...\n")
  
  pe$parameters <- petab_create_parameter_df(pe)
  pe$parameters$objectivePriorType <- NA_character_
  
  # add bestfit
  # scale
  bestfitDT <- data.table(parameterId = names(bestfit), nominalValue = bestfit)
  pe$parameters <- petab_parameters_mergeParameters(pe$parameters, bestfitDT)
  
  pe
  # cat(green(paste0("PEtab files written to Export/", modelname, "/\n")))
}



getEXgrid <- function(est.grid, fixed.grid){
  
  est.grid <- as.data.table(est.grid)
  est.grid[est.grid == "dummy"] <- NA
  fixed.grid <- as.data.table(fixed.grid)
  fixed.grid[fixed.grid == "NA"] <- NA
  
  shared_columns <- intersect(names(est.grid), names(fixed.grid))
  shared_pars <- setdiff(shared_columns, c("ID", "condition"))
  mixed.grid <- rbind(est.grid[complete.cases(est.grid),..shared_columns],
                      fixed.grid[complete.cases(fixed.grid),..shared_columns])
  mixed.grid <- mixed.grid[order(ID)]
  pe_ex <- merge(est.grid[,!..shared_pars], fixed.grid[,!..shared_pars], by = c("ID", "condition"))
  pe_ex <- merge(pe_ex, mixed.grid, by = c("ID", "condition"))
  
  # pe_ex[,ID:= paste0("condition", ID)]
  pe_ex[,ID:= condition]
  setnames(pe_ex, c("ID", "condition"), c("conditionId", "conditionName"))
  
  pe_ex
}
# 
# petab_getParameterFormulas <- function(trafo){
#   
#   trafoDF <- as.data.frame(trafo)
#   trafoDF$trafo <- as.character(trafoDF$trafo)
#   formula <- NULL
#   parscale <- NULL
#   obsParMatch <- NULL
#   for(i in 1:length(trafoDF$trafo)){
#     par_value <- trafoDF$trafo[i]
#     par_name <- rownames(trafoDF)[i]
#     if(str_detect(par_value, "exp\\(")){
#       formula <- c(formula, gsub("exp\\(", "", par_value))
#       parscale <- c(parscale, "log")
#     } else if (str_detect(obs, "log\\(")){
#       formula <- c(formula, gsub(" ", "", substr(obs,5,length(strsplit(obs, "")[[1]])-1)))
#       parscale <- c(parscale, "log")
#     } else {
#       formula <- c(formula, gsub(" ", "", obs))
#       parscale <- c(parscale, "lin")
#     }
# }
