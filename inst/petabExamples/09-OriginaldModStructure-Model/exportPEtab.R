writePE <- function(modelname="test",
                    ODEmodel=reactions,
                    obs_fun=observables,
                    errormodel = errors,
                    data=pred, # with obligatory columns name, time, value, sigma, condition
                    parameters=bestfit,
                    conditions=condition.grid,
                    est_grid = est.grid,
                    fixed_grid = fixed.grid){
  # dir.create("Export/", showWarnings = FALSE)
  # exportwd <- paste0("Export/",modelname)
  # dir.create(exportwd, showWarnings = FALSE)
  
  
  
  cat("Writing model ...\n")
  parInfo <- getParInfo(reactions) # modify as it writes random stuff in it
  speciesInfo <- getSpeciesInfo(reactions) # modify as it writes random stuff in it
  pe_mo <- petab_model(reactions,events = NULL,parInfo = parInfo, speciesInfo = speciesInfo)
  
  
  
  cat("Writing conditions ...\n")
  pe_ex <- getEXgrid(est_grid, fixed_grid)
  
  
  
  cat("Writing measurements ...\n")
  pe_me <- data.table(observableId = data$name,
                      simulationConditionId = data$condition,
                      measurement = data$value,
                      time = data$time,
                      observableParameters = NA_character_,
                      noiseParameters = data$sigma,
                      datasetId = "data1",
                      replicateId = 1,
                      preequilibrationConditionId = NA_character_,
                      datapointId = 1:nrow(data)
  )
  
  
  
  cat("Writing observables ...\n")
  obsDF <- as.data.frame(obs_fun)
  obsDF$obs_fun <- as.character(obsDF$obs_fun)
  formula <- NULL
  trafo <- NULL
  for(i in 1:length(obsDF$obs_fun)){
    obs <- obsDF$obs_fun[i]
    if(str_detect(obs, "log10\\(")){
      formula <- c(formula, gsub(" ", "", substr(obs,7,length(strsplit(obs, "")[[1]])-1)))
      trafo <- c(trafo, "log10")
    } else if (str_detect(obs, "log\\(")){
      formula <- c(formula, gsub(" ", "", substr(obs,5,length(strsplit(obs, "")[[1]])-1)))
      trafo <- c(trafo, "log")
    } else {
      formula <- c(formula, gsub(" ", "", obs))
      trafo <- c(trafo, "lin")
    }
  }

  pe_ob <- data.table(observableId = rownames(obsDF), 
                      observableName = rownames(obsDF), 
                      observableFormula = formula, 
                      observableTransformation = trafo)
  pe_ob[observableId == names(errormodel),`:=`(noiseFormula =  errormodel[names(errormodel)], 
                                               noiseDistribution = "normal")]

  # write_tsv(out, path = paste0(exportwd,"/observables_",modelname,".tsv"))
  
  
  
  

  
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
  
  pe_ex[,ID:= paste0("condition", ID)]
  setnames(pe_ex, c("ID", "condition"), c("conditionId", "conditionName"))
  
  pe_ex
}
