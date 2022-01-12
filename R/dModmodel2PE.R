#' Write PE from dMod objects
#'
#' @param ODEmodel 
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
                               errormodel,
                               data,  
                               bestfit,
                               trafo,
                               estGrid,
                               fixedGrid,
                               eventList,
                               lb = 6.14e-06, 
                               ub = 162754.8){
                    
  
  cat("Writing model ...\n")
  pfi <- petab_getParameterFormulas(trafo)

  # Create pe_mo
  if(!is.null(ODEmodel$volumes)) ODEmodel <- eqnlist_addDefaultCompartment(ODEmodel, "cytoplasm")
  
  parInfo <- getParInfo(equationList = ODEmodel, 
                        eventList = eventList, 
                        parameterFormulaList = pfi) 
  speciesInfo <- getSpeciesInfo(equationList = ODEmodel,
                                parameterFormulaList = pfi)
  
  pe_mo <- petab_model(ODEmodel,
                       events = eventList, 
                       parInfo = parInfo, 
                       speciesInfo = speciesInfo)
  
  
  cat("Writing conditions ...\n")
  pe_ex <- getEXgrid(estGrid, fixedGrid)
  
  
  cat("Writing observables ...\n")
  obsDF <- as.data.frame(obsFun)
  obsDF$obsFun <- as.character(obsDF$obsFun)
  formula <- NULL
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
    if(all(names(mycol) == mycol[[1]])) {
      pe_ex[, names(mycol) := NULL] 
      
      # exclude parameters that are fixed to the same value in all conditions
    } else if (nrow(unique(mycol))==1 & all(suppressWarnings(!is.na(as.numeric(mycol[[1]])))) ){ 
      pe_ex[, names(mycol) := NULL]
    }
  }

  
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
  # add pfi as meta
  pe$meta$parameterFormulaInjection <- pfi
  
  pe$parameters <- petab_create_parameter_df(pe)
  pe$parameters$objectivePriorType <- NA_character_
  
  # adjust bounds
  pe$parameters$lowerBound <- lb
  pe$parameters$upperBound <- ub
  
  # add bestfit
  pe$parameters$estimate <- 0
  if(attr(pfi, "generalScale")=="log") bestfit <- exp(bestfit)
  if(attr(pfi, "generalScale")=="log10") bestfit <- 10^(bestfit)
  bestfitDT <- data.table(parameterId = names(bestfit), 
                          parameterScale = attr(pfi, "generalScale"),
                          nominalValue = bestfit, 
                          estimate = 1)
  pe$parameters <- petab_parameters_mergeParameters(pe$parameters, bestfitDT)
  
  pe
}


#' Get condition specific parameters
#'
#' @param est.grid as output by dMod::getParGrids()[[1]]
#' @param fixed.grid as output by dMod::getParGrids()[[2]]
#'
#' @return merge of est.grid and fixed.grid
#' @export
#' @author Svenja Kemmer
#' @md
#' @family Parameter wrangling
#'
#' @importFrom data.table data.table 
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

#' Get parameter formulas from dMod trafo
#'
#' @param trafo as output by dMod::define()
#'
#' @return DT with columns parameterId, parameterFormula and trafoType
#' @export
#' @author Svenja Kemmer
#' @md
#' @family Parameter wrangling
#'
#' @importFrom data.table data.table
petab_getParameterFormulas <- function(trafo){

  trafoDF <- as.data.frame(trafo)
  trafoDF$trafo <- as.character(trafoDF$trafo)
  trafoDF$name <- NA
  trafoDF$scale <- NA
  parscale <- NULL
  for(i in 1:length(trafoDF$trafo)){
    par_value <- trafoDF$trafo[i]
    par_name <- rownames(trafoDF)[i]
    
    if(str_detect(par_value, "exp\\(")){
      
      par_value_spl <- strsplit2(par_value, "exp\\(", type = "before")[[1]]
      par_value_spl <- unlist(strsplit2(par_value_spl, "\\)", type = "after"))
      par_value <- NULL
      for(el in par_value_spl){
        # if (!is.na(tryCatch(eval(parse(text = el)), error = function(x) NA))) {
        #   el <- eval(parse(text = el))
        #   par_value <- paste0(par_value, el)
        #   par_scale <- "lin"
        # } else 
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
        # if (!is.na(tryCatch(eval(parse(text = el)), error = function(x) NA))) {
        #   el <- eval(parse(text = el))
        #   par_value <- paste0(par_value, el)
        #   par_scale <- "lin"
        # } else 
          if (str_detect(el, "10\\^\\(")){
          el <- gsub("10\\^\\(", "", el)
          el <- gsub("\\)", "", el)
          par_value <- paste0(par_value, el)
          par_scale <- "log10"
        } else par_value <- paste0(par_value, el)
      }
      
    } else par_scale <- "lin"
    
    trafoDF$name[i] <- par_name
    trafoDF$trafo[i] <- par_value
    trafoDF$scale[i] <- par_scale
  }
  trafoDF <- as.data.table(trafoDF)
  gscale <- setdiff(unique(trafoDF$scale), "lin")[1]
  trafoDF <- trafoDF[, list(parameterId = name, parameterFormula = trafo, trafoType = scale)]
  trafoDF <- trafoDF[parameterId!=parameterFormula]
  attr(trafoDF, "generalScale") <- gscale
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
