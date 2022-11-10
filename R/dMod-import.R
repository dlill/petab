#' Import an SBML model and corresponding PEtab objects
#'
#' @description This function imports an SBML model and corresponding PEtab files, e.g. from the Benchmark collection.
#'
#' @param modelname name of folder containing all PEtab files of the model to be imported. NULL if file paths are defined separately (see below).
#' @param path2model path to model folder
#' @param TestCases TRUE to load feature test cases
#' @param path2TestCases path to feature test case folder
#' @param compile if FALSE, g, ODEmodel and err are loaded from .RData (if present) and compilation time is saved
#' @param SBML_file SBML model as .xml
#' @param observable_file PEtab observable file as .tsv
#' @param condition_file PEtab condition file as .tsv
#' @param data_file PEtab data file as .tsv
#' @param parameter_file PEtab parameter file as .tsv
#'
#' @details Objects such as model equations, parameters or data are automatically assigned to the following standard variables and written to your current working directory (via <<-):
#' reactions, observables, errors, g, x, p0, err, obj, mydata, ODEmodel, condition.grid, trafoL, pouter, times.
#' Compiled objects (g, ODEmodel and err) are saved in .RData.
#'
#' @return name of imported model
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#'
#' @export
#'
importPEtabSBML <- function(modelname = "Boehm_JProteomeRes2014",
                            path2model = "BenchmarkModels/",
                            testCases = FALSE,
                            path2TestCases = "PEtabTests/",
                            compile = TRUE,
                            SBML_file = NULL,
                            observable_file = NULL,
                            condition_file = NULL,
                            data_file = NULL,
                            parameter_file = NULL
)
{
  ## load required packages
  
  require(libSBML)
  require(dplyr)
  require(rlang)
  require(stringr)
  require(deSolve)
  require(ggplot2)
  require(trust)
  ## Define path to SBML and PEtab files --------------------
  
  starttime <- Sys.time()
  if(testCases == FALSE){
    if(is.null(SBML_file))       SBML_file       <- paste0(path2model, modelname, "/model_", modelname, ".xml")
    if(is.null(observable_file)) observable_file <- paste0(path2model, modelname, "/observables_", modelname, ".tsv")
    if(is.null(condition_file))  condition_file  <- paste0(path2model, modelname, "/experimentalCondition_", modelname, ".tsv")
    if(is.null(data_file))       data_file       <- paste0(path2model, modelname, "/measurementData_", modelname, ".tsv")
    if(is.null(parameter_file))  parameter_file  <- paste0(path2model, modelname, "/parameters_", modelname, ".tsv")
  } else{
    SBML_file       <- paste0(path2TestCases, modelname, "/_model.xml")
    observable_file <- paste0(path2TestCases, modelname, "/_observables.tsv")
    condition_file  <- paste0(path2TestCases, modelname, "/_conditions.tsv")
    data_file       <- paste0(path2TestCases, modelname, "/_measurements.tsv")
    parameter_file  <- paste0(path2TestCases, modelname, "/_parameters.tsv")
  }
  mywd <- getwd()
  if(!file.exists(SBML_file)){       cat(paste0("The file ",mywd,SBML_file, " does not exist. Please check spelling or provide the file name via the SBML_file argument.")); return(NULL)}
  if(!file.exists(observable_file)){ cat(paste0("The file ",mywd,observable_file, " does not exist. Please check spelling or provide the file name via the observable_file argument.")); return(NULL)}
  if(!file.exists(condition_file)){  cat(paste0("The file ",mywd,condition_file, " does not exist. Please check spelling or provide the file name via the condition_file argument.")); return(NULL)}
  if(!file.exists(data_file)){       cat(paste0("The file ",mywd,data_file, " does not exist. Please check spelling or provide the file name via the data_file argument.")); return(NULL)}
  if(!file.exists(parameter_file)){  cat(paste0("The file ",mywd,parameter_file, " does not exist. Please check spelling or provide the file name via the parameter_file argument.")); return(NULL)}
  
  if(is.null(modelname)) modelname <- "mymodel"
  ## Load shared objects --------------------
  
  dir.create(paste0(mywd,"/CompiledObjects/"), showWarnings = FALSE)
  setwd(paste0(mywd,"/CompiledObjects/"))
  files_loaded <- FALSE
  if(compile == FALSE & file.exists(paste0(modelname,".RData"))){
    load(paste0(modelname,".RData"))
    files_loaded <- TRUE
  }
  setwd(mywd)
  
  ## Model Definition - Equations --------------------
  
  cat("Reading SBML file ...\n")
  mylist <- getReactionsSBML(SBML_file, condition_file)
  myreactions <- mylist$reactions
  myreactions_orig <- mylist$reactions_orig
  myevents <- mylist$events
  mypreeqEvents <- mylist$preeqEvents
  mystates <- mylist$mystates
  reactions <<- myreactions
  
  
  ## g - Model Definition - Observables --------------------
  
  cat("Reading observables ...\n")
  myobservables <- getObservablesSBML(observable_file)
  observables <<- myobservables
  
  
  cat("Compiling observable function ...\n")
  if(!files_loaded) {
    setwd(paste0(mywd,"/CompiledObjects/"))
    myg <- Y(myobservables, myreactions, compile=TRUE, modelname=paste0("g_",modelname))
    setwd(mywd)
  }
  g <<- myg
  
  
  ## data - Get Data ------------
  
  cat("Reading data file ...\n")
  mydataSBML <- getDataPEtabSBML(data_file, observable_file)
  mydata <- mydataSBML$data
  mydata <<- mydata
  
  
  ## x - Model Generation ---------------------
  
  cat("Compiling ODE model ...\n")
  
  if(!files_loaded) {
    setwd(paste0(mywd,"/CompiledObjects/"))
    myodemodel <- odemodel(myreactions, forcings = NULL, events = myevents, fixed=NULL, modelname = paste0("odemodel_", modelname), jacobian = "inz.lsodes", compile = TRUE)
    setwd(mywd)
  }
  ODEmodel <<- myodemodel
  
  
  ## err -  Check and define error model ------------
  
  cat("Check and compile error model ...\n")
  myerrors <- mydataSBML$errors
  errors <<- myerrors
  
  myerr <- NULL
  if(!files_loaded) {
    if(!is_empty(getSymbols(myerrors))){
      setwd(paste0(mywd,"/CompiledObjects/"))
      myerr <- Y(myerrors, f = c(as.eqnvec(myreactions), myobservables), states = names(myobservables), attach.input = FALSE, compile = TRUE, modelname = paste0("errfn_", modelname))
      setwd(mywd)
    }
  }
  err <<- myerr
  
  ## Define constraints, initials, parameters and compartments --------------
  
  cat("Reading parameters and initials ...\n")
  myparameters <- getParametersSBML(parameter_file, SBML_file)
  myconstraints <- myparameters$constraints
  SBMLfixedpars <- myparameters$SBMLfixedpars
  myfit_values <- myparameters$pouter
  myinitialsSBML <- getInitialsSBML(SBML_file, condition_file)
  mycompartments <- myinitialsSBML$compartments
  myinitials <- myinitialsSBML$initials
  
  ## p - Parameter transformations -----------
  
  # Generate condition.grid
  grid <- getConditionsSBML(conditions = condition_file, data = data_file, observables_file = observable_file)
  mypreeqCons <- grid$preeqCons
  mycondition.grid <- grid$condition_grid
  
  if(!is.null(SBMLfixedpars)){
    for (i in 1:length(SBMLfixedpars)) {
      if(!names(SBMLfixedpars)[i] %in% names(mycondition.grid))  mycondition.grid[names(SBMLfixedpars)[i]] <- SBMLfixedpars[i]
    }
  }
  condi_pars <- names(mycondition.grid)[!names(mycondition.grid) %in% c("conditionName","conditionId")]
  condition.grid <<- mycondition.grid
  
  cat("Generate parameter transformations ...\n")
  myinnerpars <- unique(c(getParameters(myodemodel), getParameters(myg), getSymbols(myerrors)))
  names(myinnerpars) <- myinnerpars
  trafo <- as.eqnvec(myinnerpars, names = myinnerpars)
  trafo <- replaceSymbols(names(mycompartments), mycompartments, trafo)
  # only overwrite intial if it's not defined in condition.grid
  for (i in 1:length(myinitials)) {
    if(!names(myinitials)[i] %in% names(mycondition.grid)){
      trafo <- replaceSymbols(names(myinitials)[i], myinitials[i], trafo)
    }
  }
  trafo <- replaceSymbols(names(myconstraints), myconstraints, trafo)
  
  # branch trafo for different conditions
  mytrafoL <- branch(trafo, table=mycondition.grid)
  # set preequilibration event initials to corresponding values
  if(!is.null(mypreeqEvents)){
    mypreeqEvents2replace <- filter(mypreeqEvents, !var%in%mystates)
    mytrafoL <- repar("x~y", mytrafoL , x = unique(mypreeqEvents2replace$var), y = attr(mypreeqEvents2replace, "initials"))
  }
  # set remaining event initial to 0
  mytrafoL <- repar("x~0", mytrafoL , x = setdiff(unique(myevents$var), unique(mypreeqEvents$var)))
  
  # condition-specific assignment of parameters from condition grid
  if(length(condi_pars) > 0){
    for (j in 1:length(names(mytrafoL))) {
      for (i in 1:length(condi_pars)) {
        mytrafoL[[j]] <- repar(x~y, mytrafoL[[j]], x=condi_pars[i], y=mycondition.grid[j,condi_pars[i]])
      }
    }
  }
  
  
  # transform parameters according to scale defined in the parameter PEtab file
  parscales <- attr(myfit_values,"parscale")
  mynames <- names(parscales)
  for(i in 1:length(parscales)){
    par <- parscales[i]
    par[par=="lin"] <- ""
    par[par=="log10"] <- "10**"
    par[par=="log"] <- "exp"
    parameter <- mynames[i]
    mytrafoL <- repar(paste0("x~",par,"(x)"), mytrafoL, x = parameter)
  }
  trafoL <<- mytrafoL
  
  
  ## Specify prediction functions ------
  
  # Get numeric steady state for preequilibration conditions
  if(!is.null(mypreeqCons)){
    myf <- as.eqnvec(myreactions_orig)[mystates]
    cq <- conservedQuantities(myreactions_orig$smatrix)
    if(!is.null(cq)){
      for(i in 1:nrow(cq)){
        myf[getSymbols(cq)[1]] <- paste0(as.character(conservedQuantities(myreactions_orig$smatrix)[1,]),"-1")
      }
    }
    setwd(paste0(mywd,"/CompiledObjects/"))
    pSS <- P(myf, condition = "c0", method = "implicit", compile = TRUE, modelname = paste0("preeq_", modelname))
    setwd(mywd)
  } else pSS <- NULL
  
  cat("Generate prediction function ...\n")
  tolerances <- 1e-7
  myp0 <- myx <- NULL
  for (C in names(mytrafoL)) {
    if(C%in%mypreeqCons){
      myp0 <- myp0 + pSS*P(mytrafoL[[C]], condition = C)
    } else {
      myp0 <- myp0 + P(mytrafoL[[C]], condition = C)
    }
    myx <- myx + Xs(myodemodel, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
                    optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
                    condition = C)
  }
  
  p0 <<- myp0
  x <<- myx
  
  ## Generate objective function and initial parameter set -------
  
  myouterpars <- getParameters(myp0)
  mypouter <- structure(rep(NA,length(myouterpars)), names = myouterpars)
  
  common <- intersect(names(mypouter),names(myfit_values))
  mypouter[common] <- myfit_values[common]
  attr(mypouter, "parscales") <- parscales[which(names(myfit_values)%in%names(mypouter))]
  attr(mypouter, "lowerBound") <- attr(myfit_values,"lowerBound")[which(names(myfit_values)%in%names(mypouter))]
  attr(mypouter, "upperBound") <- attr(myfit_values,"upperBound")[which(names(myfit_values)%in%names(mypouter))]
  pouter <<- mypouter
  
  
  ## Define objective function -------
  
  cat("Generate objective function ...\n")
  if(!is.null(myerrors)){
    myobj <- normL2(mydata, myg*myx*myp0, errmodel = myerr) #+ constraintL2(prior, sigma=16)
  } else myobj <- normL2(mydata, myg*myx*myp0)
  obj <<- myobj
  
  mytimes <- seq(0,max(do.call(c, lapply(1:length(mydata), function(i) max(mydata[[i]]$time)))), len=501)
  times <<- mytimes
  
  if(!files_loaded){
    setwd(paste0(mywd,"/CompiledObjects/"))
    save(list = c("myg","myodemodel","myerr"),file = paste0(modelname,".RData"))
    setwd(mywd)
  }
  
  model_name <<- modelname
  
  endtime <- Sys.time()
  mytimediff <- as.numeric(difftime(endtime, starttime, unit="secs"))
  if(mytimediff > 3600) cat(paste0(modelname, " imported in ",as.character(format(as.numeric(difftime(endtime, starttime, unit="hours")), digits=3)), " hours.\n")) else
    if(mytimediff > 60) cat(paste0(modelname, " imported in ",as.character(format(as.numeric(difftime(endtime, starttime, unit="mins")), digits=3)), " minutes.\n")) else
      cat(paste0(modelname, " imported in ",as.character(format(as.numeric(difftime(endtime, starttime, unit="secs")), digits=3)), " seconds.\n"))
  # return(modelname)
  
  
  
  # .. Collect list -----
  symbolicEquations <- list(
    reactions   = myreactions,
    events      = myevents,
    observables = myobservables,
    errors      = myerrors,
    trafo       = trafoL)
  fns <- list(
    g = myg,
    x = myx,
    p0 = myp0
  )
  prd <- (myg*myx*myp0)
  obj_data <- normL2(mydata, prd, errmodel = myerr,
                     times = seq(0,max(as.data.frame(mydata)$time), len=501))
  
  filenameParts = list(modelname = modelname, .currentFolder = mywd,
                       .compiledFolder = "CompiledObjects",type = "classic")
  pe <- readPetab(filename = file.path(path2model, modelname))
  # .. Collect final list -----
  pd <- list(
    # petab
    pe                 = pe,
    # Basic dMod elements
    dModAtoms          = list(
      # [ ] add events!
      symbolicEquations  = symbolicEquations,
      odemodel           = myodemodel,
      data               = mydata,
      gridlist           = NULL,
      e                  = myerr,
      fns                = fns
    ),
    # other components: Dump your stuff here
    filenameParts = filenameParts,
    # Parameters + Time
    pars               = dMod::unclass_parvec(myfit_values),
    times              = dMod::predtimes(pe$measurementData$time, Nobjtimes = 200),
    # High-level functions
    prd                = prd,
    obj_data           = obj_data
  )
  
  
  
  pd
}



#' Fit a model imported via importPEtabSBML
#'
#' @description A wrapper function to use \link{mstrust} with imported PEtabSBML models. Some reasonable standard arguments for mstrust are used. Results of mstrust are written to Results folder.
#'
#' @param objfun Objective function to be minimized as created by \link{importPEtabSBML}.
#' @param nrfits numeric, Number of fits to be performed
#' @param nrcores numeric, Number of cores to be used
#' @param useBounds  boolean, if TRUE, parameter bounds are taken as provided in PEtab format, if FALSE no parameter bounds are applied
#'
#' @return parframe with the parameter estimated of the multi-start optimization
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#'
#' @export
#' @importFrom dMod mstrust msParframe
#'
fitModelPEtabSBML <- function(objfun=obj, nrfits=4, nrcores=4, useBounds=TRUE){
  prior <- structure(rep(0,length(pouter)))
  names(prior) <- names(pouter)
  mywd <- getwd()
  dir.create(paste0(mywd,"/Test/mstrust/"), showWarnings = FALSE)
  if(useBounds) out <- dMod::mstrust(objfun=objfun, center=dMod::msParframe(prior, n = nrfits+1, seed=47)[-1,], studyname=model_name, rinit = 0.1, rmax = 10,
                                     fits = nrfits, cores = nrcores, samplefun = "rnorm", resultPath = "Test/mstrust/",
                                     parlower = attr(pouter, "lowerBound"), parupper=attr(pouter, "upperBound"),
                                     stats = FALSE, narrowing = NULL, iterlim=400, sd = 3)
  else out <- dMod::mstrust(objfun=objfun, center=dMod::msParframe(prior, n = nrfits, seed=47), studyname=model_name, rinit = 0.1, rmax = 10,
                            fits = nrfits, cores = nrcores, samplefun = "rnorm", resultPath = "Test/mstrust/",
                            stats = FALSE, narrowing = NULL, iterlim=400, sd = 3)
  if(any(lapply(out, function(fgh) fgh$converged)==TRUE)) return(as.parframe(out)) else {cat("No fit converged."); return(NULL)}
}



#' Test PEtabSBML import
#'
#' @description This function imports, evaluates and tests the PEtab model.
#'
#' @param models model names to test
#'
#' @return evaluation data.frame
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#'
#' @export
#'
testPEtabSBML <- function(models = c(
  #"Boehm_JProteomeRes2014"
  # "Fujita_SciSignal2010",
  # "Borghans_BiophysChem1997",
  # "Elowitz_Nature2000",
  # "Sneyd_PNAS2002",
  # "Crauste_CellSystems2017",
  # "Schwen_PONE2014",
  # "Raia_CancerResearch2011",
  # "Zheng_PNAS2012",
  # "Beer_MolBioSystems2014",
  # "Brannmark_JBC2010",
  # "Bruno_JExpBio2016",
  # "Chen_MSB2009",
  # "Fiedler_BMC2016",
  # "Weber_BMC2015",
  # "Swameye_PNAS2003"
  # "Bachmann_MSB2011"
  # "Lucarelli_CellSystems2018",
  "0001",
  "0002",
  "0003",
  "0004",
  "0005",
  "0006",
  "0007",
  "0008",
  "0009",
  "0010",
  "0011",
  "0012",
  "0013",
  "0014",
  "0015",
  "0016"
), testFit = TRUE, timelimit = 5000, testCases = TRUE) {
  try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf) {
    y <- try(
      {
        setTimeLimit(cpu, elapsed)
        expr
      },
      silent = TRUE
    )
    if (inherits(y, "try-error")) NULL else y
  }
  require(crayon)
  require(trust)
  require(deSolve)
  require(ggplot2)
  cat(green("Start test function...\n"))
  mywd <- getwd()
  teststarttime <- Sys.time()
  output <- NULL
  predictions <- NULL
  for (model in models) {
    setwd(mywd)
    importtest <- F
    plottest <- F
    bestfit <- NA
    cat(blue(paste0("Testing ", model, "\n")))
    fgh <- try_with_time_limit(
      {
        test <- try(importPEtabSBML(model, compile = T, testCases = testCases), silent = T)
        if (inherits(test, "try-error")) "import error" else test
      },
      timelimit
    )
    if (fgh == "import error") {
      cat(yellow("Import error or time limit exceeded for", model, "\n\n\n"))
      output <- rbind(output, data.frame(
        modelname = model, import = importtest,
        fitting_time = NA, plot = plottest, chi2 = NA, LL = NA, bestfit = NA, difference = NA
      ))
    } else {
      importtest <- T
      testobj <- try(obj(pouter))
      if (inherits(testobj, "try-error")) {
        cat(red("Warning: Error in calculation of objective function.\n"))
        output <- rbind(output, data.frame(
          modelname = model, import = importtest,
          fitting_time = NA, plot = plottest, chi2 = NA, LL = NA, bestfit = NA, difference = NA
        ))
      } else {
        if (is.numeric(testobj$value)) {
          cat(green("Calculation of objective function successful.\n"))
          if (testCases){
            # calculate predictions for trajectory comparison
            mysimulations <- read.csv(paste0("PEtabTests/", model, "/_simulations.tsv"), sep = "\t")
            simu_time <- unique(mysimulations$time)
            prediction <- (g*x*p0)(simu_time, pouter)
            predictions <- rbind(predictions, data.frame(
              modelname = model, pred = prediction, obs.transformation = NA
            ))
            # append observable scale to predictions
            for (i in 1:length(observables)) {
              scale <- attr(observables, "obsscales")[i]
              predictions <- predictions %>% mutate(obs.transformation = ifelse(modelname == model & pred.name == names(observables)[i], scale, obs.transformation))
            }
          }
        } else {
          cat(red("Warning: obj(pouter) is not numeric.\n"))
        }
        # objLL <- mynormL2(mydata, g * x * p0, outputLL = T)
        # testLL <- try(-0.5 * objLL(pouter)$value)
        # if (inherits(testLL, "try-error")) testLL <- NA
        if (testFit) {
          fitstarttime <- Sys.time()
          myframe <- fitModelPEtabSBML(nrfits = 20)
          fitendtime <- Sys.time()
          if (is.parframe(myframe) & !is.null(myframe)) {
            if (is.numeric(obj(myframe[1, ])$value)) {
              cat(green("Fit test successful.\n"))
              bestfit <- obj(myframe[1, ])$value
            } else {
              cat(red("Warning: obj(myframe) is not numeric.\n"))
            }
          } else {
            cat(red("Warning: Fit test not successful..\n"))
          }
          mytimediff <- as.numeric(difftime(fitendtime, fitstarttime, unit = "secs"))
          if (mytimediff > 3600) {
            cat(green(paste0("Fitting done in ", as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit = "hours")), digits = 3)), " hours.\n")))
          } else
            if (mytimediff > 60) {
              cat(green(paste0("Fitting done in ", as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit = "mins")), digits = 3)), " minutes.\n")))
            } else {
              cat(green(paste0("Fitting done in ", as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit = "secs")), digits = 3)), " seconds.\n")))
            }
        }
        pdf(file = paste0("Test/", model, "_plotAll.pdf"))
        plotPEtabSBML()
        dev.off()
        pdf(file = paste0("Test/", model, "_plotTargetsObserved.pdf"))
        plotPEtabSBML(name %in% names(observables))
        dev.off()
        pdf(file = paste0("Test/", model, "_plotConditionsObserved.pdf"))
        plotPEtabSBML(condition %in% names(mydata))
        dev.off()
        plottest <- T
        cat(green("Import and plot test for ", fgh, " successful!\n\n\n"))
        
        output <- rbind(output, data.frame(
          modelname = model, import = importtest,
          # fitting_time = format(as.numeric(difftime(fitendtime, fitstarttime, unit = "mins")), digits = 3),
          plot = plottest, chi2 = attr(testobj,"chisquare"), LL = -0.5*testobj$value
          # , bestfit = bestfit, difference = bestfit - testobj$value
        ))
      }
    }
    if (testCases){
      sharedObjects <- paste0("CompiledObjects/",
                              c(paste0("g_",model,".so"),
                                paste0("g_",model,"_deriv.so"),
                                paste0("odemodel_",model,".so"),
                                paste0("odemodel_",model,"_s.so"),
                                paste0("errfn_",model,".so"),
                                paste0("errfn_",model,"_s.so")))
      for (file in sharedObjects) if(file.exists(file)) try(dyn.unload(file))
    }
  }
  if (testCases) {
    simu_output <- output[1]
    output <- cbind(output, chi2_sol = NA, tol_chi2_sol = NA, LL_sol = NA, tol_LL_sol = NA)
    for (model in models) {
      mysolution <- read_yaml(paste0("PEtabTests/", model, "/_", model, "_solution.yaml"))
      output[which(output$modelname == model), "chi2_sol"] <- mysolution$chi2
      output[which(output$modelname == model), "tol_chi2_sol"] <- mysolution$tol_chi2
      output[which(output$modelname == model), "LL_sol"] <- mysolution$llh
      output[which(output$modelname == model), "tol_LL_sol"] <- mysolution$tol_llh
      
      # extract simulation values
      simu_output[which(simu_output$modelname == model), "tol_simus_sol"] <- mysolution$tol_simulations
      mysimulations <- read.csv(paste0("PEtabTests/", model, "/_simulations.tsv"), sep = "\t")
      simu_prediction <- subset(predictions, modelname == model)
      
      # iterate through simulation points
      for (nrow in 1:nrow(mysimulations)) {
        simu_row <- mysimulations[nrow,]
        simu_time <- simu_row$time
        simu_obs <- simu_row$observableId %>% as.character()
        simu_condi <- simu_row$simulationConditionId %>% as.character()
        simu_obspars <- simu_row$observableParameters
        
        if(!is.null(simu_obspars) & length(unique(simu_prediction$pred.condition)) > 1){
          simu_condi <- paste0(simu_condi, "_", simu_obspars)
        }
        pred_row <- subset(simu_prediction, pred.time == simu_time & pred.name == simu_obs & pred.condition == simu_condi)
        # retransform simulation value according to observable transformation
        if(pred_row$obs.transformation == "log10") pred_row$pred.value <- 10**pred_row$pred.value
        if(pred_row$obs.transformation == "log") pred_row$pred.value <- exp(pred_row$pred.value)
        simu_value <- pred_row$pred.value
        
        simu_output[which(simu_output$modelname == model), paste0("simu_", nrow)] <- simu_value
        simu_output[which(simu_output$modelname == model), paste0("simu_", nrow,"_sol")] <- mysimulations$simulation[nrow]
      }
    }
  }
  
  testendtime <- Sys.time()
  mytimediff <- as.numeric(difftime(testendtime, teststarttime, unit = "secs"))
  if (mytimediff > 3600) {
    cat(green(paste0("Test done in ", as.character(format(as.numeric(difftime(testendtime, teststarttime, unit = "hours")), digits = 3)), " hours.\n")))
  } else
    if (mytimediff > 60) {
      cat(green(paste0("Test done in ", as.character(format(as.numeric(difftime(testendtime, teststarttime, unit = "mins")), digits = 3)), " minutes.\n")))
    } else {
      cat(green(paste0("Test done in ", as.character(format(as.numeric(difftime(testendtime, teststarttime, unit = "secs")), digits = 3)), " seconds.\n")))
    }
  
  if (testCases) {
    
    # check simulations
    for (model in models) {
      correctORnot <- NULL
      modelrow <- subset(simu_output, modelname == model)
      modelrow_woNA <- modelrow[colSums(!is.na(modelrow)) > 0]
      for (ncol in seq(3,(ncol(modelrow_woNA)),2)) {
        simuCompare <- abs(modelrow_woNA[[ncol]]-modelrow_woNA[[ncol+1]]) < modelrow_woNA$tol_simus_sol
        correctORnot <- c(correctORnot, simuCompare)
      }
      if(length(unique(correctORnot)) == 1){
        SimuPassed <- unique(correctORnot)
      } else SimuPassed <- FALSE
      simu_output[which(simu_output$modelname == model), "Passed"] <- SimuPassed
    }
    
    output <- cbind(output,
                    X2Passed = (abs(output$chi2 - output$chi2_sol) < output$tol_chi2_sol),
                    LLPassed = (abs(output$LL - output$LL_sol) < output$tol_LL_sol)
    )
  }
  
  if (!testCases) simu_output <- NULL
  return(list(output = output,simu_output = simu_output))
}



#' Import condition.grid from PEtab
#'
#' @description This function imports the experimental conditions from the PEtab condition file as a gondition.grid.
#'
#' @param conditions PEtab condition file as .tsv
#'
#' @return condition.grid as data frame.
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#'
#' @export
#'
#' @importFrom dplyr inner_join filter mutate
#' @importFrom stringr str_detect
getConditionsSBML <- function(conditions,data, observables_file, FLAGnormalImport = TRUE){
  condition.grid_orig <- read.csv(file = conditions, sep = "\t")
  mydata              <- read.csv(file = data, sep = "\t")
  myobservables       <- read.csv(file = observables_file, sep = "\t")
  
  # handle preequilibration conditions
  myCons <- condition.grid_orig$conditionId
  mypreeqCons <- NULL
  for (con in myCons){if(paste0("preeq_", con)%in%myCons) mypreeqCons <- c(mypreeqCons, con)}
  condition.grid_orig <- dplyr::filter(condition.grid_orig, !conditionId%in%mypreeqCons)
  condition.grid_orig$conditionId <- sub("preeq_", "", condition.grid_orig$conditionId)
  
  # check which conditions are observed
  condis_obs <- mydata$simulationConditionId %>% unique()
  # check which observables exist
  observables <- mydata$observableId %>% unique()
  
  # replace "" by NA
  if(!is.null(mydata$observableParameters)){
    mydata$observableParameters <- mydata$observableParameters %>% as.character()
    mydata <- mydata %>% dplyr::mutate(observableParameters = ifelse(observableParameters == "",NA,observableParameters))
  }
  if(!is.null(mydata$noiseParameters)){
    mydata$noiseParameters <- mydata$noiseParameters %>% as.character()
    mydata <- mydata %>% dplyr::mutate(noiseParameters = ifelse(noiseParameters == "",NA,noiseParameters))
  }
  # generate columns for observableParameters
  if(!is.numeric(mydata$observableParameters) & !is.null(mydata$observableParameters)){
    condition.grid_obs <- data.frame(conditionId = condis_obs)
    for (obs in observables){
      data_obs <- subset(mydata, observableId == obs)
      for (condition in condis_obs){
        if(condition %in% data_obs$simulationConditionId){
          row_pars <- NULL
          obs_par <- subset(data_obs, simulationConditionId == condition)$observableParameters %>% unique() %>% as.character()
          if(length(obs_par)==1){
            if(!is.na(obs_par)){
              # one or more observable parameters?
              if(stringr::str_detect(obs_par,";")){
                myobspars <- strsplit(obs_par,";")[[1]]
                for(i in 1:length(myobspars)) {
                  row_pars <- c(row_pars, myobspars[i])
                }
              } else row_pars <- c(row_pars, obs_par)
            }
            if(!is.null(row_pars)) for (par in 1:length(row_pars)) {
              col_name <- paste0("observableParameter",par,"_",obs)
              condition.grid_obs[which(condition.grid_obs$conditionId==condition),col_name] <- row_pars[par]
            }
          } else {
            col_name <- paste0("observableParameter1_",obs)
            add <- NULL
            for(j in 2:length(obs_par)){
              add <- rbind(add, subset(condition.grid_obs, conditionId==condition))
            }
            condition.grid_obs <- rbind(condition.grid_obs, add)
            condition.grid_obs[which(condition.grid_obs$conditionId==condition),col_name] <- obs_par
            condition.grid_obs$conditionId <- as.character(condition.grid_obs$conditionId)
            condition.grid_obs$conditionId[which(condition.grid_obs$conditionId==condition)] <- paste0(condition,"_", obs_par)
            
            condition.grid_orig <- rbind(condition.grid_orig, add)
            condition.grid_orig$conditionId <- as.character(condition.grid_orig$conditionId)
            condition.grid_orig$conditionId[which(condition.grid_orig$conditionId==condition)] <- paste0(condition,"_", obs_par)
          }
        }
      }
    }
    mycondition.grid <- suppressWarnings(dplyr::inner_join(condition.grid_orig,condition.grid_obs, by = "conditionId"))
    # avoid warning if not all conditions are observed
  } else mycondition.grid <- condition.grid_orig
  
  
  # Error parameters: Need to cover the four cases
  # 1. err = noiseParameter1_obs, noiseParameter1_obs = c(1)              => No columns in condition.grid
  # 2. err = noiseParameter1_obs, noiseParameter1_obs = c("sigma1")       => Need Column in condition.grid
  # 3. err = noiseParameter1_obs * obs, noiseParameter1_obs = c(1)        => Need Column in condition.grid
  # 4. err = noiseParameter1_obs * obs, noiseParameter1_obs = c("sigma1") => Need Column in condition.grid
  # Check if noise parameters need to be generated or the error model is simply "sigma"
  err_simple <- !is.na(as.numeric(myobservables$noiseFormula)) | myobservables$noiseFormula == paste0("noiseParameter1_", myobservables$observableId)
  # Check if there are any symbolic error parameters
  err_symbols <- getSymbols(mydata$noiseParameters)
  # generate columns for noiseParameters
  if((FLAGnormalImport && is.character(mydata$noiseParameters)) | (any(!err_simple) | length(err_symbols) > 0) & !is.null(mydata$noiseParameters))
  {
    if(exists("mycondition.grid")) {condition.grid_orig <- mycondition.grid}
    condition.grid_noise <- data.frame(conditionId = condis_obs)
    for (obs in observables[!err_simple | FLAGnormalImport]) # `| FLAGnormalImport` for backwards compatibility
    {
      data_obs <- subset(mydata, observableId == obs)
      for (condition in condis_obs)
      {
        if(condition %in% data_obs$simulationConditionId){
          row_pars <- NULL
          noise_par <- subset(data_obs, simulationConditionId == condition)$noiseParameters %>% unique() %>% as.character()
          if(!is.na(noise_par)){
            # one or more observable parameters?
            if(stringr::str_detect(noise_par,";")){
              myobspars <- strsplit(noise_par,";")[[1]]
              for(i in 1:length(myobspars)) {
                row_pars <- c(row_pars, myobspars[i])
              }
            } else row_pars <- c(row_pars, noise_par)
          }
          if(!is.null(row_pars)) for (par in 1:length(row_pars)) {
            col_name <- paste0("noiseParameter",par,"_",obs)
            condition.grid_noise[which(condition.grid_noise$conditionId==condition),col_name] <- row_pars[par]
          }
        }
      }
      
    }
    mycondition.grid <- suppressWarnings(dplyr::inner_join(condition.grid_orig,condition.grid_noise, by = "conditionId"))
    # avoid warning if not all conditions are observed
  }
  
  if(!exists("mycondition.grid")) mycondition.grid <- condition.grid_orig
  rownames(mycondition.grid) <- mycondition.grid$conditionId
  # mycondition.grid$conditionId <- NULL ## we need this column in cases with just one condition!
  
  # check if all conditions are observed
  if(nrow(mycondition.grid) < nrow(condition.grid_orig)) print("There exist non-observed conditions!")
  
  for(i in 1:nrow(mycondition.grid)){
    for(j in 1:ncol(mycondition.grid)){
      if(is.na(mycondition.grid[i,j])) mycondition.grid[i,j] <- "1"
    }
  }
  
  return(list(condition_grid=mycondition.grid, preeqCons=mypreeqCons))
}



#' Import initials from SBML.
#'
#' @description This function imports initial values or equations describing the same from SBML and writes them in a named vector.
#'
#' @param model SBML file as .xml
#'
#' @return Named vector of initials.
#'
#' @author Marcus Rosenblatt, Svenja Kemmer and Frank Bergmann
#' importFrom libSBML readSBML formulaToL3String
#' @export
#'
getInitialsSBML <- function(model, conditions){
  
  condition.grid <- read.csv(file = conditions, sep = "\t")
  model = readSBML(model)$getModel()
  # model = libSBML::readSBML(model)$getModel()
  initials <- NULL
  species <- NULL
  
  # now we can go through all species
  for ( i in 0:(model$getNumSpecies()-1) ){
    current <- model$getSpecies(i)
    
    # now the species can have several cases that determine
    # their initial value
    
    # CASE 1: it could be that the species is fully determined by an assignment rule
    # (that apply at all times), so we have to check rules first
    rule <- model$getRule(current$getId())
    if (!is.null(rule))
    {
      # ok there is a rule for this species so lets figure out what its type
      # is as that determines whether it applies at t0
      rule_type <- rule$getTypeCode()
      type_name <- SBMLTypeCode_toString(rule_type, 'core')
      
      # CASE 1.1: Assignment rule
      if (type_name == "AssignmentRule"){
        # the initial value is determined by the formula
        math <- rule$getMath()
        if (!is.null(math)){
          formula <- formulaToL3String(math)
          # print(paste('Species: ', current$getId(), ' is determined at all times by formula: ', formula))
          initials <- c(initials,formula)
          species <- c(species,current$getId())
          
          # no need to look at other values so continue with next species
          next
        }
      }
      
      # CASE 1.2: Rate rule
      if (type_name == "RateRule")
      {
        math <- rule$getMath()
        if (!is.null(math))
        {
          formula <- formulaToL3String(math)
          # print(paste('Species: ', current$getId(), ' has an ode rule with formula: ', formula))
          initials <- c(initials,formula)
          species <- c(species,current$getId())
          
          # there is an ODE attached to the species, its initial value is needed
        }
      }
    }
    
    
    # CASE 2 (or subcase of 1?): Initial assignment
    
    # it could have an initial assignment
    ia <- model$getInitialAssignment(current$getId())
    if (!is.null(ia))
    {
      math <- ia$getMath()
      if (!is.null(math))
      {
        formula <- formulaToL3String(math)
        # formula <- libSBML::formulaToL3String(math)
        # print(paste("Species: ", current$getId(), " has an initial assignment with formula: ", formula))
        initials <- c(initials,formula)
        species <- c(species,current$getId())
        
        # as soon as you have that formula, no initial concentration / amount applies
        # so we don't have to look at anything else for this species
        next
      }
    }
    
    # CASE 3 (or subcase of 1?): Inital amount
    # it could have an initial amount
    if (current$isSetInitialAmount())
    {
      # print (paste("Species: ", current$getId(), "has initial amount: ", current$getInitialAmount()))
      initials <- c(initials,current$getInitialAmount())
      species <- c(species,current$getId())
    }
    
    # CASE 4 (or subcase of 1?): Inital Concentration
    # it could have an initial concentration
    if (current$isSetInitialConcentration())
    {
      # print (paste("Species: ", current$getId(), "has initial concentration: ", current$getInitialConcentration()))
      initials <- c(initials,current$getInitialConcentration())
      species <- c(species,current$getId())
    }
  }
  names(initials) <- species
  
  
  ## extract compartments
  # check if compartments exist
  if(model$getNumCompartments()>0) {
    
    #initialize vectors
    comp_name <- NULL
    comp_size <- NULL
    
    for ( i in 0:(model$getNumCompartments()-1) ) {
      
      # get compartment name and size
      
      compartmentId <- model$getCompartment(i)$getId()
      if(compartmentId %in% names(condition.grid)){
        size <- condition.grid[compartmentId] %>% as.numeric()
      } else size <- model$getCompartment(i)$getSize()
      
      comp_name <- c(comp_name,compartmentId)
      comp_size <- c(comp_size,size)
      
    }
    compartments <- comp_size
    names(compartments) <- comp_name
  }
  
  # Get compartments
  if(model$getNumCompartments()>0) {
    
    #initialize vectors
    comp_name <- NULL
    comp_size <- NULL
    
    for ( i in 0:(model$getNumCompartments()-1) ) {
      
      # get compartment name and size
      
      compartmentId <- model$getCompartment(i)$getId()
      if(compartmentId %in% names(condition.grid)){
        size <- condition.grid[compartmentId] %>% as.numeric()
      } else size <- model$getCompartment(i)$getSize()
      
      comp_name <- c(comp_name,compartmentId)
      comp_size <- c(comp_size,size)
      
    }
    compartments <- comp_size
    names(compartments) <- comp_name
  }  
  return(list(initials = initials, compartments = compartments))
}





#' Import observables from PEtab.
#'
#' @description This function imports observables from the PEtab observable file.
#'
#' @param observables PEtab observable file as .tsv
#'
#' @return Eqnvec of observables.
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#'
#' @export
#'
getObservablesSBML <- function(observables){
  ## Load observables
  myobs <- read.csv(file = observables, sep = "\t") %>% as.data.frame()
  obsNames <- myobs$observableId %>% as.character()
  
  # # rename observables with _obs
  # obsNames <- paste0(obsNames,"_obs")
  
  obsFormula <- myobs$observableFormula %>% as.character()
  obsFormula[which(myobs$observableTransformation=="log")] <- paste0("log(", obsFormula[which(myobs$observableTransformation=="log")], ")")
  obsFormula[which(myobs$observableTransformation=="log10")] <- paste0("log10(", obsFormula[which(myobs$observableTransformation=="log10")], ")")
  names(obsFormula) <- obsNames
  observables <- obsFormula %>% as.eqnvec()
  
  # collect observable transformations as attribute
  if(!is.null(myobs$observableTransformation)){
    obsscales <- myobs$observableTransformation %>% as.character()
  } else obsscales <- rep("lin", length(obsNames))
  
  attr(observables,"obsscales") <- obsscales
  return(observables)
}


#' Import Parameters from PEtab
#'
#' @description This function imports fixed and fitted parameters from the PEtab parameter file as named vectors.
#'
#' @param parameters PEtab parameter file as .tsv
#'
#' @return constraints and pouter as list of named vectros.
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#' @importFrom dplyr filter "%>%"
#' @export
#'
getParametersSBML <- function(parameters, model){
  mypars <- read.csv(file = parameters, sep = "\t")
  fixed <- mypars %>% dplyr::filter(estimate == 0)
  constraints <- NULL
  if(nrow(fixed)>0){
    for(i in 1:length(fixed$parameterScale)) {
      parscale <- fixed$parameterScale[i]
      par <- fixed$parameterId[i] %>% as.character()
      value <- fixed$nominalValue[i]
      if(parscale == "lin") constraints <- c(constraints, value)
      if(parscale == "log10") constraints <- c(constraints, log10(value))
      if(parscale == "log") constraints <- c(constraints, log(value))
      else paste("This type of parameterScale is not supported.")
      names(constraints)[i] <- par
    }
    parscales <- fixed$parameterScale %>% as.character()
    pars <- fixed$parameterId %>% as.character()
    names(parscales) <- pars
    attr(constraints,"parscale") <- parscales
  }
  estimated <- mypars %>% filter(estimate == 1)
  pouter <- NULL
  parlower <- NULL
  parupper <- NULL
  if(nrow(estimated)>0){
    for(i in 1:length(estimated$parameterScale)) {
      parscale <- estimated$parameterScale[i]
      par <- estimated$parameterId[i] %>% as.character()
      value <- estimated$nominalValue[i]
      lowervalue <- estimated$lowerBound[i]
      uppervalue <- estimated$upperBound[i]
      if(parscale == "lin"){
        pouter <- c(pouter, value)
        parlower <- c(parlower, lowervalue)
        parupper <- c(parupper, uppervalue)
      } else if(parscale == "log10"){
        pouter <- c(pouter, log10(value))
        parlower <- c(parlower, log10(lowervalue))
        parupper <- c(parupper, log10(uppervalue))
      } else if(parscale == "log"){
        pouter <- c(pouter, log(value))
        parlower <- c(parlower, log(lowervalue))
        parupper <- c(parupper, log(uppervalue))
      } else paste("This type of parameterScale is not supported.")
      names(pouter)[i] <- par
      names(parlower)[i] <- par
      names(parupper)[i] <- par
    }
    parscales <- estimated$parameterScale %>% as.character()
    pars <- estimated$parameterId %>% as.character()
    names(parscales) <- pars
    attr(pouter,"parscale") <- parscales
    attr(pouter,"lowerBound") <- parlower
    attr(pouter,"upperBound") <- parupper
  }
  
  # check if additional parameters exist in SBML file
  model = readSBML(model)$getModel()
  n_pars <- model$getNumParameters()
  SBMLfixedpars <- NULL
  count <- 1
  for (i in (seq_len(n_pars)-1)) {
    mypar <- model$getParameter(i)$getId()
    if(!mypar %in% names(pouter) & !mypar %in% names(constraints)){
      value <- model$getParameter(i)$getValue()
      SBMLfixedpars <- c(SBMLfixedpars, value)
      names(SBMLfixedpars)[count] <- mypar
      count <- count + 1
    }
  }
  
  return(list(constraints=constraints, pouter=pouter, SBMLfixedpars = SBMLfixedpars))
}

#' Import reactions from SBML.
#'
#' @description This function imports reactions from SBML. Reactions are written to an eqnlist object.
#' Assignment rules for input functions are substituted in the rate column. Time is introduced as an additionel state t.
#'
#' @param model SBML file as .xml
#'
#' @return Eqnlist of reactions.
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#' importFrom libSBML readSBML
#' @importFrom cOde replaceSymbols
#' @importFrom dMod addReaction
#' @export
#'
getReactionsSBML <- function(model, conditions, rates){
  # m = libSBML::readSBML(model)$getModel()
  m = readSBML(model)$getModel()
  
  # .. Initialization -----
  reactions <- NULL
  events <- NULL
  compartments <- NULL
  
  # .. import compartments -----
  N_species <- m$getNumSpecies()
  compartments <- do.call(c, lapply(0:(N_species-1), function(i){ m$getSpecies(i)$getCompartment()}))
  names_compartments <- do.call(c, lapply(0:(N_species-1), function(i){ m$getSpecies(i)$getId() }))
  names(compartments) <- names_compartments
  
  # if(unique(compartments)[1]=="default") compartments <- NULL
  
  # .. import reactions -----
  N_reactions <- m$getNumReactions()
  for (reaction in 0:(N_reactions-1)){
    reactionDetails <- sbmlImport_getReactionDetails(m, reaction, compartments)
    reactions <- do.call(dMod::addReaction, c(list(eqnlist = reactions), reactionDetails))
  }
  
  
  # .. import functions -----
  N_fundefs <- m$getNumFunctionDefinitions()
  if (N_fundefs > 0){
    for (fun in 0:(N_fundefs-1)){
      reactions$rates <- smblImport_replaceFunctions(m = m, fun = fun, rates = reactions$rates)
      # substitute m$getRule(0)$getVariable() by m$getRule(0)$getFormula()
      #print(formulaToL3String(mymath$getChild(mymath$getNumChildren()-1)))
    }
  }
  
  # .. import inputs -----
  N_rules <- m$getNumRules()
  if (N_rules > 0){
    for (rule in 0:(N_rules-1)){
      # substitute m$getRule(0)$getVariable() by m$getRule(0)$getFormula()
      reactions$rates <- cOde::replaceSymbols(m$getRule(rule)$getVariable(),
                                              paste0("(",m$getRule(rule)$getFormula(), ")"), reactions$rates)
    }
  }
  reactions$rates <- gsub(" ","",reactions$rates)
  
  # .. Events by piecewise functions -----
  # replace function based inputs by events (done in reactions)
  events <- sbmlImport_getEventStrings(reactions, events = events)
  if(!is.null(events)) for(i in 1:length(events)){
    replace <- gsub("\\(", "\\\\\\(", events[i])
    replace <- gsub("\\*", "\\\\\\*", replace)
    replace <- gsub("\\+", "\\\\\\+", replace)
    #replace <- gsub("\\)", "\\\\\\)", replace)
    reactions$rates <- gsub(replace, paste0("event", i), reactions$rates)
    reactions <- reactions %>% dMod::addReaction("", paste0("event", i), "0")
  }
  # replace mathematical expressions
  reactions$rates <- replaceSymbols(c("t", "TIME", "T"), "time", reactions$rates)
  events <- TransformEvents(events)
  
  # .. ## check for preequilibration conditions and handle them via events -----
  preeqEvents <- NULL
  myconditions <- read.csv(file = conditions, sep = "\t")
  myCons <- myconditions$conditionId
  mypreeqCons <- NULL
  attrib <- NULL
  for (con in myCons){if(paste0("preeq_", con)%in%myCons) mypreeqCons <- c(mypreeqCons, con)}
  if(!is.null(mypreeqCons)){
    for (con in mypreeqCons){
      mycongrid <- filter(myconditions, conditionId==con | conditionId==paste0("preeq_", con))
      if(ncol(mycongrid)>1){
        for(i in 2:ncol(mycongrid)){
          preeqEvents <- addEvent(preeqEvents,
                                  var=names(mycongrid)[i],
                                  time=0,
                                  value=mycongrid[[which(mycongrid$conditionId==con),i]],
                                  root=NA,
                                  method="replace")
          attrib <- c(attrib, mycongrid[[which(mycongrid$conditionId==paste0("preeq_",con)),i]])
        }
      }
    }
  }
  mystates <- reactions$states
  reactions_orig <- reactions
  attr(preeqEvents, "initials") <- attrib
  if(!is.null(preeqEvents)) for(i in 1:nrow(preeqEvents)){
    events <- rbind(events, preeqEvents[i,])
    reactions <- reactions %>% dMod::addReaction("", preeqEvents[[i,"var"]], "0")
  }
  
  
  # for(i in 1:length(reactions$rates)){
  #   reaction <- reactions$rates[i]
  #   if(stringr::str_detect(reaction, "pow")){
  #     reaction_new <- gsub("pow", "", reaction)
  #     # reaction_new <- gsub(", ", "**", reaction_new)
  #     reaction_new <- gsub(",", "**", reaction_new)
  #     reactions$rates[i] <- reaction_new
  #   } else reactions$rates[i] <- reaction
  # }
  
  mydata <- as.data.frame(reactions)
  reactions <- as.eqnlist(mydata, compartments)
  
  return(list(reactions=reactions, events=events, reactions_orig=reactions_orig, preeqEvents=preeqEvents, mystates=mystates))
}

#' Import Data from PEtab
#'
#' @description This function imports data from the PEtab data file as a data list and defines errors if an error model is required.
#'
#' @param data PEtab data file as .tsv
#' @param observables observables as eqnvec
#'
#' @return data as data list and errors (if required) as eqnvec.
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#'
#' @export
#'
getDataPEtabSBML <- function(data, observables){
  mydata <- read.csv(file = data, sep = "\t")
  myobs <- read.csv(file = observables, sep = "\t") %>% as.data.frame()
  obs <- myobs$observableId %>% as.character()
  errors <- NULL
  
  # if(!is.null(mydata$noiseParameters) & !is.null(myobs$noiseFormula)) cat(red("Warning: errors specified in data and observable file.\n"))
  if(!mydata$noiseParameters %>% is.numeric) {
    # define errors
    errors <- myobs$noiseFormula %>% as.character()
    which_err <- c(1:length(obs))
    if(length(errors) != length(obs)) errors <- rep(errors,length(obs))
    names(errors) <- obs[which_err]
    errors <- as.eqnvec(errors)
    # set fixed sigmas
    if(is.null(mydata$noiseParameters)){
      mydata$noiseParameters <- errors %>% as.numeric()
    } else mydata$noiseParameters <- NA
  }
  
  # # rename observables with _obs
  # obs <- mydata$observableId %>% as.character() %>% paste0("_obs")
  # mydata$observableId <- obs
  # if(!is.null(errors)) names(errors) <- paste0(names(errors), "_obs")
  mydata$simulationConditionId <- as.character(mydata$simulationConditionId)
  if(!is.null(mydata$observableParameters)){
    for (observable in unique(mydata$observableId)){
      for(condition in unique(mydata$simulationConditionId)){
        sub <- subset(mydata, simulationConditionId==condition & observableId==observable)
        if(nrow(sub) > 0){
          if(length(unique(sub$observableParameters)) > 1){
            index <- which(mydata$simulationConditionId==condition & mydata$observableId==observable)
            mydata$simulationConditionId[index] <- paste0(mydata$simulationConditionId[index], "_",mydata$observableParameters[index])
          }
        }
      }
    }
  }
  
  # select necessary data columns
  data <- data.frame(name = mydata$observableId, time = mydata$time,
                     value = mydata$measurement, sigma = mydata$noiseParameters,
                     condition = mydata$simulationConditionId,
                     lloq = if (!is.null(mydata$lloq)) mydata$lloq else -Inf)
  obs2log <- myobs$observableId[which(myobs$observableTransformation=="log")]
  data$value[which(data$name%in%obs2log)] <- log(data$value[which(data$name%in%obs2log)])
  obs2log10 <- myobs$observableId[which(myobs$observableTransformation=="log10")]
  data$value[which(data$name%in%obs2log10)] <- log10(data$value[which(data$name%in%obs2log10)])
  data <- data %>% as.datalist()
  
  return(list(data=data,errors=errors))
}


#' Plot observables and states of an SBML/PEtab model
#'
#' @description This function plots data and fits of an SBML/PEtab model after the import.
#' Note: Certain objects generated during the model import by importPEtabSBML have to be present in the global environment and are used as default variables if not specified differntly.
#'
#' @param g observation function as obsfn
#' @param x prediction function as prdfn
#' @param p parameter function as parfn
#' @param mydata data as datalist
#' @param pars parameter as vector
#' @param times times as vector
#' @param ... further arguments going to \link{plotCombined}
#'
#' @return NULL
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#' @importFrom dMod theme_dMod scale_color_dMod scale_fill_dMod
#' @export
#'
plotPEtabSBML <- function(..., g1 = g,
                          x1 = x,
                          p1 = p0,
                          mydata1 = mydata,
                          pouter1 = pouter,
                          times1 = times,
                          errfn = err){
  if(!is.null(errfn)){
    prediction <- as.data.frame((g1*x1*p1)(times1, pouter1), errfn = err)
    data <- as.data.frame(mydata1)
    P1 <- ggplot() + geom_line(data=subset(prediction, ...), aes(x=time, y=value, color=condition)) +
      geom_ribbon(data=subset(prediction, ...), aes(x=time, ymin=value-sigma, ymax=value+sigma, fill=condition), color=NA, alpha=0.25) +
      geom_point(data=subset(data, ...), aes(x=time, y=value, color=condition)) +
      theme_dMod() + scale_color_dMod() + scale_fill_dMod() +
      facet_wrap(~name, scales="free")
  } else {
    prediction <- (g1*x1*p1)(times1, pouter1)
    P1 <- plotCombined(prediction, mydata1, ...)
  }
  print(P1)
  # plotPrediction(prediction)
  
  # return(modelname)
}


# -------------------------------------------------------------------------#
# Refactoring sbml import ----
# -------------------------------------------------------------------------#

#' Sanitize pow(x,y) equations
#'
#' @param rate 
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family sbmlImport
#'
#' @examples
#' sanitizePowers("EMAX * pow(Ac / Vc + eps, hill) / (pow(Ac / Vc + eps, hill) + pow(EC50, hill)) * cytoplasm")
#' rate <- "EMAX * pow(Ac /   Vc, hill+2)"
#' sanitizePowers(rate)
sanitizePowers <- function(rate) {
  parseD <- getParseData(parse(text = rate))
  pow_ids <- parseD[which(parseD$text == "pow")-1,"id"]
  for (pow_id in pow_ids){
    child_ids <- parseD[parseD$parent == pow_id & !parseD$text %in% c("(",")",","), "id"][-1]
    texts <- getParseText(parseD, child_ids)
    rate <- gsub(paste0("pow(",texts[1], ", ", texts[2], ")"), paste0("(",texts[1],")**(", texts[2], ")"),rate, fixed = TRUE)
  }
  rate
}


#' Title
#'
#' import reactions and adjust by means of compartments
#'
#' @param m SBML model
#' @param reaction index of reaction
#' @param compartments names of compartments
#'
#' @return
#' @export
#'
#' @examples
sbmlImport_getReactionDetails <- function(m, reaction, compartments) {
  Reactantstring <- ""
  Productstring <- ""
  eq <- m$getReaction(reaction)
  Reactantnr <- eq$getNumReactants(reaction)
  if(Reactantnr > 0) Reactantstring <- paste0( eq$getReactant(0)$getStoichiometry(), "*", eq$getReactant(0)$getSpecies())
  if(Reactantnr > 1) for (s in 1:(Reactantnr-1)) {
    Reactantstring <- paste0(Reactantstring, " + ",
                             paste0(eq$getReactant(s)$getStoichiometry(), "*", eq$getReactant(s)$getSpecies()))
  }
  Productnr <- eq$getNumProducts(reaction)
  if(Productnr > 0) Productstring <- paste0( eq$getProduct(0)$getStoichiometry(), "*", eq$getProduct(0)$getSpecies())
  if(Productnr > 1) for (s in 1:(Productnr-1)) {
    Productstring <- paste0(Productstring, " + ",
                            paste0(eq$getProduct(s)$getStoichiometry(), "*", eq$getProduct(s)$getSpecies()))
  }
  rate <- eq$getKineticLaw()$getFormula()
  if (stringr::str_detect(rate, "pow")) {
    rate <- sanitizePowers(rate)
  }
  if(!is.null(compartments)){
    if(Reduce("|", stringr::str_detect(rate, unique(compartments)))){
      rate <- cOde::replaceSymbols(unique(compartments), rep("1", length(unique(compartments))), rate)
    }
  }
  
  rate <- gsub(" ", "", rate)
  
  list(from = Reactantstring, to = Productstring,rate = rate)
}

#' Replace function calls by their mathematical expressions
#'
#' @param m sbml model
#' @param fun an integer with a function id (badly refactored)
#' @param rates reaction rate formulas
#'
#' @return modified rates
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML import
#' @importFrom cOde replaceSymbols
smblImport_replaceFunctions <- function(m, fun, rates) {
  mymath <- m$getFunctionDefinition(fun)$getMath()
  fun <- m$getFunctionDefinition(fun)
  funInfo <- sbmlImport_getFunctionInfo(fun)
  
  idx <- grep(funInfo$funName, rates)
  for (i in idx) {
    # 1. Replace function arguments by supplied arguments in rate, see petab_select/testcase0001 SBML for nontrivial examples
    # https://stackoverflow.com/questions/546433/regular-expression-to-match-balanced-parentheses
    regexExtractArguments <- paste0(".*", funInfo$funName, "(", "\\((?:[^)(]*(?R)?)*+\\)", ").*")
    argsInRateEqn <- gsub(regexExtractArguments, "\\1", rates[i], perl = TRUE)
    argsInRateEqn <- gsub("^\\(|\\)$", "", argsInRateEqn)
    argsInRateEqn <- strsplit(argsInRateEqn, ",", T)[[1]]
    funBody <- cOde::replaceSymbols(funInfo$funArgs, argsInRateEqn,funInfo$funBody)
    
    # 2. Replace function call in rate by funBody
    regexMatchFunCall <- paste0( "(", funInfo$funName, "\\((?:[^)(]*(?R)?)*+\\)", ")")
    cat("Replacing ", unique(gsub(paste0(".*",regexMatchFunCall, ".*"), "\\1", rates[i], perl = TRUE)), " by ", funBody, "\n")
    rates[i] <- gsub(regexMatchFunCall, funBody, rates[i], perl = TRUE)
  }
  
  rates
}



#' Get function info
#' 
#' function adapted from Frank
#'
#' @param fun some sbml-thingy representation of a function
#'
#' @return list with infos about the funciton
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML import
#' importFrom libSBML formulaToL3String
sbmlImport_getFunctionInfo <- function(fun) {
  if (is.null(fun)) return;
  id <- fun$getId()
  
  math <- fun$getMath()
  if (is.null(math)) return;
  
  num_children <- math$getNumChildren()
  
  funInfo <- list(
    funName = id,
    funArgs = vapply(seq_len(num_children-1), function(i) formulaToL3String(math$getChild(i-1)), FUN.VALUE = "a"),
    funBody = formulaToL3String(math$getChild(num_children - 1))
  )
  funInfo
}


sbmlImport_getEventStrings <- function(reactions, events = NULL) {
  for(fun in c("piecewise")){
    for(reaction in reactions$rates){
      if(stringr::str_detect(reaction, fun)){
        split <- stringr::str_split(reaction, fun)[[1]][2]
        count_bracket <- 0
        done <- F
        for(z in 1:nchar(split)){
          if(substr(split, z, z)=="(") count_bracket <- count_bracket+1
          if(substr(split, z, z)==")") count_bracket <- count_bracket-1
          if(count_bracket==0 & !done) {done <- T; pos <- z}
        }
        #pos <- which(strsplit(split, "")[[1]]==")")[2]
        event <- paste0(fun, substr(split, 1, pos))
        events <- c(events, event)
      }
    }
  }
  unique(events)
}


TransformEvents <- function(events){
  if(!is.null(events)){
    do.call(rbind, lapply(1:length(events), function(i){
      myevent <- events[i]
      if(stringr::str_detect(myevent, "piecewise") & (stringr::str_detect(myevent, "leq") | stringr::str_detect(myevent, "lt"))){
        expr1 <- strsplit(myevent, ",")[[1]][2]
        expr1 <- gsub(paste0(strsplit(expr1, "\\(")[[1]][1],"\\("), "", expr1)
        expr2 <- strsplit(strsplit(myevent, ",")[[1]][3], ")")[[1]][1]
        if(expr1=="time") timepoint <- expr2 else
          if(stringr::str_detect(expr1, "time-")) timepoint <- gsub("time-", "", expr1) else cat("Warning: Event not yet supported.")
        first <- strsplit(strsplit(myevent, "\\(")[[1]][2], ",")[[1]][1]
        second <- strsplit(strsplit(myevent, ",")[[1]][4], ")")[[1]][1]
        if(!is.na(suppressWarnings(as.numeric(timepoint)))) timepoint <- as.numeric(timepoint) # avoid warning if variable is not numeric
        if(!is.na(suppressWarnings(as.numeric(first)))) first <- as.numeric(first)
        if(!is.na(suppressWarnings(as.numeric(second)))) second <- as.numeric(second)
        return(data.frame(var=paste0("event",i), time=c(0,timepoint), value=c(first, second),root=NA, method="replace"))
      } else {cat("Warning: Event not yet supported"); return(myevent)}
    }))
  } else return(NULL)
}

#' Get Events from "events" entries in SBML
#' 
#' The difficulty is that triggers can take arbitrary expressions which doesn't play 
#' well with dMod, which has a rather stiff input format for events.
#' Therefore, this function is tailored to my own export ONLY
#'
#' @param modelfile 
#'
#' @return
#' @export
#'
#' @examples
sbmlImport_getEventsFromTrueEvents <- function(modelfile) {
  model = readSBML(modelfile)$getModel()
  nx <- Model_getNumEvents(model)
  
  eventList <- NULL
  for (i in seq_len(nx)-1) {
    ev <- model$getEvent(i)
    eventTime <- ev$getTrigger()
    eventTime <- eventTime$getMath()
    eventTime <- formulaToL3String(eventTime)
    if (!grepl("time >= ", eventTime)) stop("Event Trigger not of form 'time >= eventtime'. \nEvent will not be imported properly. \n Event trigger has value:\n", eventTime, "\nPlease add a parser for such events to petab::sbmlImport_getEventsFromTrueEvents")
    eventTime <- gsub("time >= ", "", eventTime) # hard coded for my own export which has always the same string. error prone for other models, but I can't capture all possible cases anyway
    eventTime <- as.numeric(eventTime)
    
    ea <- ev$getEventAssignment(0)              # hard coded for my own export which assumes only one eventAssignment per trigger
    eventSpecies <- ea$getVariable()
    eventValue <- ea$getMath()
    eventValue <- formulaToL3String(eventValue)
    # eventValue <- as.numeric(eventValue)  
    
    eventList <- addEvent(event = eventList, var = eventSpecies, time = eventTime, value = eventValue,
                          root = NA, method = "replace") # Again hard coded
  }
  eventList
}

#' Extract parameters from SBML
#' 
#' Assumes the parameters are given numerically
#' 
#' @param modelfile Path to SBML file
#'
#' @return [petab_parameters()] with scale = "lin" and estimate = 0
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML import
#'
#' @examples
#' modelfile <- file.path(petab_examplePath("01", "pe"), "model_petab.xml")
#' sbmlImport_getCompartmentsFromSBML(modelfile)
sbmlImport_getParametersFromSBML <- function(modelfile) {
  model = readSBML(modelfile)$getModel()
  nx <- model$getNumParameters()
  
  pars <- NULL
  for (i in (seq_len(nx)-1)) {
    mypar <- model$getParameter(i)$getId()
    value <- model$getParameter(i)$getValue()
    pars <- c(pars, setNames(value, mypar))
  }
  
  petab_parameters(parameterId = names(pars), 
                   parameterScale = "lin",
                   lowerBound = 1e-10, upperBound = 1e10,
                   nominalValue = pars,
                   estimate = 0)
}

#' Extract compartments from SBML
#' 
#' Assumes the compartments are given numerically
#'
#' @param modelfile Path to SBML file
#'
#' @return [petab_parameters()] with scale = "lin" and estimate = 0
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML import
#'
#' @examples
#' modelfile <- file.path(petab_examplePath("01", "pe"), "model_petab.xml")
#' sbmlImport_getCompartmentsFromSBML(modelfile)
sbmlImport_getCompartmentsFromSBML <- function(modelfile) {
  
  model <- readSBML(modelfile)$getModel()
  nx <- model$getNumCompartments()
  
  comp_name <- comp_size <- NULL
  for ( i in (seq_len(nx)-1) ) {
    # get compartment name and size
    compartmentId <- model$getCompartment(i)$getId()
    size          <- model$getCompartment(i)$getSize()
    comp_name <- c(comp_name,compartmentId)
    comp_size <- c(comp_size,size)
  }
  
  petab_parameters(parameterId = comp_name, 
                   parameterScale = "lin",
                   lowerBound = 1e-10, upperBound = 1e10,
                   nominalValue = comp_size,
                   estimate = 0)
}



#' Import initials from SBML.
#'
#' @description This function imports initial values or equations describing the same from SBML and writes them in a named vector.
#'
#' @param model SBML file as .xml
#'
#' @return Named vector of initials.
#'
#' @author Marcus Rosenblatt, Svenja Kemmer and Frank Bergmann
#' importFrom libSBML readSBML formulaToL3String
#' @export
#'
sbmlImport_getInitialsFromSBML <- function(modelfile){
  
  model = readSBML(modelfile)$getModel()
  
  initials <- species <- NULL
  nx <- model$getNumSpecies()
  # now we can go through all species
  for (i in (seq_len(nx)-1) ){
    current <- model$getSpecies(i)
    species <- c(species,current$getId())
    # now the species can have several cases that determine their initial value
    
    # CASE 1: it could be that the species is fully determined by an assignment rule
    # (that apply at all times), so we have to check rules first
    rule <- model$getRule(current$getId())
    if (!is.null(rule)) {
      # ok there is a rule for this species so lets figure out what its type
      # is as that determines whether it applies at t0
      rule_type <- rule$getTypeCode()
      type_name <- SBMLTypeCode_toString(rule_type, 'core')
      # CASE 1.1: Assignment rule
      if (type_name == "AssignmentRule"){
        # the initial value is determined by the formula
        math <- rule$getMath()
        if (!is.null(math)){
          formula <- formulaToL3String(math)
          # print(paste('Species: ', current$getId(), ' is determined at all times by formula: ', formula))
          initials <- c(initials,formula)
          # no need to look at other values so continue with next species
          next
        }
      }
      
      # CASE 1.2: Rate rule
      if (type_name == "RateRule") {
        math <- rule$getMath()
        if (!is.null(math)) {
          formula <- formulaToL3String(math)
          # print(paste('Species: ', current$getId(), ' has an ode rule with formula: ', formula))
          initials <- c(initials,formula)
          # there is an ODE attached to the species, its initial value is needed
        }
      }
    }
    
    # CASE 2 (or subcase of 1?): Initial assignment
    # it could have an initial assignment
    ia <- model$getInitialAssignment(current$getId())
    if (!is.null(ia)) {
      math <- ia$getMath()
      if (!is.null(math)) {
        formula <- formulaToL3String(math)
        # formula <- libSBML::formulaToL3String(math)
        # print(paste("Species: ", current$getId(), " has an initial assignment with formula: ", formula))
        initials <- c(initials,formula)
        # as soon as you have that formula, no initial concentration / amount applies
        # so we don't have to look at anything else for this species
        next
      }
    }
    
    # CASE 3 (or subcase of 1?): Inital amount
    # it could have an initial amount
    if (current$isSetInitialAmount()) {
      # print (paste("Species: ", current$getId(), "has initial amount: ", current$getInitialAmount()))
      initials <- c(initials,current$getInitialAmount())
    }
    
    # CASE 4 (or subcase of 1?): Inital Concentration
    # it could have an initial concentration
    if (current$isSetInitialConcentration())
    {
      # print (paste("Species: ", current$getId(), "has initial concentration: ", current$getInitialConcentration()))
      initials <- c(initials,current$getInitialConcentration())
    }
  }
  names(initials) <- species
  
  
  
  
  
  # Post process initials: Divide into symbolic and numeric
  # initials <- c(a = "1", b = "2")
  # initials <- c(a = "1", b = "a+f")
  # initials <- c(a = "c+d", b = "a+f")
  
  inits_sym <- suppressWarnings(initials[ is.na(as.numeric(initials))])
  inits_sym <- data.table(parameterId = names(inits_sym),
                          parameterFormula = inits_sym)
  
  inits_num <- suppressWarnings(initials[!is.na(as.numeric(initials))])
  inits_num <- suppressWarnings(petab_parameters(parameterId = names(inits_num), 
                                                 parameterScale = "lin",
                                                 lowerBound = 1e-10, upperBound = 1e10,
                                                 nominalValue = inits_num,
                                                 estimate = 0))
  inits_num <- inits_num[!is.na(parameterId)] # To catch the case that no numeric inits are given
  
  list(inits_sym = inits_sym, inits_num = inits_num)
}




# get parameters function: Logic:
# get 

# -------------------------------------------------------------------------#
# PETabIndiv ----
# -------------------------------------------------------------------------#

#' Rename parscales to the names needed in the base trafo
#'
#' @param parscales setNames(PETABPars$parameterScale, PETABPars$parameterId)
#' @param est.grid data.table
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
#' # Example 1: k1_A and k1_B are duplicated to two inner parameters 
#' est.grid <- data.frame(ID = 1:2,
#'                        condition = c("A", "B"),
#'                        k1 = c("k1_A", "k1_B"),
#'                        k1DUPE = c("k1_A", "k1_B"),
#'                        k3 = c("k3outer", NA),
#'                        stringsAsFactors = FALSE)
#' scales_outer <- c(k1_A = "log", k1_B = "log", k3outer = "lin")
#' updateParscalesToBaseTrafo(scales_outer, est.grid)
#' 
#' # Example 2: SHOULD FAIL k1_A and k1_B map to same inner parameter, but have idfferent scales
#' est.grid <- data.frame(ID = 1:2,
#'                        condition = c("A", "B"),
#'                        k1 = c("k1_A", "k1_B"),
#'                        k1DUPE = c("k1_A", "k1_B"),
#'                        k3 = c("k3outer", NA),
#'                        stringsAsFactors = FALSE)
#' scales_outer <- c(k1_A = "log", k1_B = "log10", k3outer = "lin")
#' updateParscalesToBaseTrafo(scales_outer, est.grid)
#' 
#' # Example 3: k4 is a parameter not in est.grid. It might be a parameter in fix.grid and should be returned as is
#' est.grid <- data.frame(ID = 1:2,
#'                        condition = c("A", "B"),
#'                        k1 = c("k1_A", "k1_B"),
#'                        k1DUPE = c("k1_A", "k1_B"),
#'                        k3 = c("k3outer", NA),
#'                        stringsAsFactors = FALSE)
#' scales_outer <- c(k1_A = "log", k1_B = "log", k3outer = "lin", k4 = "log")
#' updateParscalesToBaseTrafo(scales_outer, est.grid)
updateParscalesToBaseTrafo <- function(scales_outer, est.grid) {
  # parscales = c(outername = scaleouter)
  
  # Get name mapping between est.grid pars and outer pars
  pars_inner_outer <- getEstGridParameterMapping(est.grid) # c(inner = outer)
  pars_inner_outer <- pars_inner_outer[!is.na(pars_inner_outer)]
  pars_outer_inner <- setNames(names(pars_inner_outer), pars_inner_outer) # c(outer = inner)
  
  # Get scales for inner pars
  scales_inner <- setNames(scales_outer[pars_inner_outer], names(pars_inner_outer))
  
  # Determine if there are outer parameters mapping to the same inner parameter with different scales
  dupes <- names(scales_inner)[duplicated(names(scales_inner))]
  dupes <- unique(dupes)
  for (d in dupes) {
    if (length(unique(scales_inner[names(scales_inner) == d])) > 1)
      stop("The following parameter refers to the same structural model parameter, but has different ",
           "scales in different conditions. This is not allowed. \n",
           "Parameter: ", d , "\n",
           "Outer pars: ", paste0(names(pars_outer_inner[pars_outer_inner == d]), collapse = ", "))
  }
  
  # Remove the duplicates and return updated parscales
  scales_inner <- scales_inner[!duplicated(names(scales_inner))]
  
  # Append scales which did not appear in est.grid
  scales_fix <- scales_outer[setdiff(names(scales_outer),names(pars_outer_inner))]
  
  c(scales_inner, scales_fix)
}

#' Determine the type of trafo for each element of a vector
#'
#' @param trafo_string named character
#'
#' @return character(lenth(trafo_string)):
#' TRAFO = functional trafo
#' SYMBOL = Parameter name mapping
#' NUMBER = fixed value
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
#' trafo_string <- c(STAT5A = "207.6 * ratio", STAT5B = "207.6 - 207.6 * ratio",
#'                   pApB = "alpha", pApA = "beta", pBpB = "pBpB", nucpApA = "0.1", nucpApB = "0",
#'                   nucpBpB = "0")
#' getTrafoType(trafo_string)
getTrafoType <- function(trafo_string) {
  opt.keep.source <- getOption("keep.source") # important for Rscript --vanilla
  options(keep.source = TRUE)
  
  out <- vapply(names(trafo_string), function(nm) {
    ts <- trafo_string[nm]
    pd <- getParseData(parse(text = ts))
    if (nrow(pd) > 2) return("TRAFO")
    if (pd[1,"token"] == "SYMBOL") {
      if (pd[1,"text"] == nm) return("SYMBOL")
      else return("TRAFO")
    }
    if (pd[1,"token"] == "NUM_CONST") return("NUMBER")
    stop("Unkown trafo type: ", ts)
  }, FUN.VALUE = "TYPE")
  
  options(keep.source = opt.keep.source)
  
  out
}

#' Title
#'
#' @param mycondition.grid 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML Import
#' @importFrom data.table as.data.table setnames setcolorder copy
#'
#' @examples
pdIndiv_initializeGridlist <- function(mycondition.grid) {
  # Copy condition.grid, take unique identifying column only
  cg <- copy(mycondition.grid)
  cg <- lapply(cg, as.character)
  cg <- data.table::as.data.table(cg)
  cg <- cg[,!"conditionName"]
  data.table::setnames(cg, "conditionId", "condition")
  data.table::setcolorder(cg, "condition")
  # Determine which columns contain values and/or parameter names
  is_string  <- suppressWarnings(vapply(cg[,-1], function(x) any(is.na(as.numeric(x))), FUN.VALUE = TRUE))
  is_string  <- which(is_string) + 1
  is_numeric <- suppressWarnings(vapply(cg[,-1], function(x) any(!is.na(as.numeric(x))), FUN.VALUE = TRUE))
  is_numeric <- which(is_numeric) + 1
  
  # Initialize fix.grid and est.grid
  # For mixed columns (string & numeric), need NA in the respective places
  fix.grid <- data.table::copy(cg)
  fix.grid <- fix.grid[,.SD,.SD = c(1, is_numeric)]
  suppressWarnings(fix.grid[,(names(fix.grid)[-1]) := lapply(.SD, as.numeric), .SDcols = -1])
  fix.grid[,`:=`(ID = 1:.N)]
  
  est.grid <- data.table::copy(cg)
  est.grid <- est.grid[,.SD,.SD = c(1, is_string)]
  suppressWarnings(est.grid[,(names(est.grid)[-1]) := lapply(.SD, function(x) {replace(x, !is.na(as.numeric(x)), NA)}), .SDcols = -1])
  est.grid[,`:=`(ID = 1:.N)]
  
  gl <- list(est.grid = copy(est.grid), fix.grid = copy(fix.grid))
  gl
}


#' Title
#'
#' @param x 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML Import
#'
#' @examples
pdIndiv_getParameters_conditionGrid <- function(x) {
  x <- suppressWarnings(x[,!c("conditionId", "conditionName", "ID")])
  x <- unlist(x, use.names = FALSE)
  x <- unique(x)
  x <- suppressWarnings(x[is.na(as.numeric(x))])
  x
}


#' Title
#'
#' @param cg condition.grid
#' @param scalesOuter vector of scales of outer parameters c(CS_name = scale) 
#'
#' @return vector of scales of base parameters c(BASE_name = scale) 
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family SBML Import
#'
#' @examples
#' cg <- data.table(structure(list(conditionId = c("C1", "C2"), conditionName = c("C1", "C2"), kcat = c("kcat_C1", "kcat_C2"), kon = c("kon_C1", "kon_C2"), koff = c(0.1, 0.2), observableParameter1_obsE = c("offset_E", "offset_E"), observableParameter1_obsS = c("offset_S", "offset_S_C2"), noiseParameter1_obsE = c("sigmaRel_obsE", "sigmaRel_obsE"), noiseParameter1_obsES = c("sigma_obsES", "sigma_obsES"), noiseParameter1_obsP = c("sigma_obsP", "sigma_obsP"), noiseParameter1_obsS = c("sigma_obsS", "sigma_obsS_C2"), noiseParameter2_obsE = c("sigmaAbs_obsE", "sigmaAbs_obsE"), E = c("E", "E"), ES = c("ES", "ES"), P = c("P", "P"), S = c("S", "S"), cytoplasm = c("cytoplasm", "cytoplasm")), row.names = c(NA, -2L), class = c("data.frame")))
#' scalesOuter <- c(E = "log10", ES = "log10", kcat_C1 = "log10", kcat_C2 = "log10", offset_E = "log10", offset_S = "log10", offset_S_C2 = "log10", P = "log10", S = "log10", sigma_obsES = "log10", sigma_obsP = "log10", sigma_obsS = "log10", sigma_obsS_C2 = "log10", sigmaAbs_obsE = "log10", sigmaRel_obsE = "log10", kon_C1 = "log10", kon_C2 = "log10", cytoplasm = "lin")
#' pdIndiv_getBaseScales(cg, scalesOuter)
pdIndiv_getBaseScales <- function(cg, scalesOuter) {
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
  scalesBase
}



#' Import an SBML model and corresponding PEtab objects
#'
#' @description This function imports an SBML model and corresponding PEtab files, e.g. from the Benchmark collection.
#'
#' @param filename "path/to/modelname.petab". Will look for files like
#'        "path/to/modelname/model_modelname.xml"
#' @param TestCases TRUE to load feature test cases
#' @param path2TestCases path to feature test case folder
#' @param compile if FALSE, g, ODEmodel and err are loaded from .RData (if present) and compilation time is saved
#' @param currFitName name of the current mstrust fit
#'
#' @details Objects such as model equations, parameters or data are automatically assigned to the following standard variables and written to your current working directory (via <<-):
#' reactions, observables, errors, g, x, p0, err, obj, mydata, ODEmodel, condition.grid, trafoL, pouter, times.
#' Compiled objects (g, ODEmodel and err) are saved in .RData.
#'
#' @return list of dMod model objects with entries
#' * symbolicEquations
#' * odemodel
#' * data
#' * gridlist
#' * e
#' * fns
#' * prd
#' * obj_data
#' * pars
#'
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de) building on the original function of Marcus and Svenja
#' @md
#' @importFrom dMod repar as.parframe
#' @importFrom cOde getSymbols
#' @importFrom conveniencefunctions dMod_saveMstrust
#'
#' @export
importPEtabSBML_indiv <- function(filename = "enzymeKinetics/enzymeKinetics.yaml",
                                  testCases = FALSE,
                                  path2TestCases = "PEtabTests/",
                                  .compiledFolder = file.path("CompiledObjects"),
                                  NFLAGcompile = c(Auto = 3, Recompile = 0, RebuildGrids = 1, LoadPrevious = 2)[3],
                                  SFLAGbrowser = c("0None", "1Beginning", "2BuildGrids", "3Compilation", "4CollectList",
                                                   "5Scales", "6ParameterFormulaInjection")[1],
                                  currFitName = "mstrust"
)
{
  if (grepl(SFLAGbrowser,"1Beginning")) browser()
  # .. Set exit behaviour -----
  mywd <- getwd()
  on.exit({setwd(mywd)})
  
  # .. Define path to SBML, PEtab and pd files -----
  filename      <- path.expand(filename)
  modelname     <- petab_modelname_path(filename)$modelname
  files         <- petab_files(filename)
  filenameParts <- list(modelname = modelname, .currentFolder = mywd, .compiledFolder = .compiledFolder, 
                        type = "indiv", petabYaml = if(grepl("yaml", filename)) filename else NULL,
                        .resultsFolder = file.path(dirname(.compiledFolder), "Results"),
                        .projectFolder = dirname(.compiledFolder)
  )
  rdsfile       <- pd_files(filenameParts)$rdsfile
  
  # .. Read PEtab tables -----
  pe <- readPetab(filename)
  
  # .. Read previously imported files -----
  dir.create(.compiledFolder, showWarnings = FALSE)
  if (NFLAGcompile == 3){
    if (file.exists(rdsfile)) {
      peOld <- readRDS(rdsfile)$pe
      NFLAGcompile <- as.numeric(petab_hash(pe) == petab_hash(peOld)) * 2 # Um die Ecke wegen suboptimaler NFLAGcompile Definition
    } else {
      NFLAGcompile <- 0 
    }
  }
  
  if(NFLAGcompile > 0) {
    pd <- readPd(rdsfile, currFitName)
    if (NFLAGcompile == 2) return(pd)
  }
  
  if(NFLAGcompile == 0) {
    .resultsFolder <- filenameParts$.resultsFolder
    results <- list.files(.resultsFolder, recursive = TRUE)
    if (any(grepl(paste0(currFitName, ".rds|profile|L1|"), results))) {
      cat(paste0("Deleting the following results: ", 
                 paste0(grep(paste0(currFitName, ".rds|profile|L1"), results, value = TRUE), collapse = ",\n")))
      if (readline(" Are you sure? (type yes)") != "yes") stop("import stopped")}
    unlink(.resultsFolder,recursive = TRUE)
  }
  
  # check this for some coming versions before deprecating
  if (!is.null(pe$meta$events)) {
    stop("pe$meta$events is deprecated as of version 0.4.1. Please specify your events in pe$model :)")
  }
  
  ## load required packages
  require(libSBML) # => Not very nice, better explicitly import the required functions
  
  # .. Model Definition - Equations -----
  dummy            <- getReactionsSBML(files$modelXML, files$experimentalCondition)
  myreactions      <- dummy$reactions
  myreactions_orig <- dummy$reactions_orig
  myevents         <- rbind(dummy$events, sbmlImport_getEventsFromTrueEvents(files$modelXML)) # hack because events can take many forms in sbml
  mypreeqEvents    <- dummy$preeqEvents
  myobservables    <- getObservablesSBML(files$observables)
  
  dummy    <- getDataPEtabSBML(files$measurementData, files$observables)
  mydata   <- dummy$data
  myerrors <- dummy$errors
  
  
  # [ ] Pre-Equi Events
  if(!is.null(mypreeqEvents))
    stop("Pre-Equilibration is not yet implemented")
  # [ ] Need example for preeqEvents
  
  if (!is.null(myevents)){
    # set remaining event initials to 0. <= From Original import. Don't know what this means
    inits_events <- setdiff(unique(myevents$var), unique(mypreeqEvents$var))
    inits_events <- setNames(rep(0, length(inits_events)), inits_events)
  }
  
  # .. Initialize gridlist and trafo -----
  if (grepl(SFLAGbrowser, "2BuildGrids")) browser()
  cg <- copy(pe$experimentalCondition)
  
  # Get parameter information from all possible sources
  # Estimated parameters. Possible sources
  # * pe$parameters
  parsEst <- pe$parameters[estimate == 1] # replaces myfit_values
  parsEst[,`:=`(estValue = eval(parse(text = paste0(parameterScale, "(", nominalValue, ")")))),by = 1:nrow(parsEst)]
  
  # Fixed parameters. 
  # Possible sources 
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
  
  # Possible sinks => Remove from fix
  # * pe$meta$parameterFormulaInjection
  parsFix <- parsFix[!parameterId %in% pe$meta$parameterFormulaInjection$parameterId]
  
  # .. Add measurement parameters -----
  obsParMapping <- petab_getMeasurementParsMapping(pe, column = "observableParameters")
  cg <- merge(cg, obsParMapping, by = "conditionId", all.x = TRUE)
  
  # Check if non-numeric entries are in noiseParameters of measurementData.
  # if so, an error model must be employed. useErrormodel is only FALSE, if
  # all entries in noiseParameters of measurementData can be casted as numeric without error.
  useErrormodel <- any(is.na(suppressWarnings(as.numeric(pe$measurementData$noiseParameters))))
  
  if (useErrormodel == TRUE) {
    errParMapping <- petab_getMeasurementParsMapping(pe, column = "noiseParameters")
    cg <- merge(cg, errParMapping, by = "conditionId", all.x = TRUE)
  }
  
  
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
  scalesBase <- pdIndiv_getBaseScales(cg, scalesOuter)
  
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
                                  getSymbols(as.character(myevents$value)),
                                  getSymbols(pe$meta$parameterFormulaInjection$parameterFormula))))
  trafo <- trafo[trafo != "time"]
  
  # Insert inits
  trafo_inits <- sbmlImport_getInitialsFromSBML(files$modelXML)$inits_sym
  trafo[trafo_inits$parameterId] <- trafo_inits$parameterFormula
  trafo_assignmentRules <- NULL # [] Todo, if there is anything that needs to be done...
  
  # Insert scales
  trafo <- repar("x ~ 10**(x)", trafo = trafo, x = names(which(scalesBase=="log10")))
  trafo <- repar("x ~ exp(x)" , trafo = trafo, x = names(which(scalesBase=="log")))
  
  if (grepl(SFLAGbrowser,"5InspectTrafo")) browser()
  
  
  # .. ParameterFormulaInjection -----
  if (grepl(SFLAGbrowser, "6ParameterFormulaInjection")) browser()
  trafoInjected <- NULL
  pfi <- pe$meta$parameterFormulaInjection
  if (!is.null(pfi)) {
    # Probably over-cautios: This shouldn't happen. Can probably be removed
    check_pfiEstimated <- intersect(pfi$parameterId, names(gl$est.grid))
    if (length(check_pfiEstimated)) stop("These trafoInjected parameters are in est.grid: ", paste0(check_pfiEstimated, collapse = ", "))
    # Remove injected Parameter from est-trafo
    trafo <- trafo[setdiff(names(trafo), pfi$parameterId)]
    for (pid in intersect(names(gl$fix.grid), pfi$parameterId)) gl$fix.grid[,(pid) := NULL,]
    # trafoInjected
    trafoInjected <- setNames(pfi$parameterFormula, pfi$parameterId)
  }
  
  # -------------------------------------------------------------------------#
  # .. Model Compilation -----
  # -------------------------------------------------------------------------#
  if (grepl(SFLAGbrowser, "3Compilation")) browser()
  if (NFLAGcompile == 0) {
    setwd(.compiledFolder)
    cat("Compiling g\n")
    myg <- dMod::Y(myobservables, myreactions, compile=TRUE, modelname=paste0("g_",modelname))
    
    cat("Compiling odemodel\n")
    myodemodel <- dMod::odemodel(myreactions, forcings = NULL, events = myevents, fixed=NULL,
                                 estimate = c(dMod::getParametersToEstimate(est.grid = gl$est.grid,
                                                                            trafo = c(trafo, trafoInjected),
                                                                            reactions = myreactions),
                                              pe$meta$parameterFormulaInjection$parameterId),
                                 modelname = paste0("odemodel_", modelname),
                                 jacobian = "inz.lsodes", compile = TRUE)
    
    myx <- dMod::Xs(myodemodel,
                    optionsOde = list(method = "lsoda", rtol = 1e-10, atol = 1e-10, maxsteps = 5000),
                    optionsSens = list(method = "lsodes", lrw=200000, rtol = 1e-10, atol = 1e-10))
    
    if (useErrormodel) {
      cat("Compiling errormodel\n")
      myerr <- NULL
      if(length(cOde::getSymbols(myerrors)))
        myerr <- dMod::Y(myerrors, f = c(as.eqnvec(myreactions), myobservables), states = names(myobservables),
                         attach.input = FALSE, compile = TRUE, modelname = paste0("errfn_", modelname))
    } else {
      cat("Data is passed with fixed errors. Don't use error model\n")
      myerr <- NULL
    }
    
    
    # [ ] Pre-equilibration
    # mypSS <- Id()
    
    p1 <- dMod::Id()
    if (length(trafoInjected)){
      cat("Compiling pInjected\n")
      p1 <- dMod::P(trafoInjected, compile = TRUE, modelname = paste0("PInjected_", modelname),
                    attach.input = TRUE)
    }
    
    cat("Compiling p\n")
    p0 <- dMod::P(trafo, compile = TRUE, modelname = paste0("P_", modelname))
    setwd(mywd)
    
  } else if (NFLAGcompile == 1) {
    myg        <- pd$dModAtoms$fns$g
    myx        <- pd$dModAtoms$fns$x
    p1         <- pd$dModAtoms$fns$p1
    p0         <- pd$dModAtoms$fns$p0
    myodemodel <- pd$dModAtoms$odemodel
    myerr      <- pd$dModAtoms$e
  }
  
  # -------------------------------------------------------------------------#
  # Collect list ----
  # -------------------------------------------------------------------------#
  if (grepl(SFLAGbrowser, "4CollectList")) browser()
  fns <- list(
    g = myg,
    x = myx,
    # p1 = mypSS, # [ ] Pre-Equilibration
    p1 = p1,
    p0 = p0
  )
  symbolicEquations <- list(
    reactions = myreactions,
    observables = myobservables,
    errors  = myerrors,
    trafo = trafo,
    trafoInjected = trafoInjected)
  
  pars <- setNames(parsEst$estValue, parsEst$parameterId)
  pars <- pars[setdiff(names(pars), names(trafoInjected))]
  
  times <- c(
    pe$measurementData$time,
    pe$meta$presimTimes
  )
  
  # .. Collect final list -----
  if (useErrormodel == TRUE) {
    pd <- list(
      # petab
      pe                 = pe,
      # Basic dMod elements
      dModAtoms          = list(
        # [ ] add events!
        symbolicEquations  = symbolicEquations,
        odemodel           = myodemodel,
        data               = mydata,
        gridlist           = gl,
        e                  = myerr,
        fns                = fns
      ),
      # other components: Dump your stuff here
      filenameParts = filenameParts,
      # Parameters + Time
      pars               = pars,
      times              = predtimes(times, Nobjtimes = 200)
    )
  } else {
    pd <- list(
      # petab
      pe                 = pe,
      # Basic dMod elements
      dModAtoms          = list(
        # [ ] add events!
        symbolicEquations  = symbolicEquations,
        odemodel           = myodemodel,
        data               = mydata,
        gridlist           = gl,
        e                  = dMod::Id(),
        fns                = fns
      ),
      # other components: Dump your stuff here
      filenameParts = filenameParts,
      # Parameters + Time
      pars               = pars,
      times              = predtimes(times, Nobjtimes = 200)
    )
  }
  
  
  
  # High level prediction function
  pd <- pdIndiv_rebuildPrdObj(pd = pd,Nobjtimes = 100)
  
  
  # .. Save and return -----
  saveRDS(pd, rdsfile)
  
  if (grepl(SFLAGbrowser, "7evaluateObj")) browser()
  # Try to evaluate obj for first time
  cat("Evaluating obj ... ")
  value_base <- tryCatch(pd$obj(pd$pars)$value, error = function(x) NA)
  cat("objective value: ", value_base, "\n")
  # suggestion for the future: Include parameterSetId from the beginning
  # parf_base <- dMod::parframe(data.frame(parameterSetId = "Base", value = value_base, index = 1, converged = FALSE, iterations = 1, pd$pars),parameters = names(pd$pars)) 
  # old solution
  parf_base <- dMod::as.parframe(structure(list(list(value = value_base, index = 1, converged = FALSE, iterations = 1, argument = pd$pars)), class = c("parlist", "list")))
  conveniencefunctions::dMod_saveMstrust(parf_base, dirname(dirname(rdsfile)), identifier = "base", FLAGoverwrite = TRUE)
  
  # return pd
  readPd(rdsfile, currFitName)
}

# .. pd_parameterNames -----
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @importFrom data.table data.table
#' @importFrom tibble tribble
pd_parameterNames <- data.table::data.table(tibble::tribble(
  ~STAGEID , ~STAGENAME,        
  "CS"     , "Condition specific, on est-scale"    ,
  "BASE"   , "Not condition specific, on est-scale",
  "NOMINAL", "Base parameters on nominal scale"    ,
  "INNER"  , "Base parameters after parameter formula injections"
))

