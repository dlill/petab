# -------------------------------------------------------------------------#
# 0 load packages ----
# -------------------------------------------------------------------------#
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# R packages
# library(petab)
devtools::load_all()

path <- "../inst/petab_select_examples/0002/selection_problem.yaml"
iterations <- 2

petabSelect <- function(path, iterations) {
  # load python packages
  pepy <- petab_python_setup()
  peps <- petab_python_setup(FLAGreturnpetabSelect = T)
  
  # load selection problem and model space
  selection_problem <- peps$Problem$from_yaml(path)
  model_space <- selection_problem$model_space
  
  # generate initial model (might be replaced by python function)
  mypars <- structure(rep(0, length(model_space$parameter_ids)), names = model_space$parameter_ids)
  
  initial_virtual_model <-  peps$Model(
    model_id='INITIAL_VIRTUAL_MODEL',
    petab_yaml='.',
    parameters=as.list(mypars)
  )
  
  # define candidate space and test criterion
  candidate_space <- peps$ForwardCandidateSpace(initial_virtual_model)
  mycriterion <- selection_problem$criterion
  
  for(i in 1:iterations){
    # define test_models
    test_models <- model_space$neighbors(candidate_space)
    # iteration process
    invisible(
      pd_iterate_through_models(test_models)
      )
    
    # define best model and select next candidate models
    chosen_model <- selection_problem$get_best(test_models)
    candidate_space$reset(chosen_model)
    print(paste0("Iteration ", i, " : ", chosen_model$model_id))
  }
  chosen_model
}


pd_iterate_through_models <- function(test_models){
  lapply(test_models, function(model){
    # for(model in test_models){
    id <- as.character(model$model_id)
    print(id)
    sink("/dev/null")
    # create model folder
    modelpath <- paste0("SelectionProblem/", id, "/petab")
    dir.create(modelpath, recursive = T)
    # generate petab
    model$to_petab(output_path = modelpath)
    # generate pd
    # pe<-readPetab(file.path(modelpath,"problem.yaml"))
    # petab_plotData(pe)
    pd <- importPEtabSBML_indiv(file.path(modelpath,"problem.yaml"), 
                                NFLAGcompile = "0", 
                                .compiledFolder = file.path(modelpath,"../CompiledObjects"))
    # evaluate pd (better use ms_trust!)
    # pd_predictAndPlot2(pd)
    
    # write yaml
    myreport <- pd_petabSelect_reportYaml(pd, FLAGwriteYaml = T)
    
    # include criterion in test_models object
    criterion_value <- myreport$criteria[[mycriterion]]
    model$set_criterion(mycriterion, criterion_value)
    
    sink()
  })
}

mymodel <- petabSelect(path = "../inst/petab_select_examples/0002/selection_problem.yaml",
                       iterations = 2)
