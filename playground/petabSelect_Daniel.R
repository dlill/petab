# -------------------------------------------------------------------------#
# 0 load packages ----
# -------------------------------------------------------------------------#
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
devtools::load_all()


# 
pathPSProj <- "petabSelectProject"
dir.create(pathPSProj, showWarnings = FALSE, recursive = TRUE)

pathPetabProjects <- file.path(pathPSProj, "petabProjects")
dir.create(pathPetabProjects, showWarnings = FALSE, recursive = TRUE)

path_problem <- "/home/daniel/Promotion/Promotion/Software/petab_select/test_cases/0003/petab_select_problem.yaml"


# petabSelect <- function(path) {
# load python packages
pepy <- petab_python_setup()
peps <- petab_python_setup(FLAGreturnpetabSelect = T)

# load selection problem and model space
selection_problem <- peps$Problem$from_yaml(path_problem)
candidates <- peps$candidates(problem = selection_problem)

# ..  -----
candidates <- peps$candidates(problem = selection_problem)
models <- candidates$models
m  <- (models)[[1]]
for (m  in models) {
  
  modelId <- substr(m$model_id, 1,8)
  cat(modelId, "\n")
  pathm <- file.path(pathPetabProjects, modelId, "petab")
  pathmPSYaml <- file.path(pathPetabProjects, modelId, "petab", paste0("petabSelect_", modelId, ".yaml"))
  pathmPetab <- file.path(pathPetabProjects, modelId, "petab", paste0(modelId, ".yaml"))
  pathmCompiled <- file.path(pathPetabProjects, modelId, "Compiled")
  
  if (file.exists(pathmPSYaml)) next
  
  m$to_petab(output_path = pathm)
  file.copy(file.path(pathPetabProjects, modelId, "petab", "problem.yaml"), pathmPetab)
  
  pd <- importPEtabSBML_indiv(pathmPetab, .compiledFolder = pathmCompiled, NFLAGcompile = 3)
  pd <- pd_fitMstrust(pd)
  
  reportYaml <- pd_petabSelect_reportYaml(pd)
  criterion  <- gsub("Criterion.", "", selection_problem$criterion)
  m$set_criterion(selection_problem$criterion, reportYaml$criteria[[criterion]])
  m$to_yaml(pathmPSYaml)

}

selection_problem$add_calibrated_models(models)
m_best <- selection_problem$get_best()

selection_problem$calibrated_models


# ..  -----
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

# define initial test_models
test_models <- model_space$neighbors(candidate_space)
i <- 1

while(length(test_models) > 0){
  
  # iteration process
  invisible(
    pd_iterate_through_models(test_models, mycriterion)
  )
  
  # define best model and select next candidate models
  chosen_model <- selection_problem$get_best(test_models)
  candidate_space$reset(chosen_model)
  print(paste0("Iteration ", i, " : ", chosen_model$model_id))
  
  # define new test_models
  test_models <- model_space$neighbors(candidate_space)
  i <- i+1
}
chosen_model
}



pd_iterate_through_models <- function(test_models, mycriterion) {
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
    # evaluate pd
    pd_mstrust(pd = pd, NFLAGsavePd = T, iterlim = 1000, nfits = 5, id = id)
    
    # write yaml
    myreport <- pd_petabSelect_reportYaml(pd, FLAGwriteYaml = T)
    
    # include criterion in test_models object
    criterion_value <- myreport$criteria[[mycriterion]]
    model$set_criterion(mycriterion, criterion_value)
    
    sink()
  })
}

if(F){
  mymodel <- petabSelect(path = "../inst/petab_select_examples/0002/selection_problem.yaml")
}

