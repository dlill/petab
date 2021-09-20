# -------------------------------------------------------------------------#
# 0 load packages ----
# -------------------------------------------------------------------------#
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# R packages
# library(petab)
devtools::load_all()

# petab python package
pepy <- petab_python_setup()
# petab_select python package
peps <- petab_python_setup(FLAGreturnpetabSelect = T)

# -------------------------------------------------------------------------#
# 1 select initial model  ----
# -------------------------------------------------------------------------#

selection_problem <- peps$Problem$from_yaml("../inst/petab_select_examples/0002/selection_problem.yaml")
# test <- yaml::read_yaml("inst/petab_select_examples/0002/selection_problem.yaml")
model_space <- selection_problem$model_space

# find neighbours of the model in the model space
# peps$model_space$ModelSpace$neighbors(candidate_space = model_space, limit = NULL, exclude = TRUE)

# model_space$parameter_ids

mypars <- structure(rep(0, length(model_space$parameter_ids)), names = model_space$parameter_ids)

initial_virtual_model <-  peps$Model(
  model_id='INITIAL_VIRTUAL_MODEL',
  petab_yaml='.',
  parameters=as.list(mypars)
)

candidate_space <- peps$ForwardCandidateSpace(initial_virtual_model)
test_models <-  model_space$neighbors(candidate_space)

# -------------------------------------------------------------------------#
# 2 iterate through current models and generate report.yaml(s) ----
# -------------------------------------------------------------------------#

lapply(test_models, function(model){
  id <- as.character(model$model_id)
  print(id)
  sink("/dev/null")
  # create model folder
  modelpath <- paste0("SelectionProblem/", id, "/petab")
  dir.create(modelpath, recursive = T)
  # generate petab
  test_models[[1]]$to_petab(output_path = modelpath)
  # generate pd
      # pe<-readPetab(file.path(modelpath,"problem.yaml"))
      # petab_plotData(pe)
  pd <- importPEtabSBML_indiv(file.path(modelpath,"problem.yaml"), 
                              NFLAGcompile = "0", 
                              .compiledFolder = file.path(modelpath,"../CompiledObjects"))
  # evaluate pd (better use ms_trust!)
  pd_predictAndPlot2(pd)
  # write yaml
  pd_petabSelect_reportYaml(pd, FLAGwriteYaml = T)
  sink()
})

# -------------------------------------------------------------------------#
# 3 select best model ----
# -------------------------------------------------------------------------#






# -------------------------------------------------------------------------#
# 4 select next neighbours ----
# -------------------------------------------------------------------------#
