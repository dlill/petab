# -------------------------------------------------------------------------#
# 0 load packages ----
# -------------------------------------------------------------------------#
# R packages
library(petab)

# petab python package
pepy <- petab_python_setup()
# petab_select python package
peps <- petab_python_setup(FLAGreturnpetabSelect = T)

# -------------------------------------------------------------------------#
# 1 select initial model  ----
# -------------------------------------------------------------------------#

selection_problem <- peps$Problem$from_yaml("inst/petab_select_examples/0002/selection_problem.yaml")
# test <- yaml::read_yaml("inst/petab_select_examples/0002/selection_problem.yaml")
model_space <- selection_problem$model_space

# find neighbours of the model in the model space
peps$model_space$ModelSpace$neighbors(candidate_space = model_space, limit = NULL, exclude = TRUE)

model_space$parameter_ids

mypars <- structure(rep(0, length(model_space$parameter_ids)), names = model_space$parameter_ids)

initial_virtual_model <-  peps$Model(
  model_id='INITIAL_VIRTUAL_MODEL',
  petab_yaml='.',
  parameters=as.list(mypars)
)

candidate_space <- peps$ForwardCandidateSpace(initial_virtual_model)
test_models <-  model_space$neighbors(candidate_space)


dir.create("bla")
test_models[[1]]$to_petab(output_path = "bla")

pe<-readPetab("bla")

# -------------------------------------------------------------------------#
# 2 read PEtab of initial model ----
# -------------------------------------------------------------------------#



# -------------------------------------------------------------------------#
# 3 evaluate model and write yaml ----
# -------------------------------------------------------------------------#

# -------------------------------------------------------------------------#
# 4 select next candidates ----
# -------------------------------------------------------------------------#
