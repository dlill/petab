library(petab)
library(reticulate)

petab_python_setup()

reticulate::virtualenv_list()
reticulate::use_virtualenv("petab")
reticulate::py_config()
reticulate::miniconda_path()
# ..  -----
reticulate::import("petab")
convert = TRUE
as = NULL
delay_load = FALSE
module <- "petab"
wup <- function (module, as = NULL, convert = TRUE, delay_load = FALSE)
{
  if (!is.null(as)) {
    register_class_filter(function(classes) {
      sub(paste0("^", module), as, classes)
    })
  }
  delay_load_environment <- NULL
  delay_load_priority <- 0
  delay_load_functions <- NULL
  if (is.function(delay_load)) {
    delay_load_functions <- list(on_load = delay_load)
    delay_load <- TRUE
  }
  else if (is.list(delay_load)) {
    delay_load_environment <- delay_load$environment
    delay_load_functions <- delay_load
    if (!is.null(delay_load$priority))
      delay_load_priority <- delay_load$priority
    delay_load <- TRUE
  }
  if (!delay_load || is_python_initialized()) {
    reticulate:::ensure_python_initialized(required_module = module)
    reticulate:::py_module_import(module, convert = convert)
  }
  else {
    if (is.null(.globals$delay_load_module) || (delay_load_priority >
                                                .globals$delay_load_priority)) {
      .globals$delay_load_module <- module
      .globals$delay_load_environment <- delay_load_environment
      .globals$delay_load_priority <- delay_load_priority
    }
    module_proxy <- new.env(parent = emptyenv())
    module_proxy$module <- module
    module_proxy$convert <- convert
    if (!is.null(delay_load_functions)) {
      module_proxy$get_module <- delay_load_functions$get_module
      module_proxy$before_load <- delay_load_functions$before_load
      module_proxy$on_load <- delay_load_functions$on_load
      module_proxy$on_error <- delay_load_functions$on_error
    }
    attr(module_proxy, "class") <- c("python.builtin.module",
                                     "python.builtin.object")
    module_proxy
  }
}


# ..  -----
