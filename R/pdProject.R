
filename_petab <- "bla/petab/mypetab.yaml"
pdFolder_files <- function(filename_petab) {
  
  fpe <- filename_petab
  if (basename(dirname(fpe)) != "petab") stop("petab file must lie in folder called 'petab'")
  
  fpeParts <- petab_modelname_path(fpe)
  
  # fpd <-
  
}
