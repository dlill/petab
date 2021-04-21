
#' Determine if intput file is newer than output file
#'
#' @param filename_input,filename_output file path
#'
#' @return TRUE if inputfile is newer than outputfile or outputfile does not exist
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
inputFileChanged <- function(filename_input, filename_output) {

  if (!file.exists(filename_output)) return(TRUE)

  outputFileTimeChanged <- file.info(filename_output)[,"mtime"]
  inputFileTimeChanged  <- file.info(filename_input)[,"mtime"]

  inputFileTimeChanged > outputFileTimeChanged
}
