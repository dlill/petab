#' List all available examples
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Example functions
#' @examples
#' petab_exampleList()
petab_exampleList <- function() {
  # files <- list.files(system.file("petabExamples", package = "petab"), recursive = TRUE)
  # unique(dirname(files))
  list.files(system.file("petabExamples", package = "petab"))
}

#' Get paths of an example
#'
#' @param exampleName (partial match of) example name
#' @param object
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Example functions
#'
#' @examples
#' exampleName <- "01"
#' petab_examplePath(exampleName)
#' petab_examplePath("01-EnzymeKinetics", "pd)
petab_examplePath <- function(exampleName, object = c("pe", "pd", "dir")[1]) {
  cat("Reading ", object, "...\n")
  fileEnding <- if (object == "pe") {"petab"
  } else if (object == "pd") {file.path("Compiled", "petab_indiv.rds")
  } else if (object =="dir") {NULL}

  exampleFile <- system.file(file.path("petabExamples", petab_exampleList()), package = "petab")
  exampleFile <- exampleFile[grep(exampleName, exampleFile, fixed = TRUE)]
  if (length(exampleFile) > 1) stop("Ambiguous specification. Choose one of \n* ", paste0(exampleFile, collapse = "\n* "))
  exampleFile <- file.path(exampleFile, fileEnding)
  exampleFile
}

#' Read an example
#'
#' @return pd or pe
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Example functions
#' @examples
#' pd <- petab_exampleRead("01-EnzymeKinetics", "pd")
petab_exampleRead <- function(exampleName, object = c("pe", "pd")[1]) {
  filename <- petab_examplePath(exampleName, object)
  if (object=="pe"){
    return(readPetab(filename))
  } else if (object == "pd") {
      return(readPd(filename))
  }
}



