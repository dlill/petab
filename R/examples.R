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
  examplesAvailable <- list.files(system.file("petabExamples", package = "petab"))

  for (ex in examplesAvailable) {
    cat("\n======================================================", ex, "======================================================", sep ="\n")
    
    notes <- if (file.exists(petab_examplePath(ex, "notes"))) readLines(petab_examplePath(ex, "notes")) else ""
    cat(notes, sep = "\n")
    
    pe <- readPetab(petab_examplePath(ex, "pe"))
    cat("\n\n")
    petab_overviewObsPerCond(pe)
  }
  exampleNotes <- lapply(examplesAvailable, function(ex) {


    })
}

#' Get paths of an example
#'
#' @param exampleName (partial match of) example name
#' @param object "pe", "pd" or "dir"
#'
#' @return file path to the specified object
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @family Example functions
#'
#' @examples
#' exampleName <- "01"
#' petab_examplePath(exampleName)
#' petab_examplePath("01-EnzymeKinetics", "pd")
petab_examplePath <- function(exampleName, object = c("pe", "pd", "dir")[1]) {
  cat("Reading ", object, "...\n")
  fileEnding <- if (object == "pe") {"petab"
  } else if (object == "pd") {file.path("Compiled", "petab_indiv.rds")
  } else if (object =="dir") {NULL
  } else if (object == "notes") {"exampleNotes.txt"}

  examplesAvailable <- list.files(system.file("petabExamples", package = "petab"))
  exampleFile <- system.file(file.path("petabExamples", examplesAvailable), package = "petab")
  exampleFile <- exampleFile[grep(exampleName, exampleFile, fixed = TRUE)]
  if (length(exampleFile) > 1) stop("Ambiguous specification. Choose one of \n* ", paste0(exampleFile, collapse = "\n* "))
  exampleFile <- file.path(exampleFile, fileEnding)
  exampleFile
}

#' Get paths of an select example
#'
#' @param exampleName (partial match of) example name
#' @param object "pe", "pd" or "dir"
#'
#' @return file path to the specified object
#' @export
#' @author Svenja Kemmer
#' @md
#' @family Example functions
#'
#' @examples
#' exampleName <- "01"
#' petab_examplePath(exampleName)
#' petab_examplePath("01-EnzymeKinetics", "pd")
petab_select_examplePath <- function(exampleName, object = c("pe", "pd", "dir")[1]) {
  cat("Reading ", object, "...\n")
  # fileEnding <- if (object == "pe") {"petab"
  # } else if (object == "pd") {file.path("Compiled", "petab_indiv.rds")
  # } else if (object =="dir") {NULL
  # } else if (object == "notes") {"exampleNotes.txt"}
  
  examplesAvailable <- list.files(system.file("petab_select_examples", package = "petab"))
  exampleFile <- system.file(file.path("petab_select_examples", examplesAvailable), package = "petab")
  exampleFile <- exampleFile[grep(exampleName, exampleFile, fixed = TRUE)]
  if (length(exampleFile) > 1) stop("Ambiguous specification. Choose one of \n* ", paste0(exampleFile, collapse = "\n* "))
  exampleFile <- file.path(exampleFile, fileEnding)
  exampleFile
}

#' Read an example
#'
#' @inheritParams petab_examplePath
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



