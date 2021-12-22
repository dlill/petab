context("Compilation of examples")

examplesAvailable <- list.files(system.file("petabExamples", package = "petab"))

testthat::test_that("Compilation of Examples",{
ex <- (examplesAvailable)[[3]]
  for (ex in examplesAvailable){
    cat("\n================================\n", ex,"\n================================\n")
    
    # >>>> Ugly, but exclude Bachmann for now <<<<<<<<<<< ----
    if (ex == "03-Bachmann") next
    
    # Load old examples
    # pd0 <- petab_exampleRead(ex, "pd") # Does not need to be loaded
    # sim0 <- petab_exampleRead(ex, "simulation") # Truth needs to be saved as well, with base parameters
    # obj0 <- petab_exampleRead(ex, "objval")     # Truth needs to be saved as well, with base parameters
    
    # Compile
    tempd <- "temp"
    pd <- try({importPEtabSBML_indiv(filename = petab_examplePath(ex, "pe"), .compiledFolder = tempd, NFLAGcompile = 0)})
    unlink(tempd, recursive = TRUE) # Remove immediately: Not necessary to store data
    expect_true(!inherits(x= pd, what = "try-error"))
    
    # Simulate
    sim <- pd_predictAndPlot2(pd, FLAGreturnPlotData = TRUE)
    testthat::expect_true(all(names(sim) %in% c("dplot", "pplot", "pplotRibbon")))
    # Todo compare to truth: Difficulty is that the base simulation needs to be saved first
    
    # Objective function
    objval <- pd$obj(pd$pars)
    testthat::expect_true(class(objval) == "objlist")
    testthat::expect_true(!is.na(objval$value))
    # Todo compare to truth: Difficulty is that the base simulation needs to be saved first
    
  }})

