# -------------------------------------------------------------------------#
# CF Functions ----
# -------------------------------------------------------------------------#

cf_PRD_indiv <- function(prd0, est.grid, fixed.grid) {
  prd <- function(times, pars, fixed = NULL, deriv = FALSE, conditions = est.grid$condition,
                  FLAGbrowser = FALSE,
                  FLAGverbose = FALSE,
                  FLAGrenameDerivPars = FALSE
  ) {
    out <- lapply(setNames(nm = conditions), function(cn) {
      if (FLAGbrowser) browser()
      ID <- est.grid$ID[est.grid$condition == cn]
      if (FLAGverbose) cat(ID, cn, "\n", sep = " ---- ")
      dummy <- cf_make_pars(pars, fixed, est.grid, fixed.grid, ID)
      pars_ <- dummy$pars
      fixed_ <- dummy$fixed
      if (length(setdiff(getParameters(prd0), names(c(pars_, fixed_)))))
        stop("The following parameters are missing: ", paste0(setdiff(getParameters(prd0), names(c(pars_, fixed_))), collapse = ", "))
      pred0 <-try(prd0(times, pars_, fixed = fixed_, deriv = deriv, conditions = NULL)[[1]])
      if (inherits(pred0, "try-error")) {
        browser()
        # Try this code to debug your model
        # 1 Parameters
        pinner <- p(pars_, fixed = fixed_)
        compare(names(pinner[[1]]), getParameters(x)) #setdiff(y,x) should be empty!
        # 2 ode-model
        pinner_test <- setNames(runif(length(getParameters(x))),getParameters(x))
        x(times, pinner_test, deriv = FALSE)
      }
      if (deriv && FLAGrenameDerivPars) pred0 <- cf_renameDerivPars(pred0, pars, est.grid, cn)
      pred0
    })
    dMod::as.prdlist(out)
  }
  class(prd) <- c("prdfn", "fn")
  prd
}

cf_make_pars <- function(pars, fixed = NULL, est.grid, fixed.grid, ID){
  if ("dummy" %in% names(pars))
    stop("'dummy' should not appear in est.vec (parameter vector passed to objective function)\n")
  pars        <- unclass(pars) # changed from unclass_parvec
  fixed       <- unclass(fixed) # changed from unclass_parvec
  pars_outer  <- pars
  fixed_outer <- fixed
  
  pars <- c(pars, fixed)
  pars <- c(pars, dummy = 1)
  parnames  <- unlist(est.grid[est.grid$ID == ID, setdiff(names(est.grid), c("ID", "condition"))])
  # remove dummy
  parnames <- parnames[parnames != "dummy"]
  pars <- setNames(pars[parnames], names(parnames))
  fixed <- unlist(fixed.grid[fixed.grid$ID == ID, setdiff(names(fixed.grid), c("ID", "condition"))])
  # remove NAs
  fixed <- fixed[!is.na(fixed)]
  fixed <- c(fixed, pars[parnames %in% names(fixed_outer)])
  # fixednames <- names(fixed)
  # if(logTransform){
  # doesn't work for log(0) = -Inf
  #   logfixed <- paste0("log(", fixed, ")")
  #   eval_vec <-  Vectorize(eval.parent, vectorize.args = "expr")
  #   logfixed <- eval_vec(parse(text = logfixed))
  #   names(logfixed) <- fixednames
  #   fixed <- logfixed
  # }
  pars <- pars[!parnames %in% names(fixed_outer)]
  parnames <- parnames[!parnames %in% names(fixed_outer)]
  return(list(pars = unlist(pars), fixed = unlist(fixed), parnames = parnames))
}

cf_normL2_indiv <- function (data, prd0, errmodel = NULL, est.grid, fixed.grid, times = NULL, attr.name = "data", fixed.conditions = NULL) {
  timesD <- sort(unique(c(0, do.call(c, lapply(data, function(d) d$time)))))
  if (!is.null(times))
    timesD <- sort(union(times, timesD))
  x.conditions <- est.grid$condition
  data.conditions <- names(data)
  e.conditions <- names(attr(errmodel, "mappings"))
  controls <- list(times = timesD, attr.name = attr.name, conditions = intersect(x.conditions,
                                                                                 data.conditions))
  force(errmodel)
  force(fixed.grid)
  force(est.grid)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = controls$conditions, simcores = 1,
                   FLAGbrowser = FALSE,
                   FLAGbrowser2 = FALSE,
                   FLAGverbose = FALSE,
                   FLAGNaNInfwarnings = FALSE,
                   FixedConditions = fixed.conditions) {
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    
    pars <- arglist[[1]]
    calc_objval <- function(cn) {
      
      if (FLAGbrowser) browser()
      
      ID <- est.grid$ID[est.grid$condition == cn]
      if (FLAGverbose) cat(ID, cn, "\n", sep = " ---- ")
      dummy <- cf_make_pars(pars, fixed, est.grid, fixed.grid, ID)
      pars_ <- dummy$pars
      fixed_ <- dummy$fixed
      
      timesD <- controls$times
      attr.name <- controls$attr.name
      
      myderiv <- deriv
      if (!is.null(FixedConditions) && cn %in% FixedConditions) myderiv <- FALSE
      
      prediction <- try(prd0(times = timesD, pars = pars_, fixed = fixed_, deriv = myderiv))
      
      if (inherits(prediction, "try-error"))
        stop("Prediction failed in condition = ", cn, ", ID = ", ID, ".
             Try iterating p(pars), (x*p)(pars), ... to find the problem.")
      prediction <- prediction[[1]]
      
      # [] refactor: put the following stuff into own function catch_nonproblematicNanInfs(prediciton, data, cn, FLAGNaNInfWarnings)
      whichcols <- nm <- NULL
      if (any(is.na(prediction))){
        whichcols <- unique(which(is.na(prediction), arr.ind = TRUE)[,2])
        nm <- colnames(prediction)[whichcols]
        
        if (length(intersect(data[[cn]]$name, nm)))
          stop("Prediction is.na for observables present in data in condition ", cn, "\n",
               "The following observables are affected: ", paste0(intersect(data[[cn]]$name, nm), collapse = ", "))
        
        if (FLAGNaNInfwarnings)
          warning("NaN in condition ", cn , " for the following names: ", paste0(nm, collapse = ", "))
        prediction[is.na(prediction)] <- 0
        attr(prediction, "deriv")[is.infinite(attr(prediction, "deriv"))|is.na(attr(prediction, "deriv"))] <- 0
        attr(prediction, "sensitivities")[is.infinite(attr(prediction, "sensitivities"))|is.na(attr(prediction, "sensitivities"))] <- 0
      }
      if (any(is.infinite(prediction))){
        whichcols <- unique(which(is.infinite(prediction), arr.ind = TRUE)[,2])
        nm <- colnames(prediction)[whichcols]
        
        if (length(intersect(data[[cn]]$name, nm)))
          stop("Prediction is infinite for observables present in data in condition ", cn, "\n",
               "The following observables are affected: ", paste0(intersect(data[[cn]]$name, nm), collapse = ", "))
        
        if (FLAGNaNInfwarnings)
          warning("Inf in condition ", cn , " for the following names: ", paste0(nm, collapse = ", "))
        
        prediction[is.infinite(prediction)] <- 0
        attr(prediction, "deriv")[is.infinite(attr(prediction, "deriv"))|is.na(attr(prediction, "deriv"))] <- 0
        attr(prediction, "sensitivities")[is.infinite(attr(prediction, "sensitivities"))|is.na(attr(prediction, "sensitivities"))] <- 0
      }
      
      err <- NULL
      if (any(is.na(data[[cn]]$sigma))) {
        err <- errmodel(out = prediction, pars = getParameters(prediction), conditions = cn, deriv=myderiv)
        mywrss <- nll(res(data[[cn]], prediction, err[[1]]), deriv = deriv, pars = pars)
      } else {
        mywrss <- nll(res(data[[cn]], prediction), deriv = deriv, pars = pars)
      }
      if (myderiv) {
        if (FLAGbrowser2) browser()
        mywrss$gradient <- mywrss$gradient[names(dummy$parnames)]
        names(mywrss$gradient) <- unname(dummy$parnames)
        
        mywrss$hessian <- mywrss$hessian[names(dummy$parnames),names(dummy$parnames)]
        dimnames(mywrss$hessian) <- list(unname(dummy$parnames), unname(dummy$parnames))
      }
      
      # [] catch conditions with NA value, don't include them in obj-calculation and print out warning
      return(mywrss)
    }
    
    if (simcores == 1)
      objlists <- lapply(setNames(nm = conditions), calc_objval)
    if (simcores > 1)
      objlists <- parallel::mclapply(setNames(nm = conditions), calc_objval, mc.cores = simcores)
    
    out <- objlist(value = 0,
                   gradient = structure(rep(0, length(names(pars))), names = names(pars)),
                   hessian = matrix(0, nrow = length(names(pars)), ncol = length(names(pars)), dimnames = list(names(pars), names(pars))))
    out$value <- do.call(sum, lapply(objlists, function(.x) .x$value))
    if (deriv) {
      for (gr in lapply(objlists, function(.x) .x$gradient))
        out$gradient[names(gr)] <- out$gradient[names(gr)] + gr
      for (hs in lapply(objlists, function(.x) .x$hessian))
        out$hessian[rownames(hs), colnames(hs)] <- out$hessian[rownames(hs), colnames(hs)] + hs
    }
    
    # consider fixed: return only derivs wrt pouter
    out$gradient <- out$gradient[names(pars)]
    out$hessian <- out$hessian[names(pars), names(pars)]
    
    attr(out, controls$attr.name) <- out$value
    attr(out, "condition_obj") <- vapply(objlists, function(.x) .x$value, 1)
    attr(out, "AIC") <- out$value + length(pars) * 2
    attr(out, "BIC") <- out$value + length(pars) * log(nrow(as.data.frame(data)))
    return(out)
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- data.conditions
  attr(myfn, "parameters") <- attr(prd0, "parameters")
  attr(myfn, "modelname") <- modelname(prd0, errmodel)
  return(myfn)
}

# Exit ----
