# library(petab)
devtools::load_all("~/Promotion/Promotion/Projects/petab")
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# -------------------------------------------------------------------------#
# Load model ----
# -------------------------------------------------------------------------#
pd <- importPEtabSBML_indiv("petab", NFLAGcompile = 0, .compiledFolder = "Compiled", SFLAGbrowser = "0")

pd$prd(pd$times, pd$pars)
pd_predictAndPlot2(pd,i = time > 0.5)


# -------------------------------------------------------------------------#
# L1 Workflow Bernie ----
# -------------------------------------------------------------------------#

# * For parameter estimation, all three cell types CFU-E, H838 and H838-HA-hEPOR were initially implemented separately 
# * The parameters were estimated independently for CFU-E and H838 & H838-HA-hEPOR cells, and the corresponding fold-change was L 1 -regularized.
# * The re-parameterization to decouple modules of the signaling network was performed with CFU-E as reference point for H838 & H838-HA-hEPOR. 
# * For generating theregularization path, the regularization strength was scanned in log-space from 10 −4 to 10 4 .
# * Regularization was only used for selection of parameter differences, while for model selection the un-regularization solution was used for unbiased parameter estimates

# -------------------------------------------------------------------------#
# L1 Workflow new ----
# -------------------------------------------------------------------------#

# * Fit all, no L1
# * Run L1 such that each parameter goes to zero individually
#   * Write a routine with bisection algorithm: 
#     * lambda_increment = 1.41
#     * j = i + 1
#     * if nfixed_j == nfixed_i + 1 | nfixed_j == nfixed_i
#       * lambda_j = lambda_i * lambda_increment
#     * else if nfixed_j > nfixed_i + 1
#       * lambda_j = lambda_i * (1 / lambda_increment) * lambda_increment/2
# * For each model fit the model with all available parameters
# * Model selection with LRT

# .. Fit all, no L1 -----
iterlim = 1000
printIter = TRUE
traceFile = NULL
NFLAGsavePd = 1

parametersL1 <- cOde::getSymbols(unlist(pd$dModAtoms$gridlist$est.grid[,.SD,.SDcols = (paste0("L1_", pd$pe$meta$L1$parameterId_base))]))
fit_par <- pd$pars[setdiff(names(pd$pars), parametersL1)]
fit_fix <- c(pd$fixed, pd$pars[parametersL1])

parlower <- petab_getParameterBoundaries(pd$pe, "lower")[names(fit_par)]
parupper <- petab_getParameterBoundaries(pd$pe, "upper")[names(fit_par)]

fit <- dMod::trust(pd$obj, fit_par, 1,10, iterlim = iterlim, fixed = fit_fix,
                   parlower = parlower, parupper = parupper, printIter = printIter, traceFile = traceFile)
if (!fit$converged) warning("Fit not converged, please try increasing 'iterlim' (was ", iterlim,")")

pd <- pd_updateEstPars(pd, parsEst = fit$argument, FLAGupdatePE = TRUE, FLAGsavePd = NFLAGsavePd > 0)

pd

# .. Run L1 for all conditions -----
iterlim = 1000
printIter = TRUE
traceFile = NULL
NFLAGsavePd = 1

fit_par <- pd$pars
fit_fix <- pd$fixed
fit_L1 <- grep("^L1_", names(fit_par), value = TRUE)
fit_L1 <- setNames(rep(0, length(fit_L1)), fit_L1)

parlower <- petab_getParameterBoundaries(pd$pe, "lower")[names(fit_par)]
parupper <- petab_getParameterBoundaries(pd$pe, "upper")[names(fit_par)]


node <- 1
# lambdas <- 10^(seq(log10(0.0001), log10(1000000), length.out = 16))
lambdas <- 10^(seq(log10(10^4), log10(10^6), length.out = 16))

ncores <- 8
idx <- (1:16)[[16]]
fits <- parallel::mclapply(X = 1:16, mc.cores = ncores, FUN = function(idx) {
  lambda <- lambdas[idx]
  
  fit <- dMod::trustL1(
    objfun = pd$obj, parinit = pd$pars,
    mu = fit_L1, one.sided = FALSE, lambda = lambda,
    rinit = 0.1, rmax = 10, iterlim = 500,
    parupper = parupper,parlower = parlower
  )
  
  fit <- c(list(lambdaL1 = lambda), fit[c("value", "argument", "iterations", "converged")])
  dput(fit, file = sprintf("L1-%02i-%02i.R",node , idx))
  
  fit
})

# .... Load saved fit results ------
fits <- list.files(".", "^L1.*R$", full.names = T) %>% lapply(source, local = TRUE) %>% lapply(function(x) x$value)
fits <- lapply(fits, function(f) {data.table(as.data.table(f[setdiff(names(f), "argument")]), as.data.table(as.list(f$argument))) })
fits <- rbindlist(fits)
fits <- cf_parframe(fits, metanames = cf_parf_metaNames0$l1)

# .... Determine model candidates ------
# [ ] find better name, find better data structure to save the intermediate fits in the pd
fixed_L1 <- as.matrix(fits)
fixed_L1 <- fixed_L1[,grep("^L1_", colnames(fixed_L1)),drop=FALSE]
fixed_L1 <- fixed_L1 == 0
fixed_L1 <- unique(fixed_L1)

# .. Fit models -----
iterlim = 1000
printIter = TRUE
traceFile = NULL
NFLAGsavePd = 1

idx <- (seq_len(nrow(fixed_L1)))[[1]]
refits <- lapply(seq_len(nrow(fixed_L1)[-1]), function(idx) { # Last one is already fitted
  
  parametersFixed <- fixed_L1[idx,,drop = TRUE]
  parametersFixed <- names(parametersFixed)[parametersFixed]
  
  fit_par <- pd$pars[setdiff(names(pd$pars), parametersFixed)]
  fit_fix <- c(pd$fixed, pd$pars[parametersFixed])
  
  parlower <- petab_getParameterBoundaries(pd$pe, "lower")[names(fit_par)]
  parupper <- petab_getParameterBoundaries(pd$pe, "upper")[names(fit_par)]
  
  node <- 1
  fit <- dMod::trust(
    objfun = pd$obj, parinit = fit_par, fixed = fit_fix,
    rinit = 0.1, rmax = 10, iterlim = 500,
    parupper = parupper,parlower = parlower
  )
  
  fit <- fit[c("value", "argument", "iterations", "converged")]
  dput(fit, file = sprintf("L1-Refit-%02i-%02i.R",node , idx))
  
  fit
})

refits <- c(refits, list(fit[c("value", "argument", "iterations", "converged")]))

# .... Load saved fit results ------
rf <- lapply(refits, function(f) {data.table(as.data.table(f[setdiff(names(f), "argument")]), npars = length(f$argument), as.data.table(as.list(f$argument))) })
rf <- rbindlist(rf, fill = TRUE, use.names = TRUE)
rf[, (names(pd$pars)):=(lapply(.SD, function(x) replace(x, is.na(x),0))), .SDcols = names(pd$pars)]
# .. Determine maximum model -----
rf[,`:=`(nparsdiff = max(npars)-npars)]
rf[,`:=`(deltaChi = qchisq(0.95, nparsdiff))]
rf[,`:=`(objBound = min(value)+deltaChi)]
rf[,`:=`(rejectModel = value > objBound)]

# .. Use maximum model and get profiles for supervised removal of L1 parameters -----
maxModelRow <- max(which(rf$rejectModel == FALSE))

parametersFixed <- fixed_L1[maxModelRow,,drop = TRUE]
parametersFree <- names(parametersFixed)[!parametersFixed]
parametersFixed <- names(parametersFixed)[parametersFixed]

maxModel <- rf[maxModelRow]
maxModel <- parframe(as.data.frame(maxModel), parameters = names(pd$pars))

pars_MaxModel <- unclass_parvec(as.parvec(maxModel))
prof_par <- pars_MaxModel[setdiff(names(pars_MaxModel), parametersFixed)]
prof_fix <- c(pd$fixed, pars_MaxModel[parametersFixed])

profiles <-cf_profile(obj = pd$obj, pars = prof_par, whichPar = parametersFree, fixed = prof_fix)

# .. Fix final parameters found by false positives -----
profs <- lapply(profiles, confint)
profs <- rbindlist(profs)
profs[,`:=`(includesZero = sign(upper) * sign(lower) == -1)]



# -------------------------------------------------------------------------#
# OLD ----
# -------------------------------------------------------------------------#

# -------------------------------------------------------------------------#
# L1 Workflow ----
# -------------------------------------------------------------------------#

# * Fit base condition
# * Run L1 such that each parameter goes to zero individually
#   * Write a routine with bisection algorithm: 
#     * lambda_increment = 1.41
#     * j = i + 1
#     * if nfixed_j == nfixed_i + 1 | nfixed_j == nfixed_i
#       * lambda_j = lambda_i * lambda_increment
#     * else if nfixed_j > nfixed_i + 1
#       * lambda_j = lambda_i * (1 / lambda_increment) * lambda_increment/2
# * For each model fit the model with all available parameters
# * Model selection with LRT

# .. Fit base condition -----
iterlim = 1000
printIter = TRUE
traceFile = NULL
NFLAGsavePd = 1

conditions_Reference <- pd$pe$meta$L1$L1Spec[L1Spec == pd$pe$meta$L1$conditionSpecL1_reference, conditionId]
parameters_Reference <- cOde::getSymbols(unlist(pd$dModAtoms$gridlist$est.grid[condition %in% conditions_Reference, !"condition"]))
fit_par <- pd$pars[parameters_Reference]
fit_fix <- c(pd$fixed, pd$pars[setdiff(names(pd$pars), parameters_Reference)])

parlower <- petab_getParameterBoundaries(pd$pe, "lower")[names(fit_par)]
parupper <- petab_getParameterBoundaries(pd$pe, "upper")[names(fit_par)]

fit <- dMod::trust(pd$obj, fit_par, 1,10, iterlim = iterlim, fixed = fit_fix,
                   parlower = parlower, parupper = parupper, printIter = printIter, traceFile = traceFile,
                   conditions = conditions_Reference)
if (!fit$converged) warning("Fit not converged, please try increasing 'iterlim' (was ", iterlim,")")

pd <- pd_updateEstPars(pd, parsEst = fit$argument, FLAGupdatePE = TRUE, FLAGsavePd = NFLAGsavePd > 0)

pd

pd_predictAndPlot2(pd, i = time > 0.05)

# .. Run L1 for all conditions -----
iterlim = 1000
printIter = TRUE
traceFile = NULL
NFLAGsavePd = 1

conditions_Reference <- pd$pe$meta$L1$L1Spec[L1Spec == pd$pe$meta$L1$conditionSpecL1_reference, conditionId]
parameters_Reference <- cOde::getSymbols(unlist(pd$dModAtoms$gridlist$est.grid[condition %in% conditions_Reference, !"condition"]))
fit_par <- pd$pars[parameters_Reference]
fit_fix <- c(pd$fixed, pd$pars[setdiff(names(pd$pars), parameters_Reference)])

parlower <- petab_getParameterBoundaries(pd$pe, "lower")[names(fit_par)]
parupper <- petab_getParameterBoundaries(pd$pe, "upper")[names(fit_par)]

fit <- dMod::trust(pd$obj, fit_par, 1,10, iterlim = iterlim, fixed = fit_fix,
                   parlower = parlower, parupper = parupper, printIter = printIter, traceFile = traceFile,
                   conditions = conditions_Reference)
if (!fit$converged) warning("Fit not converged, please try increasing 'iterlim' (was ", iterlim,")")

pd <- pd_updateEstPars(pd, parsEst = fit$argument, FLAGupdatePE = TRUE, FLAGsavePd = NFLAGsavePd > 0)

pd


# Exit ----
