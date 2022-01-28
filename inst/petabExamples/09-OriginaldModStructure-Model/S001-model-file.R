library(petab)
library(dMod)
library(dplyr)
library(stringr)
library(conveniencefunctions)

try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
.currentWD <- getwd()
.combiledFolder <- "CompiledObjects"
dir.create(file.path(.currentWD,.combiledFolder), showWarnings = FALSE)

# -------------------------------------------------------------------------#
# Create enzyme kinetics model and data ----
# -------------------------------------------------------------------------#
modelname <- "dMod_model"

# .. 1 Reactions -----
reactions <- NULL %>% 
  addReaction(from = "R1+R2", to = "R1R2", rate = "k1 * R1 * R2 * ligand", description = "ligand induced R1R2 formation (NISS)") %>% 
  addReaction(from = "R1R2", to = "R1+R2", rate = "k2 * R1R2", description = "R1R2 decay (NISS)") %>% 
  addReaction(from = "ERK", to = "pERK", rate = "k3 * ERK * R1R2", description = "ERK to pERK (NISS)") %>% 
  addReaction(from = "ERK", to = "pERK", rate = "k4 * ERK", description = "ERK to pERK basal") %>% 
  addReaction(from = "pERK", to = "ERK", rate = "k5 * pERK", description = "pERK to ERK")
  
# reactions <- eqnlist_addDefaultCompartment(reactions, "cytoplasm") # Need compartment information for SBML

# .. 2 Observables -----
observables <- eqnvec(
  pERK_obs = "scale_pERK * pERK + offset_pERK"
)
# log-transform observables
observables <- as.eqnvec(paste("log10(", observables, ")"), names = names(observables))
  
# .. 3 Events -----
reactions <- addReaction(reactions, "0", "ligand", rate = "0", description = "stimulation (NISS)")

eventlist <- NULL %>% 
  addEvent(var = "ligand", time = 0, value = 1, method = "replace")

# .. 4 Errors -----
which_err <- c(1:length(observables))
errors <- paste("sigma", names(observables)[which_err], sep = "_")
names(errors) <- names(observables[which_err])

# .. 5 Parameters -----
# create subset of reactions describing the SS
reactions_SS <- subset(reactions, !grepl("(N)", Description))
setwd(.combiledFolder)
steadystates <- c(steadyStates(reactions_SS), R1R2 = 0)
setwd(.currentWD)

# steady state pars
modelpars <-  getParameters(reactions)
modelparsSS <-  getParameters(reactions_SS)
noSS_pars <- setdiff(modelpars, modelparsSS)

mySS_pars_fixed <- names(which(steadystates == "0"))
mySS_eqns_est <- steadystates[setdiff(names(steadystates), mySS_pars_fixed)]

# .. 6 Simulate data -----
parInfo <- data.table(tibble::tribble(
  ~parName, ~parValue, ~parUnit,
  "k1"   ,      0.5, "per_minute" , 
  "k2"  ,     0.1, "per_minute" ,
  "k3"  ,     0.4, "per_minute",
  "k4"  ,     0.1, "per_minute",
  "k5"  ,     0.2, "per_minute",
  "cytoplasm", 1, "volume"))

speciesInfo <- data.table(tibble::tribble(
  ~speciesName, ~compName, ~initialAmount,
  "R1"         ,"cytoplasm" ,             1,          # Amount, not concentration
  "R2"         ,"cytoplasm" ,             1,
  "ERK"          ,"cytoplasm" ,             1,          # Amount, not concentration
  "pERK"         ,"cytoplasm" ,             1,
  "R1R2"       ,"cytoplasm" ,             1,
  "ligand"     ,"cytoplasm" ,             1,
  ))

# compartmentInfo is left as the default getCompartmentInfo(el)
# unitInfo is left as the default getUnitInfo(): If you need other units, you need to add them

setwd(.combiledFolder)
compiled <- odemodel(f = reactions, modelname = modelname)
x <- Xs(compiled, condition = "C1", optionsOde = list(atol = 1e-12,rtol = 1e-12))
setwd(.currentWD)
mypars <- c(setNames(parInfo$parValue, parInfo$parName),
          setNames(speciesInfo$initialAmount, speciesInfo$speciesName))

pred <- x(seq(0,50, 1), mypars)
pred <- data.table(as.data.frame(pred))
pred <- rbind(pred,copy(pred) %>% .[,condition := "C2"],copy(pred) %>% .[,condition := "C3"])
pred <- pred[time %in% c(0,2,5,10,15,20,40)]
pred <- pred[name == "pERK"]
set.seed(1)
pred[,`:=`(sigma = 0.03)]
pred[,`:=`(value = exp(log(value) + rnorm(length(value), sd = sigma)))]
pred[condition == "C1", value := 1 + rnorm(length(value), sd = 0.01)]
pred[condition == "C2" & time != 0, value := value - 0.1]
pred[,`:=`(name = paste0(name, "_obs"))]
# ggplot(pred, aes(time, value, color = condition))+ geom_point() + geom_line()
mydata <- copy(pred) %>% as.data.frame() %>% as.datalist(split.by = "condition")

plotData(mydata)+geom_line()

# .. 7 condition.grid -----
condition.grid <- attr(mydata, "condition.grid")

# .. 8 trafo -----
innerpars <- setdiff(c(unique(c(getParameters(reactions), getSymbols(observables), getSymbols(errors)))), getCompartmentInfo(reactions)$compName)
trafo <- define(NULL, "x~y", x = innerpars, y = innerpars) %>% 
  insert("x~0", x = c("ligand")) %>% 
  insert("x~1", x = "cytoplasm") %>% 
  #SS
  insert("x~y", x = names(steadystates), y = steadystates) %>%
  # fixed sigma
  insert("x~0.1", x = "sigma_R1R2_obs") %>% 
  # model reduction
  insert("x~y", x = "k4", y = "k5 * k_ratio") %>% 
  # log trafo
  insert("x ~ exp(x)", x = .currentSymbols)

# .. 9 trafo list -----
trafoL <- branch(trafo, condition.grid)

trafoL <- insert(trafoL, "x~-1000", 
                 x = c("k1", "k2"), 
                 condition == "C1") %>% 
  insert("x~x_condi", 
         x = c("k1", "k2", "scale_pERK"), 
         condi = condition)


# .. 10 est.grid -----
est.grid <- getParGrids(trafo, trafoL, condition.grid, SS_pars = "pERK")[[1]]

# .. 11 fixed.grid -----
fixed.grid <- getParGrids(trafo, trafoL, condition.grid, SS_pars = "pERK")[[2]]
# for (nm in setdiff(names(fixed.grid), c("ID", "condition"))) fixed.grid[[nm]] <- as.numeric(fixed.grid[[nm]])
# .. 12 pouter -----
outerpars <- unlist(est.grid[, getSymbols(trafo)]) %>% unique()
# remove "dummy" :)
outerpars <- outerpars[outerpars != "dummy"]
pouter <- structure(rep(1, length(outerpars)), names = outerpars)

# .. 13 compile model -----
setwd(.combiledFolder)
tolerances <- 1e-8
x <- odemodel(reactions, 
              fixed = c("ligand"), 
              modelname = "x", 
              events = eventlist,
              compile = FALSE) %>% Xs(optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances), 
                                      optionsSens = list(method = "lsodes", rtol = tolerances, atol = tolerances))

g <- Y(observables, f = reactions, condition = NULL,
       compile = F, modelname = "g")

e <- Y(errors, f = c(as.eqnvec(reactions), observables), states = names(observables), 
       compile = F, modelname = "e", attach.input = FALSE)

p <- P(trafo, modelname = "p")

compile(x, g, p, e, output = modelname, cores  = detectFreeCores())
setwd(.currentWD)
# -------------------------------------------------------------------------#
# Fit model ----
# -------------------------------------------------------------------------#
mytimes <- seq(0, 40, 1)

prd0 <- (g*x*p)
prd <- cf_PRD_indiv(prd0, est.grid, fixed.grid)

# prediction <- prd(mytimes, pouter, FLAGbrowser = F)
# plotPrediction(prediction)

obj <- cf_normL2_indiv(mydata, prd0, e, est.grid, fixed.grid, times = mytimes) + 
  constraintL2(pouter, sigma = 12) # constraintL2(pouterFit, sigma = 12) + priorL2(pouterL1, lambda = "lambda")
# obj(pouter)

# .. mstrust -----
if(FALSE){
  nrfits <- 10
  nrcores <- 10
  out <- mstrust(objfun=obj, 
                 center=dMod::msParframe(pouter, n = nrfits, seed=47), 
                 studyname=modelname, 
                 rinit = 0.1, 
                 rmax = 10,
                 fits = nrfits, 
                 cores = nrcores, 
                 samplefun = "rnorm",
                 stats = FALSE, 
                 narrowing = NULL, 
                 iterlim=400, 
                 sd = 3)
  
  myfitlist <- as.parframe(out)
  bestfit <- as.parvec(myfitlist, 1)
  
  plotCombined(prd(mytimes, bestfit, FLAGbrowser = F), mydata, name == "pERK_obs")
}

# Exit ----
