<!-- badges: start -->
[![Travis build status](https://travis-ci.com/dlill/petab.svg?branch=main)](https://travis-ci.com/dlill/petab)
[![codecov](https://codecov.io/gh/dlill/petab/branch/main/graph/badge.svg?token=HBQ0ERVEOK)](https://codecov.io/gh/dlill/petab)
<!-- badges: end -->
  
# Petab

Functionality for petab

* petab-handling itself, tool agnostic
* dMod interface

Notable functions

* `petab_mutateDCO` - consistently manipulate experimentalCondition, observables and measurementData at once in very expressive data.table language
* `petab_plotData`  - Plot your data with great liberty in customizing the plot, e.g. changing aesthetics and facets and sensible defaults. Supports pagination and asynchronous outputting so you don't wait forever until you can use your R-session again


# Version history

* New in 0.1.1
    * parameterFormulaInjection in petab$meta 
        * Allows to specify arbitrary parameter transformations after the estimation parameter trafo and before the dynamic model
        * The new prd0 looks like g*x*p1*p0 where p0 is the trafo for est-scales and p1 is the "injected trafo", e.g. a steady state trafo
        * This trafo is the same for all conditions to stay consistent in the indiv-framework

