% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dMod-import.R
\name{updateParscalesToBaseTrafo}
\alias{updateParscalesToBaseTrafo}
\title{Rename parscales to the names needed in the base trafo}
\usage{
updateParscalesToBaseTrafo(scales_outer, est.grid)
}
\arguments{
\item{est.grid}{data.table}

\item{parscales}{setNames(PETABPars$parameterScale, PETABPars$parameterId)}
}
\value{

}
\description{
Rename parscales to the names needed in the base trafo
}
\examples{
# Example 1: k1_A and k1_B are duplicated to two inner parameters 
est.grid <- data.frame(ID = 1:2,
                       condition = c("A", "B"),
                       k1 = c("k1_A", "k1_B"),
                       k1DUPE = c("k1_A", "k1_B"),
                       k3 = c("k3outer", NA),
                       stringsAsFactors = FALSE)
scales_outer <- c(k1_A = "log", k1_B = "log", k3outer = "lin")
updateParscalesToBaseTrafo(scales_outer, est.grid)

# Example 2: SHOULD FAIL k1_A and k1_B map to same inner parameter, but have idfferent scales
est.grid <- data.frame(ID = 1:2,
                       condition = c("A", "B"),
                       k1 = c("k1_A", "k1_B"),
                       k1DUPE = c("k1_A", "k1_B"),
                       k3 = c("k3outer", NA),
                       stringsAsFactors = FALSE)
scales_outer <- c(k1_A = "log", k1_B = "log10", k3outer = "lin")
updateParscalesToBaseTrafo(scales_outer, est.grid)

# Example 3: k4 is a parameter not in est.grid. It might be a parameter in fix.grid and should be returned as is
est.grid <- data.frame(ID = 1:2,
                       condition = c("A", "B"),
                       k1 = c("k1_A", "k1_B"),
                       k1DUPE = c("k1_A", "k1_B"),
                       k3 = c("k3outer", NA),
                       stringsAsFactors = FALSE)
scales_outer <- c(k1_A = "log", k1_B = "log", k3outer = "lin", k4 = "log")
updateParscalesToBaseTrafo(scales_outer, est.grid)
}
\author{
Daniel Lill (daniel.lill@physik.uni-freiburg.de)
}
