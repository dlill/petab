% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dMod-import.R
\name{fitModelPEtabSBML}
\alias{fitModelPEtabSBML}
\title{Fit a model imported via importPEtabSBML}
\usage{
fitModelPEtabSBML(objfun = obj, nrfits = 4, nrcores = 4, useBounds = TRUE)
}
\arguments{
\item{objfun}{Objective function to be minimized as created by \link{importPEtabSBML}.}

\item{nrfits}{numeric, Number of fits to be performed}

\item{nrcores}{numeric, Number of cores to be used}

\item{useBounds}{boolean, if TRUE, parameter bounds are taken as provided in PEtab format, if FALSE no parameter bounds are applied}
}
\value{
parframe with the parameter estimated of the multi-start optimization
}
\description{
A wrapper function to use \link{mstrust} with imported PEtabSBML models. Some reasonable standard arguments for mstrust are used. Results of mstrust are written to Results folder.
}
\author{
Marcus Rosenblatt and Svenja Kemmer
}
