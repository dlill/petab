% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/petab_select.R
\name{pd_petabSelect_reportYaml}
\alias{pd_petabSelect_reportYaml}
\title{Create the reportYaml from a PD}
\usage{
pd_petabSelect_reportYaml(pd, FLAGwriteYaml = TRUE)
}
\arguments{
\item{pd}{a pd}

\item{FLAGwriteYaml}{TRUE: Writes reportYaml into the original petab folder}
}
\value{
reportYaml as list, see \link{petabSelect_reportYaml}
}
\description{
Create the reportYaml from a PD
}
\examples{
# create reportYaml
pd <- petab_exampleRead("01", "pd")
reportYaml <- pd_petabSelect_reportYaml(pd, FLAGwriteYaml = FALSE)

# write it to disk and look at it
tf <- tempfile(fileext = ".yaml")
petabSelect_writeReportYaml(reportYaml, tf)
file.edit(tf)
}
\seealso{
Other petab select: 
\code{\link{as.data.table.reportYaml}()},
\code{\link{pd_petabSelect_collectEstimatedParameters}()},
\code{\link{pd_petabSelect_collectParameters}()},
\code{\link{pd_petabSelect_collectReportYamlCriteria}()},
\code{\link{petabSelect_LRT}()},
\code{\link{petabSelect_collectReportYamlCriteria}()},
\code{\link{petabSelect_compareModels}()},
\code{\link{petabSelect_readReportYaml}()},
\code{\link{petabSelect_reportYaml}()},
\code{\link{petabSelect_writeReportYaml}()}
}
\author{
Daniel Lill (daniel.lill@physik.uni-freiburg.de)
}
\concept{petab select}
