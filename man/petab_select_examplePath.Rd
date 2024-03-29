% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/examples.R
\name{petab_select_examplePath}
\alias{petab_select_examplePath}
\title{Get paths of an select example}
\usage{
petab_select_examplePath(exampleName, object = c("pe", "pd", "dir")[1])
}
\arguments{
\item{exampleName}{(partial match of) example name}

\item{object}{"pe", "pd" or "dir"}
}
\value{
file path to the specified object
}
\description{
Get paths of an select example
}
\examples{
exampleName <- "01"
petab_examplePath(exampleName)
petab_examplePath("01-EnzymeKinetics", "pd")
}
\seealso{
Other Example functions: 
\code{\link{petab_exampleList}()},
\code{\link{petab_examplePath}()},
\code{\link{petab_exampleRead}()}
}
\author{
Svenja Kemmer
}
\concept{Example functions}
