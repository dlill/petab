% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dMod-pd.R
\name{pd_plotProfilePaths}
\alias{pd_plotProfilePaths}
\title{Plot affected paths of profiles}
\usage{
pd_plotProfilePaths(
  profiles,
  tol = 0.1,
  FLAGnormalizeYParameters = TRUE,
  nrow = 3,
  ncol = 4,
  ggCallback = list(),
  ...
)
}
\arguments{
\item{profiles}{profiles to plotPaths}

\item{tol}{}

\item{page}{pagination page}
}
\value{
paginated ggplot
}
\description{
Plot affected paths of profiles
}
\seealso{
Other plotting: 
\code{\link{pd_plotParsParallelLines}()},
\code{\link{pd_plotProfile}()},
\code{\link{petab_plotData}()},
\code{\link{petab_plotHelpers_aeslist}()},
\code{\link{petab_plotHelpers_meanMeasurementsValue}()},
\code{\link{petab_plotHelpers_parameterOrder}()},
\code{\link{petab_plotHelpers_variableOrder}()},
\code{\link{stepMax}}
}
\author{
Daniel Lill (daniel.lill@physik.uni-freiburg.de)
}
\concept{plotting}
