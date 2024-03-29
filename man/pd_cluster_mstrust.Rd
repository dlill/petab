% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dMod-pd.R
\name{pd_cluster_mstrust}
\alias{pd_cluster_mstrust}
\title{Fit model on cluster}
\usage{
pd_cluster_mstrust(
  pd = NULL,
  .outputFolder,
  n_startsPerNode = 16 * 3,
  n_nodes = 10,
  identifier = "mstrust",
  FLAGforcePurge = FALSE,
  opt.parameter_startpoints = "sample",
  iterlim = 500,
  passwdEnv = NULL,
  machine = "cluster"
)
}
\arguments{
\item{pd}{}

\item{.outputFolder}{}

\item{n_startsPerNode}{}

\item{n_nodes}{}

\item{id}{}

\item{type}{}
}
\value{
Characters
}
\description{
Fit model on cluster
}
\seealso{
Other Cluster: 
\code{\link{clusterStatusMessage}()},
\code{\link{pd_cluster_L1_fitUnbiasedEachMstrust}()},
\code{\link{pd_cluster_L1_fitUnbiasedEachOnce}()},
\code{\link{pd_cluster_L1_mstrust}()},
\code{\link{pd_cluster_L1}()},
\code{\link{pd_cluster_profile}()},
\code{\link{readPd}()}
}
\author{
Daniel Lill (daniel.lill@physik.uni-freiburg.de)
}
\concept{Cluster}
