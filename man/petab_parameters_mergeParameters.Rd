% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PEtab.R
\name{petab_parameters_mergeParameters}
\alias{petab_parameters_mergeParameters}
\title{Title}
\usage{
petab_parameters_mergeParameters(
  pe_pa1,
  pe_pa2,
  mergeCols = setdiff(intersect(names(pe_pa1), names(pe_pa2)), "parameterId")
)
}
\arguments{
\item{pe_pa1}{parameters_df to merge values into}

\item{pe_pa2}{parameters_df to merge values from}

\item{mergeCols}{Columns you want to update in pe_pa1}
}
\value{
updated pe_pa1
}
\description{
Title
}
\seealso{
Other Parameter wrangling: 
\code{\link{petab_getPars_estScale}()},
\code{\link{petab_getPars_linScale}()},
\code{\link{petab_mergeParameters}()},
\code{\link{petab_setPars_estScale}()},
\code{\link{petab_setPars_linScale}()},
\code{\link{petab_transformPars_est2Lin}()},
\code{\link{petab_transformPars_lin2Est}()}
}
\author{
Daniel Lill (daniel.lill@physik.uni-freiburg.de)
}
\concept{Parameter wrangling}
