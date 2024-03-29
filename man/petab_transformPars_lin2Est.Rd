% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PEtab.R
\name{petab_transformPars_lin2Est}
\alias{petab_transformPars_lin2Est}
\title{Title}
\usage{
petab_transformPars_lin2Est(pe, parsLin)
}
\arguments{
\item{parsLin}{}
}
\value{

}
\description{
Title
}
\examples{
pe <- petab_exampleRead("01")
parsEst <- petab_getPars_estScale(pe)
parsLin <- petab_getPars_linScale(pe)
parsLin <- parsLin[sample(1:length(parsLin), 5)] + rnorm(5)
}
\seealso{
Other Parameter wrangling: 
\code{\link{petab_getPars_estScale}()},
\code{\link{petab_getPars_linScale}()},
\code{\link{petab_mergeParameters}()},
\code{\link{petab_parameters_mergeParameters}()},
\code{\link{petab_setPars_estScale}()},
\code{\link{petab_setPars_linScale}()},
\code{\link{petab_transformPars_est2Lin}()}
}
\author{
Daniel Lill (daniel.lill@physik.uni-freiburg.de)
}
\concept{Parameter wrangling}
