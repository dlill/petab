% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SBML-Export.R
\name{sbml_addOneEvent}
\alias{sbml_addOneEvent}
\title{Title}
\usage{
sbml_addOneEvent(model, eventSpecies, eventTrigger, eventFormula)
}
\arguments{
\item{model}{}

\item{eventSpecies}{}

\item{eventTrigger}{}

\item{eventFormula}{}
}
\value{

}
\description{
Rather generic already with eventTrigger and eventFormula,
but I don't know if import can handle anything else but "replace"-events at the moment
}
\examples{
eventTrigger <-  "time >= 2"
}
\seealso{
Other SBML export: 
\code{\link{eqnlist_addDefaultCompartment}()},
\code{\link{getCompartmentInfo}()},
\code{\link{getParInfo}()},
\code{\link{getReactionInfo}()},
\code{\link{getSpeciesInfo}()},
\code{\link{getUnitInfo}()},
\code{\link{sbml_addOneCompartment}()},
\code{\link{sbml_addOneParameter}()},
\code{\link{sbml_addOneReaction}()},
\code{\link{sbml_addOneSpecies}()},
\code{\link{sbml_addOneUnit}()},
\code{\link{sbml_exportEquationList}()},
\code{\link{sbml_initialize}()},
\code{\link{sbml_kineticLawAddParameters}()},
\code{\link{sbml_reactionAddKineticLaw}()},
\code{\link{sbml_reactionAddModifiers}()},
\code{\link{sbml_reactionAddProducts}()},
\code{\link{sbml_reactionAddReactants}()},
\code{\link{sbml_validateSBML}()}
}
\author{
Daniel Lill (daniel.lill@physik.uni-freiburg.de)
}
\concept{SBML export}
