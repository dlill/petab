% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dMod-import.R
\name{sbmlImport_getEventsFromTrueEvents}
\alias{sbmlImport_getEventsFromTrueEvents}
\title{Get Events from "events" entries in SBML}
\usage{
sbmlImport_getEventsFromTrueEvents(modelfile)
}
\arguments{
\item{modelfile}{}
}
\value{

}
\description{
The difficulty is that triggers can take arbitrary expressions which doesn't play 
well with dMod, which has a rather stiff input format for events.
Therefore, this function is tailored to my own export ONLY
}
