% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rw-seamove.R
\name{get.B.ts}
\alias{get.B.ts}
\title{function extracts the total biomass in the region before the movement. It is equal to the model state vector saved at previous time step and before ageing (subject to change once the outputs in SEAPODYM will be written after ageing). The biomass extracted from 'fluxes' files is extracted at quarterly time step (depends on the time unit of biomass flow rates). The variable 'sp' should be defined.}
\usage{
get.B.ts(dir, age, reg, years)
}
\arguments{
\item{dir}{is the full (or relative if working in local directory) path to the output directory.}

\item{age}{is the vector of ages, for which the extraction is to be done.}

\item{reg}{is the region, for which the biomass extraction is to be done.}

\item{years}{is the set of years, for which the extraction is to be done.}
}
\value{
Returns the list with dates and biomass time series.
}
\description{
function extracts the total biomass in the region before the movement. It is equal to the model state vector saved at previous time step and before ageing (subject to change once the outputs in SEAPODYM will be written after ageing). The biomass extracted from 'fluxes' files is extracted at quarterly time step (depends on the time unit of biomass flow rates). The variable 'sp' should be defined.
}
\examples{
dir<-"./output"; sp<-"skj"
dat<-get.B.ts(dir,3:47,1,1983:2010) # extracting the total biomass of skipjack
dat<-get.B.ts(dir,1:1,1,1983:2010)  # extracting the biomass of skipjack larvae
}
