% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_varDYM.R
\name{write.var.dym}
\alias{write.var.dym}
\title{Writing DYM file}
\usage{
write.var.dym(file.out, tvect, x, y, mask, data, verbose = TRUE)
}
\arguments{
\item{file.out}{is the name of the DYM file.}

\item{tvect}{is the date vector in decimal format.}

\item{x}{is the vector of longitude.}

\item{y}{is the vector of latitude.}

\item{mask}{is the land mask - the matrix of integers 0 - 2, with 0 for land, 1 for the sea with only one layer (nearshore), 2 for the sea with two layers and 3 for the sea with three pelagic layers though the depth column.}

\item{data}{is the 3d array of data to be written in DYM file.}

\item{verbose}{flag (default value is TRUE) controlling the prompt for the nominal function execution.}
}
\description{
Writes data to DYM file
}
\examples{
write.var.dym("sst.dym",tt,x,y,mask,var);
}
