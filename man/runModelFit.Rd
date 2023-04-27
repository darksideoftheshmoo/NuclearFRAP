\name{runModelFit}
\alias{runModelFit}
\title{runs the model fit to the data, for a single cell}
\usage{
runModelFit(db = NULL, start.kEVfrac.kI = NULL, expData = NULL,
  p = param.FRAP.pairs, fixed = p$FIXED.PARAMETERS)
}
\arguments{
  \item{db}{a data.frame with the data of a single cell
  (required if expData=NULL). 'time', 'FRAP', nuclear and
  cytoplasmic fluorescence variables required.}

  \item{start.kEVfrac.kI}{a matrix of Nx2 containing the
  initial guesses for kEVfrac and kI}

  \item{expData}{a EXPDATA list as returned by
  \code{makeEXPDATA}, which contains the data for a cell
  (required if db=NULL).}

  \item{p}{list of options for the FRAP experiment}

  \item{fixed}{a named integer vector indicating if a given
  parameter is fixed (1) or not(0). Default is
  p$FIXED.PARAMETERS}
}
\value{
a data.frame of fited parameters
}
\description{
runs the model fit to the data, for a single cell
}
\details{
wrapper function over matlabs function 'modelFit.m'
}

