\name{getSimulatedNoiseData}
\alias{getSimulatedNoiseData}
\title{gets the simulated model with added noise for a single cell and a set of parameters}
\usage{
getSimulatedNoiseData(db, p, param = NULL, sigmaFactor = 1)
}
\arguments{
  \item{db}{a data.frame with the data of a single cell
  (required if expData=NULL). 'time', 'FRAP' and nuclear
  and cytoplasmic fluorescence variables required.}

  \item{p}{list of options for the FRAP experiment}

  \item{param}{a data.frame with the parameters of the fit.
  If missing \code{runModelFit} is called.}

  \item{sigmaFactor}{multiplying factor to multiplicate the
  sigma before generating the artificial residuals}
}
\value{
a data.frame, as required by runModelFit and other
functions
}
\description{
gets the simulated model with added noise for a single cell
and a set of parameters
}
\details{
use this function to get simulated data ready to be
re-fitted
}

