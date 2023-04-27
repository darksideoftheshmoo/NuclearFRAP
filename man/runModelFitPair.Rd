\name{runModelFitPair}
\alias{runModelFitPair}
\title{runs the model fit to the data, for a cell pair}
\usage{
runModelFitPair(db, start.kEVfrac.kI = NULL, p = param.FRAP.pairs,
  fixed = p$FIXED.PARAMETERS)
}
\arguments{
  \item{db}{a data.frame with the data of a single pair.
  'time', 'FRAP', 'type', nuclear and cytoplasmic
  fluorescence variables for mother and daughter required}

  \item{p}{list of options for the FRAP experiment}

  \item{fixed}{a named integer vector indicating if a given
  parameter is fixed (1) or not(0). Default is
  p$FIXED.PARAMETERS}
}
\value{
a data.frame of fitted parameters, with a 'type' variable
}
\description{
runs the model fit to the data, for a cell pair
}
\details{
wrapper over matlabs function 'modelFit.m'
}

