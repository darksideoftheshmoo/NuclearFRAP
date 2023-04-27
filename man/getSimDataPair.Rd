\name{getSimDataPair}
\alias{getSimDataPair}
\title{gets the simulated model for a pair and a set of parameters}
\usage{
getSimDataPair(db, param, p = param.FRAP.pairs)
}
\arguments{
  \item{db}{a data.frame with the data of a cell pair
  (required if expData=NULL). 'time', 'FRAP', 'type', and
  nuclear and cytoplasmic fluorescence variables required.}

  \item{param}{a data.frame with the parameters of the fit
  (required). Must have a variable for each free parameter
  of the model, and a 'type' variable.}

  \item{p}{list of options for the FRAP experiment}
}
\value{
a data.frame with the predicted values for the nuclear and
cytoplasmic fluorescence, and a 'type' variable
}
\description{
gets the simulated model for a pair and a set of parameters
}
\details{
wrapper over matlabs function 'nuclearFrapModel.m'
}

