\name{getSimData}
\alias{getSimData}
\title{gets the simulated model for a single cell and a set of parameters}
\usage{
getSimData(db = NULL, param = NULL, expData = NULL,
  p = param.FRAP.pairs)
}
\arguments{
  \item{db}{a data.frame with the data of a single cell
  (required if expData=NULL). 'time', 'FRAP' and nuclear
  and cytoplasmic fluorescence variables required.}

  \item{param}{a data.frame with the parameters of the fit
  (required). Must have a variable for each free parameter
  of the model.}

  \item{expData}{a EXPDATA list as returned by
  \code{makeEXPDATA}, which contains the data for a cell
  (required if db=NULL).}

  \item{p}{list of options for the FRAP experiment}
}
\value{
a data.frame with the predicted values for the nuclear and
cytoplasmic fluorescence
}
\description{
gets the simulated model for a single cell and a set of
parameters
}
\details{
wrapper over matlabs function 'nuclearFrapModel.m'
}

