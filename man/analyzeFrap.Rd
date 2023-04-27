\name{analyzeFrap}
\alias{analyzeFrap}
\title{run a basic analysis of a Nuclear FRAP for a single cell}
\usage{
analyzeFrap(db, p, nStarts = 5, theta = seq(0, (2 - 1/8) * pi, pi/8),
  profileLikelihood = TRUE, saveBasename = NULL, overwrite = FALSE)
}
\arguments{
  \item{db}{a data.frame with the data of a pair (required
  if expData=NULL). 'time', 'FRAP', 'type' and nuclear and
  cytoplasmic fluorescence variables required.}

  \item{p}{list of options for the FRAP experiment}

  \item{nStarts}{number of random start points use by
  matlab's MultiStart function in order to get the global
  minimum}

  \item{theta}{paramater to be passed to
  \code{\link{runProfileLikelihoodGridPair}}}

  \item{profileLikelihood}{boolean indicanting if the
  profile likelihood for kI and kEVfrac will be run}

  \item{saveBasename}{if not NULL, a charcater with the
  basename of the RData file where the results of the
  analysis will be saved.}

  \item{overwrite}{boolean indicating what to do if a
  .RData with the same name as the one to be created is
  found.  If overwrite=FALSE, the .RData is loaded and its
  content returned by the function, therefore avoiding to
  re-run the analysis.  If overwrite=TRUE, the analysis is
  redone and the file overwritten}
}
\value{
a object of type 'frap' containing all the results of the
analysis
}
\description{
run a basic analysis of a Nuclear FRAP for a single cell
}
\details{
this funcion runs a fit of the nuclear FRAP model to the
data, a profile likelihood analysis of parameters kI and
kEVfrac, and calculates the residuals of the model, for a
single cell.
}

