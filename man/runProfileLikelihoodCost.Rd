\name{runProfileLikelihoodCost}
\alias{runProfileLikelihoodCost}
\title{get the profile likelihood cost for a point in the parameter space for a single cell}
\usage{
runProfileLikelihoodCost(db = NULL, param = NULL, expData = NULL,
  fixed = c("kEVfrac", "kI"), p = param.FRAP.pairs, nStarts = 5)
}
\arguments{
  \item{db}{a data.frame with the data of a single cell
  (required if expData=NULL). 'time', 'FRAP' and nuclear
  and cytoplasmic fluorescence variables required.}

  \item{param}{a data.frame with the parameters for which
  the profile likelihood cost is wanted.}

  \item{expData}{a EXPDATA list as returned by
  \code{makeEXPDATA}, which contains the data for a cell
  (required if db=NULL).}

  \item{fixed}{a character vector with the names of the
  model variables that are to be kepet fixed during the
  profile likelihood cost calculation}

  \item{p}{list of options for the FRAP experiment}

  \item{nStarts}{number of random start points use by
  matlab's MultiStart function in order to get the global
  minimum}
}
\value{
a data.frame with the point in the parameter space and the
profile likelihood cost
}
\description{
get the profile likelihood cost for a point in the
parameter space for a single cell
}
\details{
wrapper over matlabs function 'profileLikelihoodCost.m'
}

