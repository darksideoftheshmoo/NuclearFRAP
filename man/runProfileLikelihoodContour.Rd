\name{runProfileLikelihoodContour}
\alias{runProfileLikelihoodContour}
\title{run the profile likelihhod analysis on 'fix directions' from the fitted point}
\usage{
runProfileLikelihoodContour(db = NULL, param = NULL, expData = NULL,
  fixed = p$FIXED.PARAMETERS, scan = c("kEVfrac", "kI"),
  boundaryCost = qchisq(p = c(0.05), df = sum(!p$FIXED.PARAMETERS)),
  p = param.FRAP.pairs, nStarts = 5, theta = seq(pi/4, 2 * pi, pi/4),
  rescale = TRUE)
}
\arguments{
  \item{scan}{a character vector of length 2 with the names
  of the parameters that define the plane for the profile
  likelihood contour.}

  \item{boundaryCost}{chi square cost increment between
  best fit and boundary.}

  \item{theta}{numeric vector with angles (in radians) that
  define the directions where the points of cost = min.cost
  + boundaryCost will be found}

  \item{rescale}{boolean indicating if the coordinates
  should be rescale before setting the direction defined by
  theta. The rescaling is x'=x/x0, where x0 is the value in
  param.}

  \item{db}{a data.frame with the data of a single cell
  (required if expData=NULL). 'time', 'FRAP' and nuclear
  and cytoplasmic fluorescence variables required.}

  \item{param}{a data.frame with the parameters for which
  the profile likelihood cost is wanted.}

  \item{expData}{a EXPDATA list as returned by
  \code{makeEXPDATA}, which contains the data for a cell
  (required if db=NULL).}

  \item{fixed}{a named interger vector specifying wich
  parameters are fixed (1) and which are free (0)}

  \item{p}{list of options for the FRAP experiment}

  \item{nStarts}{number of random start points use by
  matlab's MultiStart function in order to get the global
  minimum}
}
\value{
a data.frame with the points along the directions defined
by theta where cost = min.cost + boundaryCost
}
\description{
run the profile likelihhod analysis on 'fix directions'
from the fitted point
}

