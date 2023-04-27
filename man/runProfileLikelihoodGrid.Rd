\name{runProfileLikelihoodGrid}
\alias{runProfileLikelihoodGrid}
\title{run the profile likelihhod analysis on a regular grid of points in the parameter space}
\usage{
runProfileLikelihoodGrid(db = NULL, param = NULL, expData = NULL,
  p = param.FRAP.pairs, fixed = p$FIXED.PARAMETERS, scan = c("kEVfrac",
  "kI"), absolute = FALSE, nStarts = 5, nGrid = 10,
  boundaryCost = qchisq(p = c(0.05), df = sum(!p$FIXED.PARAMETERS)),
  theta = seq(pi/4, 2 * pi, pi/4), rescale = TRUE)
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

  \item{p}{list of options for the FRAP experiment}

  \item{fixed}{a named interger vector specifying wich
  parameters are fixed (1) and which are free (0)}

  \item{scan}{either a character vector with a variables
  names to be scaned, or a named list with the values of
  the parameters with which to make a regular grid}

  \item{absolute}{boolean indicating if the values passed
  in scan should be treated as 'absolute' values, or
  increments with respect to the optimum fit passed in
  \code{param}}

  \item{nStarts}{number of random start points use by
  matlab's MultiStart function in order to get the global
  minimum}

  \item{nGrid}{number of points in the regular grid, when
  scan is a character vector.}

  \item{boundaryCost}{chi square cost increment between
  best fit and boundary. Used to estimate scan region when
  scan is a character vector}

  \item{theta}{numeric vector with angles (in radians) that
  define the directions for
  \code{runProfileLikelihoodContour}}

  \item{rescale}{boolean indicating if the coordinates
  should be rescale before setting the direction defined by
  theta. See \code{runProfileLikelihoodContour}.}
}
\value{
a data.frame with the points in the parameter space and the
profile likelihood cost
}
\description{
run the profile likelihhod analysis on a regular grid of
points in the parameter space
}
\details{
wrapper over matlabs function 'profileLikelihoodGrid.m'
}

