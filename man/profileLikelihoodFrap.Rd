\name{profileLikelihoodFrap}
\alias{profileLikelihoodFrap}
\title{run a profile likelihood analysis for FRAP data}
\usage{
profileLikelihoodFrap(db, p, vars = c("log10_kEVfrac", "log10_kI"),
  param = NULL, pairwise = FALSE, nGrid = 10, theta = seq(pi/8, 2 * pi,
  pi/8), nStarts = 5, saveBasename = NULL, overwrite = FALSE)
}
\arguments{
  \item{db}{a data.frame with the data of a pair (required
  if expData=NULL). 'time', 'FRAP' and nuclear and
  cytoplasmic fluorescence variables required.}

  \item{p}{list of options for the FRAP experiment}

  \item{param}{a data.frame with the best fits parameters}

  \item{pairwise}{boolean indicating if pairwise profile
  likelihood contours should be done}

  \item{nStarts}{number of random start points use by
  matlab's MultiStart function in order to get the global
  minimum}

  \item{theta}{paramater to be passed to
  \code{\link{runProfileLikelihoodContour}}}

  \item{nGrid}{number of points for the 1D analysis for
  each var}

  \item{saveBasename}{if not NULL, a charcater with the
  basename of the RData file where the results of the
  profile likelihood analysis will be saved.}

  \item{overwrite}{boolean indicating what to do if a
  .RData with the same name as the one to be created is
  found.  If overwrite=FALSE, the .RData is loaded and its
  content returned by the function, therefore avoiding to
  re-run the analysis.  If overwrite=TRUE, the analysis is
  redone and the file overwritten}
}
\description{
run a profile likelihood analysis for FRAP data
}

