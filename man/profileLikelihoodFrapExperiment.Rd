\name{profileLikelihoodFrapExperiment}
\alias{profileLikelihoodFrapExperiment}
\title{run profile likelihood analysis for all pairs in a experiment}
\usage{
profileLikelihoodFrapExperiment(d, vars = d$parameters$MODEL.PARAMETERS)
}
\arguments{
  \item{d}{a cell.data object with the experiment to
  analyze}

  \item{vars}{a character vector with the names of the
  model parameters that are going to be used for the
  profile likelihood analysis}
}
\description{
run profile likelihood analysis for all pairs in a
experiment
}
\details{
runs \code{\link{profileLikelihoodFrapPair}} in pairwise
mode for all pairs in the experiment
}
\note{
this function can take very long to run. In my Intel i5
quadcore a 10 cells experiment runs overnight (~12 hours)
Saves results to disk in .RData files, in the directory
named as d$experiment
}

