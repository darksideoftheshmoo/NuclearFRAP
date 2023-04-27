\name{analyzeFrapPairExperiment}
\alias{analyzeFrapPairExperiment}
\title{Analyze all positions of a Nuclear FRAP pair experiment}
\usage{
analyzeFrapPairExperiment(d, profileLikelihood = TRUE)
}
\arguments{
  \item{d}{a cell.data object with the relevant experiment}

  \item{profileLikelihood}{boolean indicanting if the
  profile likelihood for kI and kEVfrac will be run}
}
\description{
Analyze all positions of a Nuclear FRAP pair experiment
}
\details{
this high level function calls
\code{\link{analyzeFrapPair}} for all the pairs in the
experiemnt the output is saved to a new directery named as
the current experiment
}

