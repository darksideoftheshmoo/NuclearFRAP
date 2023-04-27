\name{analyzeFrapExperiment}
\alias{analyzeFrapExperiment}
\title{Analyze all positions of a Nuclear FRAP experiment}
\usage{
analyzeFrapExperiment(d, profileLikelihood = TRUE)
}
\arguments{
  \item{d}{a cell.data object with the relevant experiment}

  \item{profileLikelihood}{boolean indicanting if the
  profile likelihood for kI and kEVfrac will be run}
}
\description{
Analyze all positions of a Nuclear FRAP experiment
}
\details{
this high level function calls \code{\link{analyzeFrap}}
for all the cells in the experiemnt the output is saved to
a new directery named as the current experiment
}

