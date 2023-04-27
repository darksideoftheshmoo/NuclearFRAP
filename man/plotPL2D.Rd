\name{plotPL2D}
\alias{plotPL2D}
\title{create a 2D plot of a profile likelihood boundary}
\usage{
plotPL2D(PL2D.db, x, y, param, signifCost, color = "type", alpha = 1,
  tauBoundary = TRUE, posLabel = color == "pos")
}
\arguments{
  \item{PL2D.db}{data.frame as returned by
  \code{\link{runProfileLikelihoodGrid}} or
  \code{\link{runProfileLikelihoodContour}}}

  \item{x}{string with name of model parameter to use in
  the x axis}

  \item{y}{string with name of model parameter to use in
  the y axis}

  \item{param}{a data.frame with the best fits parameters}

  \item{signifCost}{chi square cost significant difference,
  defines boundary for
  \code{\link{runProfileLikelihoodGrid}} datasets}

  \item{color}{string indicating wich variable map to the
  aesthetic color}

  \item{alpha}{transparency level of the contours}

  \item{tauBoundary}{boolean indicating if tau boundaries
  should be drawn}

  \item{posLabel}{boolean indicating if position number
  must be used as label}
}
\value{
a ggplot object with the required plot.
}
\description{
create a 2D plot of a profile likelihood boundary
}

