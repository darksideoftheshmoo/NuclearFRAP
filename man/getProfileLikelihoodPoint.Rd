\name{getProfileLikelihoodPoint}
\alias{getProfileLikelihoodPoint}
\title{get a sample point of a profile likelihood analysis}
\usage{
getProfileLikelihoodPoint(PL.db, param, cost = NULL, point = NULL)
}
\arguments{
  \item{PL.db}{data.frame as returned by
  \code{\link{runProfileLikelihoodGrid}}}

  \item{cost}{target chi square cost. The returned point
  will be the one in PL.db, closest in cost to this
  parameter.}

  \item{point}{target point in paramter space. The returned
  point will be the one in PL.db closest to this point.}
}
\value{
a data.frame containing an element of PL.db
}
\description{
get a sample point of a profile likelihood analysis
}

