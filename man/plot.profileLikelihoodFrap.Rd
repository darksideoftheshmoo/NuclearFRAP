\name{plot.profileLikelihoodFrap}
\alias{plot.profileLikelihoodFrap}
\title{plot method for profileLikelihoodFrap objects}
\usage{
\method{plot}{profileLikelihoodFrap}(object, ncol = NULL, nrow = NULL,
  addParam = NULL)
}
\arguments{
  \item{object}{'profileLikelihoodFrapPair' object as
  returned by \code{\link{profileLikelihoodFrapPair}}}

  \item{ncol}{number of columns for facets, when only 1D
  analysis are persent}

  \item{nrow}{number of rows for facets, when only 1D
  analysis are persent}

  \item{addParam}{additional parameter points to be
  plotted}
}
\description{
plot method for profileLikelihoodFrap objects
}
\details{
plots all 1D profile likelihood analysis and, if present,
all 2D analysis
}

