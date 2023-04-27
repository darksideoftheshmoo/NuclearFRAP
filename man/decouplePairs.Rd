\name{decouplePairs}
\alias{decouplePairs}
\title{Decouple a cell.data pair object}
\usage{
decouplePairs(X, p = X$parameters)
}
\arguments{
  \item{X}{a cell.data object (as returned by
  \code{load.FRAP.pairs} to restructure}

  \item{p}{list of options for the FRAP experiment. Usually
  the default \code{X$parameters} will do.}
}
\value{
a restructured cell.data object
}
\description{
Decouple a cell.data pair object
}
\details{
Restructures from a one-pos-one-pair to a one-pos-one-cell
structure.
}

