\name{restructurePairs}
\alias{restructurePairs}
\title{Restructures a cell.data object}
\usage{
restructurePairs(X, p = X$parameters)
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
Restructures a cell.data object
}
\details{
Restructures from a one-pos-one-cell to a one-pos-one-pair
structure. A 'pair' column is required in the description
file (p$desc).
}

