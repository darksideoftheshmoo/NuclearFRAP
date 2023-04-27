\name{pairFun}
\alias{pairFun}
\title{creates a 'Pair' function from a 'cell' function}
\usage{
pairFun(FUN, db, param, p = param.FRAP.pairs, ...)
}
\arguments{
  \item{FUN}{the single cell function to be transformed to
  a pair function}

  \item{db}{either a data.frame with a 'type' column, or a
  list with 'mother' and 'daughter' elements containing the
  EXPDATA objects.}

  \item{param}{a data.frame with the parameters of the fit
  (required). Must have 'type' vatiable and a variable for
  each free parameter of the model.}

  \item{p}{list of options for the FRAP experiment}

  \item{...}{further arguments for FUN}
}
\value{
a data.frame, resulting from combining both calls to FUN,
with a 'type' column
}
\description{
creates a 'Pair' function from a 'cell' function
}
\details{
transforms a function that operates over a single cell, to
a function that operates over a mother-daughter pair
}

