\name{getResidualsPair}
\alias{getResidualsPair}
\title{gets the residuals of the data to the model, for a pair of cells}
\usage{
getResidualsPair(...)
}
\arguments{
  \item{db}{a data.frame with the data of a pair cell
  (required if expData=NULL). 'time', 'FRAP', 'type',
  nuclear and cytoplasmic fluorescence variables for the
  mother and daughter required.}

  \item{param}{a data.frame with the parameters of the fit
  (required). Must have 'type' vatiable and a variable for
  each free parameter of the model.}

  \item{a}{list of EXPDATA as returned by
  \code{makePairEXPDATA}, which contains the data for a
  pair of cells (required if db=NULL).}

  \item{p}{list of options for the FRAP experiment}

  \item{melt}{boolean specifying if the output should be
  'melted' by \code{melt}}
}
\value{
a data.frame. If melt=FASE => variables
"FRAP","time","res.nuc","res.cyt","type". If melt=TRUE =>
variables "FRAP","time","type","variable","value"
}
\description{
gets the residuals of the data to the model, for a pair of
cells
}

