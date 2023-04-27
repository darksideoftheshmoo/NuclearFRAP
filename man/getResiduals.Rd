\name{getResiduals}
\alias{getResiduals}
\title{gets the residuals of the data to the model}
\usage{
getResiduals(db = NULL, param = NULL, expData = NULL,
  p = param.FRAP.pairs, melt = FALSE)
}
\arguments{
  \item{db}{a data.frame with the data of a single cell
  (required if expData=NULL). 'time', 'FRAP', nuclear and
  cytoplasmic fluorescence variables required.}

  \item{param}{a data.frame with the parameters of the fit
  (required). Must have a variable for each free parameter
  of the model, and a sigle row.}

  \item{expData}{a EXPDATA list as returned by
  \code{makeEXPDATA}, which contains the data for a cell
  (required if db=NULL).}

  \item{p}{list of options for the FRAP experiment}

  \item{melt}{boolean specifying if the output should be
  'melted' by \code{melt}}
}
\value{
a data.frame. If melt=FASE => variables
"FRAP","time","res.nuc","res.cyt". If melt=TRUE =>
variables "FRAP","time","variable","value"
}
\description{
gets the residuals of the data to the model
}

