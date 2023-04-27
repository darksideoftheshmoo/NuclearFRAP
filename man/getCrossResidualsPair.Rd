\name{getCrossResidualsPair}
\alias{getCrossResidualsPair}
\title{gets the residuals of a local estimation model, on the variables were the FRAP is NOT taking place}
\usage{
getCrossResidualsPair(db, p = param.FRAP.pairs, melt = FALSE)
}
\arguments{
  \item{db}{a data.frame with the data of a single pair.
  'time', 'FRAP', 'type', nuclear and cytoplasmic
  fluorescence variables for mother and daughter required}

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
gets the residuals of a local estimation model, on the
variables were the FRAP is NOT taking place
}
\details{
gets the residuals of a \code{loess} model for the mother
nuc and cyt when the FRAP is done on the daughter and
vice-versa
}

