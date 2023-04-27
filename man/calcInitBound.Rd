\name{calcInitBound}
\alias{calcInitBound}
\title{calculate initial values of parameters and upper/lower bounds}
\usage{
calcInitBound(EXPDATA, p = param.FRAP.pairs)
}
\arguments{
  \item{EXPDATA}{list as returned by \code{makeEXPDATA},
  which contains the data for a cell (required if
  db=NULL).}

  \item{p}{list of options for the FRAP experiment}
}
\value{
a list with elements x0, LB, UB and names
}
\description{
This funcion has build-in estimations that can be
overwriten either by p$VALUE.PARAMETERS or by
EXPDATA$expParam. if 'tot0_mean' and 'tot0_sd' are in
p$MODEL.PARAMETERS, and have NA values in
p$VALUE.PARAMETERS, they are calculated from the
fluorescence variables of the cell and geomVfrac
}

