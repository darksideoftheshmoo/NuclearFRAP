\name{makePairEXPDATA}
\alias{makePairEXPDATA}
\title{Makes a list with a pair of EXPDATA}
\usage{
makePairEXPDATA(db, p = param.FRAP.pairs)
}
\arguments{
  \item{db}{a data.frame with the data of a single pair.
  'time', 'FRAP', 'type', nuclear and cytoplasmic
  fluorescence variables for mother and daughter required}

  \item{p}{list of options for the FRAP experiment}
}
\value{
a list with a 'mother' and a 'daughter' element, each of
which contain the correspondent EXPDATA
}
\description{
Makes a list with a pair of EXPDATA
}

