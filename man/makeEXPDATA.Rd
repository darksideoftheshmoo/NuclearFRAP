\name{makeEXPDATA}
\alias{makeEXPDATA}
\title{Makes a EXPDATA list, as required by matlab's functions}
\usage{
makeEXPDATA(db, nuc = "f.nuc*", cyt = "f.cyt*", sigma.prefix = "sd.",
  p = param.FRAP.pairs)
}
\arguments{
  \item{db}{a data.frame with the data of a single cell.
  'time', 'FRAP', nuclear and cytoplasmic fluorescence
  variables required}

  \item{nuc}{string that defines the pattern that
  identifies the nuclear fluorescence variable}

  \item{cyt}{string that defines the pattern that
  identifies the cytoplasmic fluorescence variable}

  \item{p}{list of options for the FRAP experiment}
}
\value{
a EXPDATA list with elements time, y, std, FRAP, PBindex,
timeFocusDrift and sF
}
\description{
Makes a EXPDATA list, as required by matlab's functions
}

