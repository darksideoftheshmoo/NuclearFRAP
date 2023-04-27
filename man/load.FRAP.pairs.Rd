\name{load.FRAP.pairs}
\alias{load.FRAP.pairs}
\title{Creates cell.data object from FRAP experiment. Each pos is a cell or cell pair.}
\usage{
\method{load}{FRAP.pairs}(path = getwd(), img.path = path,
  exp.name = basename(path), p = param.FRAP.pairs, ...)
}
\arguments{
  \item{path}{path to the directory where the FRAP
  processed data is located}

  \item{img.path}{alternative path where the original
  images might be}

  \item{exp.name}{name of the experiment. By dafult the
  path basename}

  \item{p}{list of options for the FRAP experiment. See
  \code{param.FRAP.pairs}}

  \item{...}{elements of p to be replaced by updated value}
}
\value{
a \code{cell.data} object, as returned by the load
functiond of \code{rcell2} package
}
\description{
Creates cell.data object from FRAP experiment. Each pos is
a cell or cell pair.
}

