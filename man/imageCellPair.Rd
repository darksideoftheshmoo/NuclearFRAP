\name{imageCellPair}
\alias{imageCellPair}
\title{makes a image montage of a cell pair with the available images}
\usage{
imageCellPair(X, pos, timeBlock = X$parameters$PRE.ACTIVATION.FRAMES,
  box.size = 128/2 - 1, ...)
}
\arguments{
  \item{X}{a cell.data object (as returned by
  \code{load.FRAP.pairs} to restructure}

  \item{pos}{number of position of interest. Numeric vector
  of length one.}

  \item{timeBlock}{number of consecutive frames to average
  in shown images}

  \item{box.size}{size in pixels of the image containing
  the cells. This specifies the 'radius', i.e. the image
  will be a square of length 2*box.size+1}

  \item{...}{further arguments for \code{\link{cimage}}}
}
\value{
The function returns an invisible 'Image' object of the
EBImage package.
}
\description{
makes a image montage of a cell pair with the available
images
}

