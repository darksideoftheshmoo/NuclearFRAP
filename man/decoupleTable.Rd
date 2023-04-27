\name{decoupleTable}
\alias{decoupleTable}
\title{Decouple a data.frame}
\usage{
decoupleTable(db, mother.posfix = "_mother", daughter.posfix = "_daughter")
}
\arguments{
  \item{db}{a data.frame to be restructured. It should have
  a pos varibale and '*_mother' and '*_daughter' coupled
  vars.}

  \item{mother.posfix}{posfix identifying mother variables}

  \item{daughter.posfix}{posfix identifying daughter
  variables}
}
\value{
a restructured data.frame
}
\description{
Decouple a data.frame
}
\details{
Restructures a data.frame, transforming '*_mother' and
'*_daughter' variables to the appropiate position to be
compatible with \code{decouplePairs}.
}

