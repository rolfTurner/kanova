\encoding{UTF-8}
\name{stomata}
\alias{stomata}
\docType{data}
\title{
   Stomata patterns
}
\description{
   Point patterns of stomata at 18 locations in each of 12 leaves of
   \emph{Michelia cavaleriei var. platypetala.}
}
\usage{stomata}
\format{
  The object \code{stomata} is a hyperframe with 216 rows.  It has
  a column named \code{patterns} containing the point patterns of
  stomata locations, a column named \code{Leaf} which is a factor
  with levels 1 to 12 identifying the leaf from which the pattern
  was obtained, a column named \code{Layer} which is a factor with
  levels 1 to 6 identifying a location within each \dQuote{position}
  (see \dQuote{\code{Pos}}), and a column named \code{Pos} which
  is a factor with levels 1 to 3, position 1 being closest to the
  central stem of the leaf and position 3 being closest to the
  outer edge of the leaf (and farthest from the central stem).
}
\details{
  Each pattern was observed in a rectangular window of dimension
  1200 \eqn{\times}{x} 900 microns.
}
\source{
  The data were kindly supplied by Prof. Peijian Shi of the College
  of Biology and the Environment, Nanjing Forestry University, Nanjing,
  P.R. China.
}
\references{
Peijian Shi, Yabing Jiao, Peter J. Diggle, Rolf Turner,
Rong Wang and \enc{Ülo}{Ulo} Niinemets 2021.  Spatial
relationship between stomata of a Magnoliaceae species at the
areole level.  \emph{Annals of Botany} \bold{128}, pp. 875--885.
DOI https://doi.org/10.1093/aob/mcab106
}

\examples{
plot(stomata[1,1,drop=TRUE])
}
\keyword{datasets}
