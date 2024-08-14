\name{ripVar}
\alias{ripVar}
\title{
The Ripley Variance of K-functions.
}
\description{
Calculates the variance of the K-function of a Poisson point pattern,
according to Ripley's formula (as taken from equation (3) in
Hahn 2012).
}
\usage{
ripVar(X, r)
}
\arguments{
  \item{X}{
  A point pattern (object of class \code{"ppp"}).  The variance
  formula is valid only if \code{X} arises from a Poisson process.
}
  \item{r}{
  A numeric vector of non-negative values at which the K-function
  for \code{X} is to be evaluated.
}
}
\details{
  The vector \code{r} would normally have entries in increasing
  order and would have a first entry equal to 0.  It may be wise
  to construct \code{r} as \code{Kest(X)$r}, but this is not
  required.
}
\value{
  A number vector, of length equal to \code{length(r)} whose
  entries are the variances of \code{K(r)} where $K(r)$ is equal to
  \code{as.function(Kest(X))}.
}
\seealso{\code{\link[spatstat.explore]{Kest}()}}

\references{
Hahn, Ute (2012) A studentized permutation test for the
comparison of spatial point patterns, \emph{Journal of the
American Statistical Association}, \bold{107}, pp. 754 -- 764,
DOI: 10.1080/01621459.2012.688463.
}
\author{Rolf Turner
  \email{rolfturner@posteo.net}
}
\examples{
if(require(spatstat.random)) {
    X    <- rpoispp(100)
    vKX1 <- ripVar(X,r=0.05*(1:5))
    if(require(spatstat.explore)) {
        r    <- spatstat.explore::Kest(X)$r
        vKX2 <- ripVar(X,r=r)
        plot(r,vKX2,type="l")
        points(0.05*(1:5),vKX1)
    } else {
        stop("Required package \"spatstat.explore\" is not available.\n")
    }
} else {
    stop("Required package \"spatstat.random\" is not available.\n")
}
}
\keyword{ utility }
