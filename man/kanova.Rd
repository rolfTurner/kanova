\name{kanova}
\alias{kanova}
\title{
    Analysis of variance of K-functions
}
\description{
    One-way and two-way analysis of variance for replicated point
    patterns, grouped by one or two classification factors, on the
    basis of one of the four standard summary functions, most often
    the K-function.
}
\usage{
    kanova(fmla,data,expo=2,rsteps=128,sumFnNm=c("Kest","Fest","Gest","Jest"),
           test=TRUE,permtype=c("resids","data"),nperm=99,
           brief=TRUE,verb=TRUE)
}
\arguments{
  \item{fmla}{A formula specifying the model to be fitted.  There can
  be at most two main effect predictors (possibly plus interaction
  between them).  If \code{data} (see below) is \emph{not}
  a hyperframe, the left hand side of \code{fmla} is arbitrary,
  and may be omitted.  If \code{data} \emph{is} a hyperframe, the
  left hand side of \code{fmla} must be supplied (and must be the
  name of the column which consists of a list of point patterns).
}
  \item{data}{A hyperframe (see \code{\link[spatstat.geom]{hyperframe}()})
  containing the data to be analysed.  It must have a column with
  a name which matches that of the response variable as determined
  by \code{fmla}).  If no such column exists, an error is thrown.
  The response column may be a list of point patterns, or it may
  be a list of numeric vectors (which must all have the same length).

  These numeric vectors are notionally values of the chosen
  summary function (see below), applied to point patterns, and then
  evaluated at the attribute \code{"r"} (see below) of \code{data}.
  However that the entries in the list need not actually \emph{be}
  summary functions.  They are simply numeric vectors and \emph{not}
  objects of class \code{"fv"} as returned by a summary function.
  In particular they may be scalars, which allows the \code{kanova()}
  function to be applied to one- or two-way analysis of variance
  of scalars.

  Argument \code{data} may also have one or two columns with names
  matching those of the main effect predictors as specified by
  \code{fmla}.  If such columns cannot be found, then objects
  with these names are sought in the \dQuote{parent frame}
  (see \code{\link{parent.frame}()}.)  If the appropriate data
  still cannot be found, then an error is thrown.  These columns
  (or objects) must be factors.

  If the response consists of numeric vectors, and if these vectors
  are `\emph{not} scalars, then \code{data} must have an attribute
  \code{"r"}  which specifies the \dQuote{argument} of the pseudo
  summary functions in the list.  This attribute should be an
  numeric vector and must be of the same length as the numeric
  vectors in the list.  If the response functions are scalars,
  then no \code{"r"} attribute is needed, and if present is ignored.

  Finally, \code{data} may have a columns named \code{"wts"}
  (\dQuote{weights}) whose entries are positive scalars.
  Notionally they may be thought of as the numbers of points in
  the corresponding point patterns.  If there is no column named
  \code{"wts"} then the weights are all taken to equal 1.  If the
  response column consists of a list of point patterns then any
  column name \code{"wts"} is ignored, and the weights are indeed
  taken to be the numbers of points in those patterns.
}
  \item{expo}{Numeric scalar.  Not used unless the response is a list
  of point patterns.  Statistics in \code{kanova} are calculated as
  \emph{weighted} means with the weights being the reciprocals of the
  counts of points in the patterns, raised to the power \code{expo}.
  To use unweighted means, set \code{expo=0}.  If \code{expo} is
  equal to \code{1} then the weights are simply the pattern counts.
  See the vignette \code{"testStat"} for a brief discussion of
  the impact of using these weights.
}
  \item{rsteps}{Integer scalar.  Not used unless the response is a list
  of point patterns. The number of (equal) steps between
  values of the vector \code{r} at which the summary function is
  evaluated.  The values of \code{r} are equispaced on the interval from
  0 to \code{rtop}, the latter being calculated internally.  The value
  of \code{rtop} depends on the observation windows of the patterns
  in the response and on their intensities.  It also depends on the
  summary function being used.
}
  \item{sumFnNm}{Character string naming the summary function to
  be used.
}
  \item{test}{Logical scalar.  Should a Monte Carlo test of the
  null hypothesis be carried out?
} 
  \item{permtype}{Character string specifying what sort of
  permutations should be done to produce the Monte Carlo test
  statistics.  Ignored if \code{test} is \code{FALSE}.

  If \code{permtype} is \code{"resid"} then the Monte Carlo data
  are formed by permuting the residuals from the fitted model and
  adding them to the fitted values from that model.

  If \code{permtype} is \code{"data"} then the Monte Carlo data
  are formed by permuting the orginal data sets.  (In the two-way
  setting, when the test is for the main effect A, the data are
  permuted \emph{within} the levels of B.)

  If \code{fmla} specifies an  interaction between the main effect
  predictors, then \code{permtype} \emph{must} be \code{resid},
  otherwise an error is thrown.
}
  \item{nperm}{
  The number of permutations to be used to determine the Monte Carlo
  \eqn{p}-value.  Ignored if \code{test} is \code{FALSE}
}
  \item{brief}{Logical scalar.  Should the object returned by this
  function be \dQuote{brief}?  See \bold{Value}.
}
  \item{verb}{
  Logical scalar.  Should rudimentary \dQuote{progress reports}
  be printed out (in the course of conducting the permutation
  test for \dQuote{significance} of the test statistic)?  Such
  \dQuote{reports} consists simply of indications of how many
  permutations have been effected so far.  Ignored if \code{test}
  is \code{FALSE}
}
}
\details{
  The value of the test statistic is obtained as a sum of
  numerical integrals of certain sums of squares.  The integrals
  are computed via a rough trapezoid rule and the sums of squares
  are \dQuote{studentised}.  I.e. they are  downweighted by the
  estimated variance of the quantity being squared.  The quantities
  that are involved in the sums of squares are formed from certain
  \dQuote{fitted values} which are (as indicate above) weighted
  means of the observations.  The variance referred to is formed
  as a weighted mean of squares of the residuals, which are of
  course equal to the observations minus the fitted values.

  The integration is carried out over the value of \code{r}, the
  argument of the summary functions that are being analysed.  If the
  reponse consists of numeric vectors of length 1, i.e. of scalars,
  then no integration is in fact performed, and the corresponding
  (downweighted) sum of squares is returned.  You may, if you like,
  think of this as integrating with respect to a measure which has
  a point mass of 1 at the conceptual single value of \code{r}.

  More detail about the test statistic, the fitted values, the
  residuals and the estimated variance, can be found in the vignette
  \code{"testStat"}.
}

\section{Note}{
  Only one-way and two-way (pseudo) analyses of variance are
  accommodated.  If you feel inclined to ask why there is not
  provision for higher order analysis of variance, just look at
  the code and the answer should be obvious.  It \emph{might} be
  possible to implement higher order pseudo analysis of variance of
  summary functions, but this is unlikely to have any practical use.
  Writing the code would, for me, be a nightmare!
}

\value{If \code{brief} is \code{TRUE}, the object returned is
  a list with components
 \item{stat}{Numerical scalar equal to the value of the test statistic
  calculated from the original data.}
  \item{pvalue}{The Monte Carlo \eqn{p}-value of the test calculated
  on the basis of the test statistics obtained from each of the
  \code{nperm} permutations.  This component is present only if
  \code{test} is \code{TRUE}.}
  If \code{brief} is \code{FALSE} the list returned has addtional
  components
  \item{fmla}{The \code{fmla} argument.}
  \item{data}{The \code{data} argument.  This may possibly have been
              augmented by any predictor values which were not
              found in the original \code{data} and were located
              in the parent frame.}
  \item{sumFnNm}{The \code{sumFnNm} argument.}
  \item{permtype}{The \code{permtype} argument.}
  \item{nperm}{The \code{nperm} argument.}
  \item{Tstar}{The vector of \code{nperm} values of the test statistic
               calculated from the simulated data sets, generated by
               permutation.}

}
\references{
Diggle, Peter J., Mateu, Jorge and Clough, Helen E. (2000) A
comparison between parametric and non-parametric approaches to the
analysis of replicated spatial point patterns, \emph{Advances in
Applied Probability} \bold{32}, pp. 331 -- 343.

Diggle, P. J., Lange, N. and Benes, F. (1991) Analysis of variance
for replicated spatial point patterns in clinical neuroanatomy,
\emph{Journal of the American Statistical Association}, \bold{86},
pp. 618 -- 625.

Hahn, Ute (2012) A studentized permutation test for the
comparison of spatial point patterns, \emph{Journal of the
American Statistical Association}, \bold{107}, pp. 754 -- 764,
DOI: 10.1080/01621459.2012.688463.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}
}

\seealso{
\code{\link[spatstat.explore]{studpermu.test}()}
}

\examples{
# The following is inappropriate since there is a second
# classification factor.
set.seed(104)
s1 <- kanova(patterns ~ Pos,data=stomata,permtype="d")

# Here we are testing for a Layer effect allowing for a Pos effect.
set.seed(7)
s2 <- kanova(patterns ~ Layer + Pos, data=stomata,permtype="d")

# Here we are testing for a Pos effect allowing for a Layer effect.
set.seed(78)
s3 <- kanova(patterns ~ Pos + Layer, data=stomata)
# permtype defaults to "resperm".
\dontrun{ # Takes too long.
    set.seed(24)
    s3a <- kanova(patterns ~ Pos + Layer, data=stomata,nperm=999)
}
# permtype defaults to "resperm".
\dontrun{ # Takes too long.
    # Get a p-value of 0.001
}

# Here, we are testing for a Layer by Pos interaction.
set.seed(2)
s4 <- kanova(patterns ~ Layer * Pos, data=stomata) # permtype must be "r"

# Articial data.
if(!requireNamespace("spatstat.geom"))
    stop("Required package \"spatstat.geom\" is not available.\n")
set.seed(3)
r    <- seq(0,25,length=128)
rsp  <- lapply(1:144,function(k){pi*r^2 + runif(128,-0.1,0.1)})
rsp  <- lapply(rsp,function(x){pmax(0,x)})
fac1 <- factor(rep(1:4,12,each=3))
fac2 <- factor(rep(1:3,48))
wts  <- sample(50:100,144,replace=TRUE)
X    <- spatstat.geom::hyperframe(rsp=rsp,fac1=fac1,fac2=fac2,wts=wts)
attr(X,"r") <- r
set.seed(118)
s5   <- kanova(rsp ~ fac1*fac2,data=X,brief=FALSE)

# Winbuilder says the total time for these examples is 12.14 > 10.00 seconds.
# I get only 8.927 seconds in my timing calculations. (???) Irrespective of
# that, I am cutting the following bit to reduce the time.

# Scalar data.
\dontrun{
if(require(Devore7)) {
    X   <- spatstat.geom::as.hyperframe(xmp11.10)
    s6a <- kanova(Tempr ~ Period*Strain,data=X,nperm=999)
    s6b <- kanova(Tempr ~ Period+Strain,data=X,nperm=999)
    s6c <- kanova(Tempr ~ Strain+Period,data=X,nperm=999)
    chk <- lm(Tempr ~ Period*Strain,data=X)
# anova(chk) reveals p-values that are at least roughly
# similar to those in s6a, s6b, and s6c.
} else {
    stop("Required package \"Devore7\" is not available.\n")
}
}
}

\keyword{ htest }
