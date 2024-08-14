\name{kanova-internal}
\alias{builds2Khat}
\alias{initPrep}
\alias{ordinal}
\alias{ordinalsuffix}
\alias{permSumFns}
\alias{permWithin}
\alias{resAndFit}
\alias{testStat}
\alias{trapint}
\title{Internal kanova functions.}
\description{
  Internal kanova functions.
}
\usage{
builds2Khat(sumFns,wts,splif,do.s2=TRUE)
initPrep(data,rspNm=NULL,Anm,Bnm=NULL,sumFnNm,type,expo,rsteps)
ordinal(k)
ordinalsuffix(k)
permSumFns(sumFns,B,rAndF,permtype)
permWithin(G)
resAndFit(iDat,type)
testStat(sumFns,A,B,AB,wts,r,type)
trapint(y,r)
}
\details{
  These functions are auxiliary and are not intended to be called
  by the user.
}
\keyword{internal}
