testStat <- function(sumFns,A,B,AB,wts,r,type) {
# Test statistic.
# type <--> "oneway", "addit", "interac".
#
w   <- wts
w.. <- sum(w)
if(type %in% c("oneway","addit")) {
    xxx   <- builds2Khat(sumFns,wts,splif=A,do.s2=TRUE)
    wi.   <- sapply(split(wts,f=A),sum)
    enns  <- table(A)
    vmlt  <- 1/wi. - 1/w..
    M     <- (do.call(cbind,xxx$Khatgp) - xxx$Khat)^2
    V     <- outer(xxx$s2,vmlt,"*")
    M     <- M/V
    Tstat <- apply(M,2,trapint,r=r)
    Tstat <- sum(enns*Tstat)
} else if(type=="interac") {
    xxx   <- builds2Khat(sumFns,wts,splif=AB,do.s2=TRUE)
    na    <- length(levels(A))
    nb    <- length(levels(B))
    nc    <- length(r)
    kAB   <- xxx$Khatgp
    mAB   <- do.call(cbind,kAB)
    aAB   <- aperm(array(t(mAB),dim=c(nb,na,nc)),c(2,1,3))
    kA    <- builds2Khat(sumFns,wts,splif=A,do.s2=FALSE)$Khatgp
    mAB   <- do.call(cbind,kAB)
    mA    <- do.call(cbind,kA)
    kB    <- builds2Khat(sumFns,wts,splif=B,do.s2=FALSE)$Khatgp
    mB    <- do.call(cbind,kB)
    sAB   <- lapply(1:nc,function(k,x,y){outer(x[k,],y[k,],"+")},x=mA,y=mB)
    tAB   <- array(unlist(sAB),dim=c(na,nb,nc))
    uAB   <- aperm(array(xxx$Khat,dim=c(nc,na,nb)),c(2,3,1))
    enns  <- table(A,B)
    M     <- try((aAB - tAB + uAB)^2)
    wij.  <- matrix(sapply(split(wts,f=AB),sum),
                    nrow=na,ncol=nb,byrow=TRUE)
    wi..  <- sapply(split(wts,f=A),sum)
    w.j.  <- sapply(split(wts,f=B),sum)
    oij   <- outer(wi..,w.j.,"+")
    orij  <- outer(1/wi..,1/w.j.,"+")
    vmlt  <- 1/wij. - orij + 2*wij./oij - 1/w..
    V     <- aperm(outer(xxx$s2,vmlt,"*"),c(2,3,1))
    M     <- M/V
    Tstat <- apply(M,c(1,2),trapint,r=r)
    Tstat <- sum(enns*Tstat)
} else {
    stop(paste0("Value of \"type\" ",type," not recognised.\n"))
}
Tstat
}
