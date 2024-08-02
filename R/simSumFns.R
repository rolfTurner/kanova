simSumFns <- function(sumFns,B,rAndF,splif,s2,Khat,infertype) {
    N <- length(sumFns)
    if(infertype=="datperm") {
        if(is.null(B)) {
            ip <- sample(1:N,N)
            sSF <- sumFns[ip]
        } else {
            ip     <- permWithin(B)
            sSF <- sumFns[ip]
        }
    } else if(infertype == "resperm") {
        fV     <- do.call(cbind,rAndF$fitVals)
        ip     <- sample(1:N,N)
        pRes   <- do.call(cbind,rAndF$resids[ip])
        sSF <- as.list(as.data.frame(fV+pRes))
    } else if(infertype == "gaussSample") {
        s1  <- sqrt(s2)
        sSF <- lapply(1:N,function(k,s1,Khat){rnorm(n=length(s1),
                                              mean=Khat,sd=s1)},
                      s1=s1,Khat=Khat)
# Got here.
    }
browser()
    sSF
}
