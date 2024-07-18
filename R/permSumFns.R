permSumFns <- function(sumFns,B,rAndF,permtype) {
    N <- length(sumFns)
    if(permtype=="data") {
        if(is.null(B)) {
            ip <- sample(1:N,N)
            sumFns <- sumFns[ip]
        } else {
            ip     <- permWithin(B)
            sumFns <- sumFns[ip]
        }
    } else { # permtype == "resids"
        fV     <- do.call(cbind,rAndF$fitVals)
        ip     <- sample(1:N,N)
        pRes   <- do.call(cbind,rAndF$resids[ip])
        sumFns <- as.list(as.data.frame(fV+pRes))
    }
    sumFns
}
