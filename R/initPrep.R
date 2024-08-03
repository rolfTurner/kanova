initPrep <- function(data,rspNm=NULL,Anm,Bnm=NULL,sumFnNm,type) {
#
# Function to prepare a list whose entries correspond to the
# model cells and are (or are interpretable as) summary functions
# of replicated point patterns.
#
# The object "data" is hyperframe.  The argument "rspNm" ("response
# name") specifies the column of the hyperframe constituting
# the response.  It defaults to the name or the *first* column.
# The column specified by "rspNm" must be either a list of point
# patterns or a list whose entries are numeric vectors.
#
# In the latter case the numeric vectors will be interpreted
# as  summary functions.  They must all be of the same length and
# all of their entries must be non-negative.
#
#
# Check that "data" is a hyperframe.
if(!is.hyperframe(data)) {
    stop("Argument \"data\" must be a hyperframe.\n")
}

# Check on the response name.
if(is.null(rspNm)) {
    rspNm <- names(data)[1]
} else {
    if(!rspNm %in% names(data)) {
        whinge <- paste0("The response name is ",rspNm,", which is not the name\n",
                         "  of any of the columns of the hyperframe \"data\".\n")
        stop(whinge)
    }
}

# Check on the predictor names and dig out the predictors.
if(!(Anm %in% names(data))) {
    whinge <- paste0("The name of the first predictor is \"",Anm,"\", which",
                     " is not the\n  name of any of the columns of",
                     " the hyperframe \"data\".\n")
    stop(whinge)
}
A <- factor(data[,Anm,drop=TRUE])
if(!(is.null(Bnm) || Bnm %in% names(data))) {
    whinge <- paste0("The name of the second predictor is \"",Bnm,"\", which",
                     " is not the\n  name of any of the columns of",
                     " the hyperframe \"data\".\n")
    stop(whinge)
}
if(is.null(Bnm)) B <- NULL else B <- factor(data[,Bnm,drop=TRUE])

# Check that the response has the appropriate structure.
if(inherits(data[,rspNm,drop=TRUE],"ppplist")) {
    bldSumFns <- TRUE
} else {
    chkClass <- unique(sapply(data[,rspNm,drop=TRUE],class))
    if(length(chkClass) != 1)
        stop(paste0("All entries of ",rspNm," must be of the same class.\n"))
    if(!(chkClass %in% c("ppp","numeric"))) {
        whinge <- paste0("The response must consist either of point",
                         " patterns or of numeric vectors.\n")
        stop(whinge)
    }
    if(chkClass=="numeric") {
        bldSumFns <- FALSE
        lnth <- unique(sapply(data[,rspNm,drop=TRUE],length))
        if(length(lnth) != 1) {
            whinge <- paste0("All of the numeric vectors in \"rspNm\"",
                             " must have the same length.\n")
            stop(whinge)
        }
        if(lnth > 1) {
            r <- attr(data,"r")
            if(is.null(r)) {
                whinge <- paste0("When the reponse consists of non-scalar",
                                 " summary functions\n \"data\" must have an",
                                 " attribute named \"r\".\n")
                stop(whinge)
            }
        } else {
            r <- 1
        }
        if("wts" %in% names(data)) {
            wts <- data[,"wts",drop=TRUE]
        } else {
            wts <- rep(1,nrow(data))
        }
    }
}

if(bldSumFns) {
    mikes <- data[,rspNm,drop=TRUE]
    if(requireNamespace("spatstat.geom")) {
        wts    <- sapply(mikes,spatstat.geom::npoints)
        if(any(wts==0)) stop("Some point patterns in \"data\" are empty.\n")
        sumFn  <- get(sumFnNm)
        sFraw  <- lapply(mikes,sumFn)
        r      <- sFraw[[1]]$r
        sumFns <- lapply(sFraw,function(x){x[[attr(x,"valu")]]})
        sumFns <- lapply(1:length(sumFns),function(k,x,w){
                                          attr(x[[k]],"weight") <- w[k]
                                          x[[k]]},x=sumFns,w=wts)
    } else {
        stop("Required package \"spatstat.geom\" is not available.\n")
    }
} else {
    sumFns <- data[,rspNm,drop=TRUE]
}

# Re-normalise the weights to have maximum 1.  (Not sure
# about this idea!!! Perhaps mean 1???  Or to sum to 1???
# Or just don't fuck around with them?)  The intent is to
# avoid overflow which might result from large weights.
# Note that re-normalising the weights does *not* (I'm pretty
# sure) screw up the assertion that s2 (see builds2Khat.R) is
# an unbiased estimate of sigma^2(r).  The value of sigma^2(r)
# being estimated changes to sigma_M^2(r) = sigma^2(r)/M, where
# M is the normalising constant by which "wts" is being divided.
# But this (I'm pretty sure) does not matter.
wts <- wts/max(wts)

# Using interaction(B,A) r.t. interaction(A,B) seems counterintuitive,
# but is necessary for making "permute within" work; at least the way
# I currently have things structured.
if(is.null(B)) {
    AB <- NULL
} else {
    AB <- interaction(B,A)
}

# Build Khat (the overall estimate of the unique K function, common
# to all groups under the null hypothesis of no group effects),
# and s2, the overall sample variance.
if(type %in% c("oneway","addit")) {
    splif <- A
} else if(type == "interac") {
    splif <- AB
} else {
    stop(paste0("Unrecognised type ",type,".\n"))
}
enns <-table(splif)
if(any(enns < 2))
    stop("All cell counts must be at least 2.\n")
list(sumFns=sumFns,A=A,B=B,AB=AB,wts=wts,r=r)
}
