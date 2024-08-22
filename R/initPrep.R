initPrep <- function(data,rspNm=NULL,Anm,Bnm=NULL,sumFnNm,type,expo,rsteps,r) {
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
    if(requireNamespace("spatstat.geom")) {
        mikes <- data[,rspNm,drop=TRUE]
        wts   <- sapply(mikes,spatstat.geom::npoints)
        if(any(wts==0)) stop("Some point patterns in \"data\" are empty.\n")
        wts   <- wts^expo
        sumFn <- get(sumFnNm)
        if(is.null(r)) {
        Let   <- if(sumFnNm=="Kest") "K" else "F"
            if(requireNamespace("spatstat.explore")) {
                rtop  <- min(sapply(mikes,function(x){
                             spatstat.explore::rmax.rule(Let,
                             spatstat.geom::Window(x),spatstat.geom::intensity(x))})) 
            } else {
                stop("Required package \"spatstat.explore\" is not available.\n")
            }
            r      <- seq(0,rtop,length=rsteps+1)
        }
        sFraw  <- lapply(mikes,sumFn,r=r)
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

# Using interaction(B,A) r.t. interaction(A,B) seems counterintuitive,
# but is necessary for making "permute within" work; at least the way
# I currently have things structured.
if(is.null(B)) {
    AB <- NULL
} else {
    AB <- interaction(B,A)
}

# Check that the cell counts are adequate.
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
