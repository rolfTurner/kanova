kanova <- function(fmla,data,expo=2,rsteps=128,sumFnNm=NULL,warnSFN=TRUE,
                   test=TRUE,permtype=c("resids","data"),nperm=99,
                   brief=TRUE,verb=TRUE) {
#
# Function to conduct one or two-way analysis of variance of
# summary functions (Kest, Fest, Gest, or Jest)  of replicated point
# patterns classified by a grouping factor A or two grouping factors
# A and B.
#

if(is.null(sumFnNm)) sumFnNm <- "Kest"
if(!(sumFnNm %in% c("Kest","Fest","Gest","Jest")) & warnSFN) {
    whinge <- paste0("Argument \"sumFnNm\" is \"",sumFnNm,"\", which is not\n",
                     "  one of the standard four.  The results may be fragile.\n")
    warning(whinge)
}
permtype <- match.arg(permtype)

if(length(fmla)==2) {
    if(inherits(data,"hyperframe")) {
        whinge <- paste0("When \"data\" is a hyperframe, the left hand side",
                         " of \"fmla\" must be specified.\n")
        stop(whinge)
    } else {
        cform <- c("y",as.character(fmla))
        fmla  <- formula(paste(cform,collapse=" "))
    }
}
rspNm  <- as.character(fmla[2])
preds  <- attr(terms(fmla),"term.labels")
npreds <- length(preds)
if(npreds > 3) {
    whinge <- paste0("The length of the vector of predictor names, (",
                     paste(preds,collapse=", "), "), is ",npreds,".\n",
                     "  It must be at most 3.\n")
     stop(whinge)
}
for(i in 1:npreds) {
    if(grepl(":",preds[i])) next
    if(!preds[i] %in% names(data)) {
       xp <- try(get(preds[i],envir=parent.frame()),silent=TRUE)
       if(inherits(xp,"try-error"))
           stop(paste0("Cannot find predictor \"",preds[i],"\".\n"))
       newcol <- data.frame(xp)
       names(newcol) <- preds[i]
       data <- cbind(data,newcol)
    }
}
switch(EXPR=npreds,
    {Anm  <- preds[1]
     Bnm  <- NULL
     type <- "oneway"
     effNm <- Anm
    },
    {Anm <- preds[1]
     Bnm <- preds[2]
     type <- "addit"
     effNm <- paste0(Anm," allowing for ",Bnm)
    },
    {Anm <- preds[1]
     Bnm <- preds[2]
     if(preds[3] != paste0(Anm,":",Bnm)) {
         stop("Argument \"fmla\" is of an incorrect form.\n")
     }
     if(permtype=="data")
        stop("Cannot use permtype=\"data\" when there is interaction in the model.\n")
     type <- "interac"
     effNm <- paste0("interaction of ",Anm," with ",Bnm)
    }
)
# Initial (real) data:
iDat <- initPrep(data,rspNm=rspNm,Anm=Anm,Bnm=Bnm,sumFnNm=sumFnNm,
            type=type,expo=expo,rsteps=rsteps)

# Calculate the statistic.
Tobs <- with(iDat,testStat(sumFns,A,B,AB,wts,r,type=type))
if(!test) {
   if(brief) {
       rslt <- list(stat=Tobs)
   } else {
       rslt <- list(fmla=fmla,data=data,sumFnNm=sumFnNm,
                    permtype=permtype,stat=Tobs)
   }
   class(rslt) <- "kanova"
   return(rslt)
}

# Testing;  carry out the Monte Carlo test.
# If permtype is "resids", create the fitted values and residuals.
if(permtype=="resids") {
   rAndF <- resAndFit(iDat,type) # List with components "resids" and "fitVals".
} else {
   rAndF <- NULL
}

Tstar <- numeric(nperm)
for(i in 1:nperm) {
    pSumFns  <- permSumFns(iDat[["sumFns"]],iDat[["B"]],rAndF,permtype)
    Tstar[i] <- with(iDat,testStat(pSumFns,A,B,AB,wts,r,type=type))
    if(verb) cat(i,"")
    if(verb & i%%10 == 0) cat("\n")
}
if(verb & i%%10 != 0) cat("\n")
m    <- sum(Tstar >= Tobs)
pv   <- (m+1)/(nperm+1)
bres <- list(EffectName=effNm,stat=Tobs,pvalue=pv)
if(brief) {
    rslt <- bres
} else {
    rslt <- c(bres,list(fmla=fmla,data=data,sumFnNm=sumFnNm,
                        permtype=permtype,Tstar=Tstar))
}
class(rslt) <- "kanova"
rslt
}
