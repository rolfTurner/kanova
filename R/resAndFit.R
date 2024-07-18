resAndFit <- function(iDat,type) {
# Note:  These are the residuals and fitted values under the
# appropriate null hypothesis.

N <- length(iDat$sumFns)
switch(EXPR=type,
    oneway = { # Fitted values from one factor model; A only.  Null
               # hypothesis is that there is no A effect, so K_{ij}-hat = Khat.
        Khat  <- with(iDat,builds2Khat(sumFns,wts,splif=A,do.s2=FALSE)$Khat)
        fitz  <- lapply(1:N,function(k,c){c},c=Khat)
    },
    addit = { # Fitted values from additive A + B model.  Null hypothesis is
              # that there is no A effect, so K_{ijk}-hat = K_{.j.}-bar
        ufitzB <- with(iDat,builds2Khat(sumFns,wts,splif=B)$Khatgp)
        fitzB  <- ufitzB[match(iDat[["B"]],names(ufitzB))]
        fitz   <- lapply(1:N,function(i,b){b[[i]]},b=fitzB)
    },
    interac = { # Fitted values from model with interaction, A * B.
                # Null hypothesisis that there is no interaction,
                # i.e. that the model is additive so
                # K_{ijk}-hat = K_{i..}-bar + K_{.j.}-bar - Khat
        ufitzA <- with(iDat,builds2Khat(sumFns,wts,splif=A)$Khatgp)
        ufitzB <- with(iDat,builds2Khat(sumFns,wts,splif=B)$Khatgp)
        Khat   <- with(iDat,builds2Khat(sumFns,wts,splif=AB)$Khat)
        fitzA  <- ufitzA[match(iDat[["A"]],names(ufitzA))]
        fitzB  <- ufitzB[match(iDat[["B"]],names(ufitzB))]
        fitz   <- lapply(1:N,function(i,a,b,c){a[[i]] + b[[i]] - c},
                          a=fitzA,b=fitzB,c=Khat)
    }
)
rez  <- lapply(1:length(fitz),function(i,obs,fit){obs[[i]] - fit[[i]]},
               obs=iDat$sumFns,fit=fitz)
names(rez) <- names(fitz)
list(resids=rez,fitVals=fitz)
}
