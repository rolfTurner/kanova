builds2Khat <- function(sumFns,wts,splif,do.s2=TRUE) {
#
# sumFns is a list of summary functions; one for each replicate in
# each group.   The groups may be singly indexed (for one factor
# anova) or doubly indexed (for two factor anova).  However we
# can, w.l.o.g., treat them as if they were singly indexed.

splitS <- split(sumFns,f=splif)
G      <- length(splitS)
# splitS is a list of length "G", G = number of groups
# Each entry of splitS is the list of summary functions from the
# corresponding group. 

splitS <- lapply(splitS,function(x){
                            do.call(cbind,as.list(x))})
# Now each entry of splitS is an (s x n) matrix, where s is the
# length of the summary function argument "r", and n is the number
# of replicates in the corresponding group.  I.e. the j-th column
# of that matrix is the j-th K-function in that group.

splitw <- split(wts,f=splif)
# Each entry of splitw is a a list whose entries are weight vectors.
# The lengths of these vectors are either all equal to the length
# of the summary function argument "r", or all equal to 1.  The
# j-entry of that list is the weight vector for the j-th member of
# the corresponding group.

# Got here.
Khatgp <- lapply(1:G,function(g,S,w){
                     S[[g]]%*%(w[[g]]/sum(w[[g]]))
                     },S=splitS,w=splitw)
Khatgp <- lapply(Khatgp,as.vector)
names(Khatgp) <- levels(splif)
# Khatgp is a list of vectors, each of length r, the j-th entry
# being the K-function estimate corresponding to the j-th group.

wsum <- sapply(splitw,sum) # vector of length G; the i-th entry is
                           # the sum of the entries of the weight
                           # vector for the i-th groupp
Khat <- as.vector(do.call(cbind,Khatgp)%*%(wsum/sum(wsum)))
# Khat is the overall estimate of K = the underlying K-function
# (common to all groups under H_0); it is a weighted mean of
# Khatgp[[1]], ..., Khatgp[[G]].
if(!do.s2) return(list(Khatgp=Khatgp,Khat=Khat))

D2   <- lapply(1:G,function(g,S,Khatgp){
                                (S[[g]] - Khatgp[[g]])^2
                           },S=splitS,Khatgp=Khatgp)
S2   <- lapply(1:G,function(g,D2,w){D2[[g]] %*% w[[g]]},D2=D2,w=splitw)
SS   <- Reduce("+",S2)
ndot <- length(sumFns)
s2   <- as.vector(SS/(ndot-G))
# s2 is (under H_0) an unbiased estimate of sigma^2.

list(Khatgp=Khatgp,Khat=Khat,s2=s2)
}
