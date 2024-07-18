permWithin <- function(G) {
newind <- integer(length(G))
for(x in levels(G)) {
   ok <- G==x
   i  <- which(ok)
   ip <- sample(i,length(i))
   newind[ok] <- ip
}
newind
}
