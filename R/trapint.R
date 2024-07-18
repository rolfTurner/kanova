trapint <- function(y, r) {
# Integrate function "y" w.r.t. the variable "r"
# by the trapezoid rule.
    if(length(y)==1) return(y) # Point mass measure.
    nonan <- is.finite(y)
    nn    <- sum(nonan)
    if(nn < 2L) return(0)
    y <- y[nonan]
    r <- r[nonan]
    0.5 * sum( (y[-1] + y[-nn]) * diff(r))
}
