ripVar <- function(X,r) {
#
# Function to compute the variance of the (estimated) K-function
# corresponding to the point pattern X (which is an object of class "ppp").
# This variance is a function of the argument "r" of the K-function.
#
# Makes use of equation (3) from page 755 of the paper
# "A Studentized Permutation Test for the Comparison of Spatial
#  Point Patterns" by Ute Hahn, JASA volume 107, pp. 754 -- 764.
#
if(requireNamespace("spatstat.geom")) {
    c1 <- 0.305
    c2 <- 0.0415
    W  <- spatstat.geom::Window(X)
    n  <- spatstat.geom::npoints(X)
    A  <- spatstat.geom::area.owin(W)
    U  <- spatstat.geom::perimeter(W)
    v  <- (2*pi*r^2*A/n^2)*(1 + c1*U*r/A + c2*n*U*r^3/A^2)
} else {
    stop("Required package spatstat.geom is not available.\n")
}
v
}
