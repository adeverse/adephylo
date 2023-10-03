#' Computes Moran's index for a variable
#' 
#' This simple function computes Moran's index of autocorrelation given a
#' variable and a matrix of proximities among observations.
#' 
#' 
#' @aliases moran.idx
#' @param x a numeric vector whose autocorrelation is computed.
#' @param prox a matrix of proximities between observations, as computed by the
#' \code{\link{proxTips}}. Off-diagonal terms must be positive or null, while
#' diagonal terms must all equal zero.
#' @param addInfo a logical indicating whether supplementary info (null value,
#' minimum and maximum values) should be returned (TRUE) or not (FALSE,
#' default); if computed, these quantities are returned as attributes.
#' @return The numeric value of Moran's index.
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso \code{\link{proxTips}} which computes phylogenetic proximities
#' between tips of a phylogeny.
#' @references Moran, P.A.P. (1948) The interpretation of statistical maps.
#' \emph{Journal of the Royal Statistical Society, B} \bold{10}, 243--251.
#' 
#' Moran, P.A.P. (1950) Notes on continuous stochastic phenomena.
#' \emph{Biometrika}, \bold{37}, 17--23.
#' 
#' de Jong, P. and Sprenger, C. and van Veen, F. (1984) On extreme values of
#' Moran's I and Geary's c. \emph{Geographical Analysis}, \bold{16}, 17--24.
#' @keywords manip
#' @examples
#' 
#' 
#' ## use maples dataset
#' if(require(ape) && require(phylobase)){
#' data(maples)
#' tre <- ape::read.tree(text=maples$tre)
#' dom <- maples$tab$Dom
#' bif <- maples$tab$Bif
#' 
#' 
#' ## get a proximity matrix between tips 
#' W <- proxTips(tre, met="Abouheif")
#' 
#' ## compute Moran's I for two traits (dom and bif)
#' moran.idx(dom, W)
#' moran.idx(bif, W)
#' moran.idx(rnorm(nTips(tre)), W)
#' 
#' ## build a simple permutation test for 'bif'
#' sim <- replicate(499, moran.idx(sample(bif), W)) # permutations
#' sim <- c(moran.idx(bif, W), sim)
#' 
#' pval <- mean(sim>=sim[1]) # right-tail p-value
#' pval
#' 
#' hist(sim, col="grey", main="Moran's I Monte Carlo test for 'bif'") # plot
#' mtext("Histogram of permutations and observation (in red)")
#' abline(v=sim[1], col="red", lwd=3)
#' }
#' 
#' @rdname moranIdx
#' @export 
moran.idx <- function(x, prox, addInfo=FALSE){

    ## handle arguments
    if(any(is.na(x))) stop("NA entries in x")
    if(!is.numeric(x)) stop("x is not numeric")

    W <- as.matrix(prox)
    if(!is.matrix(W)) stop("prox is not a matrix")
    if(ncol(W) != nrow(W)) stop("prox is not a square matrix")
    if(any(is.na(W))) stop("NA entries in prox")
    diag(W) <- 0

    n <- nrow(W)


    ## main computations
    x <- x - mean(x)
    sumW <- sum(W)
    num <- n * sum(x * (W %*% x) )
    denom <- sumW * sum(x*x)

    if(denom < 1e-14) stop("denominator equals zero")

    res <- num/denom

    if(addInfo){
        I0 <- -1/(n-1)
        matToDiag <- .5 * (t(W) + W)
        rangeI <- range(eigen(matToDiag)$values)
        attr(res, "I0") <- I0
        attr(res, "Imin") <- rangeI[1]
        attr(res, "Imax") <- rangeI[2]
    }

    return(res)

} # end moran.idx
