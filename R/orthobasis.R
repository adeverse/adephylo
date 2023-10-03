#' Computes Moran's eigenvectors from a tree or a phylogenetic proximity matrix
#' 
#' The function \code{orthobasis.phylo} (also nicknamed \code{me.phylo})
#' computes Moran's eigenvectors (ME) associated to a tree. If the tree has 'n'
#' tips, (n-1) vectors will be produced. These vectors form an orthonormal
#' basis: they are centred to mean zero, have unit variance, and are
#' uncorrelated. Each vector models a different pattern of phylogenetic
#' autocorrelation. The first vectors are those with maximum positive
#' autocorrelation, while the last vectors are those with maximum negative
#' autocorrelation. ME can be used, for instance, as regressors to remove
#' phylogenetic autocorrelation from data (see references).\cr
#' 
#' ME can be obtained from a tree, specifying the phylogenetic proximity to be
#' used. Alternatively, they can be obtained directly from a matrix of
#' phylogenetic proximities as constructed by \code{\link{proxTips}}.
#' 
#' 
#' @aliases orthobasis.phylo me.phylo
#' @param x A tree of class \code{\link[ape:read.tree]{phylo}},
#' \linkS4class{phylo4} or \linkS4class{phylo4d}.
#' @param prox a matrix of phylogenetic proximities as returned by
#' \code{\link{proxTips}}.
#' @param method a character string (full or abbreviated without ambiguity)
#' specifying the method used to compute proximities; possible values are:\cr -
#' \code{patristic}: (inversed sum of) branch lengths \cr - \code{nNodes}:
#' (inversed) number of nodes on the path between the nodes \cr -
#' \code{oriAbouheif}: original Abouheif's proximity, with diagonal (see
#' details in \code{\link{proxTips}}) \cr - \code{Abouheif}: Abouheif's
#' proximity (see details in \code{\link{proxTips}}) \cr - \code{sumDD}:
#' (inversed) sum of direct descendants of all nodes on the path (see details
#' in \code{\link{proxTips}}).
#' @param f a function to change a distance into a proximity.
#' @return An object of class \code{orthobasis}. This is a data.frame with
#' Moran's eigenvectors in column, with special attributes:\cr -
#' attr(...,"values"): Moran's index for each vector - attr(...,"weights"):
#' weights of tips; current implementation uses only uniform weights
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso - \code{\link{proxTips}} which computes phylogenetic proximities
#' between tips.\cr
#' 
#' - \code{\link{treePart}} which can compute an orthobasis based on the
#' topology of a phylogeny.\cr
#' @references Peres-Neto, P. (2006) A unified strategy for estimating and
#' controlling spatial, temporal and phylogenetic autocorrelation in ecological
#' models \emph{Oecologica Brasiliensis} \bold{10}: 105-119.\cr
#' 
#' Dray, S.; Legendre, P. and Peres-Neto, P. (2006) Spatial modelling: a
#' comprehensive framework for principal coordinate analysis of neighbours
#' matrices (PCNM) \emph{Ecological Modelling} \bold{196}: 483-493.\cr
#' 
#' Griffith, D. and Peres-Neto, P. (2006) Spatial modeling in ecology: the
#' flexibility of eigenfunction spatial analyses \emph{Ecology} \bold{87}:
#' 2603-2613.\cr
#' @keywords manip
#' @examples
#' 
#' if(require(ape) && require(phylobase)){
#' 
#' ## SIMPLE EXAMPLE ##
#' ## make a tree
#' x <- rtree(50)
#' 
#' ## compute Moran's eigenvectors
#' ME <- me.phylo(x, met="Abouheif")
#' ME
#' 
#' ## plot the 10 first vectors
#' obj <- phylo4d(x, as.data.frame(ME[,1:10]))
#' table.phylo4d(obj, cex.sym=.7, cex.lab=.7)
#' 
#' 
#' 
#' ## REMOVING PHYLOGENETIC AUTOCORRELATION IN A MODEL ##
#' ## use example in ungulates dataset
#' data(ungulates)
#' tre <- read.tree(text=ungulates$tre)
#' plot(tre)
#' 
#' ## look at two traits
#' afbw <- log(ungulates$tab[,1])
#' neonatw <- log((ungulates$tab[,2]+ungulates$tab[,3])/2)
#' names(afbw) <- tre$tip.label
#' names(neonatw) <- tre$tip.label
#' plot(afbw, neonatw) # relationship between traits
#' lm1 <- lm(neonatw~afbw)
#' abline(lm1)
#' 
#' lm1
#' resid1 <- residuals(lm1)
#' orthogram(resid1, tre) # residuals are autocorrelated
#' 
#' ## compute Moran's eigenvectors (ME)
#' myME <- me.phylo(tre, method="Abou")
#' lm2 <- lm(neonatw ~ myME[,1] + afbw) # use for ME as covariable
#' resid2 <- residuals(lm2)
#' orthogram(resid2, tre) # there is no longer phylogenetic autocorrelation
#' 
#' ## see the difference
#' table.phylo4d(phylo4d(tre, cbind.data.frame(resid1, resid2)))
#' }
#' 
#' 
#' @rdname orthobasis
#' @import phylobase ade4
#' @export
orthobasis.phylo <- function(x=NULL, prox=NULL,
                             method=c("patristic","nNodes","oriAbouheif","Abouheif","sumDD"),
                             f=function(x) {1/x} ){
    ## if(!require(phylobase)) stop("phylobase package is not installed")
    ## if(!require(ade4)) stop("ade4 package is not installed")

    ## handle arguments
    method <- match.arg(method)

    if(is.null(prox)){ # have to compute prox
        x <- as(x, "phylo4")
        if (is.character(checkval <- checkPhylo4(x))) stop(checkval)
        W <- proxTips(x, tips="all", method=method, f=f, normalize="row", symmetric=TRUE)
    } else { # prox is provided
        W <- as.matrix(prox)
        if(!is.matrix(W)) stop("W is not a matrix")
        if(ncol(W) != nrow(W)) stop("W is not a square matrix")
         diag(W) <- 0
        W <- 0.5 * (t(W) + W) # re-symmetrization
    }

    n <- nrow(W)


    ## main computation -> call to orthobasis.mat
    res <- orthobasis.mat(W, cnw=FALSE)

    ## build output
    row.names(res) <- rownames(W)
    names(res) <- paste("ME", 1:ncol(res))
    names(attr(res,"values")) <- names(res)
    attr(res,"call") <- match.call()
    attr(res,"class") <- c("orthobasis","data.frame")

    return(res)
} # end orthobasis.phylo





###########
# me.phylo
###########

#' @export
me.phylo <- orthobasis.phylo
