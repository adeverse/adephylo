###########
# distTips
###########


#' Compute some phylogenetic distance between tips
#' 
#' The function \code{distTips} computes a given distance between a set of tips
#' of a phylogeny. A vector of tips is supplied: distances between all possible
#' pairs of these tips are computed.  The distances are computed from the
#' shortest path between the tips. Several distances can be used, defaulting to
#' the sum of branch lengths (see argument \code{method}).
#' 
#' An option (enabled by default) allows computations to be run using compiled
#' C code, which is much faster than pure R code. In this case, a matrix of all
#' pairwise distances is returned (i.e., \code{tips} argument is ignored).
#' 
#' \code{Abouheif} distance refers to the phylogenetic distance underlying the
#' test of Abouheif (see references). Let P be the set of all the nodes in the
#' path going from \code{node1} to \code{node2}. Let DDP be the number of
#' direct descendants from each node in P. Then, the so-called 'Abouheif'
#' distance is the product of all terms in DDP.\cr
#' 
#' \code{sumDD} refers to a phylogenetic distance quite similar to that of
#' Abouheif. We consider the same sets P and DDP. But instead of computing the
#' product of all terms in DDP, this distance computes the sum of all terms in
#' DDP.
#' 
#' @param x a tree of class \code{\link[ape:read.tree]{phylo}},
#' \linkS4class{phylo4} or \linkS4class{phylo4d}.
#' @param tips A vector of integers identifying tips by their numbers, or a
#' vector of characters identifying tips by their names. Distances will be
#' computed between all possible pairs of tips.
#' @param method a character string (full or abbreviated without ambiguity)
#' specifying the method used to compute distances ; possible values are:\cr -
#' \code{patristic}: patristic distance, i.e. sum of branch lengths \cr -
#' \code{nNodes}: number of nodes on the path between the nodes \cr -
#' \code{Abouheif}: Abouheif's distance (see details) \cr - \code{sumDD}: sum
#' of direct descendants of all nodes on the path (see details) \cr
#' @param useC a logical indicating whether computations should be performed
#' using compiled C code (TRUE, default), or using a pure R version (FALSE). C
#' version is several orders of magnitude faster, and R version is kept for
#' backward compatibility.
#' @return An object of class \code{dist}, containing phylogenetic distances.
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso \code{\link{distTips}} which computes several phylogenetic
#' distances between tips.
#' @references Pavoine, S.; Ollier, S.; Pontier, D. & Chessel, D. (2008)
#' Testing for phylogenetic signal in life history variable: Abouheif's test
#' revisited. \emph{Theoretical Population Biology}: \bold{73}, 79-91.
#' @keywords manip
#' @examples
#' 
#' if(require(ape) & require(phylobase)){
#' ## make a tree
#' x <- as(rtree(10),"phylo4")
#' plot(x, show.node=TRUE)
#' axisPhylo()
#' ## compute different distances
#' distTips(x, 1:3)
#' distTips(x, 1:3, "nNodes")
#' distTips(x, 1:3, "Abouheif")
#' distTips(x, 1:3, "sumDD")
#' 
#' ## compare C and pure R code outputs
#' x <- rtree(10)
#' all.equal(as.matrix(distTips(x)), as.matrix(distTips(x, useC=FALSE)))
#' all.equal(as.matrix(distTips(x, meth="nNode")),
#'    as.matrix(distTips(x, meth="nNode", useC=FALSE)))
#' all.equal(as.matrix(distTips(x, meth="Abou")),
#'    as.matrix(distTips(x, meth="Abou", useC=FALSE)))
#' all.equal(as.matrix(distTips(x, meth="sumDD")),
#'    as.matrix(distTips(x, meth="sumDD", useC=FALSE)))
#' 
#' ## compare speed
#' x <- rtree(50)
#' tim1 <- system.time(distTips(x, useC=FALSE)) # old pure R version
#' tim2 <- system.time(distTips(x)) # new version using C
#' tim1[c(1,3)]/tim2[c(1,3)] # C is about a thousand time faster in this case
#' }
#' 
#' @useDynLib adephylo
#' @import phylobase
#' @export distTips
distTips <- function(x, tips="all",
                      method=c("patristic","nNodes","Abouheif","sumDD"), useC=TRUE){

    ## if(!require(phylobase)) stop("phylobase package is not installed")

    if(useC){
        tre <- as(x, "phylo")
        n <- as.integer(nTips(tre))
        resSize <- as.integer(n*(n-1)/2)
        res <- double(resSize)
        method <- match.arg(method)
        method <- match(method, c("patristic","nNodes","Abouheif","sumDD"))
        if(is.null(tre$edge.length)){
            tre$edge.length <- as.double(rep(1, nrow(tre$edge)))
        }

        temp <- .C("distalltips", as.integer(tre$edge[,1]), as.integer(tre$edge[,2]), as.double(tre$edge.length), nrow(tre$edge), n, res, resSize, as.integer(method), PACKAGE="adephylo")
        res <- temp[[6]]

        class(res) <- "dist"
        attr(res, "Size") <- nTips(tre)
        attr(res, "Diag") <- FALSE
        attr(res, "Upper") <- FALSE
        attr(res, "method") <- paste("Phylogenetic: ",method,sep="")
        attr(res, "call") <- match.call()
        attr(res, "Labels") <- tre$tip.label
    } else {

        ## handle arguments
        x <- as(x, "phylo4")
        method <- match.arg(method)
        N <- nTips(x)
        if(tips[1]=="all") { tips <- 1:N }
        tips <- getNode(x, tips)
        tips.names <- names(tips)

        ## some checks
        if (is.character(checkval <- checkPhylo4(x))) stop(checkval)
        if(any(is.na(tips))) stop("wrong tips specified")

        ## create all couples of observations
        findAllPairs <- function(vec){
            res <- list(i=NULL,j=NULL)
            k <- 0
            for(i in 1:(length(vec)-1)){
                for(j in (i+1):length(vec)){
                    k <- k+1
                    res[[1]][k] <- i
                    res[[2]][k] <- j
                }
            }
            res <- data.frame(res)
            return(res)
        }

        allPairs <- findAllPairs(tips) # this contains all possible pairs of tips

        ## get the shortest path between all pairs of tips
        if(method != "patristic") {
            allPath <- sp.tips(x, allPairs$i, allPairs$j, useTipNames=TRUE, quiet=TRUE)
        } else {
            allPath <- sp.tips(x, allPairs$i, allPairs$j, useTipNames=TRUE, quiet=TRUE,
                               include.mrca=FALSE)
        }

        ## compute distances
        if(method=="patristic"){
            if(!hasEdgeLength(x)) stop("x does not have branch length")
            ## add tip1 and tip2 to the paths, so that these edges are counted
            allPath.names <- names(allPath)
            allPath <- lapply(1:length(allPath), function(i)
                              c(allPath[[i]], allPairs[i,1], allPairs[i,2]) )
            names(allPath) <- allPath.names

            edge.idx <- lapply(allPath, function(e) getEdge(x, e) ) # list of indices of edges
            allEdgeLength <- edgeLength(x)
            res <- lapply(edge.idx, function(idx) sum(allEdgeLength[idx], na.rm=TRUE) )
        } # end patristic

        if(method=="nNodes"){
            res <- lapply(allPath, length)
        } # end nNodes

        if(method=="Abouheif"){
            E <- x@edge
            f1 <- function(onePath){ # computes product of dd for one path
                temp <- table(E[,1])[as.character(onePath)] # number of dd per node
                return(prod(temp))
            }
            res <- lapply(allPath, f1)
        } # end Abouheif

        if(method=="sumDD"){
            E <- x@edge
            f1 <- function(onePath){ # computes sum of dd for one path
                temp <- table(E[,1])[as.character(onePath)] # number of dd per node
                return(sum(temp))
            }
            res <- lapply(allPath, f1)
        } # end sumDD

        ## convert res to a dist object
        res <- unlist(res)
        class(res) <- "dist"
        attr(res, "Size") <- length(tips)
        attr(res, "Diag") <- FALSE
        attr(res, "Upper") <- FALSE
        attr(res, "method") <- paste("Phylogenetic: ",method,sep="")
        attr(res, "call") <- match.call()
        attr(res, "Labels") <- tips.names
    }
    return(res)

} # end distTips







###########
# distRoot
###########


#' Compute the distance of tips to the root
#' 
#' The function \code{distRoot} computes the distance of a set of tips to the
#' root. Several distances can be used, defaulting to the sum of branch
#' lengths.
#' 
#' \code{Abouheif} distance refers to the phylogenetic distance underlying the
#' test of Abouheif (see references). Let P be the set of all the nodes in the
#' path going from \code{node1} to \code{node2}. Let DDP be the number of
#' direct descendants from each node in P. Then, the so-called 'Abouheif'
#' distance is the product of all terms in DDP.\cr
#' 
#' \code{sumDD} refers to a phylogenetic distance quite similar to that of
#' Abouheif. We consider the same sets P and DDP. But instead of computing the
#' product of all terms in DDP, this distance computes the sum of all terms in
#' DDP.
#' 
#' @param x a tree of class \code{\link[ape:read.tree]{phylo}},
#' \linkS4class{phylo4} or \linkS4class{phylo4d}.
#' @param tips A vector of integers identifying tips by their numbers, or a
#' vector of characters identifying tips by their names.
#' @param method a character string (full or abbreviated without ambiguity)
#' specifying the method used to compute distances ; possible values are:\cr -
#' \code{patristic}: patristic distance, i.e. sum of branch lengths \cr -
#' \code{nNodes}: number of nodes on the path between the nodes \cr -
#' \code{Abouheif}: Abouheif's distance (see details) \cr - \code{sumDD}: sum
#' of direct descendants of all nodes on the path (see details) \cr
#' @return A numeric vector containing one distance value for each tip.
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso \code{\link{distTips}} which computes the same phylogenetic
#' distances, but between tips.
#' @references Pavoine, S.; Ollier, S.; Pontier, D. & Chessel, D. (2008)
#' Testing for phylogenetic signal in life history variable: Abouheif's test
#' revisited. \emph{Theoretical Population Biology}: \bold{73}, 79-91.
#' @keywords manip
#' @examples
#' 
#' if(require(ape) & require(phylobase)){
#' ## make a tree
#' x <- as(rtree(50),"phylo4")
#' ## compute 4 different distances
#' met <- c("patristic","nNodes","Abouheif","sumDD")
#' D <- lapply(met, function(e) distRoot(x, method=e) )
#' names(D) <- met
#' D <- as.data.frame(D)
#' 
#' ## plot these distances along with the tree
#' temp <- phylo4d(x, D)
#' table.phylo4d(temp, show.node=FALSE, cex.lab=.6)
#' }
#' 
#' @import phylobase
#' @export distRoot
distRoot <- function(x, tips="all", method=c("patristic","nNodes","Abouheif","sumDD") ){
    ## if(!require(phylobase)) stop("phylobase package is not installed")

    ## handle arguments
    x <- as(x, "phylo4")
    method <- match.arg(method)
    N <- nTips(x)
    if(tips[1]=="all") { tips <- 1:N }
    tips <- getNode(x, tips)
    tips.names <- names(tips)
    x <- as(x, "phylo4")
    root <- getNode(x, N+1) # so that we have a named node

    ## some checks
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    if(any(is.na(tips))) stop("wrong tips specified")


    ## main computations

    ## get path from root to tops
    allPath <- lapply(tips, function(tip) .tipToRoot(x, tip, root, include.root = TRUE))

    ## compute distances
    if(method=="patristic"){
        if(!hasEdgeLength(x)) stop("x does not have branch length")
        ## add the concerned tips to the paths, so that these edges are counted
        allPath.names <- names(allPath)
        allPath <- lapply(1:length(allPath), function(i) c(allPath[[i]], tips[i]) )
        names(allPath) <- allPath.names

        edge.idx <- lapply(allPath, function(e) getEdge(x, e) ) # list of indices of edges
        allEdgeLength <- edgeLength(x)
        res <- sapply(edge.idx, function(idx) sum(allEdgeLength[idx], na.rm=TRUE) )
    } # end patristic

    if(method=="nNodes"){
        res <- sapply(allPath, length)
    } # end nNodes

    if(method=="Abouheif"){
        E <- x@edge
        f1 <- function(onePath){ # computes product of dd for one path
            temp <- table(E[,1])[as.character(onePath)] # number of dd per node
            return(prod(temp))
        }

        res <- sapply(allPath, f1)
    } # end Abouheif

    if(method=="sumDD"){
        E <- x@edge
        f1 <- function(onePath){ # computes sum of dd for one path
            temp <- table(E[,1])[as.character(onePath)] # number of dd per node
            return(sum(temp))
        }

        res <- sapply(allPath, f1)
    } # end sumDD


    ## the output is a named numeric vector
    return(res)
} # end distRoot
