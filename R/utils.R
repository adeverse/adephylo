
#' Internal function for finding path from tips to root
#'
#' @rdname tipToRoot
#'
#' @aliases .tipToRoot
#' @param x A valid tree of class \linkS4class{phylo4}.
#' @param tip An integer identifying a tip by its numbers.
#' @param root An integer identifying the root of the tree by its number.
#' @param include.root a logical stating whether the root must be included as a
#' node of the path from tip to root (TRUE), or not (FALSE, default).
#' @return \code{.tipToRoot}: a vector of named integers identifying nodes.\cr
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @keywords manip
#'
#' @export
#'
#' @examples
#'
#' #' if(require(ape) & require(phylobase)){
#' ## make a tree
#' x <- as(rtree(20),"phylo4")
#' plot(x,show.node=TRUE)
#'
#' ## .tipToRoot
#' root <- rootNode(x)
#' .tipToRoot(x, 1, root)
#' lapply(1:nTips(x), function(i) .tipToRoot(x, i, root))
#' }
#'
.tipToRoot <- function(x, tip, root, include.root=FALSE){
    E <- x@edge
    path <- NULL
    curNode <- tip
    while(curNode != root){
        curNode <- E[(curNode==E[,2]),1] # one node <- its ancestor
        path <- c(path, curNode)
    } # end while

    if(!include.root) {
        path <- path[-length(path)] # exclude the root
    }

    return(getNode(x, path))
} # end tipToRoot





#' Find the shortest path between tips of a tree
#'
#' The function \code{sp.tips} finds the shortest path between tips of a tree,
#' identified as \code{tip1} and \code{tip2}.  This function applies to trees
#' with the class \code{\link[ape:read.tree]{phylo}}, \linkS4class{phylo4} or
#' \linkS4class{phylo4d}. Several tips can be provided at a time.
#'
#' The function checks if there are cases where tip1 and tip2 are the same.
#' These cases are deleted when detected, issuing a warning (unless
#' \code{quiet} is set to TRUE).
#'
#' @param x A tree of class \code{\link[ape:read.tree]{phylo}},
#' \linkS4class{phylo4} or \linkS4class{phylo4d}.
#' @param tip1 A vector of integers identifying tips by their numbers, or a
#' vector of characters identifying tips by their names. Recycled if needed.
#' @param tip2 A vector of integers identifying tips by their numbers, or a
#' vector of characters identifying tips by their names. Recycled if needed.
#' @param useTipNames a logical stating whether the output must be named using
#' tip names in all cases (TRUE), or not (FALSE). If not, names of \code{tip1}
#' and \code{tip2} will be used.
#' @param quiet a logical stating whether a warning must be issued when
#' tip1==tip2, or not (see details).
#' @param include.mrca a logical stating whether the most recent common
#' ancestor shall be included in the returned path (TRUE, default) or not
#' (FALSE).
#' @return A list whose components are vectors of named nodes forming the
#' shortest path between a couple of tips.
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso \code{\link[phylobase]{shortestPath}} which does the same thing as
#' \code{sp.tips}, for any node (internal or tip), but much more slowly. \cr
#' @keywords manip
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' if(require(ape) & require(phylobase)){
#' ## make a tree
#' x <- as(rtree(20),"phylo4")
#' plot(x,show.node=TRUE)
#' ## get shortest path between tip 1 and all other tips.
#' sp.tips(x, "t1", "t2")
#' sp.tips(x, 1, 2:20, TRUE)
#' }
#' }
#'
sp.tips <- function(x, tip1, tip2, useTipNames=FALSE, quiet=FALSE, include.mrca=TRUE){
    ## if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")

    ## some checks
    if (is.character(checkval <- checkPhylo4(x))) stop(checkval)
    t1 <- getNode(x, tip1)
    t2 <- getNode(x, tip2)
    if(any(is.na(c(t1,t2)))) stop("wrong tip specified")
    if(any(c(t1,t2) > nTips(x))) stop("specified nodes are internal nodes")
    if(length(t1) != length(t2)) { # recycle tip1 and tip2
        maxLength <- max(length(t1), length(t2))
        t1 <- rep(t1, length.out=maxLength)
        t2 <- rep(t2, length.out=maxLength)
    }
    toRemove <- (t1==t2)
    if(sum(toRemove)>0) {
        t1 <- t1[!toRemove]
        t2 <- t2[!toRemove]
        if(length(t1)==0) stop("tip1 and tip2 are the same vectors")
        if(!quiet) warning("tip1 and tip2 are sometimes the same; erasing these cases")
    }


    ## some global variables
    N <- nTips(x)
    root <- getNode(x, N+1)
    E <- x@edge
    allTips <- unique(c(t1,t2))


    ##   ## tipToRoot -> call to .tipToRoot
    ##     tipToRoot <- function(E, tip){
    ##         path <- NULL
    ##         curNode <- tip
    ##         while(curNode != root){
    ##             curNode <- E[(curNode==E[,2]),1] # one node <- its ancestor
    ##             path <- c(path, curNode)
    ##         } # end while

    ##         path <- getNode(x, path)
    ##         return(path)
    ##     } # end tipToRoot


    ## function pathTwoTips (takes two path-to-root as args)
    pathTwoTips <- function(path1, path2){
        cpath <- c(path1, rev(path2))
        temp <- factor(cpath, levels=unique(cpath))
        CA <- temp[table(temp)==2][1] # most recent common ancestor (MRCA)
        CA <- as.integer(as.character(CA)) # retrieve integer type
        path1 <- path1[1:(which(path1==CA))] # cut path1 after MRCA (keep MRCA)
        temp <- which(path2==CA)
        if(temp==1) return(path1)
        path2 <- path2[1:(temp-1)] # cut path2 after MRCA (erase MRCA)
        return(c(path1,path2))
    } # end pathTwoTips


    pathTwoTips.no.mrca <- function(path1, path2){
        cpath <- c(path1, rev(path2))
        temp <- intersect(path1, path2)
        res <- setdiff(cpath, temp)
        return(res)
    } # end pathTwoTips



    ## main computations
    allPathToRoot <- lapply(allTips, function(i) .tipToRoot(x, i, root, include.root=TRUE))
    names(allPathToRoot) <- allTips

    allPath1 <- allPathToRoot[as.character(t1)]
    allPath2 <- allPathToRoot[as.character(t2)]

    if(include.mrca) {
        res <- lapply(1:length(allPath1), function(i) pathTwoTips(allPath1[[i]], allPath2[[i]]) )
    } else {
        res <- lapply(1:length(allPath1), function(i) pathTwoTips.no.mrca(allPath1[[i]], allPath2[[i]]) )
        temp.names <- names(res)
        temp <- sapply(res, function(vec) length(vec)>0)
        res[temp] <- lapply(res[temp], function(vec) getNode(x, vec) ) # name the nodes
        names(res) <- temp.names
    }

    if(useTipNames) {
        names(res) <- paste(names(t1), names(t2), sep="-")
    } else {
        names(res) <- paste(t1,t2,sep="-")
    }

    return(res)
} # end sp.tips



# examples
# source("/home/master/dev/adephylo/pkg/R/utils.R")
#phy <- as(rtree(15),"phylo4")
## plot(phy,show.n=T)
## tip1 <- "t1"
## tip2 <- "t2"


## sp.tips(phy, "t1", "t2")
## sp.tips(phy, rep(1,15), 1:15)
## sp.tips(phy, rep(1, 15), 1:15, TRUE)

## heavier tree
# x <- as(rtree(1000), "phylo4")
# system.time(sp.tips(x,1,1:1000))




#' List direct descendants for all nodes of a tree
#'
#' The function \code{listDD} lists the direct descendants from each node of a
#' tree. The tree can be of class \code{\link[ape:read.tree]{phylo}},
#' \linkS4class{phylo4} or \linkS4class{phylo4d}.
#'
#'
#' @param x A tree of class \code{\link[ape:read.tree]{phylo}},
#' \linkS4class{phylo4} or \linkS4class{phylo4d}.
#' @param nameBy a character string indicating whether the returned list must
#' be named by node labels ("label") or by node numbers ("number").
#' @return A list whose components are vectors of named nodes (or tips) for a
#' given internal node.
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso \code{\link{listTips}} which lists the tips descending from each
#' node. \cr
#'
#' \code{\link{treePart}} which defines partitions of tips according to the
#' tree topology.
#' @keywords manip
#'
#' @export
#'
#' @examples
#'
#' if(require(ape) & require(phylobase)){
#' ## make a tree
#' x <- as(rtree(20),"phylo4")
#' plot(x,show.node=TRUE)
#' listDD(x)
#' }
#'
listDD <- function(x, nameBy=c("label","number")){
    ## if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")
    nameBy <- match.arg(nameBy)

    ## check phylo4 object
    if (is.character(checkval <- checkPhylo4(x))) stop(checkval)

    ## computations
    nodIdx <- nTips(x)+1
    nodIdx <- nodIdx:(nodIdx+nNodes(x)-1)
    res <- lapply(nodIdx, function(i) children(x, i))

    if(hasNodeLabels(x) & nameBy=="label") {
        names(res) <- nodeLabels(x)
    } else {
        names(res) <- nodIdx
    }

    return(res)
} # end listDD




