

#' List tips descendings from all nodes of a tree
#'
#' The function \code{listTips} lists the tips descending from each node of a
#' tree. The tree can be of class \code{\link[ape:read.tree]{phylo}},
#' \linkS4class{phylo4} or \linkS4class{phylo4d}.
#'
#'
#' @param x A tree of class \code{\link[ape:read.tree]{phylo}},
#' \linkS4class{phylo4} or \linkS4class{phylo4d}.
#' @return A list whose components are vectors of named tips for a given node.
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso \code{\link{listDD}} which lists the direct descendants for each
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
#' listTips(x)
#' }
#'
listTips <- function(x){
    ## if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")

    ## check phylo4 object
    if (is.character(checkval <- checkPhylo4(x))) stop(checkval)

    ## computations
    nodIdx <- nTips(x)+1
    nodIdx <- nodIdx:(nodIdx+nNodes(x)-1)
    res <- lapply(nodIdx, function(i) descendants(x, i))

    if(hasNodeLabels(x)) {names(res) <- nodeLabels(x)}

    return(res)
} # end listTips





###########
# treePart
###########


#' Define partitions of tips according from a tree
#'
#' The function \code{treePart} defines partitions of tips reflecting the
#' topology of a tree. There are two possible outputs (handled by the argument
#' \code{result}):\cr - \code{basis} mode: each node but the root is translated
#' into a dummy vector having one value for each tip: this value is '1' if the
#' tip descends from this node, and '0' otherwise.\cr - \code{orthobasis}: in
#' this mode, an orthonormal basis is derived from the basis previously
#' mentionned. This orthobasis was proposed in the orthogram (Ollier \emph{et
#' al.} 2006).
#'
#' Orthobasis produced by this function are identical to those stored in the
#' \$Bscores component of deprecated \link[ade4]{phylog} objects, from the ade4
#' package.
#'
#' @param x a tree of class \code{\link[ape:read.tree]{phylo}},
#' \linkS4class{phylo4} or \linkS4class{phylo4d}.
#' @param result a character string specifying the type of result: either a
#' basis of dummy vectors (\code{dummy}), or an orthobasis derived from these
#' dummy vectors (\code{orthobasis}).
#' @return A matrix of numeric vectors (in columns) having one value for each
#' tip (rows).
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso - \code{\link{listDD}} which is called by \code{treePart}.\cr -
#' \code{\link{orthogram}}, which uses by default the orthobasis produced by
#' \code{treePart}.\cr
#' @references Ollier, S., Chessel, D. and Couteron, P. (2005) Orthonormal
#' Transform to Decompose the Variance of a Life-History Trait across a
#' Phylogenetic Tree. \emph{Biometrics}, \bold{62}, 471--477.
#' @keywords manip
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' if(require(ape) & require(phylobase)){
#' ## make a tree
#' x <- as(rtree(10),"phylo4")
#' partition <- treePart(x)
#' partition
#'
#' ## plot the dummy vectors with the tree
#' temp <- phylo4d(x, partition)
#' table.phylo4d(temp, cent=FALSE, scale=FALSE)
#' }
#' }
#'
treePart <- function(x, result=c("dummy", "orthobasis")){
    ## if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")
    result <- match.arg(result)

    ## check phylo4 object
    if (is.character(checkval <- checkPhylo4(x))) stop(checkval)

    n <- nTips(x) # number of tips
    HTU.idx <- (n+1):(n+nNodes(x)) # index of internal nodes (HTU)

    if(!hasNodeLabels(x)) { # node labels will be used after
        nodeLabels(x) <- as.character(HTU.idx)
    }

    ## function coding one dummy vector
    fDum <- function(vec){ # vec is a vector of tip numbers
        dum <- integer(n)
        dum[vec] <- 1
        return(dum)
    }

    ## main computations
    temp <- listTips(x)
    res <- data.frame(lapply(temp,fDum))
    row.names(res) <- tipLabels(x)
    res <- res[,-1, drop=FALSE]

    if(result=="dummy"){
        return(res) # res is a data.frame of dummy vectors
    }



    ## If orthobasis is required ##

    ## Find values 'w' for all nodes
    ##
    ## Notations:
    ## - n: an internal node (HTU)
    ## - Dn: the set of all internal nodes descending from 'n'
    ## - En: the set 'n U Dn' (that is, Dn plus n itself)
    ## - ndd(e): the number of direct descendants from a node 'e'
    ##
    ## Then the values 'w' are computed as:
    ##
    ## w(n) = sum_{e \in En} lgamma( ndd(e) + 1)
    ##

    listDDx <- listDD(x)

    nbOfDD <- sapply(listDDx, length) # nb of DD for each node
    names(nbOfDD) <- HTU.idx # used to match the results of Dn

    findAlldHTU <- function(node){ # find all HTU descending from a node
        res <- descendants(x, node, type="all") # tips and HTU
        res <- res[res > n] # only HTU (here, just node numbers are kept
        if(length(res)==0) return(NULL)
        return(res)
    }


    listAlldHTU <- lapply(HTU.idx, function(node) c(node,findAlldHTU(node))) # ='Dn': for each HTU, list all HTU descending from it

    w <- sapply(listAlldHTU, function(e) sum(lgamma(nbOfDD[as.character(e)]+1))) # w(n)
    ## from now on, 'w' stores the w(n) values.

    ## add dummy vectors for tips
    res <- cbind(diag(1, n), root=rep(1,n), res) # sorted from first tip to last node
    colnames(res) <- 1:(nTips(x) + nNodes(x))
    valDum <- c(rep(-1, n), w) # dummy vectors of tips are given a negative value
    ## note: valDum is the w values for all nodes, sorted from first tip to last node

    ## Discard dummy vectors with lowest valDum (value of dummy vectors, w).
    ## -> for each node, a dummy vector associated to its DD is removed
    ## this one is that with the lowest valDum.

    discardOneDum <- function(node, DDnode){ # node is a node label, not a node number
        if(length(DDnode)==1) return(NULL)
        val <- valDum[DDnode]
        toRemove <- which.min(val)
        keptDD <- DDnode[-toRemove]
        return(keptDD)
    } # end discardOneDum

    dumToKeep <- lapply(1:length(listDDx), function(i) discardOneDum(i, listDDx[[i]]))
    dumToKeep <- unlist(dumToKeep) # contains indices of kept dummy vectors

    res <- res[dumToKeep] # retained dummy vectors
    res <- res[,order(valDum[dumToKeep], decreasing=TRUE)] # reorder vectors by decreasing w

    ## orthonormalization
    res <- cbind(root=rep(1,n), res) # for centring: vectors will be orthogonal to 1_n
    res <- qr.Q(qr(res)) # Gram-Schmidt orthogonalization
    res <- res[,-1] # keep only centred vectors; orthogonal for identity
    res <- res * sqrt(n) # render vectors orthogonal for 1/n

    rownames(res) <- tipLabels(x)
    colnames(res) <- paste("V",1:ncol(res))

    return(as.data.frame(res))

} # end treePart

