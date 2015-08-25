\name{Misc utils}
\alias{.tipToRoot}
\title{Low-level auxiliary functions for adephylo}
\description{
  These hidden functions are utils for adephylo, used by other
  functions. Regular users can use them as well, but no validity checks
  are performed for the arguments: speed is here favored over safety.
  Most of these functions handle trees inheriting \linkS4class{phylo4}
  class.\cr

  \code{.tipToRoot} finds the set of nodes between a tip and the root of
  a tree.\cr

}
\usage{
.tipToRoot(x, tip, root, include.root=FALSE)
}
\arguments{
  \item{x}{A valid tree of  class \linkS4class{phylo4}.}
  \item{tip}{An integer identifying a tip by its numbers.}
  \item{root}{An integer identifying the root of the tree by its
    number.}
  \item{include.root}{a logical stating whether the root must be
    included as a node of the path from tip to root (TRUE), or not
    (FALSE, default).}
}
\value{
  \code{.tipToRoot}: a vector of named integers identifying nodes.\cr
}
\author{ Thibaut Jombart \email{tjombart@imperial.ac.uk} }
\examples{
if(require(ape) & require(phylobase)){
## make a tree
x <- as(rtree(20),"phylo4")
plot(x,show.node=TRUE)

## .tipToRoot
root <- rootNode(x)
.tipToRoot(x, 1, root)
lapply(1:nTips(x), function(i) .tipToRoot(x, i, root))
}
}
\keyword{manip}
