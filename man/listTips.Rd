\name{listTips}
\alias{listTips}
\title{List tips descendings from all nodes of a tree}
\description{
  The function \code{listTips} lists the tips descending from each node
  of a tree. The tree can be of class \code{\link[ape:read.tree]{phylo}},
  \linkS4class{phylo4} or \linkS4class{phylo4d}.
}
\usage{
listTips(x)
}
\arguments{
  \item{x}{A tree of  class \code{\link[ape:read.tree]{phylo}},
    \linkS4class{phylo4} or \linkS4class{phylo4d}.}
}
\value{
  A list whose components are vectors of named tips for a given node.
}
\author{ Thibaut Jombart \email{tjombart@imperial.ac.uk} }
\seealso{\code{\link{listDD}} which lists the direct descendants for
  each node. \cr
  
  \code{\link{treePart}} which defines partitions of tips
  according to the tree topology.}
\examples{
if(require(ape) & require(phylobase)){
## make a tree
x <- as(rtree(20),"phylo4")
plot(x,show.node=TRUE)
listTips(x)
}
}
\keyword{manip}
