% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adephylo-package.R
\docType{data}
\name{maples}
\alias{maples}
\title{Phylogeny and quantitative traits of flowers}
\format{
\code{tithonia} is a list containing the 2 following objects : -
tre: a character string giving the phylogenetic tree in Newick format.\cr -
tab: a data frame with 17 species and 31 traits.\cr
}
\description{
This data set describes the phylogeny of 17 flowers as reported by Ackerly
and Donoghue (1998). It also gives 31 traits corresponding to these 17
species.
}
\note{
This dataset replaces the former version in ade4.
}
\examples{


if(require(ape) && require(phylobase)){

data(maples)

## see the tree
tre <- read.tree(text=maples$tre)
plot(tre)
axisPhylo()

## look at two variables
dom <- maples$tab$Dom
bif <- maples$tab$Bif
plot(bif,dom,pch = 20)
abline(lm(dom~bif)) # a strong negative correlation ?
summary(lm(dom~bif))
cor.test(bif,dom)

## look at the two variables onto the phylogeny
temp <- phylo4d(tre, data.frame(dom,bif, row.names=tre$tip.label))
table.phylo4d(temp) # correlation is strongly linked to phylogeny

## use ape's PIC (phylogenetic independent contrasts)
pic.bif <- pic(bif, tre)
pic.dom <- pic(dom, tre)
cor.test(pic.bif, pic.dom) # correlation is no longer significant
}


}
\references{
Ackerly, D. D. and Donoghue, M.J. (1998) Leaf size, sappling
allometry, and Corner's rules: phylogeny and correlated evolution in Maples
(Acer). \emph{American Naturalist}, \bold{152}, 767--791.
}
\keyword{datasets}
