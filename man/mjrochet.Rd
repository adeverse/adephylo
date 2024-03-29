% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adephylo-package.R
\docType{data}
\name{mjrochet}
\alias{mjrochet}
\title{Phylogeny and quantitative traits of teleos fishes}
\format{
\code{mjrochet} is a list containing the 2 following objects :
\describe{ \item{tre}{is a character string giving the phylogenetic tree in
Newick format.} \item{tab}{is a data frame with 49 rows and 7 traits.} }
}
\description{
This data set describes the phylogeny of 49 teleos fishes as reported by
Rochet et al. (2000). It also gives life-history traits corresponding to
these 49 species.
}
\details{
Variables of \code{mjrochet$tab} are the following ones : tm (age at
maturity (years)), lm (length at maturity (cm)), l05 (length at 5 per cent
survival (cm)), t05 (time to 5 per cent survival (years)), fb (slope of the
log-log fecundity-length relationship), fm (fecundity the year of maturity),
egg (volume of eggs (\eqn{mm^{3}}{mm^3})).
}
\note{
This dataset replaces the former version in ade4.
}
\examples{


if(require(ape) && require(phylobase)){

data(mjrochet)
tre <- read.tree(text=mjrochet$tre) # make a tree
traits <- log((mjrochet$tab))

## build a phylo4d
mjr <- phylo4d(tre, traits)

## see data
table.phylo4d(mjr,cex.lab=.5,show.node=FALSE,symb="square")

## perform Abouheif's test for each trait
mjr.tests <- abouheif.moran(mjr, nrep=499)
mjr.tests


}

}
\references{
Rochet, M. J., Cornillon, P-A., Sabatier, R. and Pontier, D.
(2000) Comparative analysis of phylogenic and fishing effects in life
history patterns of teleos fishes.  \emph{Oikos}, \bold{91}, 255--270.
}
\keyword{datasets}
