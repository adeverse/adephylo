\encoding{latin1}
\name{palm}
\alias{palm}
\docType{data}
\title{Phylogenetic and quantitative traits of amazonian palm trees}
\description{
This data set describes the phylogeny of 66 amazonian palm trees. It
also gives 7 traits corresponding to these 66 species.
}
\usage{data(palm)}
\format{
\code{palm} is a list containing the 2 following objects: 
\describe{ 
   \item{tre}{is a character string giving the phylogenetic tree in Newick format.} 
   \item{traits}{is a data frame with 66 species (rows) and 7 traits (columns).}
} }
\details{
Variables of \code{palm$traits} are the following ones: \cr
- rord: specific richness with five ordered levels\cr
- h: height in meter (squared transform)\cr
- dqual: diameter at breast height in centimeter with five levels
 \code{sout : subterranean}, \code{ d1(0, 5 cm)}, \code{ d2(5, 15 cm)},
 \code{ d3(15, 30 cm)} and \code{ d4(30, 100 cm)}\cr
- vfruit: fruit volume in \eqn{mm^{3}}{mm^3} (logged transform)\cr
- vgrain: seed volume in \eqn{mm^{3}}{mm^3} (logged transform)\cr
- aire: spatial distribution area (\eqn{km^{2}}{km^2})\cr
- alti: maximum altitude in meter (logged transform)\cr       
}
\source{
This data set was obtained by Clémentine Gimaret-Carpentier\cr
 \email{gimaret@biomserv.univ-lyon1.fr}.
}
\note{
  This dataset replaces the former version in ade4.
}
\examples{
\dontrun{
if(require(ape) && require(phylobase)){

## load data, make a tree and a phylo4d object
data(palm)
tre <- read.tree(text=palm$tre)
rord <- as.integer(palm$traits$rord) # just use this for plotting purpose
traits <- data.frame(rord, palm$traits[,-1])
x <- phylo4d(tre, traits)

## plot data
par(mar=rep(.1,4))
table.phylo4d(x, cex.lab=.6)

## test phylogenetic autocorrelation
if(require(ade4)){
prox <- proxTips(x, method="sumDD")
phylAutoTests <- gearymoran(prox, traits[,-3], nrep=499)
plot(phylAutoTests)
}
}
}
}
\keyword{datasets}
