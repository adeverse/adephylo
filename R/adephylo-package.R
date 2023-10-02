

#' The adephylo package
#' 
#' This package is devoted to exploratory analysis of phylogenetic comparative
#' data. It re-implements and extends phylogenetic procedures from the
#' \code{ade4} package (which are now deprecated).\cr
#' 
#' Comparative data (phylogeny+traits) are handled as \linkS4class{phylo4d}
#' objects, a canonical class implemented by the \code{phylobase} package.
#' Trees are handled as \code{\link[ape:read.tree]{phylo}} objects (from the
#' \code{ape} package) or as \linkS4class{phylo4} objects (\code{phylobase}'s
#' extension of \code{phylo} objects).\cr
#' 
#' Main functionalities of \code{adephylo} are summarized below.\cr
#' 
#' === TOPOLOGICAL INFORMATION ===\cr Several functions allow one to retrieve
#' topological information from a tree; such information can be used, for
#' instance, as a basis to compute distances or proximities between tips.\cr
#' 
#' - \code{\link{listDD}}: lists the direct descendants from each node of a
#' tree.\cr
#' 
#' - \code{\link{listTips}}: lists the tips descending from each node of a
#' tree.\cr
#' 
#' - \code{\link{.tipToRoot}}: finds the set of nodes between a tip and the
#' root of a tree.\cr
#' 
#' - \code{\link{sp.tips}}: finds the shortest path between tips of a tree.\cr
#' 
#' - \code{\link{treePart}}: defines partitions of tips reflecting the topology
#' of a tree. This function can output non-independent dummy vectors, or
#' alternatively an orthonormal basis used by the orthogram procedure.\cr
#' 
#' === PHYLOGENETIC PROXIMITIES/DISTANCES ===\cr Several phylogenetic
#' proximities and distances are implemented. Auxiliary function easing the
#' computation of other distances/proximities are also provided:\cr
#' 
#' - \code{\link{distRoot}}: computes different distances of a set of tips to
#' the root.\cr
#' 
#' - \code{\link{distTips}}: computes different pairwise distances in a set of
#' tips.\cr
#' 
#' - \code{\link{proxTips}}: computes different proximities between a set of
#' tips.\cr
#' 
#' === MEASURES/TESTS OF PHYLOGENETIC AUTOCORRELATION ===\cr Several procedures
#' allow one to measure, and/or test phylogenetic signal in biological
#' traits:\cr
#' 
#' - \code{\link{abouheif.moran}}: performs Abouheif's test, designed to detect
#' phylogenetic autocorrelation in a quantitative trait. This implementation is
#' not based on original heuristic procedure, but on the exact formulation
#' proposed by Pavoine et al. (2008), showing that the test is in fact a
#' Moran's index test. This implementation further extends the procedure by
#' allowing any measure of phylogenetic proximity (5 are proposed).\cr
#' 
#' - \code{\link{orthogram}}: performs the orthonormal decomposition of
#' variance of a quantitative variable on an orthonormal basis as in Ollier et
#' al. (2005). It also returns the results of five non parametric tests
#' associated to the variance decomposition.\cr
#' 
#' - \code{\link{moran.idx}}: computes Moran's index of autocorrelation given a
#' variable and a matrix of proximities among observations (no test).\cr
#' 
#' === MODELLING/INVESTIGATION OF PHYLOGENETIC SIGNAL ===\cr Rather than
#' testing or measuring phylogenetic autocorrelation, these procedures can be
#' used for further investigation of phylogenetic signal. Some, like
#' \code{\link{me.phylo}}, can be used to remove phylogenetic autocorrelation.
#' Others can be used to understand the nature of this autocorrelation (i.e.,
#' to ascertain which traits and tips are concerned by phylogenetic
#' non-independence).\cr
#' 
#' - \code{\link{me.phylo}}/\code{\link{orthobasis.phylo}}: these synonymous
#' functions compute Moran's eigenvectors (ME) associated to a tree. These
#' vectors model different observable phylogenetic signals. They can be used as
#' covariables to remove phylogenetic autocorrelation from data.\cr
#' 
#' - \code{\link{orthogram}}: the orthogram mentioned above also provides a
#' description of how biological variability is structured on a phylogeny.\cr
#' 
#' - \code{\link{ppca}}: performs a phylogenetic Principal Component Analysis
#' (pPCA, Jombart et al. 2010). This multivariate method investigates
#' phylogenetic patterns in a set of quantitative traits.\cr
#' 
#' === GRAPHICS ===\cr Some plotting functions are proposed, most of them being
#' devoted to representing phylogeny and a quantitative information at the same
#' time.\cr
#' 
#' - \code{\link{table.phylo4d}}: fairly customisable way of representing
#' traits onto the tips of a phylogeny. Several traits can be plotted in a
#' single graphic.\cr
#' 
#' - \code{\link{bullseye}}: an alternative to \code{\link{table.phylo4d}}
#' based on fan-like representation, better for large trees.\cr
#' 
#' - \code{\link{scatter.ppca}}, \code{\link{screeplot.ppca}},
#' \code{\link{plot.ppca}}: several plots associated to a phylogenetic
#' principal component analysis (see \code{\link{ppca}}).\cr
#' 
#' === DATASETS ===\cr Several datasets are also proposed. Some of these
#' datasets replace former version from \code{ade4}, which are now deprecated.
#' Here is a list of available datasets: \code{\link{carni19}},
#' \code{\link{carni70}}, \code{\link{lizards}}, \code{\link{maples}},
#' \code{\link{mjrochet}}, \code{\link{palm}}, \code{\link{procella}},
#' \code{\link{tithonia}}, and \code{\link{ungulates}}.\cr
#' 
#' To cite adephylo, please use the reference given by
#' \code{citation("adephylo")}.
#' 
#' \tabular{ll}{ Package: \tab adephylo\cr Type: \tab Package\cr License: \tab GPL (>=2) }
#' 
#' @name adephylo-package
#' @aliases adephylo-package adephylo
#' @docType package
#' @author Thibaut Jombart <tjombart@@imperial.ac.uk>\cr with contributions
#' Stephane Dray <stephane.dray@@univ-lyon1.fr> and Anders Ellern Bilgrau
#' <abilgrau@@math.aau.dk>. \cr Parts of former code from \code{ade4} by Daniel
#' Chessel and Sebastien Ollier.
#' @seealso The \code{ade4} package for multivariate analysis.
#' @keywords manip multivariate
NULL





#' Phylogeny and quantative trait of carnivora
#' 
#' This data set describes the phylogeny of carnivora as reported by
#' Diniz-Filho et al. (1998). It also gives the body mass of these 19 species.
#' 
#' 
#' @name carni19
#' @docType data
#' @format \code{carni19} is a list containing the 2 following objects :
#' \describe{ \item{tre}{is a character string giving the phylogenetic tree in
#' Newick format.} \item{bm}{is a numeric vector which values correspond to the
#' body mass of the 19 species (log scale).} }
#' @note This dataset replaces the former version in ade4.
#' @source Diniz-Filho, J. A. F., de Sant'Ana, C.E.R. and Bini, L.M. (1998) An
#' eigenvector method for estimating phylogenetic inertia. \emph{Evolution},
#' \bold{52}, 1247--1262.
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' if(require(ape) && require(phylobase)){
#' 
#' data(carni19)
#' tre <- read.tree(text=carni19$tre)
#' x <- phylo4d(tre, data.frame(carni19$bm))
#' table.phylo4d(x, ratio=.5, center=FALSE)
#' }
#' }
#' 
NULL





#' Phylogeny and quantitative traits of carnivora
#' 
#' This data set describes the phylogeny of 70 carnivora as reported by
#' Diniz-Filho and Torres (2002). It also gives the geographic range size and
#' body size corresponding to these 70 species.
#' 
#' 
#' @name carni70
#' @docType data
#' @format \code{carni70} is a list containing the 2 following objects:
#' \describe{ \item{tre}{is a character string giving the phylogenetic tree in
#' Newick format.  Branch lengths are expressed as divergence times (millions
#' of years)} \item{tab}{is a data frame with 70 species and two traits: size
#' (body size (kg)) ; range (geographic range size (km)).} }
#' @note This dataset replaces the former version in ade4.
#' @source Diniz-Filho, J. A. F., and N. M. Torres. (2002) Phylogenetic
#' comparative methods and the geographic range size-body size relationship in
#' new world terrestrial carnivora. \emph{Evolutionary Ecology}, \bold{16},
#' 351--367.
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' if(require(ape) && require(phylobase)){
#' 
#' data(carni70)
#' rownames(carni70$tab) <- gsub("_", ".", rownames(carni70$tab))
#' tre <- read.tree(text=carni70$tre)
#' x <- phylo4d(tre, carni70$tab)
#' table.phylo4d(x)
#' 
#' par(mar=rep(.1,4))
#' table.phylo4d(x,cex.lab=.5, show.n=FALSE, ratio=.5)
#' 
#' 
#' ## transform size in log and test for a phylogenetic signal
#' size <- log(carni70$tab)[,1]
#' names(size) <- row.names(carni70$tab)
#' orthogram(size, tre)
#' 
#' ## transform range and test for a phylogenetic signal
#' yrange <- scale(carni70$tab)[,2]
#' names(yrange) <- row.names(carni70$tab)
#' orthogram(yrange, tre)
#' }
#' }
#' 
NULL





#' Phylogeny and quantitative traits of lizards
#' 
#' This data set describes the phylogeny of 18 lizards as reported by Bauwens
#' and D\'iaz-Uriarte (1997). It also gives life-history traits corresponding
#' to these 18 species.
#' 
#' Variables of \code{lizards$traits} are the following ones : mean.L (mean
#' length (mm)), matur.L (length at maturity (mm)), max.L (maximum length
#' (mm)), hatch.L (hatchling length (mm)), hatch.m (hatchling mass (g)),
#' clutch.S (Clutch size), age.mat (age at maturity (number of months of
#' activity)), clutch.F (clutch frequency).
#' 
#' @name lizards
#' @docType data
#' @format \code{lizards} is a list containing the 3 following objects :
#' \describe{ \item{traits}{is a data frame with 18 species and 8 traits.}
#' \item{hprA}{is a character string giving the phylogenetic tree (hypothesized
#' phylogenetic relationships based on immunological distances) in Newick
#' format.} \item{hprB}{is a character string giving the phylogenetic tree
#' (hypothesized phylogenetic relationships based on morphological
#' characteristics) in Newick format.} }
#' @note This dataset replaces the former version in ade4.
#' @references Bauwens, D., and D\'iaz-Uriarte, R. (1997) Covariation of
#' life-history traits in lacertid lizards: a comparative study.
#' \emph{American Naturalist}, \bold{149}, 91--111.
#' 
#' See a data description at \url{http://pbil.univ-lyon1.fr/R/pdf/pps063.pdf}
#' (in French).
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' if(require(ape) && require(phylobase)){
#' 
#' ## see data
#' data(lizards)
#' liz.tr <- read.tree(tex=lizards$hprA) # make a tree
#' liz <- phylo4d(liz.tr, lizards$traits) # make a phylo4d object
#' table.phylo4d(liz)
#' 
#' ## compute and plot principal components
#' if(require(ade4)){
#' liz.pca1 <- dudi.pca(lizards$traits, cent=TRUE,
#'    scale=TRUE, scannf=FALSE, nf=2) # PCA of traits
#' myPC <- phylo4d(liz.tr, liz.pca1$li) # store PC in a phylo4d object
#' varlab <- paste("Principal \ncomponent", 1:2) # make labels for PCs
#' table.phylo4d(myPC, ratio=.8, var.lab=varlab) # plot the PCs
#' add.scatter.eig(liz.pca1$eig,2,1,2,posi="topleft", inset=c(0,.15))
#' title("Phylogeny and the principal components")
#' 
#' ## compute a pPCA ##
#' ## remove size effect
#' temp <- lapply(liz.pca1$tab, function(e) residuals(lm(e~-1+liz.pca1$li[,1])) )
#' temp <- data.frame(temp)
#' row.names(temp) <- tipLabels(liz)
#' 
#' ## build corresponding phylo4d object
#' liz.noSize <- phylo4d(liz.tr, temp)
#' ppca1 <- ppca(liz.noSize, method="Abouheif", scale=FALSE, scannf=FALSE)
#' plot(ppca1)
#' 
#' }
#' }
#' }
#' 
NULL





#' Phylogeny and quantitative traits of flowers
#' 
#' This data set describes the phylogeny of 17 flowers as reported by Ackerly
#' and Donoghue (1998). It also gives 31 traits corresponding to these 17
#' species.
#' 
#' 
#' @name maples
#' @docType data
#' @format \code{tithonia} is a list containing the 2 following objects : -
#' tre: a character string giving the phylogenetic tree in Newick format.\cr -
#' tab: a data frame with 17 species and 31 traits.\cr
#' @note This dataset replaces the former version in ade4.
#' @references Ackerly, D. D. and Donoghue, M.J. (1998) Leaf size, sappling
#' allometry, and Corner's rules: phylogeny and correlated evolution in Maples
#' (Acer). \emph{American Naturalist}, \bold{152}, 767--791.
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' if(require(ape) && require(phylobase)){
#' 
#' data(maples)
#' 
#' ## see the tree
#' tre <- read.tree(text=maples$tre)
#' plot(tre)
#' axisPhylo()
#' 
#' ## look at two variables
#' dom <- maples$tab$Dom
#' bif <- maples$tab$Bif
#' plot(bif,dom,pch = 20)
#' abline(lm(dom~bif)) # a strong negative correlation ?
#' summary(lm(dom~bif))
#' cor.test(bif,dom)
#' 
#' ## look at the two variables onto the phylogeny
#' temp <- phylo4d(tre, data.frame(dom,bif, row.names=tre$tip.label))
#' table.phylo4d(temp) # correlation is strongly linked to phylogeny
#' 
#' ## use ape's PIC (phylogenetic independent contrasts)
#' pic.bif <- pic(bif, tre)
#' pic.dom <- pic(dom, tre)
#' cor.test(pic.bif, pic.dom) # correlation is no longer significant
#' }
#' }
#' 
NULL










#' Phylogeny and quantitative traits of teleos fishes
#' 
#' This data set describes the phylogeny of 49 teleos fishes as reported by
#' Rochet et al. (2000). It also gives life-history traits corresponding to
#' these 49 species.
#' 
#' Variables of \code{mjrochet$tab} are the following ones : tm (age at
#' maturity (years)), lm (length at maturity (cm)), l05 (length at 5 per cent
#' survival (cm)), t05 (time to 5 per cent survival (years)), fb (slope of the
#' log-log fecundity-length relationship), fm (fecundity the year of maturity),
#' egg (volume of eggs (\eqn{mm^{3}}{mm^3})).
#' 
#' @name mjrochet
#' @docType data
#' @format \code{mjrochet} is a list containing the 2 following objects :
#' \describe{ \item{tre}{is a character string giving the phylogenetic tree in
#' Newick format.} \item{tab}{is a data frame with 49 rows and 7 traits.} }
#' @note This dataset replaces the former version in ade4.
#' @references Rochet, M. J., Cornillon, P-A., Sabatier, R. and Pontier, D.
#' (2000) Comparative analysis of phylogenic and fishing effects in life
#' history patterns of teleos fishes.  \emph{Oikos}, \bold{91}, 255--270.
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' if(require(ape) && require(phylobase)){
#' 
#' data(mjrochet)
#' tre <- read.tree(text=mjrochet$tre) # make a tree
#' traits <- log((mjrochet$tab))
#' 
#' ## build a phylo4d
#' mjr <- phylo4d(tre, traits)
#' 
#' ## see data
#' table.phylo4d(mjr,cex.lab=.5,show.node=FALSE,symb="square")
#' 
#' ## perform Abouheif's test for each trait
#' mjr.tests <- abouheif.moran(mjr, nrep=499)
#' mjr.tests
#' 
#' }
#' }
#' 
NULL















#' Phylogenetic and quantitative traits of amazonian palm trees
#' 
#' This data set describes the phylogeny of 66 amazonian palm trees. It also
#' gives 7 traits corresponding to these 66 species.
#' 
#' Variables of \code{palm$traits} are the following ones: \cr - rord: specific
#' richness with five ordered levels\cr - h: height in meter (squared
#' transform)\cr - dqual: diameter at breast height in centimeter with five
#' levels \code{sout : subterranean}, \code{ d1(0, 5 cm)}, \code{ d2(5, 15
#' cm)}, \code{ d3(15, 30 cm)} and \code{ d4(30, 100 cm)}\cr - vfruit: fruit
#' volume in \eqn{mm^{3}}{mm^3} (logged transform)\cr - vgrain: seed volume in
#' \eqn{mm^{3}}{mm^3} (logged transform)\cr - aire: spatial distribution area
#' (\eqn{km^{2}}{km^2})\cr - alti: maximum altitude in meter (logged
#' transform)\cr
#' 
#' @name palm
#' @docType data
#' @format \code{palm} is a list containing the 2 following objects: \describe{
#' \item{tre}{is a character string giving the phylogenetic tree in Newick
#' format.} \item{traits}{is a data frame with 66 species (rows) and 7 traits
#' (columns).} }
#' @note This dataset replaces the former version in ade4.
#' @source This data set was obtained by Clementine Gimaret-Carpentier.
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' if(require(ape) && require(phylobase)){
#' 
#' ## load data, make a tree and a phylo4d object
#' data(palm)
#' tre <- read.tree(text=palm$tre)
#' rord <- as.integer(palm$traits$rord) # just use this for plotting purpose
#' traits <- data.frame(rord, palm$traits[,-1])
#' x <- phylo4d(tre, traits)
#' 
#' ## plot data
#' par(mar=rep(.1,4))
#' table.phylo4d(x, cex.lab=.6)
#' 
#' ## test phylogenetic autocorrelation
#' if(require(ade4)){
#' prox <- proxTips(x, method="sumDD")
#' phylAutoTests <- gearymoran(prox, traits[,-3], nrep=499)
#' plot(phylAutoTests)
#' }
#' }
#' }
#' 
NULL





#' Phylogeny and quantitative traits of birds
#' 
#' This data set describes the phylogeny of 19 birds as reported by Bried et
#' al. (2002). It also gives 6 traits corresponding to these 19 species.
#' 
#' Variables of \code{procella$traits} are the following ones: \cr - site.fid:
#' a numeric vector that describes the percentage of site fidelity\cr -
#' mate.fid: a numeric vector that describes the percentage of mate fidelity\cr
#' - mass: an integer vector that describes the adult body weight (g)\cr - ALE:
#' a numeric vector that describes the adult life expectancy (years)\cr - BF: a
#' numeric vector that describes the breeding frequencies\cr - col.size: an
#' integer vector that describes the colony size (no nests monitored)
#' 
#' @name procella
#' @docType data
#' @format \code{procella} is a list containing the 2 following objects:
#' \describe{ \item{tre}{is a character string giving the phylogenetic tree in
#' Newick format.} \item{traits}{is a data frame with 19 species and 6 traits}
#' }
#' @note This dataset replaces the former version in ade4.
#' @references Bried, J., Pontier, D. and Jouventin, P. (2002) Mate fidelity in
#' monogamus birds: a re-examination of the Procellariiformes. \emph{Animal
#' Behaviour}, \bold{65}, 235--246.
#' 
#' See a data description at \url{http://pbil.univ-lyon1.fr/R/pdf/pps037.pdf}
#' (in French).
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' if(require(ape) && require(phylobase)){
#' 
#' ## load data, make tree and phylo4d object
#' data(procella)
#' tre <- read.tree(text=procella$tre)
#' x <- phylo4d(tre, procella$traits)
#' par(mar=rep(.1,4))
#' table.phylo4d(x,cex.lab=.7)
#' }
#' }
#' 
NULL





#' Phylogeny and quantitative traits of flowers
#' 
#' This data set describes the phylogeny of 11 flowers as reported by Morales
#' (2000). It also gives morphologic and demographic traits corresponding to
#' these 11 species.
#' 
#' Variables of \code{tithonia$tab} are the following ones : \cr morho1: is a
#' numeric vector that describes the seed size (mm)\cr morho2: is a numeric
#' vector that describes the flower size (mm)\cr morho3: is a numeric vector
#' that describes the female leaf size (cm)\cr morho4: is a numeric vector that
#' describes the head size (mm)\cr morho5: is a integer vector that describes
#' the number of flowers per head \cr morho6: is a integer vector that
#' describes the number of seeds per head \cr demo7: is a numeric vector that
#' describes the seedling height (cm)\cr demo8: is a numeric vector that
#' describes the growth rate (cm/day)\cr demo9: is a numeric vector that
#' describes the germination time\cr demo10: is a numeric vector that describes
#' the establishment (per cent)\cr demo11: is a numeric vector that describes
#' the viability (per cent)\cr demo12: is a numeric vector that describes the
#' germination (per cent)\cr demo13: is a integer vector that describes the
#' resource allocation\cr demo14: is a numeric vector that describes the adult
#' height (m)\cr
#' 
#' @name tithonia
#' @docType data
#' @format \code{tithonia} is a list containing the 2 following objects :
#' \describe{ \item{tre}{is a character string giving the phylogenetic tree in
#' Newick format.} \item{tab}{is a data frame with 11 species and 14 traits (6
#' morphologic traits and 8 demographic).} }
#' @note This dataset replaces the former version in ade4.
#' @source Data were obtained from Morales, E. (2000) Estimating phylogenetic
#' inertia in Tithonia (Asteraceae) : a comparative approach. \emph{Evolution},
#' \bold{54}, 2, 475--484.
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' if(require(ape) && require(phylobase)){
#' 
#' data(tithonia)
#' tre <- read.tree(text=tithonia$tre)
#' traits <- log(tithonia$tab + 1)
#' rownames(traits) <- gsub("_", ".", rownames(traits))
#' 
#' ## build a phylo4d object
#' x <- phylo4d(tre, traits)
#' par(mar=rep(.1,4))
#' table.phylo4d(x)
#' 
#' }
#' }
#' 
NULL





#' Phylogeny and quantitative traits of ungulates.
#' 
#' This data set describes the phylogeny of 18 ungulates as reported by
#' Pelabon et al. (1995). It also gives 4 traits corresponding to these 18
#' species.
#' 
#' Variables of \code{ungulates$tab} are the following ones : \cr
#' 
#' - afbw: is a numeric vector that describes the adult female body weight (g)
#' \cr - mnw: is a numeric vector that describes the male neonatal weight (g)
#' \cr - fnw: is a numeric vector that describes the female neonatal weight (g)
#' \cr - ls: is a numeric vector that describes the litter size \cr
#' 
#' @name ungulates
#' @docType data
#' @format \code{fission} is a list containing the 2 following objects :
#' \describe{ \item{tre}{is a character string giving the phylogenetic tree in
#' Newick format.} \item{tab}{is a data frame with 18 species and 4 traits} }
#' @note This dataset replaces the former version in ade4.
#' @source Data were obtained from Pelabon, C., Gaillard, J.M., Loison, A. and
#' Portier, A. (1995) Is sex-biased maternal care limited by total maternal
#' expenditure in polygynous ungulates?  \emph{Behavioral Ecology and
#' Sociobiology}, \bold{37}, 311--319.
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' if(require(ape) && require(phylobase)){
#' ## load data
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
#' x <- phylo4d(tre, cbind.data.frame(afbw, neonatw)) # traits on the phylogeny
#' 
#' ## test phylogenetic inertia in residuals
#' orthogram(residuals(lm1), x) 
#' }
#' }
#' 
NULL



