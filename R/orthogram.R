#' Orthonormal decomposition of variance
#' 
#' This function performs the orthonormal decomposition of variance of a
#' quantitative variable on an orthonormal basis. It also returns the results
#' of five non parametric tests associated to the variance decomposition.  It
#' thus provides tools (graphical displays and test) for analysing
#' phylogenetic, pattern in one quantitative trait. This implementation replace
#' the (deprecated) version from the \code{ade4} package.\cr
#' 
#' Several orthonormal bases can be used. By default, basis is constructed from
#' a partition of tips according to tree topology (as returned by
#' \code{\link{treePart}}); for this, the argument \code{tre} must be provided.
#' Alternatively, one can provide an orthonormal basis as returned by
#' \code{\link{orthobasis.phylo}}/\code{\link{me.phylo}} (argument
#' \code{orthobas}), or provide a proximity matrix from which an orthobasis
#' based on Moran's eigenvectors will be constructed (argument \code{prox}).
#' 
#' The function computes the variance decomposition of a quantitative vector x
#' on an orthonormal basis B. The variable is normalized given the uniform
#' weight to eliminate problem of scales.  It plots the squared correlations
#' \eqn{R^{2}}{R^2} between x and vectors of B (variance decomposition) and the
#' cumulated squared correlations \eqn{SR^{2}}{SR^2} (cumulative
#' decomposition).  The function also provides five non parametric tests to
#' test the existence of autocorrelation. The tests derive from the five
#' following statistics :
#' 
#' - R2Max=\eqn{\max(R^{2})}{max(R^2)}. It takes high value when a high part of
#' the variability is explained by one score.\cr -
#' SkR2k=\eqn{\sum_{i=1}^{n-1}(iR^{2}_i)}{sum_i^(n-1) i*(R^2)_i}. It compares
#' the part of variance explained by internal nodes to the one explained by end
#' nodes.\cr - Dmax=\eqn{\max_{m=1,...,n-1}(\sum_{j=1}^{m}R^{2}_j -
#' }{max_(m=1,...,n-1)(sum_(j=1)^m(R^2_j) - (m/n-1))}\eqn{
#' \frac{m}{n-1})}{max_(m=1,...,n-1)(sum_(j=1)^m(R^2_j) - (m/n-1))}. It
#' examines the accumulation of variance for a sequence of scores.\cr -
#' SCE=\eqn{\sum_{m=1}^{n-1} (\sum_{j=1}^{m}R^{2}_j -
#' }{sum_(m=1)^(n-1)(sum_(j=1)^m(R^2_j) - (m/n-1))^2}\eqn{
#' \frac{m}{n-1})^{2}}{sum_(m=1)^(n-1)(sum_(j=1)^m(R^2_j) - (m/n-1))^2}. It
#' examines also the accumulation of variance for a sequence of scores.\cr -
#' ratio: depends of the parameter posinega. If posinega > 0, the statistic
#' ratio exists and equals \eqn{\sum_{i=1}^{posinega}R^{2}_i}{sum_i (R^2)_i
#' with i < posinega + 1}. It compares the part of variance explained by
#' internal nodes to the one explained by end nodes when we can define how many
#' vectors correspond to internal nodes.
#' 
#' @param x a numeric vector corresponding to the quantitative variable
#' @param tre a tree of class \code{\link[ape:read.tree]{phylo}},
#' \linkS4class{phylo4} or \linkS4class{phylo4d}.
#' @param orthobas an object of class \code{'orthobasis'}
#' @param prox a matrix of phylogenetic proximities as returned by
#' \code{\link{proxTips}}.
#' @param nrepet an integer giving the number of permutations
#' @param posinega a parameter for the ratio test. If posinega > 0, the
#' function computes the ratio test.
#' @param tol a tolerance threshold for orthonormality condition
#' @param cdot a character size for points on the cumulative decomposition
#' display
#' @param cfont.main a character size for titles
#' @param lwd a character size for dash lines
#' @param nclass a single number giving the number of cells for the histogram
#' @param high.scores a single number giving the number of vectors to return.
#' If > 0, the function returns labels of vectors that explains the larger part
#' of variance.
#' @param alter a character string specifying the alternative hypothesis, must
#' be one of "greater" (default), "less" or "two-sided"
#' @return If (high.scores = 0), returns an object of class \code{'krandtest'}
#' (randomization tests) corresponding to the five non parametric tests. \cr
#' \cr If (high.scores > 0), returns a list containg : \item{w}{: an object of
#' class \code{'krandtest'} (randomization tests)} \item{scores.order}{: a
#' vector which terms give labels of vectors that explain the larger part of
#' variance}
#' @note This function replaces the former version from the ade4 package, which
#' is deprecated. Note that if ade4 is not loaded BEFORE adephylo, then the
#' version from ade4 will erase that of adephylo, which will still be available
#' from adephylo::orthogram. In practice, though, this should never happen,
#' since ade4 is loaded as a dependence by adephylo.
#' @author Original code: Sebastien Ollier and Daniel Chessel.\cr
#' 
#' Current maintainer: Stephane Dray <stephane.dray@@univ-lyon1.fr>
#' @seealso \code{\link{orthobasis.phylo}}
#' @references Ollier, S., Chessel, D. and Couteron, P. (2005) Orthonormal
#' Transform to Decompose the Variance of a Life-History Trait across a
#' Phylogenetic Tree. \emph{Biometrics}, \bold{62}, 471--477.
#' @examples
#' 
#' \dontrun{
#' if(require(ape) && require(phylobase)){
#' 
#' ## a phylogenetic example
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
#' resid <- residuals(lm1)
#' abline(lm1)
#' 
#' ## plot the two traits and the residuals of lm1
#' x <- phylo4d(tre, cbind.data.frame(afbw, neonatw, residuals=resid))
#' table.phylo4d(x) # residuals are surely not independant
#' 
#' ## default orthogram for residuals of lm1
#' orthogram(resid, tre)
#' 
#' ## using another orthonormal basis (derived from Abouheif's proximity)
#' myOrthoBasis <- orthobasis.phylo(tre, method="oriAbouheif") # Abouheif's proximities
#' orthogram(resid, ortho=myOrthoBasis) # significant phylog. signal
#' 
#' ## Abouheif's test
#' W <- proxTips(tre, method="oriAbouheif") # proximity matrix
#' abouheif.moran(resid, W)
#' }
#' }
#' 
#' @import phylobase
#' @import ade4
#' @export orthogram
orthogram <- function (x, tre=NULL, orthobas = NULL, prox = NULL,
                        nrepet = 999, posinega = 0, tol = 1e-07, cdot = 1.5,
                        cfont.main = 1.5, lwd = 2, nclass,
                        high.scores = 0,alter=c("greater", "less", "two-sided")){

    ## some checks and preliminary assignements
    ## if(!require(ade4)) stop("The ade4 package is not installed.")

    nobs <- length(x)
    alter <- match.arg(alter)

    if(is.numeric(x)&is.vector(x)){
        type <- "numeric"
        ##  } else if(is.factor(x)){
        ##     type <- "factor"
        ##   } else if (inherits(x, "dudi")){
        ##     type <- "dudi"
    } else {
        ## stop("x must be a numeric vector, a factor or a dudi object")
        stop("x must be a numeric vector")
    }
    ##  if(type == "dudi") {
    ##     nobs <- nrow(x$tab)
    ##   } else {
    ##     nobs <- length(x)
    ##   }
    ##   if (!is.null(neig)) {
    ##     orthobas <- scores.neig(neig)
    ##   } else if (!is.null(phylog)) {
    ##     if (!inherits(phylog, "phylog")) stop ("'phylog' expected with class 'phylog'")
    ##     orthobas <- phylog$Bscores
    ##   }

    ## if (is.null(orthobas)){
    ##  stop ("'orthobas','neig','phylog' all NULL")
    ## }

    ## retrieve the orthobasis from a proximity matrix
    if(is.null(orthobas)){
        if(is.null(prox)) { # both orthobas and prox are not given -> default orthobasis
            ## check that tre is provided and valid
            if(is.null(tre)) stop("tre, orthobasis or prox must be provided")
            tre <- as(tre, "phylo4")
            if (is.character(checkval <- checkPhylo4(tre))) stop(checkval)
            orthobas <- treePart(tre, result="orthobasis")
        } else { # else orthobasis from the proxi matrix.
            orthobas <- orthobasis.phylo(prox=prox)
        }
    }

    if (!inherits(orthobas, "data.frame")) stop ("'orthobas' is not a data.frame")
    if (nrow(orthobas) != nobs) stop ("non convenient dimensions")
    if (ncol(orthobas) != (nobs-1)) stop (paste("'orthobas' has",ncol(orthobas),"columns, expected:",nobs-1))
    vecpro <- as.matrix(orthobas)
    npro <- ncol(vecpro)

    w <- t(vecpro/nobs)%*%vecpro
    if (any(abs(diag(w)-1)>tol)) {

        stop("'orthobas' is not orthonormal for uniform weighting")
    }
    diag(w) <- 0
    if ( any( abs(as.numeric(w))>tol) )
        stop("'orthobas' is not orthogonal for uniform weighting")
    if (nrepet < 99) nrepet <- 99
    if (posinega !=0) {
        if (posinega >= nobs-1) stop ("Non convenient value in 'posinega'")
        if (posinega <0) stop ("Non convenient value in 'posinega'")
    }
    if(type!="dudi"){
        if (any(is.na(x)))
            stop("missing value in 'x'")
    }
    if(type == "factor"){
        dudi1 <- dudi.acm(data.frame(x), scannf = FALSE, nf = min(nobs, nlevels(x)))
    }
    if(type == "dudi") {
        if (!all.equal(x$lw, rep(1/nobs, nobs)))
            stop("not implemented for non-uniform row weights")
        dudi1 <- redo.dudi(x, newnf = x$rank)
        if(any(colMeans(dudi1$li)>tol))
            stop("not implemented for non-centered analysis")
    }

    if(type == "numeric") {
        z <- x - mean(x)
        et <- sqrt(mean(z * z))
        if ( et <= tol*(max(z)-min(z))) stop ("No variance")
        z <- z/et
        w <- .C("VarianceDecompInOrthoBasis",
                param = as.integer(c(nobs,npro,nrepet,posinega)),
                observed = as.double(z),
                vecpro = as.double(vecpro),
                phylogram = double(npro),
                phylo95 = double(npro),
                sig025 = double(npro),
                sig975 = double(npro),
                R2Max = double(nrepet+1),
                SkR2k = double(nrepet+1),
                Dmax = double(nrepet+1),
                SCE = double(nrepet+1),
                ratio = double(nrepet+1),
                PACKAGE="adephylo"
                )
    } else {
        w <- .C("MVarianceDecompInOrthoBasis",
                param = as.integer(c(nobs,npro,nrepet,posinega)),
                observed = as.double(as.matrix(dudi1$li)),
                nvar = as.integer(ncol(dudi1$li)),
                inertot = as.double(sum(dudi1$eig)),
                vecpro = as.double(vecpro),
                phylogram = double(npro),
                phylo95 = double(npro),
                sig025 = double(npro),
                sig975 = double(npro),
                R2Max = double(nrepet+1),
                SkR2k = double(nrepet+1),
                Dmax = double(nrepet+1),
                SCE = double(nrepet+1),
                ratio = double(nrepet+1),
                PACKAGE="adephylo"
                )
    }
    ##return(w$phylogram)
    ## multiple graphical window (6 graphs)
    ## 1 pgram
    ## 2 cumulated pgram
    ## 3-6 Randomization tests

    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    layout (matrix(c(1,1,2,2,1,1,2,2,3,4,5,6),4,3))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    par(usr = c(0,1,-0.05,1))


    ylim <- max(c(w$phylogram, w$phylo95))
    names(w$phylogram) <- as.character(1:npro)
    phylocum <- cumsum(w$phylogram)
    lwd0=2
    fun <- function (y, last=FALSE) {
        delta <- (mp[2]-mp[1])/3
        sel <- 1:(npro - 1)
        segments(mp[sel]-delta,y[sel],mp[sel]+delta, y[sel],lwd=lwd0)
        if(last) segments(mp[npro]-delta,y[npro],mp[npro]+delta, y[npro],lwd=lwd0)
    }
    sig50 <- (1:npro)/npro
    y0 <- phylocum - sig50
    h.obs <- max(y0)
    x0 <- min(which(y0 == h.obs))
    par(mar = c(3.1, 2.5, 2.1, 2.1))
    if(type == "numeric"){
        z0 <- apply(vecpro, 2, function(x) sum(z * x))
        mp <- barplot(w$phylogram, col = grey(1 - 0.3 * (sign(z0) > 0)), ylim = c(0, ylim * 1.05))
    } else {
        mp <- barplot(w$phylogram, ylim = c(0, ylim * 1.05))
    }
    scores.order <- (1:length(w$phylogram))[order(w$phylogram, decreasing=TRUE)[1:high.scores]]
    fun(w$phylo95,TRUE)
    abline(h = 1/npro)
    if (posinega!=0) {
        verti = (mp[posinega]+mp[posinega+1])/2
        abline (v=verti, col="red",lwd=1.5)
    }
    title(main = "Variance decomposition",font.main=1, cex.main=cfont.main)
    box()
    obs0 <- rep(0, npro)
    names(obs0) <- as.character(1:npro)
    barplot(obs0, ylim = c(-0.05, 1.05))
    abline(h=0,col="white")
    if (posinega!=0) {
        verti = (mp[posinega]+mp[posinega+1])/2
        abline (v=verti, col="red",lwd=1.5)
    }

    title(main = "Cumulative decomposition",font.main=1, cex.main=cfont.main)
    points(mp, phylocum, pch = 21, cex = cdot, type = "b")
    segments(mp[1], 1/npro, mp[npro], 1, lty = 1)
    fun(w$sig975)
    fun(w$sig025)
    arrows(mp[x0], sig50[x0], mp[x0], phylocum[x0], angle = 15, length = 0.15,
           lwd = 2)
    box()
    if (missing(nclass)) {
        nclass <- as.integer (nrepet/25)
        nclass <- min(c(nclass,40))
    }
    plot(as.randtest (w$R2Max[-1],w$R2Max[1],call=match.call()),main = "R2Max",nclass=nclass)
    if (posinega !=0) {
        plot(as.randtest (w$ratio[-1],w$ratio[1],call=match.call()),main = "Ratio",nclass=nclass)
    } else {
        plot(as.randtest (w$SkR2k[-1],w$SkR2k[1],call=match.call()),main = "SkR2k",nclass=nclass)
    }
    plot(as.randtest (w$Dmax[-1],w$Dmax[1],call=match.call()),main = "DMax",nclass=nclass)
    plot(as.randtest (w$SCE[-1],w$SCE[1],call=match.call()),main = "SCE",nclass=nclass)

    w$param <- w$observed <- w$vecpro <- NULL
    w$phylo95 <- w$sig025 <- w$sig975 <- NULL
    if (posinega==0) {
        w <- as.krandtest(obs=c(w$R2Max[1],w$SkR2k[1],w$Dmax[1],w$SCE[1]),sim=cbind(w$R2Max[-1],w$SkR2k[-1],w$Dmax[-1],w$SCE[-1]),names=c("R2Max","SkR2k","Dmax","SCE"),alter=alter,call=match.call())
    } else {
        w <- as.krandtest(obs=c(w$R2Max[1],w$SkR2k[1],w$Dmax[1],w$SCE[1],w$ratio[1]),sim=cbind(w$R2Max[-1],w$SkR2k[-1],w$Dmax[-1],w$SCE[-1],w$ratio[-1]),names=c("R2Max","SkR2k","Dmax","SCE","ratio"),alter=alter,call=match.call())
    }

    if (high.scores != 0)
        w$scores.order <- scores.order
    return(w)
} # end orthogram
