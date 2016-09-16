#' Phylogenetic principal component analysis
#' 
#' These functions are designed to perform a phylogenetic principal component
#' analysis (pPCA, Jombart et al. 2010) and to display the results.
#' 
#' \code{ppca} performs the phylogenetic component analysis. Other functions
#' are:\cr
#' 
#' - \code{print.ppca}: prints the ppca content\cr
#' 
#' - \code{summary.ppca}: provides useful information about a ppca object,
#' including the decomposition of eigenvalues of all axes\cr
#' 
#' - \code{scatter.ppca}: plot principal components using
#' \code{\link{table.phylo4d}}\cr
#' 
#' - \code{screeplot.ppca}: graphical display of the decomposition of pPCA
#' eigenvalues\cr
#' 
#' - \code{plot.ppca}: several graphics describing a ppca object\cr
#' 
#' The phylogenetic Principal Component Analysis (pPCA, Jombart et al., 2010)
#' is derived from the spatial Principal Component Analysis (spca, Jombart et
#' al. 2008), implemented in the adegenet package (see
#' \code{\link[adegenet]{spca}}).\cr
#' 
#' pPCA is designed to investigate phylogenetic patterns a set of quantitative
#' traits. The analysis returns principal components maximizing the product of
#' variance of the scores and their phylogenetic autocorrelation (Moran's I),
#' therefore reflecting life histories that are phylogenetically structured.
#' Large positive and large negative eigenvalues correspond to global and local
#' structures.\cr
#' 
#' @aliases ppca print.ppca summary.ppca scatter.ppca screeplot.ppca plot.ppca
#' @param x a \linkS4class{phylo4d} object (for \code{ppca}) or a ppca object
#' (for other methods).
#' @param prox a marix of phylogenetic proximities as returned by
#' \code{\link{proxTips}}. If not provided, this matrix will be constructed
#' using the arguments \code{method} and \code{a}.
#' @param method a character string (full or abbreviated without ambiguity)
#' specifying the method used to compute proximities; possible values are:\cr -
#' \code{patristic}: (inversed sum of) branch lengths \cr - \code{nNodes}:
#' (inversed) number of nodes on the path between the nodes \cr -
#' \code{oriAbouheif}: original Abouheif's proximity, with diagonal (see
#' details in \code{\link{proxTips}}) \cr - \code{Abouheif}: Abouheif's
#' proximity (see details in \code{\link{proxTips}}) \cr - \code{sumDD}:
#' (inversed) sum of direct descendants of all nodes on the path (see details
#' in \code{\link{proxTips}}).
#' @param f a function to change a distance into a proximity.
#' @param center a logical indicating whether traits should be centred to mean
#' zero (TRUE, default) or not (FALSE).
#' @param scale a logical indicating whether traits should be scaled to unit
#' variance (TRUE, default) or not (FALSE).
#' @param scannf a logical stating whether eigenvalues should be chosen
#' interactively (TRUE, default) or not (FALSE).
#' @param nfposi an integer giving the number of positive eigenvalues retained
#' ('global structures').
#' @param nfnega an integer giving the number of negative eigenvalues retained
#' ('local structures').
#' @param \dots further arguments passed to other methods. Can be used to
#' provide arguments to \code{\link{table.phylo4d}} in \code{plot} method.
#' @param object a \code{ppca} object.
#' @param printres a logical stating whether results should be printed on the
#' screen (TRUE, default) or not (FALSE).
#' @param axes the index of the principal components to be represented.
#' @param useLag a logical stating whether the lagged components (\code{x\$ls})
#' should be used instead of the components (\code{x\$li}).
#' @param main a title for the screeplot; if NULL, a default one is used.
#' @return The class \code{ppca} are given to lists with the following
#' components:\cr \item{eig}{a numeric vector of eigenvalues.} \item{nfposi}{an
#' integer giving the number of global structures retained.} \item{nfnega}{an
#' integer giving the number of local structures retained.} \item{c1}{a
#' data.frame of loadings of traits for each axis.} \item{li}{a data.frame of
#' coordinates of taxa onto the ppca axes (i.e., principal components).}
#' \item{ls}{a data.frame of lagged prinpal components; useful to represent of
#' global scores.} \item{as}{a data.frame giving the coordinates of the axes of
#' an 'ordinary' PCA onto the ppca axes.} \item{call}{the matched call.}
#' \item{tre}{a phylogenetic tre with class \linkS4class{phylo4}.}
#' \item{prox}{a matrix of phylogenetic proximities.}
#' 
#' Other functions have different outputs:\cr
#' 
#' - \code{scatter.ppca} returns the matched call.\cr
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso The implementation of \code{\link[adegenet]{spca}} in the adegenet
#' package (\code{\link[adegenet]{adegenet}}) \cr
#' @references Jombart, T.; Pavoine, S.; Dufour, A. & Pontier, D. (2010, in
#' press) Exploring phylogeny as a source of ecological variation: a
#' methodological approach. doi:10.1016/j.jtbi.2010.03.038
#' 
#' Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D. (2008) Revealing
#' cryptic phylogenetic patterns in genetic variability by a new multivariate
#' method. \emph{Heredity}, \bold{101}, 92--103.
#' @keywords multivariate
#' @examples
#' 
#' data(lizards)
#' 
#' if(require(ape) && require(phylobase)){
#' 
#' #### ORIGINAL EXAMPLE FROM JOMBART ET AL 2010 ####
#' 
#' 
#' ## BUILD A TREE AND A PHYLO4D OBJECT
#' liz.tre <- read.tree(tex=lizards$hprA)
#' liz.4d <- phylo4d(liz.tre, lizards$traits)
#' par(mar=rep(.1,4))
#' table.phylo4d(liz.4d,var.lab=c(names(lizards$traits),
#'    "ACP 1\n(\"size effect\")"),show.node=FALSE, cex.lab=1.2)
#' 
#' 
#' ## REMOVE DUPLICATED POPULATIONS
#' liz.4d <- prune(liz.4d, c(7,14))
#' table.phylo4d(liz.4d)
#' 
#' 
#' ## CORRECT LABELS
#' lab <- c("Pa", "Ph", "Ll", "Lmca", "Lmcy", "Phha", "Pha",
#'    "Pb", "Pm", "Ae", "Tt", "Ts", "Lviv", "La", "Ls", "Lvir")
#' tipLabels(liz.4d) <- lab
#' 
#' 
#' ## REMOVE SIZE EFFECT
#' dat <- tdata(liz.4d, type="tip")
#' dat <- log(dat)
#' newdat <- data.frame(lapply(dat, function(v) residuals(lm(v~dat$mean.L))))
#' rownames(newdat) <- rownames(dat)
#' tdata(liz.4d, type="tip") <- newdat[,-1] # replace data in the phylo4d object
#' 
#' 
#' ## pPCA
#' liz.ppca <- ppca(liz.4d,scale=FALSE,scannf=FALSE,nfposi=1,nfnega=1, method="Abouheif")
#' liz.ppca
#' tempcol <- rep("grey",7)
#' tempcol[c(1,7)] <- "black"
#' barplot(liz.ppca$eig,main='pPCA eigenvalues',cex.main=1.8,col=tempcol)
#' 
#' par(mar=rep(.1,4))
#' plot(liz.ppca,ratio.tree=.7)
#' 
#' 
#' ## CONTRIBUTIONS TO PC (LOADINGS) (viewed as dotcharts)
#' dotchart(liz.ppca$c1[,1],lab=rownames(liz.ppca$c1),main="Global principal
#' component 1")
#' abline(v=0,lty=2)
#' 
#' dotchart(liz.ppca$c1[,2],lab=rownames(liz.ppca$c1),main="Local principal
#' component 1")
#' abline(v=0,lty=2)
#' 
#' 
#' ## REPRODUCE FIGURES FROM THE PAPER
#' obj.ppca <- liz.4d
#' tdata(obj.ppca, type="tip") <- liz.ppca$li
#' myLab <- paste(" ",rownames(liz.ppca$li), sep="")
#' 
#' ## FIGURE 1
#' par(mar=c(.1,2.4,2.1,1))
#' table.phylo4d(obj.ppca, ratio=.7, var.lab=c("1st global PC", "1st local
#'    PC"), tip.label=myLab,box=FALSE,cex.lab=1.4, cex.sym=1.2, show.node.label=TRUE)
#' add.scatter.eig(liz.ppca$eig,1,1,1,csub=1.2, posi="topleft", ratio=.23)
#' 
#' 
#' ## FIGURE 2
#' s.arrow(liz.ppca$c1,xlim=c(-1,1),clab=1.3,cgrid=1.3)
#' 
#' 
#' 
#' #### ANOTHER EXAMPLE - INCLUDING NA REPLACEMENT ####
#' ## LOAD THE DATA
#' data(maples)
#' tre <- read.tree(text=maples$tre)
#' x <- phylo4d(tre, maples$tab)
#' omar <- par("mar")
#' par(mar=rep(.1,4))
#' table.phylo4d(x, cex.lab=.5, cex.sym=.6, ratio=.1) # note NAs in last trait ('x')
#' 
#' ## FUNCTION TO REPLACE NAS
#' f1 <- function(vec){
#' if(any(is.na(vec))){
#' m <- mean(vec, na.rm=TRUE)
#' vec[is.na(vec)] <- m
#' }
#' return(vec)
#' }
#' 
#' 
#' ## PERFORM THE PPCA
#' dat <- apply(maples$tab,2,f1) # replace NAs
#' x.noNA <- phylo4d(tre, as.data.frame(dat))
#' map.ppca <- ppca(x.noNA, scannf=FALSE, method="Abouheif")
#' map.ppca
#' 
#' 
#' ## SOME GRAPHICS
#' screeplot(map.ppca)
#' scatter(map.ppca, useLag=TRUE)
#' plot(map.ppca, useLag=TRUE)
#' 
#' 
#' ## MOST STRUCTURED TRAITS
#' a <- map.ppca$c1[,1] # loadings on PC 1
#' names(a) <- row.names(map.ppca$c1)
#' highContrib <- a[a< quantile(a,0.1) | a>quantile(a,0.9)]
#' datSel <- cbind.data.frame(dat[, names(highContrib)], map.ppca$li)
#' temp <- phylo4d(tre, datSel)
#' table.phylo4d(temp) # plot of most structured traits
#' 
#' 
#' ## PHYLOGENETIC AUTOCORRELATION TESTS FOR THESE TRAITS
#' prox <- proxTips(tre, method="Abouheif")
#' abouheif.moran(dat[, names(highContrib)], prox)
#' 
#' }
#' 
#' @import phylobase methods
#' @importFrom stats screeplot
#' @export ppca
ppca <- function(x, prox=NULL, method=c("patristic","nNodes","oriAbouheif","Abouheif","sumDD"),
                 f=function(x) {1/x},
                 center=TRUE, scale=TRUE, scannf=TRUE, nfposi=1, nfnega=0){

    ## handle arguments
    ## if(!require(ade4)) stop("The package ade4 is not installed.")
    if (is.character(chk <- checkPhylo4(x))) stop("bad phylo4d object: ",chk)
    ##if (is.character(chk <- checkData(x))) stop("bad phylo4d object: ",chk) : no longer needed

    tre <- as(x, "phylo4")
    method <- match.arg(method)
    NEARZERO <- 1e-10

    ## proximity matrix
    if(is.null(prox)){ # have to compute prox
        W <- proxTips(x, tips="all", method=method, f=f, normalize="row", symmetric=TRUE)
    } else { # prox is provided
        W <- as.matrix(prox)
        if(!is.matrix(W)) stop("W is not a matrix")
        if(ncol(W) != nrow(W)) stop("W is not a square matrix")
        diag(W) <- 0
        W <- 0.5 * (t(W) + W) # re-symmetrization
    }

    N <- nTips(x)

    ## data matrix X
    X <- tdata(x, type="tip")
    X.colnames <- names(X)
    X.rownames <- row.names(X)
    temp <- sapply(X, is.numeric)
    if(!all(temp)) {
        warning(paste("non-numeric data are removed:", X.colnames[!temp]))
        X <- X[,temp]
        X.colnames <- X.colnames[!temp]
        X.rownames <- X.rownames[!temp]
    }

    ## replace NAs
    f1 <- function(vec){
        m <- mean(vec,na.rm=TRUE)
        vec[is.na(vec)] <- m
        return(vec)
    }

    if(any(is.na(X))) {
        warning("Replacing missing values (NA) by mean values")
        X <- as.data.frame(apply(X, 2, f1))
    }

    X <- scalewt(X, center=center, scale=scale) # centring/scaling of traits


    ## main computation ##

    ## make a skeleton of dudi
    res <- dudi.pca(X, center=center, scale=scale, scannf=FALSE,nf=2)
    Upca <- as.matrix(res$c1)

    ## computations of the ppca
    X <- as.matrix(X)
    decomp <- eigen( ((t(X) %*% W %*% X)/N), symmetric=TRUE)
    U <- decomp$vectors # U: principal axes
    lambda <- decomp$values

    ## remove null eigenvalues and corresponding vectors
    toKeep <- (abs(lambda) > NEARZERO)
    lambda <- lambda[toKeep]
    U <- U[, toKeep]
    p <- ncol(U)

    if(scannf){ # interactive part
        barplot(lambda)
        cat("Select the number of global axes: ")
        nfposi <- as.integer(readLines(n = 1))
        cat("Select the number of local axes: ")
        nfnega <- as.integer(readLines(n = 1))
    }

    nfposi <- max(nfposi, 1)
    nfnega <- max(nfnega, 0)
    posi.idx <- 1:nfposi
    if(nfnega<1) {
        nega.idx <- NULL
    } else {
        nega.idx <- (p-nfnega+1):p
    }

    axes.idx <- unique(c(posi.idx, nega.idx)) # index of kept axes
    U <- U[, axes.idx, drop=FALSE]

    S <- X %*% U # S: scores (=princ. components)
    LS <- W %*% S # LS: lagged scores
    A <- t(Upca) %*% U # A: pca princ. axes onto ppca princ. axes.

    ## build the output
    axes.lab <- paste("PA",axes.idx, sep="")
    scores.lab <- paste("PC",axes.idx, sep="")

    res$cent <- res$norm <- res$co <- NULL # cleaning

    res$eig <- lambda # eigenvalues
    res$nf <- NULL
    res$nfposi <- nfposi
    res$nfnega <- nfnega
    res$kept.axes <- axes.idx

    res$c1 <- as.data.frame(U) # principal axes
    names(res$c1) <- axes.lab
    row.names(res$c1) <- X.colnames

    res$li <-  as.data.frame(S) # scores (princ. components)
    names(res$li) <- scores.lab
    row.names(res$li) <- X.rownames

    res$ls <-  as.data.frame(LS) # lagged scores
    names(res$ls) <- scores.lab
    row.names(res$ls) <- X.rownames

    res$as <- as.data.frame(A) # PCA axes onto pPCA axes
    names(res$as) <- axes.lab
    row.names(res$as) <- paste("PCA axis", 1:nrow(A))

    res$tre <- as(tre,"phylo4") # tree

    res$prox <- W # proximity matrix

    res$call <- match.call() # call

    class(res) <- "ppca"

    return(res)
} # end ppca





#####################
# Function scatter.ppca
#####################
#' @rdname ppca
#' @export
scatter.ppca <- function(x, axes=1:ncol(x$li), useLag=FALSE, ...){
    if(useLag){
        df <- as.data.frame(x$ls)
    } else{
        df <- as.data.frame(x$li)
    }

    if(any(axes < 1 | axes > ncol(x$li)) ) stop("Wrong axes specified.")
    df <- df[, axes, drop=FALSE]

    obj <- phylo4d(x$tre,df)
    args <- list(...)
    if(is.null(args$ratio.tree)){
        args$ratio.tree <- 0.5
    }
    args <- c(obj,args)
    do.call(table.phylo4d, args)

    return(invisible(match.call()))
} # end scatter.ppca





######################
# Function print.ppca
######################
#' @rdname ppca
#' @method print ppca
#' @export
print.ppca <- function(x, ...){
  cat("\t#############################################\n")
  cat("\t# phylogenetic Principal Component Analysis #\n")
  cat("\t#############################################\n")
  cat("class: ")
  cat(class(x))
  cat("\n$call: ")
  print(x$call)
  cat("\n$nfposi:", x$nfposi, "axes-components saved")
  cat("\n$nfnega:", x$nfnega, "axes-components saved")
  cat("\n$kept.axes: index of kept axes")

  cat("\nPositive eigenvalues: ")
  l0 <- sum(x$eig >= 0)
  cat(signif(x$eig, 4)[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else cat("\n")
  cat("Negative eigenvalues: ")
  l0 <- sum(x$eig <= 0)
  cat(sort(signif(x$eig, 4))[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else cat("\n")
  cat('\n')
  sumry <- array("", c(1, 4), list(1, c("vector", "length",
                                        "mode", "content")))
  sumry[1, ] <- c('$eig', length(x$eig), mode(x$eig), 'eigenvalues')
  class(sumry) <- "table"
  print(sumry)
  cat("\n")
  sumry <- array("", c(4, 4), list(1:4, c("data.frame", "nrow", "ncol", "content")))
  sumry[1, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "principal axes: scaled vectors of traits loadings")
  sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "principal components: coordinates of taxa ('scores')")
  sumry[3, ] <- c("$ls", nrow(x$ls), ncol(x$ls), 'lag vector of principal components')
  sumry[4, ] <- c("$as", nrow(x$as), ncol(x$as), 'pca axes onto ppca axes')

  class(sumry) <- "table"
  print(sumry)

  cat("\n$tre: a phylogeny (class phylo4)")
  cat("\n$prox: a matrix of phylogenetic proximities")

  cat("\n\nother elements: ")
  if (length(names(x)) > 16)
    cat(names(x)[17:(length(names(x)))], "\n")
  else cat("NULL\n")
} #end print.ppca





###############
# summary.ppca
###############
#' @rdname ppca
#' @method summary ppca
#' @export
summary.ppca <- function (object, ..., printres=TRUE) {

    ## some checks
    if (!inherits(object, "ppca"))stop("to be used with 'ppca' object")
    ## if(!require(ade4)) stop("The package ade4 is not installed.")


    norm.w <- function(X, w) {
        f2 <- function(v) sum(v * v * w)/sum(w)
        norm <- apply(X, 2, f2)
        return(norm)
    }

    resfin <- list()

    if(printres) {
        cat("\n### Phylogenetic Principal Component Analysis ###\n")
        cat("\nCall: ")
        print(object$call)
    }

    appel <- as.list(object$call)
    ## compute original pca
    X <- object$tab # transformed data
    W <- object$prox

    nfposi <- object$nfposi
    nfnega <- object$nfnega

    dudi <- dudi.pca(X, center=FALSE, scale=FALSE, scannf=FALSE, nf=nfposi+nfnega)
    ## end of pca

    Istat <-    data.frame(attributes(moran.idx(X[,1], W,TRUE)))
    row.names(Istat) <- ""
    resfin$Istat <- Istat

    if(printres) {
        cat("\n== Moran's I statistics ==\n")
        print(Istat)
    }

    ## pca scores
    nf <- dudi$nf
    eig <- dudi$eig[1:nf]
    cum <- cumsum(dudi$eig)[1:nf]
    ratio <- cum/sum(dudi$eig)
    moran <- apply(as.matrix(dudi$l1),2,moran.idx, W)
    res <- data.frame(var=eig,cum=cum,ratio=ratio, moran=moran)
    row.names(res) <- paste("Axis",1:nf)
    if(printres) {
        cat("\n== PCA scores ==\n")
        print(res)
    }

    resfin$pca <- res


    ## ppca scores
    ## ppca is recomputed, keeping all axes
    eig <- object$eig
    nfposimax <- sum(eig > 0)
    nfnegamax <- sum(eig < 0)

    listArgs <- appel[-1]
    listArgs$nfposi <- nfposimax
    listArgs$nfnega <- nfnegamax
    listArgs$scannf <- FALSE

    ppcaFull <- do.call(ppca, listArgs) # ppca with all axes

    ndim <- dudi$rank
    nf <- nfposi + nfnega
    toKeep <- c(1:nfposi,if (nfnega>0) (ndim-nfnega+1):ndim)
    varspa <- norm.w(ppcaFull$li,dudi$lw)
    moran <- apply(as.matrix(ppcaFull$li), 2, moran.idx, W)
    res <- data.frame(eig=eig,var=varspa,moran=moran)
    row.names(res) <- paste("Axis",1:length(eig))

    if(printres) {
        cat("\n== pPCA eigenvalues decomposition ==\n")
        print(res[toKeep,])
    }

    resfin$ppca <- res

    return(invisible(resfin))
} # end summary.ppca





#################
# screeplot.ppca
#################
#' @rdname ppca
#' @export
screeplot.ppca <- function(x,...,main=NULL){

  opar <- par("las")
  on.exit(par(las=opar))

  sumry <- summary(x,printres=FALSE)

  labels <- lapply(1:length(x$eig),function(i) bquote(lambda[.(i)]))

  par(las=1)

  xmax <- sumry$pca[1,1]*1.1
  I0 <- unlist(sumry$Istat[1])
  Imin <- unlist(sumry$Istat[2])
  Imax <- unlist(sumry$Istat[3])

  plot(x=sumry$ppca[,2],y=sumry$ppca[,3],type='n',xlab='Variance',ylab="Phylogenetic autocorrelation (I)",xlim=c(0,xmax),ylim=c(Imin*1.1,Imax*1.1),yaxt='n',...)
  text(x=sumry$ppca[,2],y=sumry$ppca[,3],do.call(expression,labels))

  ytick <- c(I0,round(seq(Imin,Imax,le=5),1))
  ytlab <- as.character(round(seq(Imin,Imax,le=5),1))
  ytlab <- c(as.character(round(I0,1)),as.character(round(Imin,1)),ytlab[2:4],as.character(round(Imax,1)))
  axis(side=2,at=ytick,labels=ytlab)

  rect(0,Imin,xmax,Imax,lty=2)
  segments(0,I0,xmax,I0,lty=2)
  abline(v=0)

  if(is.null(main)) main <- ("Decomposition of pPCA eigenvalues")
  title(main)

  return(invisible(match.call()))
} # end screeplot.ppca





############
# plot.ppca
############
#' @rdname ppca
#' @export
plot.ppca <- function(x, axes = 1:ncol(x$li), useLag=FALSE, ...){

    ## some checks
    if (!inherits(x, "ppca")) stop("Use only with 'ppca' objects.")
    if(any(axes>ncol(x$li) | axes<0)) stop("wrong axes required.")

    ## par / layout
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mar = rep(.1,4))
    layout(matrix(c(1,2,3,4,4,4,4,4,4), ncol=3))

    ## some variables
    tre <- x$tre
    n <- nrow(x$li)

    ## 1) barplot of eigenvalues
    omar <- par("mar")
    par(mar = c(0.8, 2.8, 0.8, 0.8))
    r <- length(x$eig)
    col <- rep("white", r)

    keptAxes <- c( (1:r)[1:x$nfposi], (r:1)[1:x$nfnega]) # kept axes
    if(x$nfposi==0) keptAxes <- keptAxes[-1]
    if(x$nfnega==0) keptAxes <- keptAxes[-length(keptAxes)]
    col[keptAxes] <- "grey"

    repAxes <- gsub("PC","",colnames(x$li)[axes]) # represented axes
    repAxes <- as.numeric(repAxes)
    col[repAxes] <- "black"

    barplot(x$eig, col=col)
    title("Eigenvalues", line=-1)
    par(mar=rep(.1,4))
    box()


    ## 2) decomposition of eigenvalues
    par(mar=c(4,4,2,1))
    screeplot(x,main="Eigenvalues decomposition")
    par(mar=rep(.1,4))
    box()


    ## 3) loadings
    if(length(axes)==1){ # one axis retained
        par(mar=c(2.5,4,2,1))
        dotchart(x$c1[,1], labels=row.names(x$c1), main="Loadings",
                 cex=par("cex")*.66)
        abline(v=median(x$c1[,1]), lty=2)
        par(mar=rep(.1,4))
        box()

    } else{ # at least two axes retained
        s.arrow(x$c1[,axes], sub="Loadings")
    }


    ## 4) scatter plot
    ratioTree <- .6
    cexLabel <- 1
    cexSymbol <- 1

    temp <- try(scatter(x, axes=axes, ratio.tree=ratioTree,
                        cex.lab=cexLabel, cex.sym=cexSymbol,
                        show.node=FALSE, useLag=useLag), silent=TRUE) # try default plot
    scatterOk <- !inherits(temp,"try-error")

    while(!scatterOk){
        ## clear 4th screen
        par(new=TRUE)
        plot(1, type="n",axes=FALSE)
        rect(-10,-10, 10,10,col="white")
        par(new=TRUE)
        if(ratioTree > .25 & cexSymbol <= .7) {
            ratioTree <- ratioTree - .05
        }
        if(cexLabel > .65 & cexSymbol <= .5) {
            cexLabel <- cexLabel - .05
        }
        cexSymbol <- cexSymbol - .05

        temp <- try(scatter(x, axes=axes, ratio.tree=ratioTree,
                        cex.lab=cexLabel, cex.sym=cexSymbol,
                        show.node=FALSE, useLag=useLag), silent=TRUE) # try default plot
        scatterOk <- !inherits(temp,"try-error")
    }

    return(invisible(match.call()))

} # end plot.phylo


### testing
## obj <- phylo4d(read.tree(text=mjrochet$tre),mjrochet$tab)
## x@edge.length= rep(1,length(x@edge.label))
## M = cophenetic.phylo(as(x,"phylo"))
## M = 1/M
## diag(M) <- 0


## ppca1 <- ppca(obj,scannf=FALSE,nfp=1,nfn=0)

## plot(ppca1)
