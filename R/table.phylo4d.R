#############
## table.phylo4d
#############


#' Graphical display of phylogeny and traits
#' 
#' This function represents traits onto the tips of a phylogeny. Plotted objects
#' must be valid \linkS4class{phylo4d} objects (implemented by the 
#' \code{phylobase} package). Current version allows plotting of a tree and one 
#' or more quantitative traits (possibly containing missing data, represented by
#' an 'x').\cr
#' 
#' The plot of phylogenies is performed by a call to 
#' \code{\link[ape]{plot.phylo}} from the \code{ape} package. Hence, many of the
#' arguments of \code{\link[ape]{plot.phylo}} can be passed to 
#' \code{table.phylo4d}, through the \dots{} argument, but their names must be 
#' complete.
#' 
#' For large trees, consider using \code{\link{bullseye}}.
#' 
#' The function \code{table.phylo4d} is based on former plot method for 
#' \linkS4class{phylo4d} objects from the \code{phylobase} package.  It replaces
#' the deprecated \code{ade4} functions \code{\link[ade4]{symbols.phylog}} and
#' \code{\link[ade4]{table.phylog}}.
#' 
#' @param x a \linkS4class{phylo4d} object
#' @param treetype the type of tree to be plotted ("phylogram" or "cladogram")
#' @param symbol the type of symbol used to represent data ("circles", 
#'   "squares", or "colors")
#' @param repVar the numerical index of variables to be plotted
#' @param center a logical stating whether variables should be centred (TRUE, 
#'   default) or not (FALSE)
#' @param scale a logical stating whether variables should be scaled (TRUE, 
#'   default) or not (FALSE)
#' @param legend a logical stating whether a legend should be added to the plot 
#'   (TRUE) or not (FALSE, default)
#' @param grid a logical stating whether a grid should be added to the plot 
#'   (TRUE, default) or not (FALSE)
#' @param box a logical stating whether a box should be added around the plot 
#'   (TRUE, default) or not (FALSE)
#' @param show.tip.label a logical stating whether tip labels should be printed 
#'   (TRUE, default) or not (FALSE)
#' @param show.node.label a logical stating whether node labels should be 
#'   printed (TRUE, default) or not (FALSE)
#' @param show.var.label a logical stating whether labels of variables should be
#'   printed (TRUE, default) or not (FALSE)
#' @param ratio.tree the proportion of width of the figure occupied by the tree
#' @param font an integer specifying the type of font for the labels: 1 (plain 
#'   text), 2 (bold), 3 (italic, default), or 4 (bold italic).
#' @param tip.label a character vector giving the tip labels
#' @param var.label a character vector giving the labels of variables
#' @param cex.symbol a numeric giving the factor scaling the symbols
#' @param cex.label a numeric giving the factor scaling all labels
#' @param cex.legend a numeric giving the factor scaling the legend
#' @param pch is \code{symbol} is set to 'colors', a number indicating the type 
#'   of point to be plotted (see ?points)
#' @param col is \code{symbol} is set to 'colors', a vector of colors to be used
#'   to represent the data
#' @param coord.legend an optional list with two components 'x' and 'y' 
#'   indicating the lower-left position of the legend. Can be set to 
#'   \code{locator(1) to position the legend interactively.}
#' @param \dots further arguments to be passed to plot methods from \code{ape}. 
#'   See \code{\link[ape]{plot.phylo}}.
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso The \linkS4class{phylo4d} class for storing 
#'   \code{phylogeny+data}.\cr
#'   
#'   For large trees, consider using \code{\link{bullseye}}.
#'   
#'   \code{\link[ape]{plot.phylo}} from the \code{ape} package.\cr
#'   
#'   An alternative (deprecated) representation is available from 
#'   \code{\link[ade4]{dotchart.phylog}}.
#'   
#' @returns No return value, function produces only a plot. 
#' 
#' @keywords hplot multivariate
#' @examples
#' 
#' if(require(ape) & require(phylobase) & require(ade4)){
#' 
#' ## simulated data
#' tr <- rtree(20)
#' dat <- data.frame(a = rnorm(20), b = scale(1:20), c=runif(20,-2,2) )
#' dat[3:6, 2] <- NA # introduce some NAs
#' obj <- phylo4d(tr, dat) # build a phylo4d object
#' table.phylo4d(obj) # default scatterplot
#' table.phylo4d(obj,cex.leg=.6, use.edge.length=FALSE) # customized
#' table.phylo4d(obj,treetype="clad", show.node=FALSE, cex.leg=.6,
#' use.edge.length=FALSE, edge.color="blue", edge.width=3) # more customized
#' 
#' 
#' ## teleost fishes data
#' data(mjrochet)
#' temp <- read.tree(text=mjrochet$tre) # make a tree
#' mjr <- phylo4d(x=temp,tip.data=mjrochet$tab) # male a phylo4d object
#' table.phylo4d(mjr,cex.lab=.5,show.node=FALSE,symb="square")
#' 
#' 
#' ## lizards data
#' data(lizards)
#' liz.tr <- read.tree(tex=lizards$hprA) # make a tree
#' liz <- phylo4d(liz.tr, lizards$traits) # make a phylo4d object
#' table.phylo4d(liz)
#' 
#' 
#' ## plotting principal components
#' liz.pca1 <- dudi.pca(lizards$traits, scannf=FALSE, nf=2) # PCA of traits
#' myPC <- phylo4d(liz.tr, liz.pca1$li) # store PC in a phylo4d object
#' varlab <- paste("Principal \ncomponent", 1:2) # make labels for PCs
#' table.phylo4d(myPC, ratio=.8, var.lab=varlab) # plot the PCs
#' add.scatter.eig(liz.pca1$eig,2,1,2,posi="topleft", inset=c(0,.15))
#' title("Phylogeny and the principal components")
#' 
#' }
#' 
#' @import phylobase
#' @importFrom ape plot.phylo
#' @importFrom graphics par strwidth segments symbols points text strheight
#' @importFrom grDevices heat.colors
#' @export table.phylo4d
table.phylo4d <- function(x, treetype=c("phylogram","cladogram"), symbol=c("circles", "squares", "colors"),
                          repVar=1:ncol(tdata(x, type="tip")), center=TRUE, scale=TRUE, legend=TRUE, grid=TRUE, box=TRUE,
                          show.tip.label=TRUE, show.node.label=TRUE, show.var.label=TRUE,
                          ratio.tree=1/3, font=3,
                          tip.label=tipLabels(x), var.label=colnames(tdata(x,type="tip")),
                          cex.symbol=1, cex.label=1, cex.legend=1,
                          pch=20, col=heat.colors(100), coord.legend=NULL, ...)
{

    ## preliminary stuff and checks
    if (is.character(chk <- checkPhylo4(x))) stop("bad phylo4d object: ",chk)
                                        # if (is.character(chk <- checkData(x))) stop("bad phylo4d object: ",chk) <- needed?

    ## if(!require(ape)) stop("the ape package is required")
    if(cex.label<0.1) {
        show.tip.label <- FALSE
        show.node.label <- FALSE
        show.var.label <- FALSE
    }

    cex <- par("cex")
    symbol <- match.arg(symbol)
    treetype <- match.arg(treetype)

    SYMBSCALE <- 0.2 # i.e. max size of a plotted symbol is 0.2*cex.symbol inches
    if(symbol=="colors") {
        SYMBSCALE <- 0.05
    }

    ## convert the tree into phylo
    tre <- suppressWarnings(as(x,"phylo"))
    ##tre$node.label <- x@node.label # this should be done by the as(x,"phylo")
    ## plot only tree if no tip data
    if(ncol(tdata(x,type="tip")) == 0) {
        plot(tre, type=treetype, direction="rightwards", show.tip.label=show.tip.label,
                   show.node.label=show.node.label, cex=cex.label,
                   no.margin=FALSE, x.lim=NULL, y.lim=NULL, ...)
        return(invisible())
    }

#### data handling
    ## retrieve data
    dat <- tdata(x, type="tip")
    dat <- dat[, repVar, drop=FALSE]
    clas <- lapply(dat,class)
    isNum <- sapply(clas, function(e) e %in% c("integer","numeric"))
    ## keep only numeric data
    dat <- dat[, isNum, drop=FALSE]
    var.label <- var.label[repVar]
    var.label <- var.label[isNum]
    ## order data like tips
    E <- phylobase::edges(x)
    tips.ord <- E[,2][!E[,2] %in% E[,1]]
    dat <- dat[tips.ord,,FALSE]
    tip.label <- tip.label[tips.ord] # reorder tip labels
    ## centring / scaling
    dat <- as.data.frame(scale(dat,center=center,scale=scale))

    ## compute bottom margin
    ## ! use inches as units still these won't be changed by plot.phylo
    temp <- var.label[which.max(nchar(var.label))] # longest tip label
    lab.height <- strwidth(temp, units="inches", cex=cex.label) # height required by the longest var label
    lab.height <- lab.height / par("pin")[1] # returned as a fraction of the plot region

#### define plot region
    plotreg <- plotreg0 <- par("plt")
    plotreg.width <- plotreg0[2] - plotreg0[1]
    plotreg.height <- plotreg0[4] - plotreg0[3]
    plotreg[2] <- plotreg[1] + (ratio.tree)*plotreg.width # restrain the width for phylo
    plotreg[3] <- plotreg[3] + plotreg.height*ifelse(show.var.label,lab.height+0.05,0.05) ## add internal vertical margins
    plotreg[4] <- plotreg[4] - plotreg.height*0.05 # add internal vertical margins

#### plot the tree
    par(plt = plotreg)
    plotres <- plot(tre, type=treetype, direction="rightwards", show.tip.label=FALSE,
                          show.node.label=show.node.label, cex=cex.label,
                          no.margin=FALSE, x.lim=NULL, y.lim=NULL, ...)

#### plot the data
    par(plt=plotreg0)
    cur.usr.width <- par("usr")[2] - par("usr")[1] # beware: par("usr") does not adapt to the new plot region
    usr.width <- cur.usr.width / ratio.tree
    usr.height <- par("usr")[4] - par("usr")[3]

    ## x.inset is the space between tree/data and data/tip.labels (in usr units)
    x.inset <- SYMBSCALE * cex.symbol * usr.width / par("pin")[1]
    y.inset <- SYMBSCALE * cex.symbol * usr.height / par("pin")[2]
    x.base <- plotres$x.lim[2] + x.inset # start plotting from x.base rightwards
    if(show.tip.label){
        temp <- tipLabels(x)[which.max(nchar(tipLabels(x)))] # longest tip label
        lab.width <- strwidth(temp, units="user", cex=cex.label) # compute the width to keep for tip labels
    } else{
        lab.width <- 0
    }
    xrange.data <- c(x.base , (par("usr")[1]+usr.width) - lab.width - 2*x.inset) # plot data within this range

    ##    if(diff(xrange.data) < (x.inset*ncol(dat))) ("No room left to plot data; please try reducing ratio.tree or cex.label.")
    if(diff(xrange.data) < (x.inset*ncol(dat))) warning("There may not be enough room left to plot data; you may consider reducing ratio.tree or cex.label.")

    ## define x and y coordinates
    x.grid <- seq(xrange.data[1],xrange.data[2], length=ncol(dat))
    if(ncol(dat)==1) {x.grid <- mean(c(xrange.data[1],xrange.data[2]))}
    y.grid <- seq(plotres$y.lim[1],plotres$y.lim[2],length=plotres$Ntip)
    temp <- expand.grid(y.grid, x.grid) # here are coordinates for data
    xy.data <- data.frame(x=temp[,2],y=temp[,1])

    ## merge data and their coordinates
    alldat <- cbind.data.frame(xy.data, unlist(dat))
    ##    fac <- factor(rep(1:ncol(dat), rep(nrow(dat),ncol(dat))))
    ##     alldat <- split(alldat, fac)

    ## need to "reboot" the plot region without changing coordinates
    ## seems that box does the job.
    if(box) {box()} else {box(col="transparent")}
    if(grid){
        ## vertical segments
        segments(x0=x.grid, y0=rep(min(y.grid),plotres$Ntip),
                 x1=x.grid, y1=rep(max(y.grid),plotres$Ntip), col="grey")
        ## horizontal segments
        segments(x0=rep(min(x.grid),plotres$Ntip), y0=y.grid,
                 x1=rep(max(x.grid),plotres$Ntip), y1=y.grid, col="grey")
    }


    ## auxiliary function to translate a variable into colors
    makeColors <- function(x, col){ # x is a numeric vector, col is a vector of colors
        if(length(x)==1) return(col[1])
        nCol <- length(col)
        res <- x - min(x)
        res <- res / max(res)
        res <- res * (nCol-1) + 1
        res <- round(res)
        res[res>nCol] <- nCol
        res[res<1] <- 1
        return(col[res])
    }


    ## auxiliary function to plot a single variable
    ## max size of a symbol is set to SYMBSCALE*cex inches
    plotaux <- function(x,y,var,symbol,cex){
        if(any(var[!is.na(var)]<0)) {
            usebw <- TRUE
        } else {
            usebw <- FALSE
        }

        if(usebw){
            ispos <- var>0
            fg.col <- rep("black",length(var))
            fg.col[ispos] <- "white"
            bg.col <- rep("white",length(var))
            bg.col[ispos] <- "black"

            if(symbol == "squares"){
                symbols(x=x, y=y, squares=abs(var), inches=SYMBSCALE*cex, fg=fg.col, bg=bg.col, add=TRUE)
            } # end squares

            if(symbol == "circles"){
                symbols(x=x, y=y, circles=abs(var), inches=SYMBSCALE*cex, fg=fg.col, bg=bg.col, add=TRUE)
            } # end circles

            if(symbol == "colors"){
                myCol <- makeColors(var, col)
                points(x=x, y=y, pch=pch, cex=cex, col=myCol)
            } # end colors

        } else {

            if(symbol == "squares"){
                symbols(x=x, y=y, squares=var, inches=SYMBSCALE*cex, fg="white", bg="black", add=TRUE)
            } # end squares

            if(symbol == "circles"){
                symbols(x=x, y=y, circles=var, inches=SYMBSCALE*cex, fg="white", bg="black", add=TRUE)
            } # end circles

            if(symbol == "colors"){
                myCol <- makeColors(var, col)
                points(x=x, y=y, pch=pch, cex=cex, col=myCol)
            } # end colors

        } # end else

        if(any(is.na(var))){
            isNA <- is.na(var)
            points(x[isNA],y[isNA],pch=4,cex=cex.symbol)
        }
    } # end plotaux


    ## finally plot the data
    ## carefull : all variables must be plotted in as a single vector, so that
    ## scaling is the same for all variables
    ## lapply(alldat, function(X) plotaux(X[,1],X[,2],X[,3],symbol,cex.symbol))
    plotaux(alldat[,1],alldat[,2],alldat[,3],symbol,cex.symbol)

#### plot labels for variables
    if(show.var.label) text(x=x.grid, y=rep(min(y.grid)-1.5*y.inset, ncol(dat)), lab=var.label,
                            adj=1, srt=90, cex=cex.label)

#### plot tip labels
    if(show.tip.label){
        x.base <- xrange.data[2] + x.inset
        text(x=rep(x.base,plotres$Ntip), y=1:plotres$Ntip, lab=tip.label, font=font, cex=cex.label, pos=4)
    }

#### add a legend for symbols
    if(legend){

        ## Auxiliary function to add the legend
        ## (x,y): coordinates of the lower-left annotation
        ## z: a numeric vector whose values are being legended
        addLegend <- function(x,y,z,cex.legend,cex.label,cex.symbol){
            z <- z*cex.legend
            leg.values <- pretty(z,n=4, min.n=1)
            temp <- length(leg.values)
            ## make sure to get maximum 4 symbols
            if(temp>4) {
                leg.values <- leg.values[c(1,2,temp-1,temp)]
            }

            leg.txt <- as.character(leg.values)

            ## compute the maximum size taken by symbols in usr coordinates
            if(symbol=="colors") {
                sym.w <-strwidth("o",units="user",cex=cex.symbol)
                sym.w <- rep(sym.w, length(leg.values))
                sym.h <- strheight("o",units="user",cex=cex.symbol)
                sym.h <- rep(sym.h, length(leg.values))
            } else {
                usr.w <- (par("usr")[2]-par("usr")[1]) / ratio.tree # because par("usr") is the one of plot.phylo
                usr.h <- par("usr")[4]-par("usr")[3]
                sym.w <- usr.w *
                    ((abs(leg.values)/max(abs(leg.values))) * SYMBSCALE * cex.symbol * cex.legend) / par("pin")[1]
                sym.h <- usr.h * (SYMBSCALE * cex.symbol * cex.legend) / par("pin")[2]
            }

            ## compute the maximum size taken by annotations in usr coordinates
            ann.w <- strwidth(leg.txt,units="user",cex=cex.label*cex.legend)
            ann.h <- strheight(leg.txt,units="user",cex=cex.label*cex.legend)

            ## retain relevant spaces between symbols / annotations
            space.w.sym <- sapply(1:(length(sym.w)-1),function(i)  sum(sym.w[c(i,i+1)]))
            space.w.ann <- sapply(1:(length(ann.w)-1),function(i)  sum(ann.w[c(i,i+1)])) / 2
            temp <- cbind(space.w.sym, space.w.ann)
            space.w <- apply(temp,1,max)
            if(symbol=="colors"){
                space.h <- sym.h + ann.h
            } else {
                space.w <- space.w + 0.01*usr.w
                space.h <- sym.h + ann.h + 0.01*usr.h
            }

            ## define coordinates of annotations and symbols
            ann.coordX <- c(x, x + cumsum(space.w)) + max(sym.w[1],ann.w[1]) + 0.01*usr.w
            ann.coordY <- y
            sym.coordX <- ann.coordX
            sym.coordY <- y + space.h

            ## plot annotations
            text(ann.coordX, ann.coordY, leg.txt, cex=cex.label*cex.legend)

            ## plot symbols
            plotaux(sym.coordX, sym.coordY, leg.values, symbol, cex.symbol*cex.legend)
        } # end addLegend

        if(!is.null(coord.legend)){
            x.leg <- coord.legend$x
            y.leg <- coord.legend$y
        } else {
            usr.w <- (par("usr")[2]-par("usr")[1]) / ratio.tree
            usr.h <- par("usr")[4]-par("usr")[3]

            temp <- lab.height * usr.height / (1 - lab.height) ## need to substract temp from par("usr")[3]
            y.base <- par("usr")[3] - temp - y.inset ## to get closer the actual par("usr")[3] !

            x.leg <- par("usr")[1] + 0.01 * usr.w
            y.leg <- y.base ## remember to use y.base instead of par("usr3")[3], which is wrong
        }

        addLegend(x=x.leg, y=y.leg, z=alldat[,3],
                  cex.legend=cex.legend, cex.label=cex.label, cex.symbol=cex.symbol)
        ## FIXME ##
        ## draw a rectangle around the legend
                                        #rect.size <- c(diff(range(leg.x)) , diff(c(y.base, max(leg.y))) )
                                        #rect(min(leg.x)- rect.size[1]*0.05,
                                        #     min(y.base) - rect.size[2]*0.05,
                                        #     max(leg.x) + rect.size[1]*0.05,
                                        #     max(y.base) + rect.size[2]*0.05)
    } ## end legend

    return(invisible())
} # end table.phylo4d
