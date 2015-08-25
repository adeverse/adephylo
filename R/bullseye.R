


#' Fan-like phylogeny with possible representation of traits on tips
#'
#' This function represents a phylogeny as a fan, using circles to provide a
#' legend for distances and optionally colored symbols to represent traits
#' associated to the tips of the tree. This function uses and is compatible
#' with ape's \code{\link[ape]{plot.phylo}}.
#'
#'
#' @param phy a tree in \code{phylo}, \linkS4class{phylo4} or
#' \linkS4class{phylo4d} format.
#' @param traits an optional data.frame of traits.
#' @param col.tips.by an optional vector used to define colors for tip labels;
#' if unamed, must be ordered in the same order as \code{phy$tip.label}.
#' @param col.pal a function generating colors according to a given palette;
#' several palettes can be provided as a list, in the case of several traits;
#' the first palette is always reserved for the tip colors; this argument is
#' recycled.
#' @param circ.n the number of circles for the distance annotations.
#' @param circ.bg the color of the circles.
#' @param circ.unit the unit of the circles; if NULL, determined automatically
#' from the data.
#' @param legend a logical specifying whether a legend should be plotted; only
#' one legend is displayed, with priority to tip colors first, and then to the
#' first trait.
#' @param leg.posi,leg.title,leg.bg position, title and background for the
#' legend.
#' @param traits.inset inset for positioning the traits; 1 corresponds to the
#' circle crossing the furthest tip, 0 to the center of the plot.
#' @param traits.space a coefficient indicating the spacing between traits.
#' @param traits.pch,traits.cex type and size of the symbols used for the
#' traits; recycled if needed.
#' @param alpha alpha value to be used for the color transparency, between 0
#' (invisible) and 1 (plain).
#' @param axis a logical indicating whether an axis should be displayed.
#' @param \dots further arguments to be passed to plot methods from \code{ape}.
#' See \code{\link[ape]{plot.phylo}}.
#' @author Thibaut Jombart \email{tjombart@@imperial.ac.uk}
#' @seealso \code{\link{table.phylo4d}} for non-radial plots.\cr
#'
#' The \linkS4class{phylo4d} class for storing \code{phylogeny+data}.\cr
#'
#' \code{\link[ape]{plot.phylo}} from the \code{ape} package.\cr
#'
#' \code{\link[ade4]{dotchart.phylog}}.
#' @keywords hplot multivariate
#'
#' @export
#'
#' @examples
#'
#' if(require(ape) && require(phylobase) && require(adegenet)){
#'
#' data(lizards)
#' tre <- read.tree(text=lizards$hprA) # make a tree
#'
#' ## basic plots
#' bullseye(tre)
#' bullseye(tre, lizards$traits)
#'
#' ## customized
#' par(mar=c(6,6,6,6))
#' bullseye(tre, lizards$traits, traits.cex=sqrt(1:7), alpha=.7,
#'          legend=FALSE, circ.unit=10, circ.bg=transp("black",.1),
#'          edge.width=2)
#'
#' }
#'
bullseye <- function(phy, traits=NULL, col.tips.by=NULL, col.pal=spectral,
                     circ.n=6, circ.bg=transp("royalblue",.1), circ.unit=NULL,
                     legend=TRUE, leg.posi="bottomleft", leg.title="", leg.bg="white",
                     traits.inset=1.1, traits.space=0.05, traits.pch=19, traits.cex=1,
                     alpha=1, axis=TRUE, ...){
    ## CHECKS ##
    if(inherits(phy, c("phylo4","phylo4d"))) phy <- as(phy, "phylo")
    if(!is.list(col.pal)) col.pal <- c(col.pal)
    leg.info <- NULL

    ## REORDER DATA BY TIP LABEL ##
    ## make sure traits is a data.frame
    if(!is.null(traits)) traits <- as.data.frame(traits)
    if(!is.null(traits) && !is.null(row.names(traits))){
        if(!all(phy$tip.label %in% row.names(traits))){
            warning("tip labels and names of the traits matrix do not match")
        } else {
            traits <- traits[phy$tip.label,,drop=FALSE]
        }
    }

    ## col.tips.by
    if(!is.null(col.tips.by) && is.data.frame(col.tips.by)){
        old.names <- row.names(col.tips.by)
        col.tips.by <- unlist(col.tips.by)
        names(col.tips.by) <- old.names
    }
    if(!is.null(col.tips.by) && !is.null(names(col.tips.by))){
        col.tips.by <- col.tips.by[phy$tip.label]
    }

    ## recycle col.pal
    pal.length <- 0
    if(!is.null(traits)) pal.length <- pal.length + ncol(traits)
    if(!is.null(col.tips.by)) pal.length <- pal.length + 1
    col.pal <- rep(col.pal, length=pal.length)


    ## PLOT THE PHYLOGENY
    ## window setting
    oxpd <- par("xpd")
    par(xpd=TRUE)
    on.exit(par(oxpd))

    ## handle color info
    if(!is.null(col.tips.by)){
        tip.col.info <- any2col(col.tips.by, col.pal=col.pal[[1]])
        plot(phy, type="fan", tip.col=transp(tip.col.info$col,alpha), ...)
    } else{
        plot(phy, type="fan", ...)
    }


    ## HANDLE THE 'BULLSEYE' ##
    ## annot info
    if(is.null(circ.unit)){
        annot.max <- 0.5*diff(par("usr")[1:2])
        annot.dist <- seq(from=0, to=annot.max, length=circ.n)
    } else {
        annot.dist <- seq(from=0, by=circ.unit, length=circ.n)
        annot.max <- max(annot.dist)
    }

    ## trace the disks
    symbols(rep(0,circ.n), rep(0,circ.n), circles=annot.dist, inches=FALSE,
        bg=circ.bg, fg=NA, add=TRUE)

    ## axis annotation
    if(axis){
        segments(-annot.dist[2],0,-annot.dist[3],0)
        text(-mean(annot.dist[2:3]),-annot.dist[2]/5,
             label=format(annot.dist[2], scientific=TRUE, digits=3),cex=.7)
    }


    ## PLOT TRAITS ##
    if(!is.null(traits)){
        ## recycle pch and cex
        traits.pch <- rep(traits.pch, length=ncol(traits))
        traits.cex <- rep(traits.cex, length=ncol(traits))

        ## get tips coordinates
        tips.x <- get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[1:length(phy$tip.label)]
        tips.y <- get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[1:length(phy$tip.label)]

        ## use furthest tip from the root to define new base coords
        vec.length <- sqrt(tips.x^2 + tips.y^2)

        x.base <- (tips.x/vec.length) * max(vec.length) * traits.inset
        y.base <- (tips.y/vec.length) * max(vec.length) * traits.inset

        ## plot traits
        for(i in 1:ncol(traits)){
            col.info <- any2col(traits[,i], col.pal=col.pal[[i]])
            temp.x <- x.base * (traits.inset + i*traits.space)
            temp.y <- y.base * (traits.inset + i*traits.space)
            points(temp.x, temp.y, pch=traits.pch[i], col=transp(col.info$col,alpha), cex=traits.cex[i])

            ## save info for legend if needed
            if(is.null(col.tips.by) && i==1){
                leg.info <- list(col=transp(col.info$leg.col,alpha), txt=col.info$leg.txt)
            }
        }
    }


    ## ADD LEGEND ##
    ## legend info
    if(!is.null(legend)){
        ## legend for tip colors
        if(!is.null(col.tips.by)){
            leg.col <- transp(tip.col.info$leg.col,alpha)
            leg.txt <- tip.col.info$leg.txt
            leg.info <- list(col=transp(tip.col.info$leg.col,alpha), txt=tip.col.info$leg.txt)
        }

        ## plot legend
        if(!is.null(leg.info) && legend){
            leg.info$posi <- leg.posi
            legend(x=leg.info$posi, legend=leg.info$txt, fill=leg.info$col, title=leg.title, bg=leg.bg)
            return(invisible(leg.info))
        }
    }

    return(invisible())
} # end bullseye
