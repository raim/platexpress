
#' Plot a heatmap of \code{platexpress} data.
#' 
#' @param data \code{platexpress} data object
#' @param yid ID of the \code{platexpress} data to use
#' @param wells subset or order of wells
#' @param well.norm normalize each well
#' @param q.cut cut values higher then this quantile for color scale
#' @param ylab y-axis label
#' @param xid x-axis data ID in \code{platexpress} object
#' @param xlab x-axis label
#' @param wcol well colors (named vector) for axis labels,
#' eg. from plate layout map
#' @param col color palette to use in heatmap
#' @param ... arguments passed to \code{\link{image}}
#' @export
image_plate <- function(data, yid, wells, well.norm=FALSE, q.cut=.99,
                        ylab="", xid, xlab, wcol,
                        col=gray.colors(100,start=1,end=0), ...) {

    if ( missing(yid) )
        stop("y-axis data ID required, see data$dataIDs for options")
    dat <- data[[yid]]$data[,wells]
    
    if ( missing(wells) ) wells <- colnames(dat)
    dat <- dat[,wells]

    ## transform data such that "wells" are plotted from
    ## top to bottom and time (xid) from left to right
    dat <- dat[,ncol(dat):1]

    ## smooth data
    ##dat <- t(apply(dat, 1, ma, 10))
    
    ## normalize each well
    if ( well.norm )
        dat <- t(t(dat)/apply(dat,2,max,na.rm=T))
    ## cut upper quantiles to rm outliers
    dat[dat>quantile(dat,q.cut,na.rm=TRUE)] <- quantile(dat,q.cut,na.rm=TRUE)

    if ( missing(xid) ) xid <- data$xids[1]
    if ( missing(xlab) ) xlab <- xid
    x <- data[[xid]]

    ## TODO: transform as in image_matrix, or use segmenTools?
    image(z=dat, x=x, y=1:ncol(dat), col=col,
          ylab=ylab, xlab=xlab, axes=FALSE, ...)
    ## colored well axis labels
    if ( missing(wcol) )
        axis(2, at=ncol(dat):1, labels=wells, las=2)
    else 
        for ( j in ncol(dat):1 ) {
            well <- wells[ncol(dat)-j+1]
            axis(2, at=j, labels=well, las=2,
                 col.axis=wcol[well], col=wcol[well])
        }
    axis(1)
}


#' Draw arrows of growth phases.
#'
#' draw arrows of linear segments calculated by
#' \code{\link{dpseg_plate}} or \code{\link{segmented_plate}}
#' with local slope (growth rate) as scaling factor, optionally in
#' well colors.
#' @param segments a list of segmentations for each well as returned
#' by \code{\link{dpseg_plate}} or \code{\link{segmented_plate}}
#' @param wells subset and plot order of wells
#' @param col default color for arrows
#' @param wcol well colors (named vector) for arrows and axis labels,
#' eg. from plate layout map
#' @param xlim correct left and right borders required for
#' results from \code{segmented_plate}
#' @param maxmu maximal slope (growth rate) for line width scaling
#' @param minmu minimal slope (growth rate) for line width scaling
#' @param lwd.max maximal line width
#' @param head.length arrow head length, argument \code{length} to
#' \code{\link{arrows}}
#' @param add add to existing plot
#' @param axis2 add wells as y-axis tick labels
#' @param ... arguments passed to \code{\link{arrows}}
#' @export
arrows_plate <- function(segments, wells, col="#000000", wcol, xlim,
                         maxmu, minmu, lwd.max=10, head.length=0.02, 
                         add=TRUE, axis2=!add, ...) {

    ## take list of tables of growth phases
    ## as returned by dpseg and adapt segmented output

    ## TODO: adapt segmented output
    
    ## breakpoints for growth phases plots
    if ( class(segments[[1]])[1]=="segmented" ) {
        ## estimate xlim from all breakpoints
        ## TODO: betterer!
        if ( missing(xlim) )
            xlim <- range(unlist(lapply(segments, function(o) o$model$x)))
        ## get breakpoints
        segs <- lapply(segments, function(o) c(xlim[1],
                                               confint(o)$x[,1],
                                               xlim[2]))
        #xlim=c(range(unlist(segs))
        ## get slopes
        slps <- lapply(segments, function(o) segmented::slope(o)$x[,1])

    } else if ( class(segments[[1]])[1]=="dpseg" ) {
        ## get breakpoints
        ## TODO: smarter recognition of length-segments without slopes
        segs <-  lapply(segments, function(x) unique(c(x$segments[1,"x1"],
                                                       x$segments[,"x2"])))
        ## get slopes
        slps <- lapply(segments, function(x)
            unlist(na.omit(x$segments[1:nrow(x$segments),"slope"])))
    }

    if ( missing(maxmu) )
        maxmu <- max(unlist(slps),na.rm=TRUE)
    if ( missing(minmu) )
        minmu <- min(unlist(slps),na.rm=TRUE)

    if ( missing(xlim) ) xlim <- c(range(unlist(segs)))
    if ( !add )
        plot(0,col=NA,ylim=c(0,1+length(wells)),xlim=xlim,axes=FALSE)

    for ( i in 1:length(wells) ) {
        sg <- segs[[wells[i]]]
        sl <- slps[[wells[i]]]
        ## dashed lines for negative mu
        lt <- as.numeric(sl<0) +1
        ## colors and lwd scaled with mu, but cutoff
        ## TODO: better control of scales and or alpha?
        sl[sl>maxmu] <- maxmu
        sl[sl<minmu] <- minmu
        ## line width
        lw <- lwd.max*abs(sl)
        ## color
        cl <- rep(col, length(sg)-1)
        if ( !missing(wcol) )
            cl <- rep(wcol[wells[i]], length(sg)-1)
        ##if ( alpha.scale )
        ##alphas <- (sl-minmu.col)/(maxmu-minmu.col) # scale between max and min mu
        ##cl <- add_alphas(cl,alpha=alphas)
        ## well coordinate
        arrows(y0=rep(length(wells)-i+1,length(sg)-1), #rep(i,length(sg)-1),
               x0=sg[2:length(sg)-1],
               x1=sg[2:length(sg)], length=head.length,
               col=cl, lty=lt, lwd=lw, ...)
       # points(sg[3:length(sg)-1], rep(length(wells)-i+1,length(sg)-2),
                                        #rep(i,length(sg)-2),
       #        col=wcol[wells[i]],pch=3)
    }
    if ( axis2 ) 
        if ( missing(wcol) )
            axis(2, at=length(wells):1, labels=wells, las=2)
        else 
            for ( j in length(wells):1 ) {
                well <- wells[length(wells)-j+1]
                axis(2, at=j, labels=well, las=2,
                     col.axis=wcol[well], col=wcol[well])
            }
}

#' Call \code{\link[dpseg:dpseg]{dpseg}} for selected wells.
#'
#' The algorithm in \code{\link[dpseg:dpseg]{dpseg}} splits curves
#' into linear segments.  If the passed data is a biomass measure
#' (eg. OD), and option \code{log=TRUE} the slopes correspond to local
#' growth rate (per time unit).
#' @param data \code{platexpress} data object
#' @param yid ID of the \code{platexpress} data to use
#' @param wells subset and plot order of wells
#' @param xid x-axis data ID in \code{platexpress} data
#' @param log use ln of the data
#' @param verb progress messages
#' @param ... arguments passed to \code{\link[dpseg:dpseg]{dpseg}}
#' @export
dpseg_plate <- function(data, yid="OD", wells, log=TRUE, xid, verb=0, ...) {

    if ( missing(xid) ) xid <- data$xids[1]
    x <- data[[xid]]

    Y <- data[[yid]]$data
    if ( log ) Y <- log(Y)

    if ( missing(wells) )
        wells <- colnames(Y)

    ## call dpseg - TODO: use parallel?
    segments <- sapply(wells, function(well) {
        if ( verb>0 )
            cat(paste("calculating segmentation of", yid,
                      "from well", well, "\n"))
        list(dpseg::dpseg(x=x, y=Y[,well], verb=verb, ...))
    })
    names(segments) <- wells
    invisible(segments)
}

#' 

#' Call \code{\link[segmented]{segmented}} for selected wells.
#'
#' The package \code{\link[segmented]{segmented}} splits curves
#' into linear segments, given a list of pre-defined breakpoints.
#' If the passed data is a biomass measure
#' (eg. OD), and option \code{log=TRUE} the slopes correspond to local
#' growth rate (per time unit).
#' @param data \code{platexpress} data object
#' @param yid ID of the \code{platexpress} data to use
#' @param wells subset and plot order of wells
#' @param log use ln of the data
#' @param xid x-axis data ID in \code{platexpress} data
#' @param man apply a moving average \code{\link{ma}} with \code{n=man}
#' @param psis named list of breakpoints for wells, generated from
#' argument \code{psi} if missing
#' @param psi vector of breakpoints (x-values) to be passed to
#' \code{\link[segmented]{segmented}}, generated as equal spaced breakpoints
#' from argument \code{npsi} if missing
#' @param npsi number of equally spaced breakpoints, used if both
#' \code{psis} and \code{psi} are missing
#' @param plot plot each fit
#' @param verb progress messages
#' @param ... arguments passed to \code{\link[segmented]{segmented}}
#' @export
segmented_plate <- function(data, yid="OD", wells, log=TRUE, xid,
                            man=1, psis, psi, npsi=5, plot=FALSE, verb=0, ...) {

    if ( missing(xid) ) xid <- data$xids[1]
    X <- data[[xid]]

    Y <- data[[yid]]$data
    if ( man>1 )
        Y <- apply(Y, 2, ma, man)
    if ( log ) Y <- log(Y)

    if ( missing(wells) )
        wells <- colnames(Y)

    if ( missing(psis) ) {
        if ( missing(psi) ) 
            psi <-  seq(min(x),max(x),length.out=npsi)[2:(npsi-1)]
        psis <- rep(list(psi), length(wells))
        names(psis) <- wells
    }

    segments <- list()
    for ( well in wells ) {
        if ( verb>0 )
            cat(paste("calculating segments for well", well, "\t"))

        y <- Y[,well]
        x <- X[!is.na(y)]
        y <- y[!is.na(y)]
        x <- x[is.finite(y)]
        y <- y[is.finite(y)]

        if ( plot ) {
            plot(x,y, type="b",cex=.5,main=paste(well))
        }

        ## linear model
        out.lm <- lm(y~x)
        ## TODO: count and warn
        test <- NA
        class(test) <- "try-error"
        maxtry <- 5
        while( class(test)[1]=="try-error" & maxtry>0 ) {
            test <- try(o <- segmented::segmented(out.lm,psi=psis[[well]]))
            if ( maxtry<5 )
                warning(well, ": segmented failed trying again: ",maxtry)
            maxtry <- maxtry -1
        }
        

        if ( verb>0 ) cat("\n")
        if ( maxtry==0 ) {
            warning(well, ": segmented failed")
            next
        }
        segments[[well]] <- o
        if ( plot ) {
            segmented::plot.segmented(o,add=TRUE, col=2, lwd=3, shade=TRUE, rug=FALSE)
            segmented::lines.segmented(o,col=1, lwd=1)
            abline(v=confint(o)$x[,1])
            Sys.sleep(.4)
          }
    }
    segments
}

#' Add results from  \code{\link[dpseg:dpseg]{dpseg}}
#' to \code{platexpress} object.
#'
#' Calls the \code{\link[dpseg:predict.dpseg]{predict.dpseg}} method
#' to reconstruct the piecewise linear growth curve, and adds it to the
#' \code{platexpress} data object. Option \code{add.slopes}
#' allows to add a time-course of growth rates (slopes
#' of the linear segments) instead.
#' @param data \code{platexpress} data object
#' @param segments \code{\link[dpseg:dpseg]{dpseg}} object
#' @param ID data ID for the new object
#' @param add.slopes add slopes instead of reconstructed data
#' @param ... arguments passed to \code{\link{addModel}}
#'@export
addModel_dpseg <- function(data, segments, ID="y", add.slopes=FALSE, ...) {

    if ( add.slopes ) {
        ## get breakpoints and slopes
        segs <- lapply(segments, function(x) x$segments[,c("x1","x2")])
        slps <- lapply(segments, function(x) x$segments[1:nrow(x$segments),
                                                        "slope"])

        x <- data[[data$xids[1]]]

        ## add mus as xy data
        mus <- sapply(1:length(segs), function(i) {
            sg <- segs[[i]]
            sl <- slps[[i]]
            yo <- rep(sl, each=2)
            xo <- c(sapply(1:nrow(sg), function(j)
                seq(sg[j,1], sg[j,2], length.out=2)))
            
            approx(xo,yo, xout=x)$y
        })
        colnames(mus) <- names(segments)
        ## add to platexpress data to plot!
        invisible(addData(data, ID=paste(ID,"mu",sep="_"),
                          dat=mus, ...))
    } else {
        ## add modelled data
        OD_segs <- lapply(segments, function(o) predict(o)$y)
        OD_segs <- do.call(cbind, OD_segs)
        ## add to platexpress data to plot!
        invisible(addData(data, ID=paste(ID,"dpseg",sep="_"),
                          dat=exp(OD_segs), ...))
    }
}
