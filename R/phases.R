
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
#' @param wcol color to use for well axis labels, eg. from plate layout map
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


## arrows of growth phases
## option add
## draw arrows with local mu as scaling factor, optionally in colors
## or in opposite scale of heatmap
#arrows_plate <- function(segments, wells) {
#
#    ## TODO: take list of tables of growth phases
#    ## as returned by dpseg and adapt segmented output
#
#    for ( i in 1:length(wells) ) {
#        sg <- segs[[wells[i]]]
#        sl <- slps[[wells[i]]]
#        ## dashed lines for negative mu
#        lt <- as.numeric(sl<0) +1
#        ## colors and lwd scaled with mu, but cutoff 
#        sl[sl>maxmu] <- maxmu
#        sl[sl<mina] <- mina
#        alphas <- (sl-minmu)/(maxmu-minmu)
#        lw <- 3*abs(sl) #alphas
#        cl <- rep(cols[wells[i]], length(sg)-1)
#        #cl <- add_alphas(cl,alpha=alphas)
#        arrows(y0=rep(length(wells)-i+1,length(sg)-1), #rep(i,length(sg)-1),
#               x0=sg[2:length(sg)-1],
#               x1=sg[2:length(sg)], length=0.1,
#               col=cl, lty=lt, lwd=lw)
#       # points(sg[3:length(sg)-1], rep(length(wells)-i+1,length(sg)-2),
#                                        #rep(i,length(sg)-2),
#       #        col=cols[wells[i]],pch=3)
#    }
#
#    
#}

#' Call  \code{\link[dpseg:dpseg]{dpseg}} for selected wells.
#'
#' The algorithm in \code{\link[dpseg:dpseg]{dpseg}} splits
#' curves into linear segments.
#' If the passed data is a biomass measure (eg. OD),
#' and option \code[log=TRUE} the slopes correspond to local
#' growth rate (per time unit).
#' @param data \code{platexpress} data object
#' @param yid ID of the \code{platexpress} data to use
#' @param wells subset or order of wells
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
        slps <- lapply(segments, function(x) x$segments[1:nrow(x$segments),"slope"])

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

    ## breakpoints for growth phases plots
    ## TODO: smarter recognition of length-segments without slopes
    ###segs <-  lapply(segments, function(x) unique(c(x$segments[1,"x1"],
    ###                                               x$segments[,"x2"])))
    ###slps <- lapply(segments, function(x) unlist(na.omit(x$segments[1:nrow(x$segments),"slope"])))
    ## -> use this for above heatmap plots

    ##sgdat[["seg_dp"]] <- segments
    #sgdat[["segIDs"]] <- c(sgdat[["segIDs"]], "segments_dp")
    
    #sgdat
}
