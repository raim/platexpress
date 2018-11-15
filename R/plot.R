
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
        sl <- abs(sl)
        sl[sl>maxmu] <- maxmu
        sl[sl<minmu] <- minmu
        ## line width
        lw <- lwd.max*sl
        ## color
        cl <- rep(col, length(sg)-1)
        if ( !missing(wcol) )
            cl <- rep(wcol[wells[i]], length(sg)-1)
        ##if ( alpha.scale )
        ##alphas <- (sl-minmu.col)/(maxmu-minmu.col) 
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


## TODO
## wrapper for viewGroups taking list of plates
## and a merged layout file
viewGroupl <- function(data, groups, groups2, ...) {}

#' View Group Averages
#'
#' plot grouped wells as summary plots, incl. means and confidence intervals.
#' This function gives a first overview of the reproducability between
#' replicates.
#' @param data the list of measurement data as provided by
#' \code{\link{readPlateData}}
#' @param groups a list of well grouped wells, as produced by
#' \code{\link{getGroups}}(platemap, by=c("media")); cf. \code{groups2}
#' @param groups2 sub-groups of \code{groups}, group2 must be constructed as
#' \code{groups}, but with one additional grouping, e.g.
#' \code{\link{getGroups}}(platemap, by=c("media","strain")) following the
#' example for parameter see \code{groups}
#' @param xid ID of a data-set in the input data that can be used as x-axis
#' instead of the default Time vector
#' @param yids IDs of the data to be plotted; if missing, all data will
#' be plotted
#' @param dtype type of the data to be plotted, default is the main 'data', but
#' e.g. 'orig' allows to plot the original data w/o processing (e.g.
#' interpolation, blank correction, etc.)
#' @param xscale use a global range for the x-axes; only relevant if xid
#' specifies a subset of the data as x-axis
#' @param xlim plot range of the x-axis
#' @param xlab custom x-axis label
#' @param ylab custom y-axis label
#' @param pcols a named list of RGB colors to be used the plotted data types;
#' the color vector must have names according to the data IDs
#' @param group2.col a named list of RGB colors for group2 (which are shown
#' in one plot)
#' @param yscale if TRUE (default) global y-axis limits will be calculated from
#' all plotted wells; if FALSE each well be locally scaled
#' @param ylims a named list of y-axis ranges pairs for each data ID
#' @param ylim one y-axis limit range that will be used for all plotted data
#' @param log plot logarithmic axis, use equivalent to normal plot 'log', i.e.,
#' log="y" for a log y-axis, log="x" for x-axis and log="yx" for both axes
#' @param show.ci95 show the 95\% condifence intervals of groups
#' @param show.mean show the mean of groups as lines
#' @param emphasize.mean show the means in black instead of group colors,
#' @param lwd.mean line width for the mean lines, see \code{(par("lwd")}
#' @param lty.mean line style for the mean lines, see \code{(par("lty")}
#' this can help to emphasize the means in noisy data (broad overlapping
#' confidence intervals
#' @param g1.legend show the main legend, giving plot colors
#  of the plotted data types (\code{yids})
#' @param g1.legpos position of the \code{groups} legend, see \code{g1.legend}
#' @param g2.legend plot a legend for groups in argument \code{groups2}
#' @param g2.legpos position of the \code{groups2} legend, see \code{g2.legend}
#' @param lwd.orig line-width of the original single data, set to 0 to
#' supress plotting of all original data
#' @param lty.orig line type of the original single data, set to 0 to supress
#' plotting of all original data
#' @param nrow number of plot rows, number of columns will be selected
#' automatically; NOTE: to change the ordering of the plots
#' you can change the ordering of the input \code{groups}/
#' @param mai set the outer margins around plot areas, see ?par
#' @param xmgp x-axis \code{mgp} settings: margin line in \code{mex} units
#' of axis title, tick labels and axis line; used only if \code{no.par==FALSE}
#' @param ymgp y-axis \code{mgp} settings: margin line in \code{mex} units
#' of axis title, tick labels and axis line
#' @param ytcl y-axis \code{tcl} settings: tick length and direction
#' @param xaxis plot x-axis if TRUE
#' @param yaxis the data types for which axes are to be plotted, corresponds
#' to the order in argument \code{yids} and as plotted in the legend
#' (see argument \code{g1.legend}
#' @param embed setting TRUE allows to embed plots of single groups within
#' in layouted plots, see ?layout and par("mfcol"); also see argument
#' \code{no.par}
#' @param no.par setting TRUE supresses all internal plot defaults (e.g.,
#' mai, mgp), useful to style your own plots, also see argument \code{embed}
#' @seealso \code{\link{viewPlate}}, \code{\link{getGroups}},
#' \code{\link{readPlateMap}}
#' @param verb print messages if true
#' @examples
#' data(ap12)
#' groups <- getGroups(plate=ap12plate, by=c("strain"))
#' vg <- viewGroups(ap12data,groups=groups,lwd.orig=0.1,nrow=1)
#' vg <- viewGroups(ap12data,groups2=groups,lwd.orig=0.1,nrow=1)
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
viewGroups <- function(data, groups, groups2,
                       xid, xscale=FALSE, xlim, xlab,
                       yids, dtype="data", ylab, pcols,group2.col,
                       yscale=TRUE,ylims,ylim, log="",
                       show.ci95=TRUE,show.mean=TRUE,emphasize.mean=FALSE,
                       lty.orig=1,lwd.orig=0.1,lty.mean=1,lwd.mean=2,
                       g2.legpos="topleft", g2.legend=TRUE,
                       embed=FALSE, no.par=FALSE,
                       mai=c(0.5,0,0,0), xmgp=c(1.5,.5,0),
                       ymgp=c(0,-1,-.05), ytcl=.25,
                       nrow=1, xaxis=TRUE, yaxis=c(1,2),
                       g1.legpos="topright", g1.legend=TRUE, verb=TRUE) {


    if ( missing(groups) ) {
        groups <- list(unlist(groups2))
        names(groups) <- "*"
    }
    wells <- unique(unlist(groups))

    ## filter for present wells
    pwells <- unique(c(sapply(data$dataIDs,
                              function(id) colnames(data[[id]][[dtype]]))))
    if ( sum(!wells%in%pwells)>0 ) {
        warning("wells ", wells[!wells%in%pwells]," not present, skipped!")
        wells <- wells[wells%in%pwells]
    }

    ## get x-axis data: time and temperature or another data set
    if ( missing(xid) )
        xid <- data$xids[1]

    ## TODO: interpolate data here on the fly, if not done
    ## upon parsing data or subsequently
    global.x <- xid %in% data$xids
    if ( global.x )
        time <- data[[xid]]
    else if ( xid %in% data$dataIDs )
        xdat <- data[[xid]][[dtype]][,wells,drop=FALSE]
    else
      stop("x-axis data: \"", xid, "\" not found")
    if ( verb )
      cat(paste("x-axis:", xid, "\n"))

    ## COLOR SELECTION default
    ## get plot params - colors
    ## TODO: also allow individual group2 colors
    ## (if only one data set is plotted)
    if ( missing(pcols) ) # if not missing: overrides set&default colors
        if ( !is.null(data$colors) )
            pcols <- data$colors # set colors
        else
            pcols <- getColors(data$dataIDs) # default colors

    ## reduce data to plotted data
    data <- data[data$dataIDs]
    ptypes <- names(data)
    if ( !missing(yids) ) # only use requested data for plots
      ptypes <- ptypes[ptypes%in%yids]
    if ( length(ptypes)==0 ) {
        stop("no data to plot")
    }else if ( verb )
      cat(paste("y-axis:", paste(ptypes,collapse=";"),"\n"))

    ## plot params
    if ( missing(xlim) )
      if ( !global.x )
        xlim <- range(c(xdat),na.rm=TRUE)
      else
        xlim <- range(time,na.rm=TRUE)
    ## set local ylims
    ## TODO: generalize and align with code in viewPlate
    if ( missing(ylim) ) {
        ylim <- list()
        for ( k in 1:length(ptypes) ) {
            dat <- data[[ptypes[k]]][[dtype]][,wells,drop=FALSE] # get plotted wells
            if ( global.x ) {
                dat <- dat[time>=xlim[1]&time<=xlim[2],,drop=FALSE]
            } else {
                dt <- NULL
                for ( i in 1:ncol(dat) )
                    dt <- c(dt,dat[xdat[,i]>=xlim[1]&xdat[,i]<=xlim[2],i])
                dat <- dt
            }
            ylm <- range(c(dat[is.finite(dat)]),na.rm=TRUE)
            ylim <- append(ylim,list(ylm))
        }
        names(ylim) <- ptypes
        if ( !missing(ylims) ) # update by argument ylims
            ylim[names(ylims)] <- ylims[names(ylims)]
        ylims <- ylim
    } else if ( missing(ylims) & !missing(ylim) )  { # expand single ylim
        ylims <- rep(list(ylim),length(ptypes))
        names(ylims) <- ptypes
    }

    ## colors
    #if ( missing(pcols) )
    #    pcols <- getColors(ptypes)

    ## colors - TODO - add only missing colors
    #if ( !"colors" %in% names(data) )
    #    ## as color palette 1:n, but in RGB to allow alpha
    #    pcols <- getRGB(length(ptypes))
    #    names(pcols) <- ptypes
    #}
    orig.par <- NULL
    if ( length(groups)>1 | !embed ) {
        ncol <- ceiling(length(groups)/nrow)
        orig.par <- list(mfcol=par("mfcol"))
        par(mfcol=c(nrow,ncol))
    }
    if ( !no.par ) {
        orig.par <- append(orig.par, par(c("mai","mgp")))
        par(mai=mai,mgp=xmgp)
    } else { ## THIS OVERRIDES options ymgp and ytcl with external
        ytcl <- par("tcl")
        ymgp <- par("mgp")
    }
    for ( g in 1:length(groups) ) {

        wells <- groups[[g]]
        wells <- wells[wells%in%pwells] # filter for present wells

        id <- names(groups)[g]
        ## x data other then time
        if ( !global.x )
          x <- xdat[,wells]
        else x <- time

        ## get subgroups:
        ## each subgroup must be fully contained in
        ## the group!
        ## TODO: better check of conflicts!
        if ( !missing(groups2) ) {
            gidx <- NULL
            if ( id=="*" )
              gidx <- 1:length(groups2)
            else
              for ( i in 1:length(groups2) )
                if ( sum(wells%in%groups2[[i]])==length(groups2[[i]]) )
                  gidx <- c(gidx,i)
                else if ( any(wells%in%groups2[[i]]) )
                  stop("Some but not all wells from groups group ",id,
                       " in groups2 group ",names(groups2)[i],
                       ". Groups in groups2 MUST be sub-sets of groups groups")
            sgroups <- groups2[gidx]
        } else
          sgroups <- groups[g]

        parnew <- FALSE
        for ( i in 1:length(ptypes) ) {
            ptyp <- ptypes[i]

            ## COLOR SELECTION group 1
            col.orig <- pcols[ptyp]
            ## TODO: unify !global.x and group2.col, they seem similar
            ## override colors if
            orig.cols <- NA
            if ( !global.x ) # there is no common x-axis
                orig.cols <- getColors(names(sgroups))
            if ( !missing(group2.col) ) {
                ## colors were explicitly provided
                orig.cols <- group2.col[names(sgroups)]
            }
            for ( sg in 1:length(sgroups) ) {
                wells <- sgroups[[sg]]
                wells <- wells[wells%in%pwells] # filter for present wells
                sid <- names(sgroups)[sg]
                ## x data other then time
                if ( !global.x )
                  x <- xdat[,wells]
                else x <- time


                ## get data for selected wells
                dat <- data[[ptyp]][[dtype]][,wells]
                ## calculate stats only for common x!
                ## TODO: instead bin data on x and calculate ci there
                ## or interpolate data to common x (on the fly)?
                ## TODO: do we get NAs or empty vals from ci?
                if ( is.null(dim(x)) & length(wells)>1 ) {
                    mn <- apply(dat,1,function(x) mean(x,na.rm=TRUE))
                    if ( show.ci95 )
                        ci <- apply(dat,1,function(x) ci95(x,na.rm=TRUE))
                }
                ## PLOT
                par(new=parnew) #i!=1)
                parnew <- TRUE
                ## override lty.orig=0 and lwd.orig=0 if x is data-specific
                if ( !global.x )
                    col.orig <- orig.cols[sid]
                if ( !is.null(dim(x)) | length(wells)==1 ) {
                    if ( lwd.orig==0 ) lwd.orig <- 0.1
                    lty.orig <- ifelse(g2.legend,sg,
                                ifelse(lty.orig==0,1,lty.orig))
                    #col.orig <- ifelse(g
                }
                ## COLOR SELECTION group 2
                g2.col <- col.orig
                ## group2 colors
                ## TODO: line-types by data, if group2.col is used!
                if ( !missing(group2.col) ) {
                    g2.col <- group2.col[sid]
                    ##lty.mean <- i
                    ##cat(paste(sg, g2.col, "\n"))
                }

                ## override color to allow lwd.orig=0 to work for PDF as well
                tmp <- ifelse(lwd.orig==0,NA, col.orig)

                if ( log=="y"&ylims[[ptyp]][1]<=0 ) {
                    ylims[[ptyp]][1] <- .01
                    warning("ylim <= 0 for ", ptyp, ". Raised to ",
                            ylims[[ptyp]][1], " for log='y'")
                }
                
                matplot(x,dat,type="l",lty=lty.orig,lwd=lwd.orig,axes=FALSE,
                        ylab=NA,xlab=NA,
                        ylim=ylims[[ptyp]],col=tmp,xlim=xlim,log=log)
                ## plot mean and confidence intervals
                if ( is.null(dim(x)) & length(wells)>1 ) { # only for common x!
                    if ( show.ci95 ) {

                        px <- c(x,rev(x))
                        py <- c(mn-ci,rev(mn+ci))
                        px <- px[!is.na(py)]
                        py <- py[!is.na(py)]
                        px <- c(px,px[1])
                        py <- c(py,py[1])
                        if ( log=="y" & any(py<0) ) {
                            py[py<0] <- ylims[[ptyp]][1]
                            warning("mean - 95% c.i. is below 0;",
                                    "raised to ylim for log. y-axis polygon: ",
                                    min(py))
                        }
                        polygon(x=px,y=py,border=NA,
                                col=paste(g2.col,"55",sep=""))
                    }
                    if ( show.mean )
                        lines(x=x,mn,col=ifelse(emphasize.mean,1,g2.col),
                              lwd=lwd.mean,
                              lty=ifelse(!missing(group2.col),i,sg))
                }

                ## add axes for first two values
                acol <- ifelse(global.x & missing(group2.col),col.orig,1)
                if ( yaxis[1]==i & sg==1 ) 
                    axis(2, tcl=ytcl, mgp=ymgp, col=acol, col.axis=acol)
                if ( yaxis[2]==i & sg==1 ) 
                    axis(4, tcl=ytcl, mgp=ymgp, col=acol, col.axis=acol)
            }
        }
        if ( length(sgroups)>1 & g2.legend )
            if ( global.x & missing(group2.col) )
                legend(g2.legpos,names(sgroups),lty=1:length(sgroups),
                       col=1,bty="n")
            else {
                lty <- ifelse(!missing(group2.col),
                              1:length(ptypes),
                              1:length(sgroups))
                legend(g2.legpos,names(sgroups),lty=lty,
                       col=orig.cols,bg="#FFFFFFAA",box.lwd=NA) # TODO: use g2cols
            }
        else if ( g2.legend )
            legend(g2.legpos,id, bty="n")
        if ( xaxis ) axis(1)
        if ( missing(xlab) )
          xlab <- xid
        mtext(xlab, 1, par("mgp")[1])
        if ( !missing(ylab) )
            mtext(ylab, 2, par("mgp")[1])
  
    }
    ## add legend to last plot
    if ( g1.legend ) {
      if ( !missing(group2.col) ) {
        lty <- 1:length(ptypes)
        col <- 1
      } else {
        lty <- 1
        col <- pcols[ptypes]
      }
      legend(g1.legpos,ptypes,lty=lty,col=col,bg="#FFFFFFAA")
    }
    ## reset par!
    par(orig.par)

    ## TODO: return meaningful and/or non-plotted information
    ## assigning it makes it silent!
    #if ( global.x ) xid <- "Time"
    plotparams <- list(ylims=ylims, xid=xid, xlim=xlim,  colors=pcols, orig.cols=orig.cols)
}


#' View Plate
#'
#' plots all data in plate format for a quick inspection of
#' the current data set. This function should be called immediately
#' after loading a new platereader dataset to inspect the data
#' and decide on further processing steps.
#' @param data the list of measurement data as provided by
#' \code{\link{readPlateData}}
#' @param wells a list of wells to plot, overrules \code{rows} and \code{cols}
#' @param wcols named color vector for wells, names should correspond to
#' column IDs in \code{data}, ie. the column 'well' in the plate layout map
#' @param rows a list of strings/characters used as row ID in the composite
#' row:col well description in the plate layout (map) and plate data
#' @param cols as rows but plate column IDs
#' @param xid ID of a data-set in the input data that can be used as x-axis
#' instead of the default Time vector
#' @param yids IDs of the data to be plotted; if missing, all data will
#' be plotted
#' @param dtype type of the data to be plotted, default is the main 'data', but
#' e.g. 'orig' allows to plot the original data w/o processing (e.g.
#' interpolation, blank correction, etc.)
#' @param xscale use a global range for the x-axes; only relevant if xid specifies a subset of the data as x-axis
#' @param xlim plot range of the x-axis
#' @param pcols a named list of RGB colors to be used the plotted data types; the color vector must have names according to the data IDs
#' @param yscale if TRUE (default) global y-axis limits will be calculated from
#' all plotted wells; if FALSE each well be locally scaled
#' @param ylims a named list of y-axis ranges pairs for each data ID
#' @param ylim one y-axis limit range that will be used for all plotted data
#' @param log plot logarithmic axis, use equivalent to normal plot 'log', i.e.,
#' log="y" for a log y-axis, log="x" for x-axis and log="yx" for both axes
#' @param axes axes numbers to draw, integers from 1 to 4 indicate the side
#' to plot (argument \code{side} in function \code{\link[graphics:axis]{axis}});
#' multiple y-axis values (multiple yids) are handled in argument\code{yaxis}
#' @param yaxis the ID of maximally two data types for which y-axes are to be plotted
#' @param legpos position of the well IDs on the plots
#' @param add.legend add a legend for the plotted data types (see
#' argument \code{yids}) in the last plotted well
#' @examples
#' data(ap12)
#' # view all data on the plate
#' viewPlate(ap12data)
#' # inspect natural logarithm of OD_600 values, ie. log(X(t))
#' viewPlate(ap12data, yids="600", log="y")
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
viewPlate <- function(data, wells, wcols,
                      rows=toupper(letters[1:8]),cols=1:12,
                      xid, xscale=FALSE,xlim,
                      yids, dtype="data", pcols, yscale=TRUE, ylims, ylim,
                      log="", axes, yaxis="", legpos="topleft", add.legend=TRUE) {

    ## which wells to plot?
    if ( missing(wells) ) {
        wells <-  paste(rep(rows,each=length(cols)),cols,sep="")
    } else {
        ## TODO: allow multiple cols, by groups
        ##       allow to plot single rows or columns
        rows <- wells
        cols <- ""
    }
    if ( missing(wcols) ) {
        wcols <- rep(1,length(wells))
        names(wcols) <- wells
    }

    ## filter for present wells
    pwells <- unique(unlist(sapply(data$dataIDs,
                                   function(id)
                                       colnames(data[[id]][[dtype]]))))
    if ( sum(!wells%in%pwells)>0 ) {
        #warning("wells ", wells[!wells%in%pwells]," not present, skipped!")
        wells <- wells[wells%in%pwells]
    }

    ## get x-axis data: time and temperature or another data set
    if ( missing(xid) )
        xid <- data$xids[1]

    ## TODO: interpolate data here on the fly, if not done
    ## upon parsing data or subsequently
    global.x <- xid %in% data$xids
    if ( global.x ) {
        time <- data[[xid]]
    } else if ( xid %in% data$dataIDs )
        xdat <- data[[xid]][[dtype]][,wells,drop=FALSE]
    else
        stop("x-axis data: \"", xid, "\" not found")
    cat(paste("x-axis:", xid, "\n"))
    ## get x-data: global data (time, temperature)
    ## or a data set specified via xid
    ## TODO: implement absence of master time, if interpolate=FALSE
    ## upon reading of data
    ## TODO: adapt to new data$xids,
    ## if xid is not present in data$xids, get it from dataIDs
    #time <- data[[data$xids[1]]]
    #if ( !missing(xid) ) { # or use other data-set as x
    #    xdat <- data[[xid]][[dtype]][,wells,drop=FALSE]
    #} else xid <- NULL

    ## get plot params - colors
    if ( missing(pcols) ) # if not missing: overrides set&default colors
        if ( !is.null(data$colors) )
            pcols <- data$colors # set colors
        else
            pcols <- getColors(data$dataIDs) # default colors

    ## reduce matrix to plot data
    data <- data[data$dataIDs]
    ptypes <- names(data)
    if ( !missing(yids) ) # only use requested data for plots
      ptypes <- ptypes[ptypes%in%yids]
    if ( length(ptypes)==0 ) {
        cat(paste("no data to plot\n"))
        return()
    } else
        cat(paste("plotting", paste(ptypes,collapse=";"),"\n"))

    ## PLOT PARAMS
    ## xlim
    if ( missing(xlim) )
      if ( !global.x )
        xlim <- range(c(xdat),na.rm=TRUE)
      else
        xlim <- range(time)
    else
        xscale <- TRUE

    ## set local ylims
    ## TODO: generalize and align with code in viewGroup
    if ( missing(ylim) ) {
        pidx <- 1:length(ptypes)
        ylim <- rep(list(c(NA,NA)),length(ptypes)) # initiliaze
        for ( k in pidx ) {
            dat <- data[[ptypes[k]]][[dtype]][,wells,drop=FALSE] # get plotted wells
            if ( global.x ) # if x is time, limit to plot xlim
              dat <- dat[time>=xlim[1]&time<=xlim[2],,drop=FALSE]
            ylm <- range(c(dat[is.finite(dat)]),na.rm=TRUE) # get range
            ylim[[k]] <- ylm
        }
        names(ylim) <- ptypes
        if ( !missing(ylims) ) # update by argument ylims
            ylim[names(ylims)] <- ylims[names(ylims)]
        ylims <- ylim
    } else  { # .. or expand single ylim argument c(y1,y2)
        ylims <- rep(list(ylim),length(ptypes))
        names(ylims) <- ptypes
    }

    # y-axis: take first two data sets as y-axis it was requested
    if ( !missing(axes) )
      if ( any(c(2,4) %in% axes) )
        if ( missing(yaxis) )
          yaxis <- na.omit(ptypes[1:sum(c(2,4) %in% axes)])

    ## plot plate
    orig.par <- par(c("mfcol","mai")) # store parameters
    par(mfcol=c(length(rows),length(cols)),mai=rep(0,4))
    for ( j in cols )
      for ( i in rows ) {
          well <- paste(i,j,sep="")
          if ( !well%in%wells ) { ## skipped wells
              plot(1,axes=FALSE,col="gray",pch=4,cex=5)
              next
          }
          ## x data
          if ( !global.x )
            x <- xdat[,well]
          else x <- time
          for ( k in 1:length(ptypes) ) {
              # y data
              ptyp <- ptypes[k]
              y <- data[[ptyp]][[dtype]][,well]
              if ( k>1 ) par(new=TRUE)
              #cat(paste("hallo",k))
              ## TODO: obsolete? since we filter for wells above and
              ## skipWells now removes columns
              if ( !any(!is.na(y)) ) ## skip well - all NA?
                plot(1,axes=FALSE,col="gray",pch=4,cex=5)
              else if ( yscale & xscale )
                plot(x,y,axes=FALSE,xlim=xlim,ylim=ylims[[ptyp]],
                     type="l",log=log,col=pcols[ptyp])
              else if ( yscale & !xscale )
                plot(x,y,axes=FALSE,ylim=ylims[[ptyp]],
                     type="l",log=log,col=pcols[ptyp])
              else if ( !yscale & xscale )
                plot(x,y,axes=FALSE,xlim=xlim,
                     type="l",log=log,col=pcols[ptyp])
              else if ( !yscale & !xscale )
                plot(x,y,axes=FALSE,
                     type="l",log=log,col=pcols[ptyp])
              ## plot axes
              if ( !missing(axes) )
                for ( ax in axes ) {
                  if ( ax%in%c(1,3))
                    axis(ax)
                  else if ( ax%in%c(2,4) )
                      if ( ptyp%in%yaxis ) {
                          ## TODO: clean this mess!
                          if ( sum(c(2,4)%in%axes)==1 )
                              axs <- c(2,4)[which(c(2,4)%in%axes)]
                          else
                              axs <- c(2,4)[which(yaxis==ptyp)]
                          axis(axs,
                               col.ticks = pcols[ptyp], col.axis=pcols[ptyp],
                               tcl=-par("tcl"), mgp=c(0,-1.5,0))
                      }
                }
              
              if ( k==1 ) {
                  box()
                  legend(legpos,well,text.col=wcols[well],bty="n")
              }
          }
      }
    ## add legend to last plot
    if ( add.legend )
      legend("topright",ptypes,lty=1,col=pcols[ptypes],bg="#FFFFFFAA")

    ## reset par
    par(orig.par)

    ## TODO: return meaningful and/or non-plotted information
    ## assigning it makes it silent!
    plotparams <- list(ylims=ylims, xid=xid, xlim=xlim,  colors=pcols)
}

