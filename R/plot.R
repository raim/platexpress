
#' Plot a representation of microreader "plate"
#' 
#' @param nrow number of rows of the plate, with letter index from top to bottom
#' @param ncol number of columns of the plate
#' @param map the plate map, see  \code{\link{readPlateMap}}
#' @param color color column in the \code{map}
#' @param text text column in the \code{map}
#' @param title plot title (argument \code{main} to plot)
#' @param add add to an existing plot
#' @param ... arguments passed to \code{\link{image}}
#' @export
viewMap <- function(map, nrow=8, ncol=12, color="color", text="amount", title="plate layout", add=FALSE) {

    bright <- "#FFFFFF"
    dark <- "#000000"

    ## set-up rows and columns, using standard plate design,
    ## with letters for rows from top to bottom, and
    ## numbers for columns
    cols <- 1:ncol
    names(cols) <- cols
    rows <- 1:nrow
    names(rows) <- rev(toupper(letters[1:nrow]))
    
    ## plot canvas
    if ( !add ) {
        plot(rep(cols,nrow),rep(rows, ncol), col=NA, axes=FALSE, main=title,
             xlab=NA, ylab=NA)
        axis(1, at=cols)
        axis(2, at=rows, labels=names(rows), las=2)
    }

    ## plot wells
    for ( i in rows )
        for ( j in cols ) {
            well <- paste0(names(rows)[i],names(cols)[j])
            widx <- which(as.character(map$well)==well)
            col <- map[widx, color]

            ## skip missing wells
            if ( length(widx)==0 ) next 

            ## select text color: dark or bright based on RGB luminescence 
            crgb <- col2rgb(col)/255
            L <- 0.2126 * crgb[1,1] + 0.7152 * crgb[2,1] + 0.0722 * crgb[3,1]
            col.txt <- ifelse(L>.5, dark, bright)

            ## plot well
            points(cols[j], rows[i], col=map[widx, color], cex=3, pch=19)
            text(cols[j], rows[i], map[widx, text], cex=1, col=col.txt)
        }
}


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
#' @param col default color for arrows from positive slope segments
#' @param ncol default color for arrows from negative slope segments
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
arrows_plate <- function(segments, wells,
                         col="#000000", ncol="#F81894", wcol, xlim,
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
        ## growth phases: slope, direction as colors,lty,lwd
        ## dashed lines for negative mu
        lt <- as.numeric(sl<0) +1
        ## color
        cl <- rep(col, length(sg)-1)
        cl[sl<0] <- ncol
        if ( !missing(wcol) )
            cl <- rep(wcol[wells[i]], length(sg)-1)
        ##if ( alpha.scale )
        ##alphas <- (sl-minmu.col)/(maxmu-minmu.col) 
        ##cl <- add_alphas(cl,alpha=alphas)
        ## TODO: better control of scales and or alpha?
        lwsl <- abs(sl)
        lwsl[lwsl>maxmu] <- maxmu
        lwsl[lwsl<minmu] <- minmu
        ## line width
        lw <- lwd.max*lwsl
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


## TODO: wrapper for viewGroups taking list of plates
## and a merged layout file
## TODO: rename groups and groups2, eg. into sets and replicates
viewGroupl <- function(data, groups, groups2, ...) {}

#' View Group Averages
#'
#' plot grouped wells as summary plots, incl. means and confidence intervals.
#' This function gives a first overview of the reproducability between
#' replicates.
#' @param data the platereader data object as provided by
#' \code{\link{readPlateData}}
#' @param groups a list of well groups, wells from each
#' group will be plotted together in one plot, and statistics plotted
#' for sub-groups (replicates) in \code{groups2}. Groupings can also
#' be present as a list item in the \code{data} object.
#' @param groups2 replicate well groups for which statistics are calculated
#' and plotted  (see argument \code{stats}), each group must be a subset
#' of exactly on well group in argument \code{groups}. 
#' @param xid ID of a data-set in the input data that can be used as x-axis
#' instead of the default Time vector
#' @param yids IDs of the data to be plotted; if missing, all data will
#' be plotted
#' @param stats type of replicate statistics to plot, where replicates
#' are defined in option \code{groups2} or defined in the plate data object.
#' Options are either only a line for each replicate group ("mean" or
#' "median") or a line with transparent ranges. "CI": mean and 95\%
#' confidence interval (default), "SD": mean and standard deviation,
#' "SE": mean and standard error, "boxplot": median and 25\% and 75\%
#' quantiles, "range": median and full data range (min/max)
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
#' @param show.ci95 deprecated, but kept for backward compatibility, see
#' option \code{stats="CI"}: show the 95\% condifence intervals of groups;
#' set \code{stats="mean"} to reproduce previous \code{show.ci95=FALSE}
#' behaviour
#' @param show.mean show the mean of groups as lines
#' @param emphasize.mean show the means in black instead of group colors,
#' @param lwd.mean line width for the mean lines, see \code{(par("lwd")}
#' @param lty.mean line style for the mean lines, see \code{(par("lty")}
#' this can help to emphasize the means in noisy data (broad overlapping
#' confidence intervals
#' @param alpha hexadecimal alpha value for range plot opaqueness
#' @param lwd.orig line-width of the original single data, set to 0 to
#' supress plotting of all original data
#' @param lty.orig line type of the original single data, set to 0 to supress
#' plotting of all original data
#' @param g1.legend show the main legend, giving plot colors
#  of the plotted data types (\code{yids})
#' @param g1.legpos position of the \code{groups} legend, see \code{g1.legend}
#' @param g2.legend plot a legend for groups in argument \code{groups2}
#' @param g2.legpos position of the \code{groups2} legend, see \code{g2.legend}
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
#' @param verb print progress messages
#' @examples
#' data(ap12)
#' groups <- getGroups(plate=ap12plate, by=c("strain"))
#' vg <- viewGroups(ap12data,groups=groups,lwd.orig=0.1,nrow=1)
#' vg <- viewGroups(ap12data,groups2=groups,lwd.orig=0.1,nrow=1)
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
viewGroups <- function(data, yids, stats="CI", groups, groups2, 
                       xid, xscale=FALSE, xlim, xlab,
                       dtype="data", ylab, pcols,group2.col,
                       yscale=TRUE,ylims,ylim, log="",
                       lty.mean=1,lwd.mean=2,alpha=55,lty.orig=1,lwd.orig=0.1,
                       show.ci95=FALSE,show.mean=TRUE,emphasize.mean=FALSE,
                       g2.legpos="topleft", g2.legend=TRUE,
                       embed=FALSE, no.par=FALSE,
                       mai=c(0.5,0,0,0), xmgp=c(1.5,.5,0),
                       ymgp=c(0,-1,-.05), ytcl=.25,
                       nrow=1, xaxis=TRUE, yaxis=c(1,2),
                       g1.legpos="topright", g1.legend=TRUE, verb=FALSE) {

    ## STATISTICS: deprecated argument
    if ( show.ci95 ) stats <- "CI"
    
    ## PARSE GROUPINGS
    if ( "groups" %in% names(data) ) {
        if ( missing(groups2) & "group2" %in% names(data$groups) )
            groups2 <- data$groups$group2
        if ( missing(groups) & "group1" %in% names(data$groups) )
            groups <- data$groups$group1
        if ( missing(group2.col) & "group2.color" %in% names(data$groups) )
            group2.col <- data$groups$group2.color
    }

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
    if ( (length(groups)>1 | !embed) & !no.par ) {
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
                cat(paste("USE GROUP COLORS\n"))
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


                ## STATISTICS
                ## get data for selected wells
                dat <- data[[ptyp]][[dtype]][,wells, drop=FALSE]
                ## calculate stats only for common x!
                ## TODO: instead bin data on x and calculate ci there
                ## or interpolate data to common x (on the fly)?
                ## TODO: do we get NAs or empty vals from ci?
                if ( is.null(dim(x)) & length(wells)>=1 ) {
                    linefun <- mean
                    if ( stats%in%c("median","range") )
                        linefun <- median
                    if ( stats!="boxplot" ) # via boxplot.stats below
                        mn <- apply(dat,1,function(x) linefun(x,na.rm=TRUE))
                    ci1 <- NULL
                    if ( stats=="CI" ) {
                        ci <- apply(dat,1,function(x) ci95(x,na.rm=TRUE))
                        ci1 <- mn+ci
                        ci2 <- mn-ci
                    } else if ( stats=="SD" ) {
                        ci <- apply(dat,1,function(x) sd(x,na.rm=TRUE))
                        ci1 <- mn+ci
                        ci2 <- mn-ci
                    } else if ( stats=="SE" ) {
                        ci <- apply(dat,1,function(x) se(x,na.rm=TRUE))
                        ci1 <- mn+ci
                        ci2 <- mn-ci
                    } else if ( stats=="boxplot" ) {
                        ## NOTE: overrules median
                        ci <- apply(dat,1,function(x)
                            boxplot.stats(x)$stats[2:4])
                        ci1 <- ci[1,]
                        mn  <- ci[2,]
                        ci2 <- ci[3,]
                    } else if ( stats=="range" ) {
                        ci1 <-  apply(dat,1,function(x) min(x,na.rm=TRUE))
                        ci2 <-  apply(dat,1,function(x) max(x,na.rm=TRUE))
                        ci1[!is.finite(ci1)] <- NA
                        ci2[!is.finite(ci2)] <- NA
                    }
                }
                ## PLOT
                par(new=parnew) #i!=1)
                parnew <- TRUE
                ## override lty.orig=0 and lwd.orig=0 if x is data-specific
                ## 20190619: the if was commented out - why?
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
                    ## plot statistics polygon
                    if ( !is.null(ci1) ) {

                        px <- c(x,rev(x))
                        py <- c(ci1,rev(ci2))
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
                                col=paste(g2.col,alpha,sep=""))
                    }
                }
                if ( is.null(dim(x)) & length(wells)>=1 ) { # only for common x!
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
#' @param verb print progress messages
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
                      log="", axes, yaxis="", legpos="topleft", add.legend=TRUE,
                      verb=FALSE) {
    ## TODO: ADAPT viewPlate to use rows and cols from data$wells!!?
    ## currently for BioLector only!!
    if ( "wells"%in%names(data) ) {
        rows <- data$wells$rows
        cols <- data$wells$cols
    }
    
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
        if ( verb )
            cat(paste("wells ", wells[!wells%in%pwells],
                      " not present, skipped!\n"))
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
    if ( verb ) cat(paste("x-axis:", xid, "\n"))
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
        if ( verb ) cat(paste("no data to plot\n"))
        return()
    } else
        if ( verb ) cat(paste("plotting", paste(ptypes,collapse=";"),"\n"))

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

#' Dose-Response Box-Plots
#'
#' Plots a continuous value for each well as a function
#' of an \code{dose} of a \code{substance} and make boxplots
#' for all replicates at a given dose. Note that this is deprecated,
#' and plots with error bars can be generated instead by
#' function \code{\link{doseResponse}}.
#' @param map a plate layout map with columns \code{dose} and some
#' calculated value in column \code{response},
#' eg. results from \pkg{grofit} or \pkg{growthrates}
#' @param wells a list of well IDs to be used in the plot
#' @param response the name of a column in \code{map} containing numeric values
#' that should be plotted on y-axis
#' @param dose the name of a column in \code{map} providing numeric values
#' that should be plotted on x-axis, typically a substance added to wells
#' in multiple replicates; the default value "amount" is automatically
#' generated from appropriate plate layout files by \code{link{readPlateMap}}
#' @param substance the name of a column in \code{map} providing the names of
#' the substance in \code{dose}, used as x-axis label; if it doesn't
#' correspond to a column, its value is directly used as label
#' @param color the name of a column in \code{map} providing colors
#' for each well; NOTE, that wells with the same \code{dose} should have
#' the same color, only the first color for a given \code{dose} is used;
#' Colors are automatically assigned from a color ramp mapped to the
#' numerical range of \code{dose} by \code{link{readPlateMap}}
#' @param na.y value to be used to plot replicates with \code{NA} in
#' column \code{response}; set to NA to supress plotting
#' @param ylim limits of the y-axis
#' @param xnum use numerical x-axis instead of default categorical
#' @param xlab alternative label for the x-axis, default is to use argument
#' \code{substance}, or if this is a column in \code{map}, the substance
#' indicated there
#' @param ylab alternative label for the y-axis, default is to use argument \code{response}
#' @export
doseResponse.box <- function(map, wells, response, dose="amount",
                             substance="substance",
                             color="color", na.y=0, ylim, xnum=FALSE,
                             xlab, ylab) {

    if( missing(wells) ) wells <- map[,"well"]

    if ( missing(xlab) )
        xlab <- substance
    if ( missing(ylab) )
        ylab <- response
    
    wells <- match(wells,map[,"well"])

    y <- map[wells,response]
    x <- map[wells,dose]

    ## get unique colors for unique sorted x, as it will appear in
    ## boxplots
    ## UGLY, TODO: less ugly?
    if ( color %in% colnames(map) ) {
        cl <- map[wells,color]
        cl <- cl[order(x)]
        cl <- cl[which(!duplicated(sort(x)))]
        names(cl) <- sort(x[!duplicated(x)])
    } else cl <- 1

    if ( missing(ylim) )
        ylim <- range(y,na.rm=TRUE)

    if ( substance %in%  colnames(map) )
        subid <- map[wells,substance][1]
    else subid <- substance

    #xlim <- range(x)
    #if ( is.na(na.y) ) # if na.y are not plotted limit xaxis
    #  xlim <- range(x[!is.na(y)])

    y1 <- y
    x1 <- x

    if ( is.na(na.y) ) {
      y1 <- y1[!is.na(y)]
      x1 <- x1[!is.na(y)]
    }

    xaxt="s"
    if ( xnum ) {# attempt numerical x-axis
      x1 <- factor(x1, levels = do.call(seq, as.list(range(x1))))
      xaxt <- "n"
    }
    cl <- cl[as.character(unique(sort(x1)))]
    graphics::boxplot(y1~x1, border=cl,
                      xlab=xlab,ylab=ylab,ylim=ylim, #xlim=xlim,
                      na.action="na.pass",xaxt=xaxt, outline=FALSE)
    graphics::stripchart(y1~x1,add=T, vertical=TRUE,col=cl,
                         method="jitter", pch=1,cex=1,
                         na.action="na.pass")
    if ( xnum ) # numerical xaxis
      axis(1)

    if ( xnum )
      x <- factor(x, levels=do.call(seq, as.list(range(x))))

    nay <- y
    nay[is.na(y)] <- na.y
    nay[!is.na(y)] <- NA
    graphics::stripchart(nay~x,add=T, vertical=TRUE,col="red",
                         method="jitter", pch=4,cex=1/2,
                         na.action="na.pass")
}

#' Dose-Response Plots with Error Bars
#'
#' plots a continuous value for each well as a function
#' of an \code{dose} of a \code{substance} and generates error bars
#' for all replicates at a given dose; returns the plotted mean values
#' and error bar ranges.
#' @param map a plate layout map with columns \code{dose} and some
#' calculated value in column \code{response},
#' eg. results from \pkg{grofit} or \pkg{growthrates}
#' @param wells a list of well IDs to be used in the plot
#' @param response the name of a column in \code{map} containing numeric values
#' that should be plotted on y-axis
#' @param dose the name of a column in \code{map} providing numeric values
#' that should be plotted on x-axis, typically a substance added to wells
#' in multiple replicates; the default value "amount" is automatically
#' generated from appropriate plate layout files by \code{link{readPlateMap}}
#' @param substance the name of a column in \code{map} providing the names of
#' the substance in \code{dose}, used as x-axis label; if it doesn't
#' correspond to a column, its value is directly used as label
#' @param col either a valid color representation or the name of a column in
#' \code{map} providing colors for each well or a, if \code{wells} is a named
#' list of wells, a named vector of colors with matching names. NOTE for the
#' second case, that wells with the same \code{dose} should have the same color, only
#' the first color for a given \code{dose} is used; such colors are automatically
#' assigned from a color ramp mapped to the numerical range of \code{dose} by
#' \code{\link{readPlateMap}}, see \code{\link{amountColors}};
#' @param pch pch of the mean value points
#' @param bartype type of the error bars, "range" for the full range
#' of the data, where the average value will be the
#' \code{\link[stats:median]{median}}; or "ci95" for the 95\% confidence
#' interval, "sd" for the standard deviation, or "se" for the standard error
#' where the average value will be the \code{\link[base:mean]{mean}};
#' TODO: implement boxplot-like quantiles
#' @param barl length of the horizontal lines of error bars (argument \code{length}
#' of function \code{\link{arrows}})
#' @param all logical indicating whether to plot all data points additionally to
#' error bars
#' @param line logical indicating whether to plot a line connecting the mean
#' values
#' @param na.y value to be used to plot replicates with \code{NA} in
#' column \code{response}; set to NA to supress plotting
#' @param add add data to existing plot
#' @param ylim limits of the y-axis
#' @param xlim limits of the x-axis
#' @param ylab alternative label for the y-axis, default is to use argument \code{response}
#' @param xlab alternative label for the x-axis, default is to use argument
#' \code{substance}, or if this is a column in \code{map}, the substance
#' indicated there
#' @param verb print progress messages
#' @param ... arguments passed on to the main setup \code{\link{plot}}
#' @export
doseResponse <- function(map, wells, response,
                         dose="amount", substance="substance",
                         col="black", pch=1, bartype="range", barl=.05,
                         all=FALSE, line=TRUE, na.y=0, add=FALSE,
                         ylim, xlim, ylab, xlab, verb=TRUE, ...) {

  if( missing(wells) ) wells <- map[,"well"]

  ## call recursively with add option if wells is a grouping list!
  ## NOTE: interpreting col as a group color list!
  ## TODO: use classes
  if ( typeof(wells)=="list" ) {

    ## record call
    mc <- match.call(expand.dots = FALSE)

    ## set common xlims/ylims
    if ( missing(xlim) )
      mc[["xlim"]] <- range(map[map$well%in%unlist(wells),dose],na.rm=TRUE)
    if ( missing(ylim) )
      mc[["ylim"]] <- range(map[map$well%in%unlist(wells),response],na.rm=TRUE)

    ymats <- list()
    for ( i in 1:length(wells) ) {
      mc[["wells"]] <- wells[[i]]
      mc[["add"]] <- i>1
      if ( names(wells)[i]%in%names(col) )
        mc[["col"]] <- col[[names(wells)[i]]]
      if ( verb )
        cat(paste("plotting group",i, names(wells)[i],"in color",
                  mc[["col"]],"\n"))
      #scan()
      ymat <- eval(mc)
      ymats <- append(ymats,list(ymat))
    }
    names(ymats) <- names(wells)
    ## TODO: legend?
    return(invisible(ymats))
  }

  wells <- match(wells,map[,"well"])

    ## TODO: remove blanks!

  y <- map[wells,response]
  x <- map[wells,dose]

  if ( substance %in%  colnames(map) )
    subid <- map[wells,substance][1]
  else subid <- substance
  if ( missing(xlab) )
    xlab <- subid
  if ( missing(ylab) )
    ylab <- response

  y1 <- y
  x1 <- x

  y1 <- y1[!is.na(y)]
  x1 <- x1[!is.na(y)]

  ## AVERAGE FUNCTION
  avgf <- base::mean # for ci95%, SD, SE
  if ( bartype=="range" ) avgf <- stats::median

  ## calculate ranges for duplicates at each x
  xlevels <- sort(unique(x1))
  ymat <- matrix(NA,nrow=length(xlevels),ncol=4)

  ## COLOR SELECTION
  ## TODO: organize coloring schemes
  ## default, first color from plotted xlevel or direct coloring
  ## get unique colors for unique sorted x, as it will appear in
  ## boxplots
  ## UGLY, TODO: less ugly?
  if ( col %in% colnames(map) ) {
    cl <- sapply(xlevels,function(x) map[map[,dose]==x,col][1])
    names(cl) <- xlevels
    linecol <- 1
  } else {
    cl <- rep(col,length(xlevels))
    linecol <- col
  }
  ## TODO: why does this appear in recursive call at end?
  if ( length(xlevels)==0 ) {
    warning("NO XLEVELS FOUND")
    return(NULL)
  }
  for ( i in 1:length(xlevels) ) {
    xi <- xlevels[i]
    yi <- y1[x1==xi]
    ymat[i,2] <- avgf(yi,na.rm=TRUE) # average function (mean or median)
    ## TODO: allow "quantiles" for boxplot-like plot
    if ( bartype=="range" )
      ymat[i,3:4] <- range(yi,na.rm=TRUE)
    else if ( bartype=="sd" )
      ymat[i,3:4] <- c(ymat[i,2] + c(-1,1)*sd(yi,na.rm=TRUE))
    else if ( bartype=="se" )
      ymat[i,3:4] <- c(ymat[i,2] + c(-1,1)*se(yi,na.rm=TRUE))
    else if ( bartype=="ci95")
      ymat[i,3:4] <- c(ymat[i,2] + c(-1,1)*ci95(yi,na.rm=TRUE))
    else stop("bar type ", bartype, " unknown\n")
    ymat[i,1] <- xi
  }
  colnames(ymat) <- c(dose,response,paste0(bartype,c(".min",".max")))
  if ( missing(ylim) )
    ylim <- range(c(ymat[,2:4],ifelse(all,y,NA)),na.rm=TRUE)
  if ( missing(xlim) ) {
    if ( !is.na(na.y))
      xlim <- range(x,na.rm=TRUE)
    else xlim <- range(x1,na.rm=TRUE)
  }

  ## plot
  if ( !add )
    plot(x1, y1, col=NA, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, ...)
  if ( line )
    lines(ymat[,1],ymat[,2],col=linecol)
  for ( i in 1:nrow(ymat) ) {
    if ( ymat[i,3]<ymat[i,4])
      arrows(ymat[i,1], ymat[i,3], ymat[i,1], ymat[i,4], length=barl, angle=90, code=3, col=cl[i])
    else
      points(ymat[i,1], ymat[i,3],pch=3, col=cl[i])
    if ( !line )
      points(ymat[i,1],ymat[i,2], col=cl[i], pch=pch)
    if ( all )
      points(x1[x1==ymat[i,1]], y1[x1==ymat[i,1]], pch=20, cex=.5, , col=cl[i])
  }
  if ( !is.na(na.y) & sum(is.na(y)) )
    points(x[is.na(y)], rep(na.y,sum(is.na(y))), col="red", pch=4, cex=.5)
  invisible(cbind.data.frame(ymat,color=cl,stringsAsFactors=FALSE))
}

#' Box-Plots of Data Ranges
#'
#' returns data for group of wells in a given range of the x-axis;
#' if only one value is given, the values closest to this x-axis points
#' are returned; see "Value" for details.
#' @param data \code{\link{platexpress}} data, see \code{\link{readPlateData}}
#' @param groups a grouping of wells, see \code{\link{getGroups}}
#' @param rng x-axis range or a single point, for the latter the closest
#' point will be selected or, if option \code{interpolate=TRUE} the
#' value at \code{rng} will be interpolated using R base function
#' \code{\link[stats:approxfun]{approx}}
#' @param xid the x-axis ID, if multiple x-axes are present
#' @param yid the y-axis data to be grouped
#' @param interpolate interpolate data at \code{rng} (if a single value,
#' see argument \code{rng})
#' @param plot if TRUE a box-plot or bar-plot will be plotted
#' @param type either "box" (default) or "bar" for box-plot or bar-plot
#' @param etype type of statistics to be used for error bars in the bar-plot,
#' either "ci" (default) for the 95%-confidence interval or "se" for
#' the standard error
#' @param ylim y-axis limits
#' @param ylab y-axis label
#' @param ... further arguments to boxplot/barplot
#' @return Returns an annotated \code{data.frame} of the values, with well
#' and group IDs in the first two columns. The values in the third column are
#' raw values (if argument \code{rng} was a single point
#' on the x-axis) or mean values (if argument\code{rng} was a range).
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
boxData <- function(data, rng, groups, xid, yid="OD", interpolate=TRUE, plot=TRUE, type="box", etype="se", ylim, ylab, ...) {

    if ( missing(xid) )
        xid <- data$xids[1]
    ## cut data to selected range (closest point if length(rng)==1)
    if ( length(rng)==2 ) {
      cdat <- cutData(data, xid=xid, xrng=rng)
    } else {
      if ( !interpolate ) {
        cdat <- cutData(data, xid=xid, xrng=rng)
        ## get the actual point if only one points was chosen
        if ( length(rng)==1) rng <- signif(unique(range(cdat[[xid]])),4)
      } else { # interpolate using getData!
        cdat <- data
        cdat[[yid]]$data <- t(as.matrix(getData(data, ID=yid, xrng=rng, xid=xid, verb=TRUE)))
      }
    }
    bdat <- rep(list(NA),length(groups))
    names(bdat) <- names(groups)
    for ( sg in 1:length(groups) )
        bdat[[sg]] <- cdat[[yid]]$data[,groups[[sg]],drop=FALSE]
    ## get means for all groups
    pdat <- lapply(bdat, function(x) apply(x,2,mean,na.rm=TRUE))

    if ( plot ) {
        if ( missing(ylab) ) ylab <- yid
        ##par(mai=c(1,.75,.1,.1))
        if ( type=="box" ) {
            if ( missing(ylim) )
                ylim <- range(pdat, na.rm=TRUE)
            boxplot(pdat,ylab=ylab,las=2, ylim=ylim, ...)
        } else if ( type=="bar" ) {
            mn <- unlist(lapply(pdat, mean,na.rm=TRUE))
            if ( etype=="ci" ) # 95% confidence interval
                ci <- unlist(lapply(pdat, ci95,na.rm=TRUE))
            else if ( etype=="se" ) { # standard error sd/n^2
                sd <- unlist(lapply(pdat, sd, na.rm=TRUE))
                n2 <- sqrt(unlist(lapply(pdat, length)))
                ci <-  sd / n2
            }
            if ( missing(ylim) )
                ylim <- c(0,max(mn+ci,na.rm=TRUE))
            x <- barplot(mn,ylim=ylim,ylab=ylab,las=2, ...)
            arrows(x0=x,x1=x,y0=mn-ci,y1=mn+ci,code=3,angle=90,
                   length=.05,lwd=1.5)
        }
        legend("topright",paste("at",xid, "=",paste(rng,collapse="-")),
               bty="n",box.lwd=0)
    }

    ## if a range was chosen and not a single point,
    ## report mean values
    if ( length(rng)>1 ) {
        bdat <- lapply(pdat,function(x) {
            y<-matrix(x,nrow=1)
            colnames(y) <- names(x)
            y})
    }
    ## summarize results if only one value was requested!
    tmp <- data.frame(well=unlist(lapply(bdat, function(x) colnames(x))),
                      group=unlist(sapply(1:length(bdat),
                        function(x) rep(names(bdat)[x],length(bdat[[x]])),
                        simplify = FALSE)),
                      data=unlist(bdat))
    colnames(tmp)[3] <- paste(yid,paste(rng,collapse="-"),sep="@")
    bdat <- tmp
    rownames(bdat) <- NULL

    result <- bdat
}

