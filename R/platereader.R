### ANALYZING PLATE-READER GROWTH & EXPRESSION CURVES
#' platexpress: A package for analysing microbial growth & expression.
#'
#' The platexpress package provides a quick&easy interface to 
#' microbial growth & gene expression data as measured in typical
#' microplate-readers or other parallel growth systems.
#' 
#'@author Rainer Machne
#'@docType package
#'@name platexpress
#'@section Dependencies: The package uses mostly functionality from R base,
#' (graphics, grDevices, stats) but more functionality is available when
#' \code{\link[grofit:grofit]{grofit}} is installed.
#'@importFrom stats median sd qt approx spline filter
#'@importFrom graphics plot matplot boxplot barplot legend arrows locator
#' abline lines points polygon box axis par text title mtext
#'@importFrom grDevices rainbow rgb col2rgb png pdf postscript graphics.off
#'@importFrom utils read.csv read.table
NULL


### UTILS

#' plots to png, eps or pdf, taking the same arguments
#' @param file the name of the file without the extension (.pdf)
#' @param type type of the figure, either png, pdf or eps
#' @param width the figure width in inches
#' @param height the figure height in inches
#' @param res the figure resolution in ppi (pixels per inch), only used
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' png
#' @export
plotdev <- function(file="test", type="png", width=5, height=5, res=100) {
  file <- paste(file, type, sep=".")
  if ( type == "png" )
    png(file, width=width, height=height, units="in", res=res)
  if ( type == "eps" )
    postscript(file, width=width, height=height, paper="special")
  if ( type == "pdf" )
    pdf(file, width=width, height=height)
}


## return R colors as RGB, to allow setting alpha
getRGB <- function(n) {
    cols <- col2rgb(1:n)
    apply(cols,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))
}


#' \code{\link{showSpectrum}}:
#' shows the color spectrum of visible light along wavelengths in nm.
#' See \code{\link{wavelength2RGB}} for details.
#' @param wavelengths vector of wavelengths to be plotted
#' @param alpha the alpha factor (opacity) for plot symbols, min=1, max=255
#' @param pch the plot symbol type, see ?par("pch")
#' @param cex the plot symbol size, see ?par("cex")
#' @param ylab ylab, see ?plot
#' @param xlab xlab, see ?plot
#' @param main main, see ?plot
#' @param ... further arguments to plot
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @details Plots the color spectrum over the selected wavelengths in nm.
#' The utility functions  \code{\link{plotWavelength}} and
#' \code{\link{findWavelength}} can be used to plot and interactively (click
#' on plot) find wavelengths. The function \code{\link{wavelength2RGB}} can
#' then be used to convert this wavelength to RGB colors.
#' @seealso \code{\link{wavelength2RGB}}, \code{\link{findWavelength}},
#' \code{\link{plotWavelength}}
#' @export
showSpectrum <- function(wavelengths=380:780, alpha=99, pch=19, cex=3,
                         ylab="approximate color",
                         main="use findWavelength(3)",
                         xlab="wavelength, nm", ...) {
    cols <- wavelength2RGB(wavelengths)
    plot(wavelengths,rep(1,length(wavelengths)),
         axes=FALSE,ylab=ylab,xlab=xlab,main=main,
         col=paste(cols,alpha,sep=""),ylim=c(.9,1.1),pch=pch,cex=cex,...)
    axis(1,at=pretty(wavelengths))
}
#' \code{\link{plotWavelength}}:
#' a color selector for \code{\link{showSpectrum}}; draws a vertical line,
#' the wavelength, and a filled circle at a given wavelength in nm.
#' @param x wavelength in nm
#' @param y position of the text, from 0.9 to 1.1
#' @param ... further arguments to ?text
#' @param ych y position of the symbol
#' @param pch the plot symbol type, see ?par("pch")
#' @param cex the plot symbol size, see ?par("cex")
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{showSpectrum}}, \code{\link{findWavelength}},
#' \code{\link{wavelength2RGB}}
#' @export
plotWavelength <- function(x=534, y=1.09, ych=0.95, cex=5, pch=19, ...) {
    col <- wavelength2RGB(x)
    abline(v=x,col=col)
    text(x, y, x,col=col, ...)
    points(x, ych, col=col, cex=cex, pch=pch)
    axis(4)
}

#' \code{\link{findWavelength}}:
#' a color selector for \code{\link{showSpectrum}}; expects the user to
#' click on the spectrum, then draws a vertical line at the clicked
#' wavelength using \code{\link{plotWavelength}} and records the wavelength
#' in nm.
#' @param n number of clicks/wavelengths to record
#' @param ... further arguments to \code{\link{plotWavelength}}, for
#' selecting plot symbol and text size, and positions
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{showSpectrum}}, \code{\link{plotWavelength}}, \code{\link{wavelength2RGB}}
#' @export
findWavelength <- function(n=1, ...) {
    cat(paste("PLEASE CLICK ON THE PLOT ... "))
    rgb <- matrix(NA,ncol=2,nrow=n)
    colnames(rgb) <- c("nm","RGB")
    for ( i in 1:n ) {
        xy <- locator(1)
        rgb[i,] <-  c(round(xy$x,1), wavelength2RGB(xy$x))
        plotWavelength(round(xy$x), ych=xy$y, ...)
        cat(paste("\n\t", round(xy$x,1), "nm, in RGB:", rgb[i,2],"\n",
                  ifelse(i==n,"",paste(n-i, "LEFT, CLICK AGAIN ..."))))
    }
    cat("\n")
    rgb
}

#' \code{\link{wavelength2RGB}}:
#' converts wavelength in nm (visible light: 380:780 nm) to RGB.
#' @details implemented following the java code posted at
#' http://stackoverflow.com/questions/1472514/convert-light-frequency-to-rgb .
#' Also see http://www.fourmilab.ch/documents/specrend/ for original code and
#' why not all wavelengths can be converted to RGB.
#' @param wavelength the wavelength in nm to convert to RGB
#' @examples
#' wavelengths <- seq(380,780,1)
#' cols <- sapply(wavelengths, wavelength2RGB)
#' bars <- rep(1,length(wavelengths)); names(bars) <- wavelengths
#' barplot(bars,border=cols,col=cols,las=2)
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{showSpectrum}}
#' @export
wavelength2RGB <- function(wavelength)
    cols <- sapply(wavelength, lambda2RGB)
## translated from
## http://stackoverflow.com/questions/1472514/convert-light-frequency-to-rgb
lambda2RGB <- function(wavelength) {

    ## intensity scaling near visibility
    gamma <- 0.80 
    intensityMax <- 255

    Red <- 0.0
    Green <- 0.0
    Blue <- 0.0

    if ( (wavelength >= 380) & (wavelength<440) ) {
        Red <- -(wavelength - 440) / (440 - 380)
        Green <- 0.0
        Blue <- 1.0
    } else if ( (wavelength >= 440) & (wavelength<490) ) {
        Red <- 0.0
        Green <- (wavelength - 440) / (490 - 440)
        Blue <- 1.0
    }else if ( (wavelength >= 490) & (wavelength<510) ) {
        Red <- 0.0;
        Green <- 1.0;
        Blue <- -(wavelength - 510) / (510 - 490);
    } else if ( (wavelength >= 510) & (wavelength<580) ) {
        Red <- (wavelength - 510) / (580 - 510)
        Green <- 1.0
        Blue <- 0.0
    } else if ((wavelength >= 580) & (wavelength<645) ) {
        Red <- 1.0
        Green <- -(wavelength - 645) / (645 - 580)
        Blue <- 0.0
    } else if ( (wavelength >= 645) & (wavelength<781) ) {
        Red <- 1.0
        Green <- 0.0
        Blue <- 0.0
    }

    ## intensity scaled near the vision limits
    factor <- 0.0
    if ( (wavelength >= 380) & (wavelength<420) ) {
        factor <- 0.3 + 0.7*(wavelength - 380) / (420 - 380)
    } else if ( (wavelength >= 420) & (wavelength<701) ) { 
        factor <- 1.0
    } else if ( (wavelength >= 701) & (wavelength<781) ) {
        factor <- 0.3 + 0.7*(780 - wavelength) / (780 - 700)
    }
    
    rgb <- rep(NA,3)
    ## Don't want 0^x = 1 for x <> 0
    rgb[1] <- ifelse(Red==0,  0, round(intensityMax * (Red   * factor)^gamma))
    rgb[2] <- ifelse(Green==0,0, round(intensityMax * (Green * factor)^gamma))
    rgb[3] <- ifelse(Blue==0, 0, round(intensityMax * (Blue  * factor)^gamma))
    
    rgb(rgb[1],rgb[2],rgb[3],maxColorValue=intensityMax)
}


## trim leading and trailing white-space from parsed strings
trim <- function(str)
  gsub(" *$", "", gsub("^ *", "", str) )


### STATS

## moving average
ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

# calculate 95% confidence intervals for the given
# data vector using a t-distribution
ci95 <- function(data,na.rm=FALSE) {
    if ( na.rm ) data <- data[!is.na(data)]
    n <- length(data)
    if ( n<2 ) return(NA)
    error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
    return(error)
}


## calculcates standard error
se <- function(data,na.rm=TRUE) {
    if ( na.rm ) data <- data[!is.na(data)]
    n <- length(data)
    if ( n<1 ) return(NA)
    error <- sd(data)/sqrt(n)
    error
}



#' \code{\link{prettyData}} : set colors, rename, order or filter the data set
#' @param data \code{\link{platexpress}} data, see \code{\link{readPlateData}}
#' @param dids a vector of data IDs, data will be filtered and sorted by this list; if the vector is named the IDs will be replaced by these names
#' @param colors a vector of plot colors as RGB strings, optionally already named by dataIDs 
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{readPlateData}}
#' @export
prettyData <- function(data, dids, colors) {

    ## dids: data sorting, if missing, take all
    if ( missing(dids) )
        dids <- data$dataIDs

    ## get mids (global x-axis data)
    mids <- data$mids
    
    ## re-order: mids first and IDs last
    data$dataIDs <- dids
    ## store colors
    origcols <- data$colors
    data <- data[c(match(mids,names(data)),
                   match("mids",names(data)),
                   match(dids,names(data)),
                   match("dataIDs",names(data)))]
    data$colors <- origcols

    ## rename!
    if ( !is.null(names(dids)) ) {
        names(data)[match(dids,names(data))] <- names(dids)
        data$dataIDs[match(dids,data$dataIDs)] <- names(dids)
    }
    
    ## generate colors, if none are present
    if ( missing(colors) & !"colors"%in%names(data) ) {
        data$colors <- getColors(data$dataIDs)
    }
    ## set new or replace existing 
    if ( !missing(colors) ) {
        if ( is.null(names(colors)) )
            names(colors) <- data$dataIDs
        data$colors[names(colors)] <- colors
    }
    data
}

## TODO: implement other types? 
getColors <- function(ptypes,type="R") {
    ## colors
    if ( type=="R" )
        pcols <- getRGB(length(ptypes))
    if ( type=="rainbow" )
        pcols <- sub("FF$","",rainbow(length(ptypes)))
    names(pcols) <- ptypes
    pcols
}

## add data
#' \code{\link{addData}} : add a data set, e.g., calculated ratios
#' @param data the current platexpress data set
#' @param ID the ID of the new data 
#' @param dat the new data, must be a matrix akin to other data in the set
#' @param col plot color for the new data, an RGB string w/o alpha suffix
#' @param processing optional processing note
#' @param replace replace existing data, default: FALSE
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{readPlateData}}
#' @export
addData <- function(data, ID, dat, col, processing,replace=FALSE) {
    ## replace existing
    if ( ID%in%data$dataIDs ) {
        if ( replace ) {
            data[[ID]]$data <- dat
            if ( !missing(col) )
                if ( "colors" %in% names(data) )
                    data$colors[ID] <- col
            if ( missing(processing) )
                processing <- paste("replaced",date())
            data[[ID]]$processing <- processing
        } else
            stop("\"",ID, "\" already present, use replace=TRUE to replace")
    } else {
        ## add new data
        data$dataIDs <- c(data$dataIDs, ID) # add to ID list
        if ( "colors" %in% names(data) ) { # add color
            if ( missing(col) ) 
                col <- getRGB(length(data$colors)+1)
            data$colors <- c(data$colors, col[length(col)])
            names(data$colors) <-
                c(names(data$colors)[2:length(data$colors)-1],ID)
        } else if ( !missing(col) ) {
            data$colors <- getColors(data$dataIDs)
            data$colors[ID] <- col
        }
        
        if ( missing(processing) ) # add date as processing note
            processing <- date()
        data <- append(data, list(list(data=dat, processing=processing)))
        names(data) <- c(names(data)[2:length(data)-1],ID)
    }
    data
}


#' \code{\link{rmData}} : remove a data set
#' @param data a \code{\link{platexpress}} data set
#' @param ID a vector of IDs of the data to be removed
#' @return Returns the new \\code{\link{platexpress}} data, with 
#' data set <\code{ID}> removed
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{addData}}
#' @export
rmData <- function(data, ID) {
  
  for ( id in ID ) {
    if ( !id%in%data$dataIDs )
      warning("\"",id,"\" not found")
    data$dataIDs <- data$dataIDs[!data$dataIDs%in%id] # rm ID
    if ( "colors" %in% names(data) ) # rm color
        data$colors <- data$colors[!names(data$colors)%in%id]
    data <- data[-which(names(data)%in%id)] # rm data
  }
  data
}


#' \code{\link{getData}} : get a specific data set, returns a data matrix
#' @param data the current platexpress data set
#' @param ID the ID of the data to be obtained
#' @param type the type of the data to be returned (default "data"),
#' the data list also contains original values for some of the processing
#' steps
#' @export
getData <- function(data, ID, type="data") {
    data[[ID]][[type]] # just return the current data or old versions
}


#' \code{\link{shiftData}} : shift x-axis by a lag-phase to align growth curves
#' @param data \code{\link{platexpress}} data set
#' @param lag a named vector given the lag-phase to be removed; the names
#' correspond to the wells in data
#' @param mid the x-axis to be used, defaults to the first available
#' (usually "Time")
#' @export
shiftData <- function(data, lag, mid) {

    ## NOTE: not optional, since data don't fit anymore
    ## if only some are shifted for a given well list (names of lag)
    dids <- data$dataIDs

    if ( missing(mid) )
        mid <- data$mids[1]
    xdat <- data[[mid]]
    
    for ( i in 1:length(lag) ) {
        well <- names(lag)[i]
        idx <- which(xdat >= lag[i])[1]
        end <- length(xdat)
        for ( did in dids ) {
            data[[did]]$data[,well] <-
                          c(data[[did]]$data[idx:end, well],
                            rep(NA,idx-1))
        }
        data[[did]]$processing <- c(data[[did]]$processing,
                                    paste("well",well,"shifted by lag", lag[i]))
    }
    data
}

## TODO: cut data either by time, or by a chosen data set
#' \code{\link{cutData}} : cut data in a range of the x-axis
#' @param data \code{\link{platexpress}} data, see \code{link{readPlateData}}
#' @param rng a single value or a data range
#' @param mid ID of the x-axis data to be used for cutting
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @details Cuts the passed \code{\link{platexpress}} data to ranges of
#' of the x-axis (time or other, see \code{data$mids}). If paramter \code{rng}
#' is a single value, data for the closest x value will be returned
#' if rng is a single value
#' @export
cutData <- function(data, rng, mid) {

    if ( missing(mid) )
        mid <- data$mids[1]
    xdat <- data[[mid]]
    if ( length(rng)==2 )
        filter <- xdat >= rng[1] & xdat <= rng[2]
    else if ( length(rng)==1 )
        filter <- abs(rng-xdat) == min(abs(rng-xdat))
    for ( did in data$dataIDs ) {
        data[[did]]$data <- data[[did]]$data[filter,,drop=FALSE]
        data[[did]]$processing <- c(data[[did]]$processing,
                                    paste("cut at", paste(rng,collapse="-")))
    }
    data[[mid]] <- xdat[filter]
    data
}

#' dose-response box-plots
#'
#' plots a continuous value for each well as a function
#' of an \code{amount} of a \code{substance} and make boxplots
#' for all replicates at a given amount
#' @param map a plate layout map with columns \code{amount} and some 
#' calculated value in column \code{val}, 
#' eg. results from \code{\link[grofit:grofit]{grofit}}
#' @param wells a list of well IDs to be used in the plot
#' @param val the name of a column in \code{map} containing numeric values
#' that should be plotted on y-axis
#' @param amount the name of a column in \code{map} providing numeric values
#' that should be plotted on x-axis, typically a substance added to wells
#' in multiple replicates; the default value "amount" is automatically
#' generated from appropriate plate layout files by \code{link{readPlateMap}} 
#' @param substance the name of a column in \code{map} providing the names of 
#' the substance in \code{amount}, used as x-axis label; if it doesn't 
#' correspond to a column, its value is directly used as label
#' @param color the name of a column in \code{map} providing colors
#' for each well; NOTE, that wells with the same \code{amount} should have
#' the same color, only the first color for a given \code{amount} is used;
#' Colors are automatically assigned from a color ramp mapped to the
#' numerical range of \code{amount} by \code{link{readPlateMap}} 
#' @param na.y value to be used to plot replicates with \code{NA} in 
#' column \code{val}; set to NA to supress plotting
#' @param ylim limits of the y-axis
#' @param xnum use numerical x-axis instead of default categorical
#' @export
doseResponse.box <- function(map, wells, val, amount="amount", substance="substance", 
                             color="color", na.y=0, ylim, xnum=FALSE) {

  if( missing(wells) ) wells <- map[,"well"]
  
    wells <- match(wells,map[,"well"])
  
    y <- map[wells,val]
    x <- map[wells,amount]
    
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
    boxplot(y1~x1, border=cl,
            xlab=subid,ylab=val,ylim=ylim, #xlim=xlim,
            na.action="na.pass",xaxt=xaxt)
    stripchart(y1~x1,add=T, vertical=TRUE,col=cl,
               method="jitter", pch=1,cex=1,
               na.action="na.pass")
    if ( xnum ) # numerical xaxis
      axis(1)    
    
    if ( xnum )
      x <- factor(x, levels=do.call(seq, as.list(range(x))))
    
    nay <- y
    nay[is.na(y)] <- na.y
    nay[!is.na(y)] <- NA
    stripchart(nay~x,add=T, vertical=TRUE,col="red",
               method="jitter", pch=4,cex=1/2,
               na.action="na.pass")
}

#' dose-response plots with error bars for replicates
#'
#' plots a continuous value for each well as a function
#' of an \code{amount} of a \code{substance} and generates error bars
#' for all replicates at a given amount; returns the plotted mean values
#' and error bar ranges.
#' @param map a plate layout map with columns \code{amount} and some 
#' calculated value in column \code{val}, 
#' eg. results from \code{\link[grofit:grofit]{grofit}}
#' @param wells a list of well IDs to be used in the plot
#' @param val the name of a column in \code{map} containing numeric values
#' that should be plotted on y-axis
#' @param amount the name of a column in \code{map} providing numeric values
#' that should be plotted on x-axis, typically a substance added to wells
#' in multiple replicates; the default value "amount" is automatically
#' generated from appropriate plate layout files by \code{link{readPlateMap}} 
#' @param substance the name of a column in \code{map} providing the names of 
#' the substance in \code{amount}, used as x-axis label; if it doesn't 
#' correspond to a column, its value is directly used as label
#' @param col either a valid color representation or the name of a column in 
#' \code{map} providing colors for each well; NOTE, that wells with the same 
#' \code{amount} should have the same color, only the first color for a given 
#' \code{amount} is used; such colors are automatically assigned from a color ramp 
#' mapped to the numerical range of \code{amount} by \code{\link{readPlateMap}},
#' see \code{\link{amountColors}} 
#' @param pch pch of the mean value points
#' @param bartype type of the error bar range, "range" for the full range
#' of the data, "ci95" for the 95% confidence interval, "sd" for the standard
#' deviation, "se" for the standard error
#' @param barl length of the horizontal lines of error bars (argument \code{length}
#' of function \code{\link{arrows}})
#' @param all logical indicating whether to plot all data points additionally to 
#' error bars
#' @param line logical indicating whether to plot a line connecting the mean
#' values
#' @param na.y value to be used to plot replicates with \code{NA} in 
#' column \code{val}; set to NA to supress plotting
#' @param add add data to existing plot
#' @param ylim limits of the y-axis
#' @param xlim limits of the x-axis
#' @param ylab alternative label for the y-axis, default is to use argument \code{val}
#' @param xlab alternative label for the x-axis, default is to use argument 
#' \code{substance}, or if this is a column in \code{map}, the substance
#' indicated there
#' @param ... arguments passed on to the main setup \code{\link{plot}}
#' @export
doseResponse <- function(map, wells, val, amount="amount", substance="substance", 
                         col="color", pch=1, bartype="range", barl=.05, 
                         all=FALSE, line=TRUE, na.y=0, add=FALSE, ylim, xlim, ylab, xlab, ...) {
  
  if( missing(wells) ) wells <- map[,"well"]
  
  wells <- match(wells,map[,"well"])
  
  y <- map[wells,val]
  x <- map[wells,amount]

  ## get unique colors for unique sorted x, as it will appear in
  ## boxplots
  ## UGLY, TODO: less ugly?
  if ( col %in% colnames(map) ) {
    cl <- map[wells,col]
    cl <- cl[order(x)]
    cl <- cl[which(!duplicated(sort(x)))]
    names(cl) <- sort(x[!duplicated(x)])
    linecol <- 1
  } else { 
    cl <- rep(col,length(unique(x)))
    linecol <- col
  }
  
  if ( substance %in%  colnames(map) ) 
    subid <- map[wells,substance][1]
  else subid <- substance
  if ( missing(xlab) )
    xlab <- subid
  if ( missing(ylab) )
    ylab <- val
  
  y1 <- y
  x1 <- x
    
  y1 <- y1[!is.na(y)]
  x1 <- x1[!is.na(y)]
  

  ## calculate ranges for duplicates at each x
  xlevels <- sort(unique(x1))
  ymat <- matrix(NA,nrow=length(xlevels),ncol=4)
  for ( i in 1:length(xlevels) ) {
    xi <- xlevels[i]
    yi <- y1[x1==xi]
    ymat[i,2] <- mean(yi,na.rm=TRUE)
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
  colnames(ymat) <- c(amount,val,bartype,"")
  if ( missing(ylim) )
    ylim <- range(c(ymat[,2:4],ifelse(all,y,NA)),na.rm=TRUE)
  if ( missing(xlim) ) {
    if ( !is.na(na.y))
      xlim <- range(x,na.rm=TRUE)
    else xlim <- range(x1,na.rm=TRUE)
  }
  
  ## plot
  if ( !add )
    plot(x1, y1, col=NA, xlab=subid, ylab=ylab, ylim=ylim, xlim=xlim, ...)
  if ( line )
    lines(ymat[,1],ymat[,2],col=linecol)
  for ( i in 1:nrow(ymat) ) {
    if ( !line ) 
      points(ymat[i,1],ymat[i,2], col=cl[i], pch=pch)
    arrows(ymat[i,1], ymat[i,3], ymat[i,1], ymat[i,4], length=barl, angle=90, code=3, col=cl[i])
    if ( all )
      points(x1[x1==ymat[i,1]], y1[x1==ymat[i,1]], pch=20, cex=.5, , col=cl[i])
  }
  if ( !is.na(na.y) & sum(is.na(y)) )
    points(x[is.na(y)], rep(na.y,sum(is.na(y))), col="red", pch=4, cex=.5)
  invisible(ymat)
}

#' returns data for group of wells in a given range of the x-axis
#' @param data \code{\link{platexpress}} data, see \code{\link{readPlateData}}
#' @param groups a grouping of wells, see \code{\link{getGroups}}
#' @param rng x-axis range or a single point, for the latter the closest
#' point will be selected
#' @param mid the x-axis ID, if multiple x-axes are present
#' @param did the y-axis data to be grouped
#' @param plot if TRUE a box-plot or bar-plot will be plotted
#' @param type either "box" (default) or "bar" for box-plot or bar-plot
#' @param etype type of statistics to be used for error bars in the bar-plot,
#' either "ci" (default) for the 95%-confidence interval or "se" for
#' the standard error
#' @return Returns an annotated \code{data.frame} of the values, with well
#' and group IDs in the first two columns. The values in the third column are
#' raw values (if argument \code{rng} was a single point
#' on the x-axis) or mean values (if argument\code{rng} was a range).
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
boxData <- function(data, rng, groups, mid, did="OD", plot=TRUE, type="box", etype="ci") {

    if ( missing(mid) )
        mid <- data$mids[1]
    ## cut data to selected range (closest point if length(rng)==1)
    cdat <- cutData(data, rng, mid)

    ## get the actual point if only one points was chosen
    if ( length(rng)==1) rng <- signif(unique(range(cdat[[mid]])),4)

    bdat <- rep(list(NA),length(groups))
    names(bdat) <- names(groups)
    for ( sg in 1:length(groups) ) 
        bdat[[sg]] <- cdat[[did]]$data[,groups[[sg]],drop=FALSE]
    ## get means for all groups
    pdat <- lapply(bdat, function(x) apply(x,2,mean,na.rm=TRUE))

    if ( plot ) {
        #par(mai=c(1,.75,.1,.1))
        if ( type=="box" )
            boxplot(pdat,ylab=did,las=2)
        else if ( type=="bar" ) {
            mn <- unlist(lapply(pdat, mean,na.rm=TRUE))
            if ( etype=="ci" ) # 95% confidence interval
                ci <- unlist(lapply(pdat, ci95,na.rm=TRUE))
            else if ( etype=="se" ) { # standard error sd/n^2
                sd <- unlist(lapply(pdat, sd, na.rm=TRUE))
                n2 <- sqrt(unlist(lapply(pdat, length)))
                ci <-  sd / n2
            }
            
            x <- barplot(mn,ylim=c(0,max(mn+ci,na.rm=TRUE)),ylab=did,las=2)
            arrows(x0=x,x1=x,y0=mn-ci,y1=mn+ci,code=3,angle=90,
                   length=.05,lwd=1.5)            
        }
        legend("topright",paste("at",mid, "=",paste(rng,collapse="-")),
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
    colnames(tmp)[3] <- paste(did,paste(rng,collapse="-"),sep="@")
    bdat <- tmp
    rownames(bdat) <- NULL

    result <- bdat
}



#' \code{\link{skipWells}} rm wells from both data, plate maps and groupings
#' @param data data structures from \code{\link{platexpress}}; either data
#' (\code{\link{readPlateData}}), a plate layout map
#' (\code{\link{readPlateMap}}) or a well grouping (\code{\link{getGroups}})
#' @param skip a list of strings identifiying the wells to be skipped,
#' e.g. "B3" to skip the well in row B/column 3
#' @details Removes specific wells from \code{\link{platexpress}} data,
#' groupings and plate layout maps. If the first argument is
#' \code{\link{platexpress}} data, the specified wells will be set to NA.
#' If the first argument is a \code{platexpress} well grouping or plate
#' layout map, the specified wells will be removed.
#' @examples
#' data(ap12)
#' raw <- skipWells(ap12data, skip="A9") # rm data from well "A9"
#' plate <- skipWells(ap12plate, skip="A9") # rm well "A9" from the plate layout
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
skipWells <- function(data, skip) {

    if ( "dataIDs" %in% names(data) ) ## rm from data
      for ( id in data$dataIDs ) {
          wells <- colnames(data[[id]]$data)
          data[[id]]$data <- data[[id]]$data[,-which(wells%in%skip),drop=FALSE]
      }
    else if ( !is.null(dim(data)) ) ## rm from plate layout map
        ##data[match(skip,data[,"well"]),2:ncol(data)] <- NA
        data <- data[-match(skip,data[,"well"]),]
    else ## rm from grouping
      for ( g in 1:length(data) )
        data[[g]] <- data[[g]][!data[[g]]%in%skip]
    data
}

#' get a filtered list of wells
#' @param plate the plate layout map, see \code{\link{readPlateMap}}
#' @param blanks if set to \code{TRUE} (default) the argument \code{values}
#' is optional, and only blank values will be returned.
#' @param values a named list of strings, key:value pairs, where the
#' list names (keys) correspond to column names in the plate layout map
#' and the values are entries in these columns
#' @return Returns the list of wells according to argument \code{values},
#' or a list of blanks if \code{blanks==TRUE}
#' @export
getWells <- function(plate, blanks=FALSE, values) {
    if ( !missing(values) ) {
        for ( i in 1:length(values) ) {
            key <- names(values)[i]
            val <- values[[i]]
            plate <- plate[as.character(plate[,key])%in%val,]
        }
    } else blanks <- TRUE # return blank as default
    if ( blanks ) 
      plate <- plate[plate[,"blank"]==TRUE,]
    res <- plate[,"well"]
    return(as.character(res[!is.na(res)]))
}

#' \code{\link{correctBlanks}} correct for blanks
#' @param data the \code{\link{platexpress}} data list to be blank-corrected
#' @param plate the plate layout where column "blanks" indicates which wells
#' are to be treated as blanks
#' @param dids IDs of the data which should be blank-corrected, all will be
#' blanked if missing
#' @param by a list of column IDs of the plate layout; separate blank
#' correction will be attempted for groups in these columns; each group
#' must have at least one specified blank associated
#' @param type TODO
#' @param mid ID of the x-axis data to be used, if blanked along x-axis, set
#' by \code{mbins}>1
#' @param mbins the number of bins the x-axis is to be divided, if blanked
#' along the x-axis, see \code{mid}
#' @param max.mid the maximal x-axis value where blanks should be used
#' @param verb issued progress messages and info
#' @seealso \code{\link{adjustBase}}
#' @examples
#' data(ap12)
#' data <- correctBlanks(data=ap12data, plate=ap12plate, by="strain")
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
correctBlanks <- function(data, plate, type="median", by, dids, mid, max.mid, mbins=1, verb=TRUE) {

### TODO: correct by time point, eg. for fluorescence in ecoli_rfp_iptg_20160908

    if ( missing(mid) )
        mid <- data$mids[1]
    
    ## start new data list
    corr <- data
    time <- data[[mid]] ## TODO: take from data mids
  
    ## reduce matrix to requested data
    data <- data[data$dataIDs]
    ptypes <- names(data)
    if ( !missing(dids) ) # only use requested data 
      ptypes <- ptypes[ptypes%in%dids]
    if ( length(ptypes)==0 )
        stop("no data to blank")
    else if ( verb )
      cat(paste("blanking", paste(ptypes,collapse=";"),"\n"))

    ## get present wells
    pwells <- unique(c(sapply(corr$dataIDs,
                              function(id) colnames(data[[id]]$data))))
    
    ## blank wells
    blanks <- plate[,"blank"]
    blanks[is.na(blanks)] <- FALSE
    
    ## get different types to be blanked separately
    ## convert platemap to char
    types <- rep(TRUE,nrow(plate))
    if ( !missing(by) ) {
        lpl <- matrix(unlist(lapply(plate,function(x) as.character(x))),
                      ncol=ncol(plate))
        colnames(lpl) <- colnames(plate)
        ## collapse requested combinations into new type
        types <- rep("",nrow(plate))
        for ( b in by ) 
          types <- paste(types,lpl[,b],sep="_")
    }
    btypes <- unique(types)

    ## CORRECT BY TIME
    ## sometimes blanks show trends, e.g., higher fluorescence
    ## in the beginning

    
    ## blank each type separately
    for ( i in 1:length(btypes) ) {
        btyp <- btypes[i]
        ## data and blank wells of the correct type
        dwells <- as.character(plate[types==btyp & !blanks,"well"])
        bwells <- as.character(plate[types==btyp &  blanks,"well"])
        ## filter for actually present wells
        dwells <- dwells[dwells%in%pwells]
        bwells <- bwells[bwells%in%pwells]
        
        if ( verb )
          cat(paste("blanking", btyp, ":", length(dwells), "wells, using",
                    length(bwells), "blank wells\n"))
        for ( k in 1:length(ptypes) ) {
            ptyp <- ptypes[k]
            dat <- data[[ptyp]]$data
            ## TODO: do this in time bins
            
            timebins <- unique(c(seq(1,nrow(dat),nrow(dat)/mbins),nrow(dat)))
            nbin <- length(timebins)
            timebins <- cbind(ceiling(timebins[1:(nbin-1)]),
                              floor(timebins[2:nbin]))
            if ( verb ) cat(paste(ptyp, "\n"))
            ## calculate and subtract blanks for time bins (default: all)
            for ( t in 1:nrow(timebins) ) {
                bin <- timebins[t,1]:timebins[t,2]
                #if ( nrow(timebins)>1 )
                if ( verb )
                  cat(paste("\ttime bin:",t,timebins[t,1],"-",timebins[t,2]))

                ## cut maximal time for blanking
                bbin <- bin
                if ( !missing(max.mid) ) {
                    if ( verb ) 
                      cat(paste("\tskipping",sum(time[bbin]>max.mid),
                                "bins at",max.mid,"\n"))
                    bbin <- bbin[time[bbin]<=max.mid]
                }

                ## TODO: this should only happen if time is
                if ( length(bbin)==0 ) {
                    #cat(paste("skipping time bin", t, "at", max.mid, "\n"))
                    warning("no blank data at time bin", t, "at", max.mid)
                    blank <- 0
                    #next # TODO: warning?
                } else {
                    ## calculate blank!
                    if ( type=="median" )
                        blank <- median(c(dat[bbin,bwells]),na.rm=TRUE)
                    else if ( type=="mean" )
                        blank <- mean(c(dat[bbin,bwells]),na.rm=TRUE)
                    else if ( type=="ci95" ) # TODO: why is ci95 subtracted?
                        blank <- mean(c(dat[bbin,bwells]),na.rm=TRUE) -  ci95(c(dat[bbin,bwells]),na.rm=TRUE)
                }
                ## subtract blank
                corr[[ptyp]]$data[bin,c(dwells,bwells)] <-
                               dat[bin,c(dwells,bwells)] - blank
                if ( verb ) cat(paste("\tblank:",blank,"\n"))
                ##cat(paste("\tdata wells:",paste(dwells,collapse=";"),"\n",
                ##          "\tblank wells:",paste(bwells,collapse=";"),"\n"))
            }
            corr[[ptyp]]$processing <- c(corr[[ptyp]]$processing,
                                         paste("blank-corrected by",btyp))
        }
    }
    corr
}

#' \code{\link{adjustBase}} adjust data to a minimal base
#' @details Adjusts data to a new mininum, this is useful for adjustment
#' of negative values after blank corrections
#' @param data \code{\link{platexpress}} data, see \code{\link{readPlateData}}
#' @param dids vector of ID strings for which base correction should be
#' executed
#' @param base the new minimum for the data, default is 0, but it could
#' e.g. be the OD used for inoculation
#' @param xlim min and max row number of the data to be adjusted
#' @param add.fraction a fraction of the whole data range, added to base
#' @param each add base for each well separately!
#' @param verb print messages if true
#' @return Returns `data' where all data sets or only those selected by option
#' dids where raised to a minimum level in 
#' @seealso \code{\link{correctBlanks}}
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
adjustBase <- function(data, base=0, dids, add.fraction, xlim, each=FALSE, verb=TRUE) {

    if ( missing(dids) ) # only use requested data 
        dids <- data$dataIDs # use all
    
    for ( did in dids ) {

        ## each well separately?
        if ( each )
            bins <- as.list(1:ncol(data[[did]]$data))
        else
            bins <- list(1:ncol(data[[did]]$data))

        for ( bin in bins ) {
            
            dat <- data[[did]]$data[,bin,drop=FALSE]

            if ( missing(xlim) )
                xlim <- c(1,nrow(dat))
            xrng <- xlim[1]:xlim[2]

            if ( verb )
              cat(paste(paste(bin,collapse=";"),
                        "adding", min(dat[xrng,],na.rm=TRUE), "\n"))

            ## TODO: smarter? only if any value is <0?
            dat <- dat - min(dat[xrng,],na.rm=TRUE) + base
            
            if ( !missing(add.fraction) ) 
                dat <- dat + diff(range(dat[xrng,],na.rm=TRUE))*add.fraction
            
            data[[did]]$data[,bin] <- dat
        }
        data[[did]]$processing <- c(data[[did]]$processing,
                                    paste("corrected to base",base))
    }
    data
}

## helper function to calculate an average value
## for a variable given in several list items,
## used for average (master) time and temperatures
## in interpolatePlateTimes()
## TODO: allow different lengths?
listAverage <- function(lst, id) {

    ## reduce to entries that have <id>
    lst <- lst[unlist(lapply(lst, function(x) id%in%names(x)))]
    if ( length(lst)==0 ) return(NULL)
    ## collect values with the same ID for different data sets
    vals <- lapply(lst, function(x) x[[id]])
    ## check length of lists and find missing time-points
    if ( any(diff(unlist(lapply(vals, length)))!=0) ) {
        stop("data have different lengths: ", id)
        ## TODO: simple solution: take longest vector
        
        ## TODO: find point of discrepancy
        all <- unlist(vals) # concat all values
        # get lengths
        len <- unlist(sapply(1:length(vals),
                             function(x) rep(x, length(vals[[x]]))))
        ## for each set there should be one measurement per time point
        idx <- which(diff(len[order(all)])==0)
        mn <- min(unlist(lapply(vals, length)))
        tmp <- lapply(vals, function(x) x[1:mn])
        tmp <- matrix(unlist(tmp), ncol = length(tmp), byrow = FALSE)
    }
    vals <- matrix(unlist(vals), ncol = length(vals), byrow = FALSE)
    ## take average of each time-point
    avg <- apply(vals,1,mean)

    ## if this is used form "time" check difference between data
    if ( id %in% paste(c("time","Time"),rep(c("","s"),each=2),sep="") ) {
        tsd <- apply(vals,1,sd) ## standard deviation within timepoint
        tdf <- apply(vals,2,diff) ## difference between timepoints
        if ( max(tsd) > 0.01*median(c(tdf)) ) {
            td <- max(tsd)/median(c(tdf))
            warning(paste(id, ": max. SD within timepoint is",
                          round(td,3)*100, "% of median difference between",
                          "time points.\n"))
        }
    }    
    avg
}

#' \code{\link{interpolatePlateTimes}} interpolate all data to an average
#' master time: calculates average time for each measurement point
#' and interpolates all values to this time; this is also used for
#' well temperatures
#' @param data TODO
#' @param verb TODO
#' @param xid TODO
#' @return returns a copy of the full data list with a master time and
#' temperature added at the top level
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
interpolatePlateTimes <- function(data, verb=TRUE, xid) {


    ## catch single data case
    if ( length(data$dataIDs)==1 ) {
        data$Time <- data[[data$dataIDs[1]]]$time
        data$mids <- c(data$mids, "Time")
        return(data)
    }
    
    if ( verb )
      cat(paste("Interpolating all data to a single master time.\n"))
    
    ## 0) TODO: check whether all data items have the same
    ## number of time-points, and cut end - this can stem from
    ## termination of the measurement before programmed end
    ## 0a: first, cut data at non-increasing time steps, 00:00:00 in Synergy
    ## 0b: check length
    
    ## 1) calculate average (MASTER) time
    mtime <- listAverage(data, "time")
    ## TODO: add back temperature
    #mtemp <- listAverage(data, "temp") 
    
    ## 2) interpolate all data to MASTER time
    for ( id in data$dataIDs ) {
        data[[id]]$orig <- data[[id]]$data
        mdat <- data[[id]]$data
        for ( j in 1:ncol(data[[id]]$data) ) {
            x <- data[[id]]$time
            y <- data[[id]]$data[,j]
            ## interpolate data, NOTE that rule=2 will fill the end points
            #mdat[,j] <- approx(x=x,y=y,xout=mtime,rule=2)$y
            ## TODO: is this OK/betterer then simple approx?
            mdat[,j] <- spline(x=x,y=y,xout=mtime,method="fmm")$y
                
        }
        ## replace data
        data[[id]]$data <- mdat
        ## indicate interpolation
        data[[id]]$processing <- c(data[[id]]$processing,
                                   "interpolated")
    }
    

    ## add master time and calculate master temperature as well
    ## TODO: separate temperature:time pairs could be used to
    ## interpolate temperatures to master time as well
    data$Time <- mtime
    data$mids <- c(data$mids, "Time")
    #data$Temperature <- mtemp
    data
}

#' interpolate one dataset to common points of another
#' data set, e.g., fluorescence to OD
#' @param data the list of measurement data as provided by
#' \code{\link{readPlateData}}
#' @param xid ID of the data set which should serve as the new x-axis
#' all data will be interpolated to equi-spaced points along the range
#' of measured values
#' @param dids restrict interpolation to these data IDs
#' @param n specify the number of interpolation points, if missing the
#' original number of rows will be used
#' @param xout specify the interpolation points directly
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
interpolatePlateData <- function(data, xid, dids, n, xout) {

    if ( missing(dids) )
        dids <- data$dataIDs
    dids <- dids[dids!=xid] # rm target value

    ## store original master data
    orig.id <- data$mids[1]
    orig <- data[[orig.id]]
    
    ## get new master data
    xdat <- data[[xid]]$data
    if ( missing(n) ) n <- nrow(xdat)
    if ( missing(xout) ) {
        xout <- range(c(xdat),na.rm=TRUE)
        xout <- seq(xout[1], xout[2], length.out=n)
    }
    
    ## 2) interpolate all data to MASTER time
    for ( id in dids ) {
        data[[id]]$orig <- data[[id]]$data
        mdat <- data[[id]]$data
        tmp <- mdat ## reverse interpolation of original master data
        tmp[] <- NA
        for ( j in 1:ncol(data[[id]]$data) ) {
            x <- xdat[,j]
            y <- data[[id]]$data[,j]
            ## TODO: split into non-NA ranges of data
            if ( sum(!is.na(y))<2 ) next
            ## interpolate data, NOTE that rule=2 will fill the end points
            if ( length(unique(x)) == 1 )
                mdat[,j] <- NA
            else {
                ## TODO: replace by cubic spline fit with smart end handling!
                #mdat[,j] <- spline(x=x,y=y,xout=xout,method="natural")$y
                mdat[,j] <- approx(x=x,y=y,xout=xout,rule=1)$y
                ## TODO: reverse interpolation - store x data!
                tmp[,j] <- approx(x=y,y=orig,xout=mdat[,j])$y
            }
        }
        ## replace data
        data[[id]]$data <- mdat
        ## indicate interpolation
        data[[id]]$processing <- c(data[[id]]$processing,
                                   paste("interpolated to", xid))
        data[[id]]$reverse <- tmp
        names(data[[id]])[length(data[[id]])] <- orig.id
    }

    ## keep original master data to check
    ##old.id <- paste("original_",xid,sep="")
    ##names(data)[which(names(data)==xid)] <- 
    ##    data$dataIDs[data$dataIDs==xid] <- 
    ##   names(data$colors)[which(names(data$colors)==xid)] <- old.id
    
    data <- append(data, list(xout), after=0)
    names(data)[1] <- xid
    ## rm old master x-axes
    data <- data[-which(names(data)%in%data$mids)]
    ## and set new master
    data$mids <- xid
    ## set new data IDs
    data$dataIDs <- dids
    data
}


#' \code{\link{viewPlate}} plots all data in plate format
#' @param data the list of measurement data as provided by
#' \code{\link{readPlateData}}
#' @param wells a list of wells to plot, overrules \code{rows} and \code{cols}
#' @param wcols named color vector for wells, used for the well ID
#' @param rows a list of strings/characters used as row ID in the composite
#' row:col well description in the plate layout (map) and plate data
#' @param cols as rows but plate column IDs
#' @param xid ID of a data-set in the input data that can be used as x-axis
#' instead of the default Time vector
#' @param dids IDs of the data to be plotted; if missing, all data will
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
#' @param legpos position of the well IDs on the plots
#' @param add.legend add a legend for the plotted data types (see
#' argument \code{dids}) in the last plotted well
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
viewPlate <- function(data, wells, wcols,
                      rows=toupper(letters[1:8]),cols=1:12,
                      xid, xscale=FALSE,xlim,
                      dids, dtype="data", pcols, yscale=TRUE,ylims,ylim,log="",
                      legpos="topleft",add.legend=TRUE) {

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
    pwells <- unique(c(sapply(data$dataIDs,
                              function(id) colnames(data[[id]][[dtype]]))))
    if ( sum(!wells%in%pwells)>0 ) {
        #warning("wells ", wells[!wells%in%pwells]," not present, skipped!")
        wells <- wells[wells%in%pwells]
    }
    
    ## get x-axis data: time and temperature or another data set
    if ( missing(xid) )
        xid <- data$mids[1]
    
    ## TODO: interpolate data here on the fly, if not done
    ## upon parsing data or subsequently
    global.x <- xid %in% data$mids
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
    ## TODO: adapt to new data$mids,
    ## if xid is not present in data$mids, get it from dataIDs
    #time <- data[[data$mids[1]]]
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
    if ( !missing(dids) ) # only use requested data for plots
      ptypes <- ptypes[ptypes%in%dids]
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
## @example
## data(ap12)
## groups <- getGroups(plate=ap12plate, by=c("strain"))


#' group wells by experiment annotations (in plate map file)
#' @param plate the plate layout map, see \code{\link{readPlateMap}}
#' @param by a list of column IDs of the plate layout
#' @param order if TRUE groups will be alphabetically ordered
#' @param verb if TRUE report messages are more detailed
#' @details Calculates the distinct groups from the plate layout by the selected
#' experimental parameters.
#' @return Returns a list of well IDs for the identified grouping. This list
#' can be used, e.g., in viewGroups(data,groups) or \code{link{groupStats}}
#' to summarize data for these groups.
#' @seealso \code{\link{readPlateMap}}, \code{\link{viewGroups}}
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
getGroups <- function(plate, by="medium", order=FALSE, verb=TRUE) {

    ## get different types to be grouped as replicates
                                        #types <- rep(TRUE,nrow(plate))
    lpl <- matrix(unlist(lapply(plate,function(x) as.character(x))),
                  ncol=ncol(plate))
    colnames(lpl) <- colnames(plate)
    ## collapse requested combinations into new type
    types <- rep("",nrow(plate))
    for ( b in by ) 
      types <- paste(types,lpl[,b],sep="_")
    types <- sub("^_","",types) # rm leading _
    btypes <- unique(types)

    ## rm NA
    btypes <- btypes[btypes!=paste(rep("NA",length(by)),collapse="_")]
    
    ## blank wells
    blanks <- plate[,"blank"]
    blanks[is.na(blanks)] <- FALSE
    
    ## collect wells of each group
    groups <- rep(list(NA),length(btypes))
    names(groups) <- btypes

    for ( i in 1:length(btypes) ) {
        btyp <- btypes[i]
        ## data and blank wells of the correct type
        dwells <- as.character(plate[types==btyp & !blanks,"well"])
        bwells <- as.character(plate[types==btyp &  blanks,"well"])
        dwells <- dwells[!is.na(dwells)]
        bwells <- bwells[!is.na(bwells)]
        groups[[btyp]] <- dwells
        if ( verb ) 
          cat(paste("\tgroup", btyp, ":", length(dwells), "wells, skipping",
                    length(bwells), "blank wells\n"))
    }
    ## return non-empty groups
    groups <- groups[lapply(groups,length)>0]
    ## order groups by "by"
    if ( order )
        groups <- groups[order(names(groups))]
    groups
}

## TODO
## wrapper for viewGroups taking list of plates
## and a merged layout file
viewGroupl <- function(data, groups, groups2, ...) {}

## TODO: use this in viewGroups as well?
#' calculates simple statistics for groups as plotted in
#' \code{\link{viewGroups}}
#' @param data \code{\link{platexpress}} data, see \code{\link{readPlateData}}
#' @param groups a list of well grouped wells, as produced by
#' \code{\link{getGroups}}
#' @param dids data IDs for which statistics should be reported,
#' if missing stats for all data will be reported
#' @details Calculates the simple statistics over grouped wells
#' (means, 95% confidence intervals, stdandard errors) along the x-axis
#' (usually time).
#' @return Returns a data structure similar to the input data, but
#' where actual data is replaced  by statistics over grouped wells.
#' @seealso \code{\link{readPlateMap}}, \code{\link{viewGroups}}
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
groupStats <- function(data, groups, dids) {


    
    if ( missing(dids) )
        dids <- data$dataIDs
        
    for ( did in dids ) {
        SE <- matrix(NA,nrow=nrow(data[[did]]$data),ncol=length(groups))
        colnames(SE) <- names(groups)
        MN<-CI<-SE
        for ( sg in 1:length(groups) ) {
            wells <- groups[[sg]]
            #wells <- wells[wells%in%pwells] # filter for present wells
            sid <- names(groups)[sg]
            
            ## get data for selected wells
            dat <- data[[did]]$data[,wells]
            ## calculate stats only for common x!
            MN[,sg] <- apply(dat,1,function(x) mean(x,na.rm=TRUE))
            SE[,sg] <- apply(dat,1,function(x) se(x,na.rm=TRUE))
            CI[,sg] <- apply(dat,1,function(x) ci95(x,na.rm=TRUE))
        }
        data[[did]]$stats <- list(mean=MN,ci05=CI,se=SE)
    }
    data
}

#' get colors for all wells in a group
#'
#' returns a named vector of colors for a grouping if all
#' colors in that group are unique
#' @param map the plate layout map, see \code{\link{readPlateMap}},
#' for which a coloring scheme has been specified, eg. by
#' \code{\link{amountColors}}
#' @param group a well grouping, see \code{\link{getGroups}}
#' @param color the name of the column in \code{map} providing colors
#' @export
groupColors <- function(map, group, color="color") {
    grcols <- rep(NA,length(group))
    names(grcols) <- names(group)
    #wells <- rownames(map)
    #if ( is.null(wells) ) wells <- map[,"well"]
    for ( i in 1:length(group) ) {
        grp <- group[[i]]
        cols <- unique(map[match(grp,map$well),color])
        if ( length(cols)!=1 )
            stop("different colors observed in group", i, names(group)[i])
        grcols[i] <- cols
    }
    grcols
}


#' plot grouped wells as summary plots, incl. confidence intervals and means
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
#' @param dids IDs of the data to be plotted; if missing, all data will
#' be plotted
#' @param dtype type of the data to be plotted, default is the main 'data', but
#' e.g. 'orig' allows to plot the original data w/o processing (e.g.
#' interpolation, blank correction, etc.)
#' @param xscale use a global range for the x-axes; only relevant if xid
#' specifies a subset of the data as x-axis
#' @param xlim plot range of the x-axis
#' @param xlab custom x-axis label
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
#  of the plotted data types (\code{dids})
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
#' @param mgp set the position of axis title, tick marks and tick lengths
#' @param xaxis plot x-axis if TRUE
#' @param yaxis the data types for which axes are to be plotted, corresponds
#' to the order in argument \code{dids} and as plotted in the legend 
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
                       dids, dtype="data",pcols,group2.col,
                       yscale=TRUE,ylims,ylim, log="",
                       show.ci95=TRUE,show.mean=TRUE,emphasize.mean=FALSE,
                       lty.orig=1,lwd.orig=0.1,lty.mean=1,lwd.mean=2,
                       g2.legpos="topleft", g2.legend=TRUE,
                       embed=FALSE, no.par=FALSE,
                       mai=c(0.5,0,0,0), mgp=c(1.5,.5,0),
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
        xid <- data$mids[1]
    
    ## TODO: interpolate data here on the fly, if not done
    ## upon parsing data or subsequently
    global.x <- xid %in% data$mids
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
    if ( !missing(dids) ) # only use requested data for plots
      ptypes <- ptypes[ptypes%in%dids]
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
        par(mai=mai,mgp=mgp)
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
            if ( !missing(group2.col) ) # colors were explicitly provided
                orig.cols <- group2.col[names(sgroups)]
          
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
                if ( !missing(group2.col) ) 
                    g2.col <- group2.col[sid]
                ##cat(paste(sg, g2.col, "\n"))
                
                ## override color to allow lwd.orig=0 to work for PDF as well
                tmp <- ifelse(lwd.orig==0,NA, col.orig)

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
                        if ( log=="y") {                            
                            py[py<0] <- ylims[[ptyp]][1]
                            warning("mean - 95% c.i. is below 0;",
                                    "raised to ylim for log. y-axis polygon")
                        }
                        polygon(x=px,y=py,border=NA,
                                col=paste(g2.col,"55",sep=""))
                    }
                    if ( show.mean )
                        lines(x=x,mn,col=ifelse(emphasize.mean,1,g2.col),
                              lwd=lwd.mean,lty=ifelse(g2.legend,sg,lty.mean))
                }

                ## add axes for first two values
                if ( yaxis[1]==i & sg==1 ) 
                    axis(2, tcl=.25, mgp=c(0,-1,-.05),
                         col=ifelse(global.x,col.orig,1),
                         col.axis=ifelse(global.x,col.orig,1))
                if ( yaxis[2]==i & sg==1 ) axis(4, tcl=.25, mgp=c(0,-1,-.05),
                            col=ifelse(global.x,col.orig,1),
                            col.axis=ifelse(global.x,col.orig,1))
            }
        }
        if ( length(sgroups)>1 & g2.legend )
            if ( global.x & missing(group2.col) ) 
                legend(g2.legpos,names(sgroups),lty=1:length(sgroups),
                       col=1,bty="n")
            else
                legend(g2.legpos,names(sgroups),lty=1:length(sgroups),
                       col=orig.cols,bg="#FFFFFFAA",box.lwd=NA) # TODO: use g2cols
        else
            legend(g2.legpos,id, bty="n")
        if ( xaxis ) axis(1)
        if ( missing(xlab) ) 
          xlab <- xid
        mtext(xlab, 1, par("mgp")[1])
    }
    ## add legend to last plot
    if ( g1.legend )
        legend(g1.legpos,ptypes,lty=1,col=pcols[ptypes],bg="#FFFFFFAA")

    ## reset par!
    par(orig.par)
    
    ## TODO: return meaningful and/or non-plotted information
    ## assigning it makes it silent!
    #if ( global.x ) xid <- "Time"
    plotparams <- list(ylims=ylims, xid=xid, xlim=xlim,  colors=pcols, orig.cols=orig.cols)
}


### COMMENTS FOR EXAMPLE DATA
  
 
#' ap12plate. example plate layout map by Dennis Dienst and Alice Pawloski,
#' fitting to the experimental data in \code{\link{ap12data}}
#' The plate layout table indicates the different strains tested, biological
#' replicates (B1 to B3), and blank wells (containing only growth medium) 
#'
#' @name ap12plate
#' @docType data
NULL
## @keywords datasets
## @format the plate layou map as a a data.frame, as produced by
## readPlateMap("AP12_layout.csv", fields=c("strain","samples"))
## @seealso \code{\link{ap12data}}, \code{\link{readPlateData}} and
## \code{\link{readPlateMap}}


#' ap12data. data by Dennis Dienst and Alice Pawloski, incl. the
#' plate reader measurements of E.coli growth, expressing a fluorescent
#' proteins, in a Synergy Mx platereader; the corresponding plate layout
#' map is in \code{\link{ap12plate}}.
#' 
#' \itemize{
#'   \item Data:
#'   \item Time: the interpolated 'master' time
#'   \item Temperature: the temperature time-course
#'   \item Data matrix '600': well absorbance at 600 nm, i.e., the OD,
#'   \item Data matrix 'YFP_50:500,535': the YFP fluorescence measured by excitation at 500 nm and emission at 535 nm 
#' }
#'
#' @name ap12data
#' @docType data
NULL
## @keywords datasets
## @format a list of time-courses of absorbance and fluorescence data, read
## in by readPlateData("AP12.csv", type="Synergy", data.ids=c("600","YFP_50:500,535"), time.format="%H:%M:%S", time.conversion=1/3600)
## @seealso \code{link{ap12plate}}, \code{\link{readPlateData}} and
## \code{\link{readPlateMap}} 
