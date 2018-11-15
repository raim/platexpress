#' platexpress: quick inspection and analysis of microbial growth & expression.
#'
#' The platexpress package provides a quick & easy interface to
#' inspect microbial growth & gene expression data as measured in typical
#' platereaders ("96-well plates") or other parallel growth systems.
#' Interfaces to existing packages allow growth model fitting
#' (currently \pkg{grofit}, \pkg{growthrates}, and
#' partially/untested \code{\link[cellGrowth:fitCellGrowth]{cellGrowth}}.)
#'
#'@author Rainer Machne
#'@docType package
#'@name platexpress
#'@section Dependencies: The package uses mostly functionality from R base,
#' (graphics, grDevices, stats) but more functionality is available when
#' \pkg{grofit} or \pkg{growthrates} are installed.
#'@importFrom stats median sd qt approx spline filter predict coef na.omit quantile confint lm
#'@importFrom graphics plot matplot boxplot barplot legend arrows locator
#' abline lines points polygon box axis par text title mtext stripchart image
#'@importFrom grDevices rainbow rgb col2rgb png pdf svg tiff jpeg postscript graphics.off gray.colors
#'@importFrom tidyr separate
#'@importFrom readxl read_excel
#'@importFrom utils read.csv read.table
NULL


### UTILS

#' Select Plot Device
#'
#' plots to png, eps, pdf, tiff, svg or jpeg devices taking the same
#' width, height and resolution arguments for all
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
    grDevices::png(file, width=width, height=height, units="in", res=res)
  if ( type == "eps" )
    grDevices::postscript(file, width=height, height=width,paper="special",
                          horizontal = FALSE, onefile = FALSE)
  if ( type == "pdf" )
    grDevices::pdf(file, width=width, height=height)
  if ( type == "tiff" )
    grDevices::tiff(file, width=width, height=height, units="in", res=res)
  if ( type == "svg" )
    grDevices::svg(file, width=width, height=height)
  if ( type == "jpeg" )
    grDevices::jpeg(file, width=width, height=height, units="in", res=res)
}


## return R colors as RGB, to allow setting alpha
getRGB <- function(n) {
    cols <- col2rgb(1:n)
    apply(cols,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))
}


#' Show Visual Ligh Spectrum
#'
#' shows the color spectrum of visible light along wavelengths in nm.
#' NOTE that this is not a fully correct spectrum,
#' see \code{\link{wavelength2RGB}} for details.
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

#' Plot Wavelength in Spectrum
#'
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

#' Find Wavelength in Spectrum
#'
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

#' Convert Wavelength to RGB Colors
#'
#' Converts wavelength in nm (visible light: 380:780 nm) to RGB.
#' NOTE that this is not a fully correct spectrum, see Details.
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



### STATS

## moving average
#' Moving Average
#'
#' calculate a moving average using \code{\link[stats]{filter}}
#' @param x data vector along which a moving average will be calculated
#' @param n moving average window size
#' @param circular logical see help of function \code{\link[stats]{filter}}
#' @export
ma <- function(x, n=5, circular=FALSE) {
  filter(x,rep(1/n,n), sides=2, circular=circular)
}

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



#' Set Data IDs and colors
#'
#' Sets colors, rename, order or filter the data set
#' @param data \code{\link{platexpress}} data, see \code{\link{readPlateData}}
#' @param yids a vector of data IDs, data will be filtered and sorted by this list; if the vector is named the IDs will be replaced by these names
#' @param colors a vector of plot colors as RGB strings, optionally already named by dataIDs
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{readPlateData}}
#' @export
prettyData <- function(data, yids, colors) {

    ## yids: data sorting, if missing, take all
    if ( missing(yids) )
        yids <- data$dataIDs

    ## get xids (global x-axis data)
    xids <- data$xids

    if ( any(is.na(match(yids,names(data)))))
      stop("IDs (yids) not found in data: ", paste(yids,sep=";"))

    ## re-order: xids first and IDs last
    data$dataIDs <- yids
    ## store colors
    origcols <- data$colors
    data <- data[c(match(xids,names(data)),
                   match("xids",names(data)),
                   match(yids,names(data)),
                   match("dataIDs",names(data)))]
    data$colors <- origcols

    ## rename!
    if ( !is.null(names(yids)) ) {
        names(data)[match(yids,names(data))] <- names(yids)
        data$dataIDs[match(yids,data$dataIDs)] <- names(yids)
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
#' Add Data
#'
#' Adds a data set, e.g., calculated ratios, smoothed data, etc.
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
    
    ## copy into existing first data frame
    ## only if colnames are available!
    datf <- data[[data$dataIDs[1]]]$data
    if ( all(colnames(dat)%in%colnames(datf))) {
        datf[] <- NA
        datf[,colnames(dat)] <- dat
        dat <- datf
    } 
    
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
            processing <- paste("added",date())
        data <- append(data, list(list(data=dat, processing=processing)))
        names(data) <- c(names(data)[2:length(data)-1],ID)
    }
    data
}


#' Remove Data
#'
#' removes a data set from plate data
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


#' Get Data
#'
#' Returns a specific data set as data matrix
#' @param data the current platexpress data set
#' @param ID the ID of the data to be obtained
#' @param type the type of the data to be returned (default "data"),
#' the data list also contains original values for some of the processing
#' steps
#' @param xrng x-axis range or single point; similar to cutData, but
#' use interpolation for single values!
#' @param xid ID of the x-axis data to use
#' @param verb print messages
#' @export
getData <- function(data, ID, type="data", xrng, xid, verb=TRUE) {
  if ( missing(xrng) )
    return(data[[ID]][[type]]) # just return the current data or old versions
  else {
    #stop("x-axis range not implemented yet")
    if ( missing(xid) )
      xid <- data$xids[1]
    xdat <- data[[xid]]
    if ( length(xrng)==2 ) {
      if ( verb ) cat(paste("cutting data to x within", paste(xrng,collapse=":"), "\n"))
      filter <- xdat >= xrng[1] & xdat <= xrng[2]
      return(data[[ID]][[type]][filter,,drop=FALSE])
    } else if ( length(xrng)==1 ) {
      if ( verb ) cat(paste("interpolating to", xrng, "\n"))
      if ( xrng > max(xdat, na.rm=TRUE) )
        stop("requested value: ",xrng,
             " is outside of data range (",paste(range(xdat),collapse=":"),
             "); sorry, can't extrapolate")
      return(apply(data[[ID]][[type]], 2, function(y)
        ifelse(sum(!is.na(y))<2, NA, approx(x=xdat, y=y, xout=xrng)$y)))
    }
  }
}


#' Shift Data
#'
#' Shifts the x-axis, eg. by calculated lag-phases, to align growth curves
#' @param data \code{\link{platexpress}} data set
#' @param lag a named vector providing the lag-phase to be removed; the names
#' correspond to the wells in data
#' @param xid the x-axis to be used, defaults to the first available
#' (usually "Time")
#' @export
shiftData <- function(data, lag, xid) {

    ## NOTE: not optional, since data don't fit anymore
    ## if only some are shifted for a given well list (names of lag)
    yids <- data$dataIDs

    if ( missing(xid) )
        xid <- data$xids[1]
    xdat <- data[[xid]]

    for ( i in 1:length(lag) ) {
        well <- names(lag)[i]
        idx <- which(xdat >= lag[i])[1]
        end <- length(xdat)
        for ( yid in yids ) {
            data[[yid]]$data[,well] <-
                          c(data[[yid]]$data[idx:end, well],
                            rep(NA,idx-1))
        }
        data[[yid]]$processing <- c(data[[yid]]$processing,
                                    paste("well",well,"shifted by lag", lag[i]))
    }
    data
}

## TODO: cut data either by time, or by a chosen data set
#' Cut Data Range
#'
#' Cuts data to a range or single point of the x-axis, and/or cuts
#' all y-axis values within a range
#' @param data \code{\link{platexpress}} data, see \code{link{readPlateData}}
#' @param xrng a single or two value(s) for x-axis cuts, only within the range
#' of two values will be retained; if only one value is passed only the
#' data closest to this point will be retained!
#' @param xid ID of the x-axis data to cut, default is the main x-axis ('Time')
#' @param yid ID of the y-axis data to cut
#' @param yrng a single or two value(s) for x-axis cuts, only data smaller
#' then single value, or within the range of two values will be retained
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @details Cuts the passed \code{\link{platexpress}} data to ranges of
#' of the x-axis (time or other, see \code{data$xids}) and/or y-axis.
#' Note that the behaviour is different for passing single values to
#' \code{xrng} or \code{yrng}:
#' If \code{xrng} is a single value, data for the closest x value will be
#' returned. If \code{yrng} is a single value, all data smaller than this
#' value will be returned. Note that data outside \code{yrng} are simply
#' set to NA, which may cause problems downstream.
#' @export
cutData <- function(data, xrng, xid, yid, yrng) {

    ## default x-axis cut
    if ( missing(xid) & missing(yid) )
        xid <- data$xids[1]

    ## TODO: allow cutting both x and y data
    ## default no xid and yid: cut xid
    ## only xid: cut xid
    ## only yid: cut yid
    ## both yid and xid: cut both, first xid

    ## first, cut x-axis
    if ( !missing(xid) ) {
        xdat <- data[[xid]]
        if ( length(xrng)==2 )
            filter <- xdat >= xrng[1] & xdat <= xrng[2]
        else if ( length(xrng)==1 )
            filter <- abs(xrng-xdat) == min(abs(xrng-xdat))
        for ( yid in data$dataIDs ) {
            data[[yid]]$data <- data[[yid]]$data[filter,,drop=FALSE]
            data[[yid]]$processing <- c(data[[yid]]$processing,
                                        paste("cut at", paste(xrng,collapse="-")))
        }
        data[[xid]] <- xdat[filter]
    }
    ## second, cut y-axis
    ## TODO: cut ALL Y-DATA at these points?
    ## TODO: cut time at points where all are NA!!
    if ( !missing(yid) & !missing(yrng) ) {
        ## get indices where to cut
        ydat <- data[[yid]]$data

        ## expand x range to minimum - maximum
        if ( length(yrng)==1 ) ## max value
            yrng <- c(min(apply(ydat, 2, min, na.rm=T),na.rm=T),
                      yrng)
        ## just set to NA
        ydat <- apply(ydat, 2, function(y) {
            y[y < yrng[1] | y > yrng[2]]<-NA;y})
        data[[yid]]$data <- ydat
    }
    data
}



#' Skip Wells
#'
#' Removes wells from plate data, plate layout maps and well groupings
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

    ## TODO: use classes
    
    if ( "dataIDs" %in% names(data) ) ## rm from data
      for ( id in data$dataIDs ) {
          wells <- colnames(data[[id]]$data)
          if ( any(wells%in%skip) )
              data[[id]]$data <- data[[id]]$data[,-which(wells%in%skip),
                                                 drop=FALSE]
      }
    else if ( !is.null(dim(data)) ) ## rm from plate layout map
        ##data[match(skip,data[,"well"]),2:ncol(data)] <- NA
        data <- data[-match(skip,data[,"well"]),]
    else ## rm from grouping
      for ( g in 1:length(data) )
        data[[g]] <- data[[g]][!data[[g]]%in%skip]
    data
}

#' Get List of Wells
#'
#' get a filtered list of wells that match argument \code{values}
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

#' Subtract Blank Values
#'
#' The function subtracts values from "blank" wells. Optionally this can
#' be done in bins over time (or the current x-axis value) to account
#' for time-dependent blanks. E.g. fluorescence blanks from LB medium
#' sometimes show time-dependence, perhaps due to light-dependent
#' degradation of LB fluorescence. Separate blanks for each condition can
#' be used via the \code{by} option, to be used the same way as
#' in \code{\link{getGroups}}.
#' @param data the plate data list to be blank-corrected
#' @param plate the plate layout where column "blanks" indicates which wells
#' are to be treated as blanks
#' @param yids IDs of the data which should be blank-corrected, all will be
#' blanked if missing
#' @param by a list of column IDs of the plate layout; separate blank
#' correction will be attempted for groups in these columns; each group
#' must have at least one specified blank associated
#' @param type calculation of blank values from multiple time-points and
#' wells; "median", "mean" or "ci95", where the latter subtracts the mean
#' minus  the 95\% confidence interval to avoid blanked values below 0
#' @param xid ID of the x-axis data to be used, if blanked along x-axis, set
#' by \code{mbins}>1
#' @param max.xid the maximal x-axis value where blanks should be used
#' @param mbins the number of bins the x-axis is to be divided, if blanked
#' along the x-axis, see \code{xid}
#' @param base optional minimal value; all values will be raised by
#' the same amount using the function \code{\link{adjustBase}}
#' @param verb issued progress messages and info
#' @param ... further arguments to \code{\link{adjustBase}}
#' @seealso \code{\link{adjustBase}}
#' @examples
#' data(ap12)
#' data <- correctBlanks(data=ap12data, plate=ap12plate, by="strain")
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
correctBlanks <- function(data, plate, type="median", by, yids,
                          xid, max.xid, mbins=1, base, verb=TRUE, ...) {

### TODO: correct by time point, eg. for fluorescence in ecoli_rfp_iptg_20160908

    if ( missing(xid) )
        xid <- data$xids[1]

    ## start new data list
    corr <- data
    time <- data[[xid]] ## TODO: take from data xids

    ## reduce matrix to requested data
    data <- data[data$dataIDs]
    ptypes <- names(data)
    if ( !missing(yids) ) # only use requested data
      ptypes <- ptypes[ptypes%in%yids]
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
                if ( !missing(max.xid) ) {
                    if ( verb )
                      cat(paste("\tskipping",sum(time[bbin]>max.xid),
                                "bins at",max.xid,"\n"))
                    bbin <- bbin[time[bbin]<=max.xid]
                }

                ## TODO: this should only happen if time is
                if ( length(bbin)==0 ) {
                    #cat(paste("skipping time bin", t, "at", max.xid, "\n"))
                    warning("no blank data at time bin", t, "at", max.xid)
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
    if ( !missing(base) )
      corr <- adjustBase(corr, base=base, yids=yids, verb=verb, ...)
    corr
}

#' Adjusts to Minimal Base
#'
#' The function raises all data to a "base" level, default 0, to avoid
#' negative values that sometimes occur after blank correction
#' in \code{\link{correctBlanks}}. The function can be optionally
#' called directly in \code{\link{correctBlanks}} by option \code{base}.
#' @details Adjusts data to a new mininum, this is useful for adjustment
#' of negative values after blank corrections.
#' @param data \code{\link{platexpress}} data, see \code{\link{readPlateData}}
#' @param yids vector of ID strings for which base correction should be
#' executed
#' @param base the new minimum for the data, default is 0, but it could
#' e.g. be the OD used for inoculation
#' @param xlim min and max row number of the data to be adjusted
#' @param add.fraction a fraction of the whole data range, added to base
#' @param each add base for each well separately!
#' @param verb print messages if true
#' @return Returns `data' where all data sets (or only those specified in
#' option \code{yids}) where raised to a minimum level in
#' @seealso \code{\link{correctBlanks}}
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
adjustBase <- function(data, base=0, yids, add.fraction, xlim, each=FALSE, verb=TRUE) {

    if ( missing(yids) ) # only use requested data
        yids <- data$dataIDs # use all

    if ( verb )
      cat(paste("adjusting base", paste(yids, collapse=";"),"\n"))

    for ( yid in yids ) {

        ## each well separately?
        if ( each )
            bins <- as.list(1:ncol(data[[yid]]$data))
        else
            bins <- list(1:ncol(data[[yid]]$data))

        for ( bin in bins ) {

            dat <- data[[yid]]$data[,bin,drop=FALSE]

            if ( missing(xlim) )
                xlim <- c(1,nrow(dat))
            xrng <- xlim[1]:xlim[2]

            if ( verb )
              cat(paste("\t\t",length(bin), #paste(bin,collapse=";"),
                        "wells, adding", min(dat[xrng,],na.rm=TRUE), "\n"))

            ## TODO: smarter? only if any value is <0?
            dat <- dat - min(dat[xrng,],na.rm=TRUE) + base

            if ( !missing(add.fraction) )
                dat <- dat + diff(range(dat[xrng,],na.rm=TRUE))*add.fraction

            data[[yid]]$data[,bin] <- dat
        }
        data[[yid]]$processing <- c(data[[yid]]$processing,
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

#' Interpolate to Common Timepoints
#'
#' interpolates all time-series of a plate to an average
#' master time, using the R base function \code{\link[stats:splinefun]{spline}}.
#' An average time is calculated for each measurement point and all values
#' are interpolated to these new time points. This is automatically done when
#' parsing the raw data with \code{\link{readPlateData}}, unless
#' explicitly suppressed. The same is also done for well temperatures.
#' @param data  \code{\link{platexpress}} data, see \code{\link{readPlateData}}
#' @param verb  print messages if true
#' @return returns a copy of the full data list with a master time and
#' temperature added at the top level
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
interpolatePlateTimes <- function(data, verb=TRUE) {


    ## catch single data case
    if ( length(data$dataIDs)==1 ) {
        data$Time <- data[[data$dataIDs[1]]]$time
        data$xids <- c(data$xids, "Time")
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
            mdat[,j] <- stats::spline(x=x,y=y,xout=mtime,method="fmm")$y

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
    data$xids <- c(data$xids, "Time")
    #data$Temperature <- mtemp
    data
}

#' Interpolate Plate Data
#'
#' interpolate one dataset to common points of another
#' data set. This is used to change the x-axis of a data set, e.g.,
#' the measured OD values can become the new x-axis and all fluorescence
#' values will be interpolated to common OD vlues, using the
#' R base function \code{\link[stats:splinefun]{spline}}, the same way
#' as in \code{\link{interpolatePlateTimes}}.
#' @param data the list of measurement data as provided by
#' \code{\link{readPlateData}}
#' @param xid ID of the data set which should serve as the new x-axis
#' all data will be interpolated to equi-spaced points along the range
#' of measured values
#' @param yids restrict interpolation to these data IDs
#' @param n specify the number of interpolation points, if missing the
#' original number of rows will be used
#' @param xout specify the interpolation points directly
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
interpolatePlateData <- function(data, xid, yids, n, xout) {

    if ( missing(yids) )
        yids <- data$dataIDs
    yids <- yids[yids!=xid] # rm target value

    ## store original master data
    orig.id <- data$xids[1]
    orig <- data[[orig.id]]

    ## get new master data
    xdat <- data[[xid]]$data
    if ( missing(n) ) n <- nrow(xdat)
    if ( missing(xout) ) {
        xout <- range(c(xdat),na.rm=TRUE)
        xout <- seq(xout[1], xout[2], length.out=n)
    }

    ## 2) interpolate all data to MASTER time
    for ( id in yids ) {
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
    data <- data[-which(names(data)%in%data$xids)]
    ## and set new master
    data$xids <- xid
    ## set new data IDs
    data$dataIDs <- yids
    data
}


## @example
## data(ap12)
## groups <- getGroups(plate=ap12plate, by=c("strain"))


#' Get Groups of Replicates
#'
#' groups wells by experiment annotations in the plate layout map
#' (selected by option \code{by}), and returns a list of well IDs
#' that are all replicates for these groups.
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


## TODO: use this in viewGroups as well?
#' Group Statistics
#'
#' calculates simple statistics for groups as plotted in
#' \code{\link{viewGroups}}. TODO: actually use these stats
#' in viewGroups, optionally add group stats to original data structure
#' @param data \code{\link{platexpress}} data, see \code{\link{readPlateData}}
#' @param groups a list of well grouped wells, as produced by
#' \code{\link{getGroups}}
#' @param yids data IDs for which statistics should be reported,
#' if missing stats for all data will be reported
#' @details Calculates the simple statistics over grouped wells
#' (means, 95% confidence intervals, stdandard errors) along the x-axis
#' (usually time).
#' @return Returns a data structure similar to the input data, but
#' where actual data is replaced  by statistics over grouped wells.
#' @seealso \code{\link{readPlateMap}}, \code{\link{viewGroups}}
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
groupStats <- function(data, groups, yids) {

    if ( missing(yids) )
        yids <- data$dataIDs

    for ( yid in yids ) {
        SE <- matrix(NA,nrow=nrow(data[[yid]]$data),ncol=length(groups))
        colnames(SE) <- names(groups)
        MN<-CI<-SD<-SE
        for ( sg in 1:length(groups) ) {
            wells <- groups[[sg]]
            #wells <- wells[wells%in%pwells] # filter for present wells
            sid <- names(groups)[sg]

            ## get data for selected wells
            dat <- data[[yid]]$data[,wells]
            ## calculate stats only for common x!
            MN[,sg] <- apply(dat,1,function(x) mean(x,na.rm=TRUE))
            SD[,sg] <- apply(dat,1,function(x) sd(x,na.rm=TRUE))
            SE[,sg] <- apply(dat,1,function(x) se(x,na.rm=TRUE))
            CI[,sg] <- apply(dat,1,function(x) ci95(x,na.rm=TRUE))
        }
        data[[yid]]$stats <- list(mean=MN,sd=SD,se=SE,ci05=CI)
    }
    data
}

#' Group Colors
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
