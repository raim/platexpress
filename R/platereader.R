
### ANALYZING PLATE-READER GROWTH & EXPRESSION CURVES
#' platexpress: A package for analysing microbial growth & expression.
#'
#' The platexpress package provides a quick&easy interface to 
#' microbial growth & gene expression data as measured in typical
#' microplate-readers or other parallel growth systems.
#' 
#' @section Platexpress Workflow
#' @examples
#' ### A TYPICAL WORKFLOW
#' ## 1) parse the plate layout map
#' 
#' plate.file <- system.file("extdata", "AP12_layout.csv", package = "platexpress")
#' plate <- readPlateMap(file=plate.file, blank.id="blank",fsep="\n", fields=c("strain","samples"))
#' 
#' ## 2) parse the data, exported from platereader software
#' 
#' data.file <- system.file("extdata", "AP12.csv", package = "platexpress")
#' raw <- readPlateData(file=data.file, type="Synergy", data.ids=c("600","YFP_50:500,535"), dec=",")
#' 
#' ## 3) inspect the raw data
#' 
#' vp <- viewPlate(raw)
#' 
#' ## 4) Note that there is no growth in A9, so let's skip it
#' 
#' raw <- skipWells(raw, skip="A9")
#' 
#' ## 5) Now correct for blank well measurements, and view only present
#' ## rows/cols
#' 
#' data <- correctBlanks(data=raw, plate=plate)
#' vp <- viewPlate(data, rows=c("A","B","C"),cols=1:9)
#' 
#' ## 6) group replicates and view summarized growth/exprssion curves
#' 
#' groups <- getGroups(plate, by=c("strain","samples"))
#' vg <- viewGroups(data,groups=groups,lwd.orig=0.5,nrow=3)
#' 
#' @docType package
#' @name platexpress
NULL


### HELPERS

## return R colors as RGB, to allow setting alpha
getRGB <- function(n) {
    cols <- col2rgb(1:n)
    apply(cols,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))
}


#' \code{\link{viewSpectrum}}:
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
#' @details Plots the color spectrum over the selected wavelengths in nm.
#' The utility functions  \code{\link{plotWavelength}} and
#' \code{\link{findWavelength}} can be used to plot and interactively (click
#' on plot) find wavelengths. The function \code{\link{wavelength2RGB}} can
#' then be used to convert this wavelength to RGB colors.
#' @seealso \code{\link{wavelength2RGB}}, \code{\link{findWavelength}},
#' \code{\link{plotWavelength}}
#' @export
viewSpectrum <- function(wavelengths=380:780, alpha=99, pch=19, cex=3,
                         ylab="approximate color",
                         main="use findWavelength(353)",
                         xlab="wavelength, nm", ...) {
    cols <- wavelength2RGB(wavelengths)
    plot(wavelengths,rep(1,length(wavelengths)),
         axes=FALSE,ylab=ylab,xlab=xlab,main=main,
         col=paste(cols,alpha,sep=""),ylim=c(.9,1.1),pch=pch,cex=cex,...)
    axis(1,at=pretty(wavelengths))
}
#' \code{\link{plotWavelength}}:
#' a color selector for \code{\link{viewSpectrum}}; draws a vertical line,
#' the wavelength, and a filled circle at a given wavelength in nm.
#' @param x wavelength in nm
#' @param y position of the text, from 0.9 to 1.1
#' @param ... further arguments to ?text
#' @param ych y position of the symbol
#' @param pch the plot symbol type, see ?par("pch")
#' @param cex the plot symbol size, see ?par("cex")
#' @seealso \code{\link{viewSpectrum}}, \code{\link{findWavelength}},
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
#' a color selector for \code{\link{viewSpectrum}}; expects the user to
#' click on the spectrum, then draws a vertical line at the clicked
#' wavelength using \code{\link{plotWavelength}} and records the wavelength
#' in nm.
#' @param n number of clicks to record
#' @param ... further arguments to \code{\link{plotWavelength}}, for
#' selecting plot symbol and text size, and positions
#' @seealso \code{\link{viewSpectrum}}, \code{\link{plotWavelength}}, \code{\link{wavelength2RGB}}
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
#' @examples
#' wavelengths <- seq(380,780,1)
#' cols <- sapply(wavelengths, wavelength2RGB)
#' bars <- rep(1,length(wavelengths)); names(bars) <- wavelengths
#' barplot(bars,border=cols,col=cols,las=2)
#' @seealso \code{\link{viewSpectrum}}
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

# calculate 95% confidence intervals for the given
# data vector using a t-distribution
ci95 <- function(data,na.rm=FALSE) {
    if ( na.rm ) data <- data[!is.na(data)]
    n <- length(data)
    error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
    return(error)
}


### PLATE DESIGN & DATA

## Plate Design: parses a plate design file in CSV. Rows and columns
## should be named as in the datafile. Each field can
## have multiple descriptors, separated by a specified separator (e.g.
## "newline"); blanks are specified by a keyword (default: "blank"),
## and if separate blanks are used for different conditions, the
## blank field must have the same format as measurement fields, except
## one parameter replaced by the keyword. Names for the values
## in the fields can be passed via argument "fields" as a vector of
## strings.
## TODO: repair this in example:
#' \code{\link{readPlateMap}} parses a plate design file in CSV. Rows and 
#' columns should be named as in the corresponding data files.
#' @param file text file containing the plate layout.
#' @param sep column separator, as in read.table
#' @param fsep within-field separator, separating the well-specific descriptors
#' @param blank.id keyword that indicates blank wells. Blank wells can be
#'                 combined with other well descriptors for separate blanking
#' @param fields names for the field descriptors
#' @param formatted indicates whether the file is already in the required
#'                  format; all other paramaters but 'sep' will be ignored
#' @return a table of well content descriptors, where the first column 'wells'
#'         maps the plate map to the data files.
#' @seealso \code{\link{readPlateData}}
#' @examples
#' plate.file <- system.file("extdata", "AP12_layout.csv", package = "platexpress")
#' plate <- readPlateMap(file=plate.file, blank.id="blank",fsep="\n", fields=c("strain","samples"))
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
readPlateMap <- function(file, sep="\t", fsep="\n", blank.id="blank",
                         fields, nrows=8, formatted=FALSE) {


    ## already in well format?
    if ( formatted ) {
        dat <- read.table(file, sep=sep, header=TRUE, stringsAsFactors=FALSE)
        return(dat)
    }
    
    ## plate map by columns and rows
    dat <- read.table(file, sep=sep, header=TRUE, row.names=1, nrows=nrows,
                      stringsAsFactors=FALSE)

    ## generate well names from row and column names
    plate <- paste(rep(rownames(dat),ncol(dat)),
                   rep(sub("X","",colnames(dat)),each=nrow(dat)),sep="")

    ## parse field values
    vals <- strsplit(unlist(dat),fsep)
    nvals <- max(unlist(lapply(vals,length)))

    values <- matrix("",nrow=length(vals),ncol=nvals)
    for ( i in 1:nvals ) 
      values[,i] <- as.character(unlist(lapply(vals, function(x) trim(x[i]))))

    if ( missing(fields) )
      colnames(values) <- paste("X",1:nvals,sep="")
    else colnames(values) <- fields

    plate <- cbind(well=plate,values)
    
    blanks <- apply(plate,1, function(x) any(x==blank.id))
    plate[which(plate==blank.id)] <- NA

    plate <- cbind(data.frame(plate),blank=blanks)
    
    return(plate)
}

### PLATE DATA

## TODO: repair this in example
#
#' \code{\link{readPlateData}} parses data files in CSV, as exported by
#' the plate reader software. Header IDs in the data file should match with 
#' IDs in the plate map, see \code{link{readPlateMap}}. Pre-defined read-in
#' functions exist for a couple of plate-readers.
#' @param files list of one or more data files
#' @param sep column separator, as in read.table
#' @param type pre-defined formats, as exported from platereaders; currently
#' for BMG Optima/Mars, ('BMG') and Synergy Mx ('Synergy').
#' @param interpolate if true a master time, the average time between distinct
#' measurements of one timepoint, is calculated and all values are interpolated
#' to this mastertime. This is currently obligatory for further processing.
#' @param data.ids an optional sub-selection of data types in the input file,
#' as a list of strings
#' @param skip lines to skip from the data file
#' @param dec decimal operator used in the data file
#' @param verb print messages if true
#' @note The original data is all interpolated to a common/average 'master' time
#' @return a list of distinct measurement time-courses from one plate
#' @seealso \code{\link{readPlateMap}}, \code{\link{viewPlate}}
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @examples
#' data.file <- system.file("extdata", "AP12.csv", package = "platexpress")
#' raw <- readPlateData(file=data.file, type="Synergy", data.ids=c("600","YFP_50:500,535"), dec=",")
#' @export
readPlateData <- function(files, type, data.ids, interpolate=TRUE,
                          skip=0, sep="\t", dec=".", verb=TRUE, ...) {

    if ( type=="BMG" )
      readBMGPlate(files=files, data.ids=data.ids, interpolate=interpolate,
                   verb=verb, skip=5, sep=";", dec=".", ...)
    else if ( type=="Synergy" )
      readSynergyPlate(file=files, data.ids=data.ids, interpolate=interpolate,
                       verb=verb, skip=58, sep=";", dec=".", ...)
} 

# Read Synergy Mx-exported files
#' @inheritParams readPlateData
#' @seealso \code{\link{readPlateData}}
#' @export
readSynergyPlate <- function(file, data.ids, interpolate=TRUE,
                             skip=58, sep=";", dec=".", skiplastcol=FALSE,
                             pcols, verb=TRUE) {

    indat <- read.csv(file, header=FALSE,stringsAsFactors=FALSE,
                     sep=sep, dec=dec, skip=skip)

    ## data IDs are in column 1, followed by data matrices starting
    ## in column 2; skip internal result calculation
    didx <- c(which(indat[,1]!="" & indat[,1]!="Results"))
    ## filter only requested data
    dataIDs <- indat[didx,1]
    didx <- c(didx, which(indat[,1]=="Results")) # add "Results" index
    names(didx) <- c(dataIDs,"end")
    
    if ( !missing(data.ids) )
      dataIDs <- dataIDs[dataIDs%in%data.ids]

    data <- rep(list(NA),length(dataIDs))
    names(data) <- dataIDs
    for ( dataID in dataIDs  ) {

        ## get DATA:
        ## rows: the header is 2 rows after the ID, data starts 3 rows after
        ## and ends 2 rows before next
        i <- which(names(didx)==dataID)
        hidx <- didx[i]+2
        sidx <- didx[i]+3
        eidx <- didx[i+1]-2
        ## columns:
        ## time: column 2; temperature: column 3;
        ## data: columns 4 to last-1; last data set ends with nrow
        lcol <- ncol(indat)
        if ( skiplastcol  ) lcol <- lcol -1
        cols <- 4:lcol
        time <- as.numeric(strptime(indat[sidx:eidx,2],format <- "%H:%M:%S"))
        temp <-  as.numeric(sub(",",".",indat[sidx:eidx,3]))
        dat <- matrix(as.numeric(sub(",",".",unlist(indat[sidx:eidx,cols]))),
                      ncol=length(cols),nrow=length(sidx:eidx))
        ## get columns and change Temperature ID
        colnames(dat) <- as.character(indat[hidx,cols])
        ## check last columns (sometimes empty, sometimes not)
        emptycols <- which(colnames(dat)=="NA")
        if ( length(emptycols)>0 )
          dat <- dat[,-emptycols]
        data[[dataID]] <- list(time=time, temp=temp, data=dat)
    }

    ## since time here comes formatted, the current data is
    ## added -> subtract minimal time
    t0 <- min(unlist(lapply(data, function(x) x$time)))
    for ( i in 1:length(data) ) 
      data[[i]]$time <- data[[i]]$time - t0
    
    ## add colors
    ## TODO: use these in plots
    ## TODO: check passed pcols
    #if ( missing(pcols) ) 
    #    data$colors <- getColors(dataIDs)
    #else data$colors <- pcols

    ## SET DATA ID 
    data$dataIDs <- names(data)
    
    ## SET GLOBAL TIME & TEMPERATURE by INTERPOLATION:
    ## interpolate data: this adds a master time and temperature
    ## and interpolates all data to this time; if this step
    ## is omitted, there will be no global master time!
    if ( interpolate )
      data <- interpolatePlateTimes(data, verb=verb)
    data

}


## interpolates all data to a master time
## and reduced redundant temperature information
## TODO: correct times for reading delay in platereader
##       - perhaps newer versions can give exact time for each
#' Read BMG Optima/MARS-exported files
#' @inheritParams readPlateData
#' @seealso \code{\link{readPlateData}}
#' @export
readBMGPlate <- function(files, data.ids, interpolate=TRUE,
                         skip=5, sep=";", dec=".", verb=TRUE, pcols) {

    data <- list()
    ## 1) PARSE ALL DATA FILES and collect the individual measurements
    for ( i in 1:length(files) ) {
        file <- files[i]
        #id <- names(files)[i]
        if ( verb )
          cat(paste("Parsing file", file , "\n"))
        dat <- read.csv(file,header=FALSE,stringsAsFactors=FALSE,
                        skip=skip,sep=sep)

        ## last row in BMG files is usually empty, remove
        if ( sum(dat[nrow(dat),]=="")==ncol(dat) )
          dat <- dat[2:nrow(dat)-1,]

        ## convert well/row info from rows 1:2 to column names
        colnames(dat) <- paste(dat[1,],dat[2,],sep="")

        ## columns 1 and 2 are data type and measurement time
        ## and ID is in row 3
        colnames(dat)[1:2] <- as.character(dat[3,1:2])

        ## get sample/blank information from row 3, substituting "Sample "
        ## also converts to character
        samples <- sub("Blank ","",sub("Sample ","",dat[3,3:ncol(dat)]))

        ## remove first three rows from which info has been parsed
        dat <- dat[4:nrow(dat),]
        
        ## GET ALL DATA

        ## data type identified by first column
        ## substituting common stuff "Raw Data (XXX \d)"
        dtype <- sub("\\)","",sub(" .*", "",sub("Raw Data \\(", "", dat[,1])))
        ## the rest is numeric information, convert now
        dat <- data.matrix(dat[,2:ncol(dat)])

        ## collect all data
        types <- unique(dtype)
        dlst <- rep(list(NA),length(types))
        names(dlst) <- types
        for ( dtyp in types ) {
            if ( verb )
              cat(paste("\tfound data", dtyp, "\n"))
            tdat <- dat[dtype==dtyp,]
            dlst[[dtyp]] <- list(time=tdat[,1],
                                 data=tdat[,2:ncol(tdat)])            
        }

        ## handle temperature
        ## BMG gives T for each well, but actually it's the same for all
        ## add this temperature time-course to each other data
        ## and remove from list
        ## NOTE: temperature is not exported after fusing wells in BMG Mars!
        ## TODO: find out which time is used when fusing wells in BMG Mars!
        tidx <- which(names(dlst)=="Temperature")
        if ( length(tidx)>1 )
          warning(paste("BMG Temperature format has changed, multiple entries",
                        "in one data file", file, "\n\t",
                        "Perhaps check validity of code\n"))
        for ( i in tidx ) {
            temp <- dlst[[i]]$data
            if ( any(apply(temp,1,sd)>0 ) ) {
                warning(paste("BMG format changed; different temperatures for",
                              "different wells noted. UPDATE R CODE\n"))
            }
            temp <- c(matrix(apply(temp,1,mean),ncol=1))
            ## TODO: check whether times are ok, but so far, BMG/Mars
            ## don't give different times
        }
        
        ## rm from data list
        dlst <- dlst[names(dlst)!="Temperature"]
        ## .. and add to each other data (if it was present
        if ( length(tidx)>1 )
            dlst <- lapply(dlst, function(x) {x$temp=temp; x} )

        ## append data
        data <- append(data,dlst)
    }

    ## add colors
    ## TODO: use these in plots
    ## TODO: check passed pcols
    #if ( missing(pcols) ) 
    #    data$colors <- getColors(ptypes)
    #else data$colors <- pcols

    ## NOTE: at this stage, data between different plate-readers
    ## should already look similar; each entry containing separate
    ## time and temperature vectors

    ## SET DATA ID 
    data$dataIDs <- names(data)
    
    ## SET GLOBAL TIME & TEMPERATURE by INTERPOLATION:
    ## interpolate data: this adds a master time and temperature
    ## and interpolates all data to this time; if this step
    ## is omitted, there will be no global master time!

    if ( interpolate )
      data <- interpolatePlateTimes(data, verb=verb)
    data
}


#' \code{\link{prettyData}} : set colors and order or filter the data set
#' @param dids a vector of data IDs, data will be filtered and sorted by this list; of the vector is named the IDs will be replaced by these names
#' @param colors a vector of plot colors as RGB strings, optionally already named by dataIDs 
#' @seealso \code{\link{readPlateData}}
#' @export
prettyData <- function(data, dids, colors,
                       mids=c(time="Time",temp="Temperature")) {

    ## dids: data sorting, if missing, take all
    if ( missing(dids) )
        dids <- data$dataIDs

    ## reduce mids
    mids <- mids[mids%in%names(data)]
    
    ## re-order: mids first and IDs last
    data$dataIDs <- dids
    data <- data[c(match(mids,names(data)),
                   match(dids,names(data)),
                   match("dataIDs",names(data)))]


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
        data$colors <- colors
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
        }
        if ( missing(processing) ) # add date as processing note
            processing <- date()
        data <- append(data, list(list(data=dat, processing=processing)))
        names(data) <- c(names(data)[2:length(data)-1],ID)
    }
    data
}


#' \code{\link{rmData}} : remove a data set
#' @param data the current platexpress data set
#' @param ID a vector of IDs of the data to be removed
#' @param dat the new data, must be a matrix akin to other data in the set
#' @seealso \code{\link{addData}}
#' @export
rmData <- function(data, ID) {
    if ( !ID%in%data$dataIDs )
        warning("\"",ID,"\" not found")
    data$dataIDs <- data$dataIDs[!data$dataIDs%in%ID] # rm ID
    if ( "colors" %in% names(data) ) # rm color
        data$colors <- data$colors[!names(data$colors)%in%ID]
    data[-which(names(data)%in%ID)] # rm data
}


#' \code{\link{getData}} : get a specific data set, returns a data matrix
#' @param data the current platexpress data set
#' @param ID the ID of the data to be obtained
#' @export
getData <- function(data, ID, type="data") {
    data[[ID]][[type]] # just return the current data or old versions
}

## TODO: cut data either by time, or by a chosen data set
###' @export
cutData <- function(data, rng, ID) {
    
}

## TODO: convert to groFit data format
#' @export
data2grofit <- function(data,dids,dose) {
    if ( missing(dose) )
        dose <- as.factor(rep(1,ncol(odn)))
    grdat <- data.frame(cbind(colnames(odn),
                              as.factor(rep("test",ncol(odn))),
                              dose),
                        t(data[[ID]]$data))
    grdat
}

#' \code{\link{skipWells}} set wells that should be skipped from all
#' analyses and plots to NA
#' @param skip a list of strings identifiying the wells to be skipped,
#' e.g. "B3" to skip the well in row B/column 3
#' @examples
#' raw <- skipWells(raw, skip="A9")
#' @export
skipWells <- function(data, skip) {
    for ( i in 1:length(data) ) {
        if ( !"data"%in%names(data[[i]]) ) next
        sk <- skip[skip%in%colnames(data[[i]]$data)]
        if ( length(sk)>0 )
          for ( j in 1:length(sk)) 
              data[[i]]$data[,sk[j]] <- NA
    }
    data
}

#' \code{\link{correctBlanks}} correct for blanks
#' @param data the data list to be blank-corrected
#' @param plate the plate layout where column "blanks" indicates which wells
#' are to be treated as blanks
#' @param dids IDs of the data which should be blank-corrected, all will be
#' blanked if missing
#' @param by a list of column IDs of the plate layout; separate blank
#' correction will be attempted for groups in these columns; each group
#' must have at least one specified blank associated
#' @examples
#' data(ap12)
#' data <- correctBlanks(data=ap12data, plate=ap12plate, by="strain")
#' @export
correctBlanks <- function(data, plate, dids, by,
                          mids=c(time="Time",temp="Temperature")) {

    ## TODO: correct by time point, eg. for fluorescence in ecoli_rfp_iptg_20160908

    ## start new data list
    corr <- data
  
    ## reduce matrix to requested data
    data <- data[data$dataIDs]
    ptypes <- names(data)
    if ( !missing(dids) ) # only use requested data 
      ptypes <- ptypes[ptypes%in%dids]
    if ( length(ptypes)==0 ) {
        cat(paste("no data to blank\n"))
        return()
    }else
      cat(paste("blanking", paste(ptypes,collapse=";"),"\n"))

    ## reducing result matrix to requested data
    #corr <- corr[names(corr)%in% c(mids,ptypes)]

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
        cat(paste("blanking", btyp, ":", length(dwells), "wells, using",
                  length(bwells), "blank wells\n"))
        for ( k in 1:length(ptypes) ) {
            ptyp <- ptypes[k]
            dat <- data[[ptyp]]$data
            blank <- median(c(dat[,bwells]),na.rm=TRUE)
            cat(paste("\t",ptyp,"- blank:",blank,"\n"))
            #cat(paste("\tdata wells:",paste(dwells,collapse=";"),"\n",
            #          "\tblank wells:",paste(bwells,collapse=";"),"\n"))
            corr[[ptyp]]$data[,c(dwells,bwells)] <-
                           dat[,c(dwells,bwells)] - blank
            corr[[ptyp]]$processing <- c(corr[[ptyp]]$processing,
                                         paste("blank-corrected by",btyp))
        }
    }
    corr
}

#' \code{\link{adjustBase}} adjust data to a minimal base
#' @details Adjusts data to a new mininum, this is useful for adjustment
#' of negative values after blank corrections
#' @param dids vector of ID strings for which base correction should be
#' executed
#' @param base the new minimum for the data, default is 0, but it could
#' e.g. be the OD used for inoculation
#' @param by TODO: choose specific groups via plate-designs
#' @return Returns `data' where all data sets or only those selected by option
#' dids where raised to a minimum level in 
#' @export
adjustBase <- function(data, dids, base=0) {

    if ( missing(dids) ) # only use requested data 
        dids <- data$dataIDs # use all
    
    for ( did in dids ) {
        dat <- data[[did]]$data
        for ( j in 1:ncol(dat) ) 
            dat[,j] <- dat[,j] - min(dat[,j],na.rm=TRUE) + base
        data[[did]]$data <- dat
        data[[did]]$processing <- c(data[[did]]$processing,
                                    paste("corrected to base",base))
    }
    data
}

## helper function to calculate an average value
## for a variable given in several list items,
## used for average (master) time and temperatures
## in interpolatePlateTimes()
listAverage <- function(lst, id) {

    ## reduce to entries that have <id>
    lst <- lst[unlist(lapply(lst, function(x) id%in%names(x)))]
    if ( length(lst)==0 ) return(NULL)
    ## collect values with the same ID for different data sets
    vals <- lapply(lst, function(x) x[[id]])
    ## check length of lists and find missing time-points
    if ( any(diff(unlist(lapply(vals, length)))!=0) ) {
        stop("data have different lengths: ", id)
        ## find point of discrepancy
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
#' @return returns a copy of the full data list with a master time and
#' temperature added at the top level
#' @export
interpolatePlateTimes <- function(data, verb=TRUE) {

    if ( verb )
      cat(paste("Interpolating all data to a single master time.\n"))
    
    ## 1) calculate average (MASTER) time
    mtime <- listAverage(data, "time")
    mtemp <- listAverage(data, "temp") 
    
    ## 2) interpolate all data to MASTER time
    for ( id in data$dataIDs ) {
        data[[id]]$orig <- data[[id]]$data
        mdat <- data[[id]]$data
        for ( j in 1:ncol(data[[id]]$data) ) {
            x <- data[[id]]$time
            y <- data[[id]]$data[,j]
            ## interpolate data, NOTE that rule=2 will fill the end points
            mdat[,j] <- approx(x=x,y=y,xout=mtime,rule=2)$y
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
    data$Temperature <- mtemp
    data
}


#' \code{\link{viewPlate}} plots all data in plate format
#' @param data the list of measurement data as provided by \code{\link{readPlateData}}
#' @param rows a list of strings/characters used as row ID in the composite
#' row:col well description in the plate layout (map) and plate data
#' @param cols as rows but plate column IDs
#' @param mids vector of named strings, indicating the IDs of the master
#' time and temperature vectors in the data
#' @param xid ID of a data-set in the input data that can be used as x-axis
#' instead of the default Time vector
#' @param dids IDs of the data to be plotted; if missing, all data will
#' be plotted
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
#' @examples
#' data(ap12data)
#' vp <- viewPlate(ap12data)
#' @export
viewPlate <- function(data,rows=toupper(letters[1:8]),cols=1:12,
                      mids=c(time="Time",temp="Temperature"), 
                      xid, xscale=FALSE,xlim,
                      dids, pcols, yscale=TRUE,ylims,ylim,log="",
                      legpos="topleft") {

    ## which wells to plot?
    wells <-  paste(rep(rows,each=length(cols)),cols,sep="")
    
    ## get x-data: global data (time, temperature)
    ## or a data set specified via xid
    ## TODO: implement absence of master time, if interpolate=FALSE
    ## upon reading of data
    time <- data[[mids["time"]]]
    if ( mids["temp"] %in% names(data) )
      temp <- data[[mids["temp"]]]
    if ( !missing(xid) ) { # or use other data-set as x
        xdat <- data[[xid]]$data[,wells,drop=FALSE]
    } else xid <- NULL
    mids <- c(mids,xid) # will all be removed from plot data

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
      if ( !is.null(xid) ) 
        xlim <- range(c(xdat),na.rm=TRUE)
      else
        xlim <- range(time)
    ## set local ylims
    ## TODO: generalize and align with code in viewGroup
    if ( missing(ylim) ) {
        pidx <- 1:length(ptypes)
        ylim <- rep(list(c(NA,NA)),length(ptypes)) # initiliaze
        for ( k in pidx ) {
            dat <- data[[ptypes[k]]]$data[,wells,drop=FALSE] # get plotted wells
            if ( is.null(xid) ) # if x is time, limit to plot xlim
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
    par(mfcol=c(length(rows),length(cols)),mai=rep(0,4))
    for ( j in cols ) 
      for ( i in rows ) {
          well <- paste(i,j,sep="")
          ## x data
          if ( !is.null(xid) ) 
            x <- xdat[,well]
          else x <- time
          for ( k in 1:length(ptypes) ) {
              # y data
              ptyp <- ptypes[k]
              y <- data[[ptyp]]$data[,well]
              if ( k>1 ) par(new=TRUE)
              #cat(paste("hallo",k))
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
                  legend(legpos,well,bty="n")
              }
          }
      }
    ## TODO: return meaningful and/or non-plotted information
    ## assigning it makes it silent!
    if ( is.null(xid) ) xid <- "Time"
    plotparams <- list(ylims=ylims, xid=xid, xlim=xlim,  colors=pcols)
}


## TODO - repair example code
## @example
## data(ap12)
## groups <- getGroups(plate=ap12plate, by=c("strain"))
#' group wells by experiment annotations (in plate map file)
#' @param by a list of column IDs of the plate layout
#' @details Calculates the distinct groups from the plate layout by the selected
#' experimental parameters.
#' @return Returns a list of well IDs for the identified grouping. This list
#' can be used in viewGroups(data,groups) to summarize data for these groups.
#' @seealso \code{\link{readPlateMap}}, \code{\link{viewGroups}}
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

## TODO - repair example
## @example
## groups <- getGroups(plate=plate, by=c("strain"))
## vg <- viewGroups(data,groups=groups,lwd.orig=0.1,nrow=3)

#' plot grouped wells as summary plots, incl. confidence intervals and means
#' @param data the list of measurement data as provided by \code{\link{readPlateData}}
#' @param groups a list of well grouped wells, as produced by \code{\link{getGroups}}
#' @param nrows number of plot rows
#' @param mids vector of named strings, indicating the IDs of the master
#' time and temperature vectors in the data
#' @param xid ID of a data-set in the input data that can be used as x-axis
#' instead of the default Time vector
#' @param dids IDs of the data to be plotted; if missing, all data will
#' be plotted
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
#' @param lwd.orig line-width of the original single data, set to 0 to supress plotting of all original data
#' @param lty.orig line type of the original single data, set to 0 to supress plotting of all original data
#' @param mai set the outer margins around plot areas, see ?par
#' @param mgp set the position of axis title, tick marks and tick lengths
#' @param xaxis plot x-axis if TRUE
#' @param yaxis plot y-axis if TRUE
#' @seealso \code{\link{viewPlate}}, \code{\link{getGroups}}, \code{\link{readPlateMap}}
#' @examples
#' data(ap12)
#' groups <- getGroups(plate=ap12plate, by=c("strain"))
#' vg <- viewGroups(ap12data,groups=groups,lwd.orig=0.1,nrow=1)
#' @export
viewGroups <- function(data, groups,
                       mids=c(time="Time", temp="Temperature"), 
                       xid, xscale=FALSE, xlim,
                       dids, pcols, yscale=TRUE, ylims, ylim, log="",
                       legpos="topleft", lty.orig=1,lwd.orig=0.1,
                       mai=c(0.5,0,0,0), mgp=c(1.5,.5,0),
                       nrow=1, xaxis=TRUE, yaxis=c(1,2)) {
    

    wells <- unique(unlist(groups))
    
    ## get master data: time and temperature
    ## TODO: interpolate data here on the fly, if not done
    ## upon parsing data or subsequently
    time <- data[[mids["time"]]]
    xlab <- mids["time"]
    if ( mids["temp"] %in% names(data) )
      temp <- data[[mids["temp"]]]
    if ( !missing(xid) ) {
        xdat <- data[[xid]]$data[,wells,drop=FALSE]
        xlab <- xid
    } else xid <- NULL
    mids <- c(mids,xid) # will all be removed from plot data

    ## get plot params - colors
    if ( missing(pcols) ) # if not missing: overrides set&default colors
        if ( !is.null(data$colors) )
            pcols <- data$colors # set colors
        else
            pcols <- getColors(data$dataIDs) # default colors

    ## reduce data to plotted data
    data <- data[data$dataIDs]
    #data <- data[!names(data)%in%mids]
    ptypes <- names(data)
    if ( !missing(dids) ) # only use requested data for plots
      ptypes <- ptypes[ptypes%in%dids]
    if ( length(ptypes)==0 ) {
        cat(paste("no data to plot\n"))
        return()
    }else
      cat(paste("plotting", paste(ptypes,collapse=";"),"\n"))
    
    ## plot params
    if ( missing(xlim) ) 
      if ( !is.null(xid) ) 
        xlim <- range(c(xdat),na.rm=TRUE)
      else
        xlim <- range(time)
    ## set local ylims
    ## TODO: generalize and align with code in viewPlate
    if ( missing(ylims) & missing(ylim) ) {
        ylims <- list()
        for ( k in 1:length(ptypes) ) {
            dat <- data[[ptypes[k]]]$data[,wells,drop=FALSE] # get plotted wells
            if ( is.null(xid) ) { 
                dat <- dat[time>=xlim[1]&time<=xlim[2],,drop=FALSE]
            } else {
                dt <- NULL
                for ( i in 1:ncol(dat) ) 
                    dt <- c(dt,dat[xdat[,i]>=xlim[1]&xdat[,i]<=xlim[2],i])
                dat <- dt
            }
            ylim <- range(c(dat[is.finite(dat)]),na.rm=TRUE)
            ylims <- append(ylims,list(ylim))
        }
        names(ylims) <- ptypes
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

    if ( length(groups)>1 ) {
        ncol <- ceiling(length(groups)/nrow)
        par(mfcol=c(nrow,ncol))
    }
    par(mai=mai,mgp=mgp)
    for ( g in 1:length(groups) ) {
        wells <- groups[[g]]
        id <- names(groups)[g]
        ## x data other then time
        if ( !is.null(xid) ) 
          x <- xdat[,wells]
        else x <- time
        for ( i in 1:length(ptypes) ) {
            ptyp <- ptypes[i]
            ## get data for selected wells
            dat <- data[[ptyp]]$data[,wells]
            ## calculate stats only for common x!
            ## TODO: instead bin data on x and calculate ci there
            ## or interpolate data to common x (on the fly)?
            ## TODO: do we get NAs or empty vals from ci?
            if ( is.null(dim(x)) ) { 
                mn <- apply(dat,1,function(x) mean(x,na.rm=TRUE))
                ci <- apply(dat,1,function(x) ci95(x,na.rm=TRUE))
            }
            ## PLOT
            par(new=i!=1)
            ## override lty.orig=0 and lwd.orig=0 if x is data-specific
            if ( !is.null(dim(x)) ) {
                if ( lwd.orig==0 ) lwd.orig <- 0.1
                if ( lty.orig==0 ) lty.orig <- 1
            }
            ## override color to allow lwd.orig=0 to work for PDF as well
            tmp <- ifelse(lwd.orig==0,NA, pcols[ptyp])
            matplot(x,dat,type="l",lty=lty.orig,lwd=lwd.orig,axes=FALSE,
                    ylim=ylims[[ptyp]],col=tmp,xlim=xlim,xlab=xlab,log=log)
            
            ## plot mean and confidence intervals
            if ( is.null(dim(x)) ) { # only for common x!
                polygon(x=c(x,rev(x)),y=c(mn-ci,rev(mn+ci)),border=NA,
                    col=paste(pcols[ptyp],"55",sep=""))
                lines(x=x,mn,col=pcols[ptyp],lwd=2)
            }

            ## add axes for first two values
            if ( yaxis[1]==i ) axis(2, tcl=.25, mgp=c(0,-1,-.05),
                                    col=pcols[ptyp],
                                    col.axis=pcols[ptyp])
            if ( yaxis[2]==i ) axis(4, tcl=.25, mgp=c(0,-1,-.05),
                                     col=pcols[ptyp],
                                     col.axis=pcols[ptyp])

        }
        legend(legpos,id,bty="n")
        if ( xaxis ) axis(1)
    }
    ## TODO: return meaningful and/or non-plotted information
    ## assigning it makes it silent!
    if ( is.null(xid) ) xid <- "Time"
    plotparams <- list(ylims=ylims, xid=xid, xlim=xlim,  colors=pcols)
}


### COMMENTS FOR EXAMPLE DATA
  
#' ap12data: example data by Dennis Dienst and Alice Pawloski, incl. the
#' plate reader measurements of E.coli growth, expressing a fluorescent
#' proteins, in a Synergy Mx platereader  
#' 
#' \itemize{
#'   \item Data:
#'   \item Time: the interpolated 'master' time
#'   \item Temperature: the temperature time-course
#'   \item Data matrix '600': well absorbance at 600 nm, i.e., the OD,
#'   \item Data matrix 'YFP_50:500,535': the YFP fluorescence measured by excitation at 500 nm and emission at 535 nm 
#'   \item Plate Layout:
#'   \item The plate layout table indicates the different strains tested, biological replicates (B1 to B3), and blank wells (containing only growth medium) 
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ap12data
#' @usage data(ap12data)
#' @format a list of time-courses of absorbance and fluorescence data, read in by readPlateData("AP12.csv", type="Synergy", data.ids=c("600","YFP_50:500,535")) and readPlateMap("AP12_layout.csv", fields=c("strain","samples"))
#' @seealso \code{\link{readPlateData}} and \code{\link{readPlateMap}} 
NULL
