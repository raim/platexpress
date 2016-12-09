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
#' Parses a plate design file in CSV. Rows and 
#' columns should be named as in the corresponding data files.
#' @param file text file containing the plate layout.
#' @param sep column separator, as in read.table
#' @param fsep within-field separator, separating the well-specific descriptors
#' @param blank.id keyword that indicates blank wells. Blank wells can be
#'                 combined with other well descriptors for separate blanking
#' @param fields names for the field descriptors
#' @param formatted indicates whether the file is already in the required
#'                  format; all other paramaters but 'sep' will be ignored
#' @param nrows number of rows to expect, defaults to 8 for rows A to H in
#' a typical 96 well plate
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

    ## rm blank column
    plate <- plate[,colnames(plate)!=blank.id]
    ## and add new blank
    plate <- cbind(data.frame(plate),blank=blanks)

    #class(plate) <- "platemap"
    return(plate)
}

### PLATE DATA

## TODO: repair this in example
#
#' Parses data files in CSV, as exported by
#' the plate reader software. Header IDs in the data file should match with 
#' IDs in the plate map, see \code{\link{readPlateMap}}. Pre-defined read-in
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
#' @param time.conversion conversion factor for the plate time, e.g., 1/3600
#' to convert from hours to seconds
#' @param ... further parameters to plate-reader specific parsing functions
#' @note The original data is all interpolated to a common/average 'master' time
#' @return a list of distinct measurement time-courses from one plate
#' @seealso \code\link{readBMGPlate}}, \code{\link{readPlateMap}}, \code{\link{viewPlate}}
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @examples
#' data.file <- system.file("extdata", "AP12.csv", package = "platexpress")
#' raw <- readPlateData(files=data.file, type="Synergy", data.ids=c("600","YFP_50:500,535"), dec=",",time.format="%H:%M:%S")
#' @export
readPlateData <- function(files, type, data.ids, 
                          verb=TRUE,
                          interpolate=TRUE, time.conversion, ...) {

    if ( type=="BMG" )
      data <- readBMGPlate(files=files, data.ids=data.ids,
                           verb=verb, ...)
    else if ( type=="Synergy" )
      data <- readSynergyPlate(files=files, data.ids=data.ids, 
                               verb=verb, ...)
    else if ( type=="simple" ) # single data item in simple spreadsheet 
      data <- readSimplePlate(files=files, data.ids=data.ids,
                              verb=verb, ...)
      
    ## NOW PREPARE DATA
    ## SET GLOBAL TIME & TEMPERATURE by INTERPOLATION:
    ## interpolate data: this adds a master time and temperature
    ## and interpolates all data to this time; if this step
    ## is omitted, there will be no global master time!
    if ( interpolate )
      data <- interpolatePlateTimes(data, verb=verb)

    if ( !missing(time.conversion) )
        data$Time <- data$Time * time.conversion

    class(data) <- "platedata"
    data
}

#' Read in simple platereader data with only one data type, and time
#' in the first column
#' @param skip lines to skip before parsing data
#' @param sep field separator used in the input tabular file
#' @param dec decimal number symbol (e.g., ',' if data export was with
#' german language settings)
#' @inheritParams readPlateData
readSimplePlate <- function(files, data.ids, skip,
                            time.format="%M:%S", verb=TRUE) {

    if ( verb )
      cat(paste("Parsing file", files, "\n"))

    ## defaults
    if ( missing(skip) )
      skip <- 0
    if ( missing(data.ids) )
      data.ids <- "data"
 
    if ( verb  )
      cat(paste("\tloading data", data.ids, "\n"))
    
    dat <- read.csv(files, sep="\t", skip=skip)
    time <- as.numeric(strptime(dat[,1],format=time.format)) # time in seconds
    time <- time - time[1]
    dat <- dat[, -1]

    ## generate local data
    data <-  list()
    data[[data.ids]] <- list()
    data[[data.ids]]$time <- time
    data[[data.ids]]$data <- data.matrix(dat)

    ## set dataID
    data$dataIDs <- data.ids

    data
}

#' Read Synergy Mx-exported files
#' @param skip lines to skip before parsing data
#' @param time.format format of the time, e.g., "%H:%M"%S" (see
#' \code{strptime}), or "numeric" if the time is provided in normal numbers
#' @param sep field separator used in the input tabular file
#' @param dec decimal number symbol (e.g., ',' if data export was with
#' german language settings)
#' @param skiplastcol the last column is often empty, set to TRUE if
#' this is the case for the current data set
#' @inheritParams readPlateData
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{readPlateData}}
#' @export
readSynergyPlate <- function(files, data.ids,
                             skip, sep=";", dec=".", skiplastcol=FALSE,
                             time.format="%H:%M:%S",
                             verb=TRUE) {

    if ( verb )
      cat(paste("Parsing file", files, "\n"))

    ## defaults
    if ( missing(skip) )
      skip <- 58


    indat <- read.csv(files, header=FALSE,stringsAsFactors=FALSE,
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

        if ( verb  )
          cat(paste("\tloading data", dataID, "\n"))
        
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
        time <- indat[sidx:eidx,2]
        #cat(paste("time format", time.format, time[1], "\n"))
        if ( time.format!="numeric" )
          time <- as.numeric(strptime(time,format=time.format))
        else
          time <- as.numeric(sub(",",".",time))
         #cat(paste("time format", time.format, time[1], "\n"))

        temp <-  as.numeric(sub(",",".",indat[sidx:eidx,3]))
        ## convert german decimal ',' to '.'
        ## NOTE: align this with command-line argument 'dec'
        dat <- matrix(as.numeric(sub(",",".",unlist(indat[sidx:eidx,cols]))),
                      ncol=length(cols),nrow=length(sidx:eidx))
        ## get columns and change Temperature ID
        colnames(dat) <- as.character(indat[hidx,cols])
        ## check last columns (sometimes empty, sometimes not)
        emptycols <- which(colnames(dat)=="NA")
        if ( length(emptycols)>0 )
          dat <- dat[,-emptycols]
        ## check last rows (sometimes empty, sometimes not)
        ## -> either NA values in time or in all wells
        present <- (!is.na(time) &
                    ncol(dat) > apply(dat,1,function(x) sum(is.na(x))))
        
        dat  <- dat[ present,]
        temp <- temp[present]
        time <- time[present]
        data[[dataID]] <- list(time=time, temp=temp, data=dat)
    }

    ## since time here comes formatted, the current data is
    ## added -> subtract minimal time
    t0 <- min(unlist(lapply(data, function(x) x$time)),na.rm=TRUE)
    for ( i in 1:length(data) ) 
      data[[i]]$time <- data[[i]]$time - t0
    
    ## SET DATA ID 
    data$dataIDs <- names(data)

    data

}


## interpolates all data to a master time
## and reduced redundant temperature information
## TODO: correct times for reading delay in platereader
##       - perhaps newer versions can give exact time for each
#' Read BMG Optima/MARS-exported files
#' @param skip lines to skip before parsing data
#' @param sep field separator used in the input tabular file
#' @param dec decimal number symbol (e.g., ',' if data export was with
#' german language settings)
#' @inheritParams readPlateData
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{readPlateData}}
#' @export
readBMGPlate <- function(files, data.ids, 
                         skip, sep=";", dec=".", verb=TRUE) {

    if ( missing(skip) )
      skip <- 5
    
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
              cat(paste("\tloading data", dtyp, "\n"))
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

    ## NOTE: at this stage, data between different plate-readers
    ## should already look similar; each entry containing separate
    ## time and temperature vectors

    ## SET DATA ID 
    data$dataIDs <- names(data)

    data
}
