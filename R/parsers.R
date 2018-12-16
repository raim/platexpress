### PLATE DESIGN & DATA

### UTILS

## get the square-most matrix for a given number
tosquare <- function(x) {
    fct <-factorize(x)
    mid <- length(fct)/2
    if ( length(fct)%%2==0 )
        dim <- c(nrow=fct[mid], ncol=fct[mid+1])
    else
        dim <- c(nrow=fct[mid-.5],ncol=fct[mid+.5])
    dim
}

## https://stackoverflow.com/questions/6424856/r-function-for-returning-all-factors
factorize <- function(x) {
    x <- as.integer(x)
    div <- seq_len(abs(x))
    factors <- div[x %% div == 0L]
    return(factors)
}

### HIGH-LEVEL WRAPPER
#' High-level interface of \code{\link{platexpress}}
#'
#' This high-level wrapper allows to load a platereader experiment
#' and its layout file, do blank corrections, and replicate groupings.
#'
#' Calls, in this order the  \code{\link{platexpress}} functions
#'  \code{\link{readPlateData}},  \code{\link{readPlateMap}};
#'  optionally: \code{\link{skipWells}}, and \code{\link{correctBlanks}},
#' and \code{\link{getGroups}}. Please see the help files of these
#' functions for details on parameters.
#' @inheritParams readPlateData
#' @inheritParams readPlateMap
#' @inheritParams skipWells
#' @inheritParams skipWells
#' @inheritParams correctBlanks
#' @param layout the text file containing the plate layout information, see
#' argument \code{file} in \code{\link{readPlateMap}}
#' @param skip.wells skip these wells from all analyses, see argument
#' \code{skip} in \code{\link{skipWells}}
#' @param blank call \code{\link{correctBlanks}} to blank data
#' @param blank.data subset of data to blank, see argument \code{yids} in
#' function \code{\link{correctBlanks}}
#' @param base optional minimal value; all values will be raised by
#' the same amount using the function \code{\link{adjustBase}} via
#' function \code{\link{correctBlanks}} 
#' @param group1 plate layout column to be used for coarse grouping, via
#' function \code{\link{getGroups}}, eg. "strains"
#' @param group2 plate layout column to be used for fine grouping, via
#' function \code{\link{getGroups}}, eg. "inducer amount")
#' @param group2.color numeric plate layout column to be used for group2
#' coloring, via function \code{\link{groupColors}}, eg. "amount"
#' @export
readExperiment <- function(files, type,
                           time.range=c("full"), time.conversion=1/3600,
                           layout, sep=";",fsep=";",asep=":",blank.id="blank",
                           fields, afields,
                           skip.wells, blank=FALSE, blank.data, base=0,
                           group1, group2, group2.color) {

    ## parse raw data
    raw <- readPlateData(files, type=type,
                         time.range=time.range, time.conversion=time.conversion,
                         verb=TRUE)


    ## read layout file
    map <- readPlateMap(layout, 
                        fields=fields, blank.id=blank.id,
                        sep=sep, fsep=fsep, asep=asep,
                        afields=afields)

    ## skip wells
    ## TODO: implement in skipWells and do after map is attached
    if ( !missing(skip.wells) ) {
        map <- skipWells(map, skip.wells)
        raw <- skipWells(raw, skip.wells)
    }
    
    
    ## correct blanks
    if ( blank | !missing(blank.data) ) {
        if ( missing(blank.data) )
            blank.data <- raw$dataIDs
        dat <- correctBlanks(raw, map, yids = blank.data, base = base)
    } else dat <- raw
    rm(raw)
    
    ## get groupings
    dat$groups <- NULL
    if ( !missing(group1) )
        dat$groups$group1 <- getGroups(map, by = group1, verb=FALSE)
    if ( !missing(group2) ) {
        dat$groups$group2 <- getGroups(map, by = group2, verb=FALSE)
        if ( !missing(group2.color) )
            dat$groups$group2.color <- groupColors(map, dat$groups$group2,
                                                   color = group2.color)
    }
    
    ## attach all layout
    dat$layout <- map
    
    dat
    
}


### TODO FUNCTIONS

## TODO: wrapper to read both layout and data
readPlates <- function(layout, data.path, type, settings) {
    ## read plate layout, potentially multiple plates in one file
    ## automatically read in each plate data
    ## merge data from multiple plates?
}
## TODO: merge several plates and layouts
## just generate lists of plate objects, and make viewGroups/viewPlate
## etc. take such lists as defaults, and convert singles to lists themselves
mergePlates <- function(data=list(), layouts=list()) {
    #just generates lists of plates, and one big layout table
}

### LAYOUT PARSERS

## Read Plate Layout Map
## parses a plate design file in CSV. Rows and columns
## should be named as in the datafile. Each field can
## have multiple descriptors, separated by a specified separator (e.g.
## "newline"); blanks are specified by a keyword (default: "blank"),
## and if separate blanks are used for different conditions, the
## blank field must have the same format as measurement fields, except
## one parameter replaced by the keyword. Names for the values
## in the fields can be passed via argument "fields" as a vector of
## strings.
## TODO: repair this in example:
#' Read Plate Layout Map
#'
#' Parses a plate design file in CSV format. Rows and
#' columns should be named as in the corresponding plate reader data files.
#' @param file text file containing the plate layout.
#' @param sep column separator, as in \code{\link[utils:read.table]{read.table}}
#' @param fsep within-field separator, separating the well-specific descriptors within
#' well fields
#' @param blank.id keyword that indicates blank wells. Blank wells can be
#'                 combined with other well descriptors for separate blanking
#' @param fields names for the field descriptors
#' @param asep a separator for substance:amount pairs, eg. to indicate amount
#' of an inducer or a nutrient, can only be used together with
#' argument \code{afields}, eg. use \code{asep=":", afields="inducer"}
#' @param afields field names which hold substance:amount pair information,
#' see argument \code{asep}
#' @param formatted indicates whether the file is already in the required
#'                  format; all other paramaters but 'sep' will be ignored
#' @param nrows number of rows to expect, defaults to 8 for rows A to H in
#' a typical 96 well plate; TODO: get rid of this and instead check for
#' rownames?
#' @param header logical argument indicating the presence/absence of a header
#' row in the layout file
#' @return a table of well content descriptors, where the first column 'wells'
#'         maps the plate map to the data files.
#' @seealso \code{\link{readPlateData}}
#' @examples
#' plate.file <- system.file("extdata", "AP12_layout.csv", package = "platexpress")
#' plate <- readPlateMap(file=plate.file, blank.id="blank",fsep="\n", fields=c("strain","samples"))
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
readPlateMap.old <- function(file, sep="\t",
                             fsep="\n", fields, asep, afields,
                             blank.id="blank",
                             nrows=8, formatted=FALSE, header=TRUE) {


    ## already in well format?
    if ( formatted ) {
        if ( length(grep("\\.xlsx?$", file)) ) {
            tmp <- data.frame(readxl::read_excel(file, col_names=header),
                              stringsAsFactors=FALSE)
        } else {
            dat <- read.table(file, sep=sep, header=header,
                              stringsAsFactors=FALSE)
        }
        ## split inducer:amount columns
        if ( !missing(asep) & !missing(afields) ) {
            dat <- replaceAmounts(dat, col=afields, sep=asep, replace=TRUE)
            dat <- amountColors(dat)
        }
        return(dat)
    }

    ## plate map by columns and rows
    ## TODO: get rid of nrows and check fields instead?
    if ( length(grep("\\.xlsx?$", file)) ) {
        tmp <- data.frame(readxl::read_excel(file, col_names=header),
                          stringsAsFactors=FALSE)
        rownames(tmp) <- tmp[,1]
        dat <- tmp[1:nrows,-1]
    } else {
        dat <- read.table(file, sep=sep,
                          header=header, row.names=1, nrows=nrows,
                          stringsAsFactors=FALSE)
    }
    if ( !header )
      colnames(dat) <- 1:ncol(dat)

    ## generate well names from row and column names
    plate <- paste(rep(rownames(dat),ncol(dat)),
                   rep(sub("X","",colnames(dat)),each=nrow(dat)),sep="")

    ## parse field values
    vals <- strsplit(unlist(dat),fsep)
    nvals <- max(unlist(lapply(vals,length)))

    values <- matrix("",nrow=length(vals),ncol=nvals)
    for ( i in 1:nvals )
      values[,i] <- as.character(unlist(lapply(vals, function(x) trimws(x[i]))))

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

    ## split inducer:amount columns
    if ( !missing(asep) & !missing(afields) ) {
        plate <- replaceAmounts(plate, col=afields, sep=asep, replace=TRUE)
        plate <- amountColors(plate)
    }

    ## use wells as rownames!
    rownames(plate) <- plate$well

    #class(plate) <- "platemap"
    return(plate)
}

## version 2 of readplatemap
## 1) parse xls/odt files
## 2) use tidyr/dplyr to parse fields - DONE
## Read Plate Layout Map
## parses a plate design file in CSV. Rows and columns
## should be named as in the datafile. Each field can
## have multiple descriptors, separated by a specified separator (e.g.
## "newline"); blanks are specified by a keyword (default: "blank"),
## and if separate blanks are used for different conditions, the
## blank field must have the same format as measurement fields, except
## one parameter replaced by the keyword. Names for the values
## in the fields can be passed via argument "fields" as a vector of
## strings.
## TODO: repair this in example:
#' Read Plate Layout Map
#'
#' Parses a plate design file in CSV format. Rows and
#' columns should be named as in the corresponding plate reader data files.
#' TODO: causes empty groups in getGroups; because separate causes empty
#' strings in first field instead of NA
#' @param file text file containing the plate layout.
#' @param sep column separator, as in \code{\link[utils:read.table]{read.table}}
#' @param fsep within-field separator, separating the well-specific descriptors within
#' well fields
#' @param blank.id keyword that indicates blank wells. Blank wells can be
#'                 combined with other well descriptors for separate blanking
#' @param fields names for the field descriptors
#' @param asep a separator for substance:amount pairs, eg. to indicate amount
#' of an inducer or a nutrient, can only be used together with
#' argument \code{afields}, eg. use \code{asep=":", afields="inducer"}
#' @param afields field names which hold substance:amount pair information,
#' see argument \code{asep}
#' @param formatted indicates whether the file is already in the required
#'                  format; all other paramaters but 'sep' will be ignored
#' @param nrows number of rows to expect, defaults to 8 for rows A to H in
#' a typical 96 well plate; TODO: get rid of this and instead check for
#' rownames?
#' @param header logical argument indicating the presence/absence of a header
#' row in the layout file
#' @return a table of well content descriptors, where the first column 'wells'
#'         maps the plate map to the data files.
#' @seealso \code{\link{readPlateData}}
#' @examples
#' plate.file <- system.file("extdata", "AP12_layout.csv", package = "platexpress")
#' plate <- readPlateMap(file=plate.file, blank.id="blank",fsep="\n", fields=c("strain","samples"))
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
readPlateMap <- function(file, sep="\t",
                          fsep="\n", fields, asep, afields,
                          blank.id="blank",
                          nrows=8, formatted=FALSE, header=TRUE) {

    ## already in well format?
    if ( formatted ) {
        if ( length(grep("\\.xlsx?$", file)) ) {
            tmp <- data.frame(readxl::read_excel(file, col_names=header),
                              stringsAsFactors=FALSE)
        } else {
            dat <- read.table(file, sep=sep, header=header,
                              stringsAsFactors=FALSE)
        }
        ## split inducer:amount columns
        if ( !missing(asep) & !missing(afields) ) {
          #if ( length(afields)>1 )
          #  warning("multiple substance:amount pairs not fully implemented;",
          #          "please add colors manually")
          for ( i in 1:length(afields) ) {
            afield <- afields[i]
            amnt.name <- ifelse(i>1,paste0(afield,".amount"),"amount")
            vdat <- tidyr::separate(vdat, col = afield,
                                    into = c(afield, amnt.name),
                                    sep = asep)
            vdat[[amnt.name]] <- as.numeric(vdat[[amnt.name]])
            vdat <- amountColors(vdat,
                                 substance=afield,
                                 amount=amnt.name,
                                 color=sub("amount","color",amnt.name))
          }
        }
    }

    ## plate map by columns and rows
    ## TODO: get rid of nrows and check fields instead?
    if ( length(grep("\\.xlsx?$", file)) ) {
        tmp <- data.frame(readxl::read_excel(file, col_names=header),
                          stringsAsFactors=FALSE)
        rownames(tmp) <- tmp[,1]
        dat <- tmp[1:nrows,-1]
    } else {
        dat <- read.table(file, sep=sep, header=header,
                          row.names=1, nrows=nrows,
                          stringsAsFactors=FALSE)
    }
    if ( !header )
        colnames(dat) <- 1:ncol(dat)

    ## generate well names from row and column names
    wells <- paste(rep(rownames(dat),ncol(dat)),
                   rep(sub("X","",colnames(dat)),each=nrow(dat)),sep="")

    ## convert to 1D
    vdat <- data.frame(well=wells,field=unlist(dat))
    ## get blanks, to add later
    blank.idx <- grep(blank.id, vdat$field)
    ## split fields - using tidyr separate
    vdat <- tidyr::separate(vdat, col = "field", into = fields, sep = fsep)
    if ( !missing(asep) & !missing(afields) ) {
      #if ( length(afields)>1 )
      #  warning("multiple substance:amount pairs not fully implemented;",
      #          "please add colors manually")
      for ( i in 1:length(afields) ) {
        afield <- afields[i]
        amnt.name <- ifelse(i>1,paste0(afield,".amount"),"amount")
        vdat <- tidyr::separate(vdat, col = afield,
                                into = c(afield, amnt.name),
                                sep = asep)
        vdat[[amnt.name]] <- as.numeric(vdat[[amnt.name]])
        vdat <- amountColors(vdat,
                             substance=afield,
                             amount=amnt.name,
                             color=sub("amount","color",amnt.name))
        ## TODO: split different substances in one afield
        ## into different columns!
        ## and replace column name
      }
    }
    ## replace all empty strings by NA
    vdat[vdat==""] <- NA
    ## add blanks
    vdat$blank <- FALSE
    vdat$blank[blank.idx] <- TRUE

    ## don't - use well column throughout!
    ##rownames(vdat) <- wells
    #class(vdat) <- "platemap"
    vdat
}

#' Numeric Amounts
#'
#' replaces amount strings of the form \code{substance:amount} in a column of
#' the plate layout map by two columns, separating substance and the amount
#' into a string and a numeric column.
#' @param map the plate map, see  \code{\link{readPlateMap}}
#' @param col the column ID or index to be split
#' @param sep the separator of the substance:amount pair strings
#' @param replace if TRUE the original columns in \code{col} will be removed
#' @seealso \code{\link{parseAmounts}}, \code{\link{readPlateMap}}
#' @export
replaceAmounts <- function(map, col, sep=":", replace=TRUE) {
    for ( cl in col )
        map <- cbind.data.frame(map, parseAmounts(map[,cl], sep=sep))
    if ( replace ) {
        if ( typeof(col)=="character" )
            col <- match(col,colnames(map))
        map <- map[,-col]
    }
    map
}

#' Parse Amount Strings
#'
#' splits a vector of strings of substance:amount pairs into
#' a matrix with 2 columns, of substances and amounts
#' @param str a vector of strings, providing substance:amount pairs
#' @param sep the separator of the substance:amount pair strings
#' @seealso \code{\link{replaceAmounts}}
#' @export
parseAmounts <- function(str, sep=":") {
    amount <- as.numeric(sub(paste0(".*",sep),"",str))
    inducer <- sub(paste0(sep,".*"),"",str)
    cbind.data.frame(substance=inducer,amount=amount)
}

#' Add Amount Colors
#'
#' add a color palette to the plate layout map, with colors along the range
#' of added amounts of a given substance. Substance and amount columns are eg.
#' auto-generated by \code{\link{readPlateMap}} with options \code{asep}
#' and \code{afields}.
#' @param map the plate map, see  \code{\link{readPlateMap}}
#' @param substance name of the substance column in \code{map}
#' @param amount name of the amount column in \code{map}
#' @param color name of the color column to be added to \code{map}
#' @param colf color palette function, default is colorRamp's
#' \code{\link[colorRamps:matlab.like2]{matlab.like2}}
#' @seealso \code{\link{readPlateMap}}, \code{\link{replaceAmounts}}
#' @export
amountColors <- function(map, substance="substance",amount="amount", color="color", colf=colorRamps::matlab.like2) {
    substances <- unique(map[,substance])
    colors <- rep(NA, nrow(map))
    for ( subst in substances ) {
        idx <- which(map[,substance]==subst)
        amnt <- as.numeric(map[idx,amount])
        amnt <- round(100*amnt/max(amnt,na.rm=TRUE))+1
        colors[idx] <- colf(101)[amnt]
    }
    if ( color %in% colnames(map) )
      map[,color] <- as.character(colors)
    else {
        map <- cbind.data.frame(map,
                                as.character(colors),
                                stringsAsFactors=FALSE)
        colnames(map)[ncol(map)] <- color
    }
    map
}

### PLATE DATA PARSERS

## TODO: repair this in example
#' Read Plate Data
#'
#' Parses platereader data files in CSV format, as exported by
#' the plate reader software. Header IDs in the data file should match with
#' IDs in the plate map, see \code{\link{readPlateMap}}. Pre-defined parsing
#' functions exist for a couple of plate-readers. If your format is
#' not implement, you can manually create simple data tables and use
#' the function \code{\link{readSimplePlate}}
#' @param files list of one or more data files
#' @param type pre-defined formats, as exported from platereaders; currently
#' for BMG Optima and Mars v3.01 R2, ('BMG'), BMG Clariostar and Mars
#' vXXX ('BMG2') and Biotek Synergy Mx ('Synergy'), BioLector ('BioLector'), BioLector Pro ('BioLectorPro');  TODO: define export
#' protocols!
#' @param data.ids an optional sub-selection of data types in the input file,
#' as a list of strings
#' @param interpolate if true a master time, the average time between distinct
#' measurements of one timepoint, is calculated and all values are interpolated
#' to this mastertime. This is currently obligatory for further processing.
#' See function \code{\link{interpolatePlateTimes}} for details.
#' @param time.range "common" requires that all data must have
#' the same number of time-points, "full" interpolates
#' the maximal available time range, using the mean
#' time step of the mean of data-specific time vectors
#' @param time.conversion conversion factor for the plate time, e.g., 1/3600
#' to convert from hours to seconds
#' @param verb print messages if true
#' @param ... further parameters to plate-reader specific parsing functions
#' @note The original data is all interpolated to a common/average 'master' time
#' @return a list of distinct measurement time-courses from one plate
#' @seealso \code{\link{readSimplePlate}}, \code{\link{readBMGPlate}}, \code{\link{readSynergyPlate}}, \code{\link{readPlateMap}}, \code{\link{viewPlate}}
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @examples
#' data.file <- system.file("extdata", "AP12.csv", package = "platexpress")
#' raw <- readPlateData(files=data.file, type="Synergy",
#'                      data.ids=c("600","YFP_50:500,535"),
#'                      dec=",", time.format="%H:%M:%S")
#' @export
readPlateData <- function(files, type, data.ids, 
                          interpolate=TRUE, time.range=c("common","full"),
                          time.conversion, verb=TRUE, ...) {

    if ( type=="BMG" )
      data <- readBMGPlate(files=files, data.ids=data.ids,
                           verb=verb, ...)
    else if ( type=="BMG2" )
      data <- readBMG2Plate(files=files, data.ids=data.ids,
                            verb=verb, ...)
    else if ( type=="Synergy" )
      data <- readSynergyPlate(files=files, data.ids=data.ids,
                               verb=verb, ...)
    else if ( type=="BioLector" ) {
      data <- readBioLectorPlate(files=files, data.ids=data.ids,
                                 verb=verb, ...)
      interpolate <- FALSE
    } else if ( type=="BioLectorPro" ) {
      data <- readBioLectorProPlate(files=files, data.ids=data.ids,
                                    verb=verb, ...)
      interpolate <- FALSE
    } else if ( type=="simple" ) {
        ## single data item in simple spreadsheet
        data <- readSimplePlate(files=files, data.ids=data.ids,
                                verb=verb, ...)
        interpolate <- FALSE
    }
    
    ## NOW PREPARE DATA
    ## SET GLOBAL TIME & TEMPERATURE by INTERPOLATION:
    ## interpolate data: this adds a master time and temperature
    ## and interpolates all data to this time; if this step
    ## is omitted, there will be no global master time!
    if ( interpolate )
      data <- interpolatePlateTimes(data, time.range=time.range, verb=verb)

    if ( !missing(time.conversion) )
        data$Time <- data$Time * time.conversion

    if ( verb )
        cat(paste0("Successfully parsed data:\n",
                  "\ttype: ", type, "\n",
                  ifelse(interpolate,
                         paste0("\ttime.range rule:", time.range, "\n"),""),
                  "\n"))

    class(data) <- "platedata"
    data
}

#' Read Simple Plate Data Tables
#'
#' Read in simple platereader data with only one data type, time
#' in the first column and data for all wells in the following columns,
#' where the column header must correspond to wells
#' in \code{\link{readPlateMap}}.
#' @param nrow organize data into a matrix with this number of rows
#' @param time.format format of the time, see \code{\link[base:strptime]{strptime}}); only used if time column is not numeric!
#' @param ... arguments passed to \code{\link{read.csv}}
#' @inheritParams readPlateData
readSimplePlate <- function(files, data.ids, 
                            nrow, time.format="%M:%S", verb=TRUE, ...) {

    if ( verb )
      cat(paste("Parsing file:", files, "\n"))

    ## defaults
    if ( missing(data.ids) )
      data.ids <- "data"

    if ( verb  )
      cat(paste0("\tloading data: \"", data.ids, "\"\n"))

    ## parse file and get data
    dat <- read.csv(files, ...)

    ## split into time (first column) and well data
    time <- dat[,1]
    ## attempt format time, if not numeric
    if (  !typeof(dat[,1])%in%c("double","integer") )
        time <- as.numeric(strptime(time,format=time.format)) # time in seconds
    time <- time - time[1]
    dat <- dat[, -1]

    ## for convenience replace leading 0, e.g. A07 -> A7
    ## in column namnes
    ## TODO: generalize somehow? 
    colnames(dat) <- gsub("([A-Z])0(.*)", "\\1\\2", colnames(dat))
    
    ## generate local data
    data <-  list()
    data[[data.ids]] <- list()
    data[[data.ids]]$data <- data.matrix(dat)

    ## set dataID
    data$dataIDs <- data.ids

    ## generate wells: factorize and take middle factors
    N <- ncol(dat)
    if ( missing(nrow) ) {
        dim <- tosquare(N)
        nrow <- dim[1]
        ncol <- dim[2]
    }
    if ( !exists("ncol", mode="numeric" ) ) ncol <- N/nrow
    cat(paste(nrow, N, "\n"))
    
    data$wells$plate <- t(matrix(colnames(data[[data$dataIDs[1]]]$data),
                                 nrow=ncol, ncol=nrow))

    ## ROWS: are used in viewPlate
    ## TODO: avoid to force letters for rows and numbers for columns?
    data$wells$rows <- toupper(letters[1:nrow])
    data$wells$cols <- as.character(1:ncol)

    ## or alternatively from well names, assuming rows in
    ## first letter, and columns in rest!?
    if ( FALSE ) {
        wells <- colnames(data[[data$dataIDs[1]]]$data)
        rows <- do.call(rbind,strsplit(wells,""))[,1]
        cols <- unlist(lapply(strsplit(wells,""),
                              function(x) paste0(x[2:length(x)])))
    }
    
    ## add global time and xids
    data$Time <- time
    data$xids <- "Time"
    
    class(data) <- "platedata"
    data
}


## TODO: get rid of german format
## TODO: what's wrong with GFP?
## TODO: get temperature and humidity
## headerline <- 23 ## line of data header in BioLector result file
readBioLectorPlate <- function(files, data.ids,
                               headerline=23, hnrw=20, sep=";", dec=",",
                               replace.zero=TRUE,
                               verb=TRUE) {


    if ( verb )
      cat(paste("Parsing file:", files, "\n"))

    ## parse header to get experiment information
    filt <- read.csv(files, nrow=hnrw, sep = sep, dec = dec,
                     stringsAsFactors=FALSE, header=FALSE, fill=TRUE)

    nrw <- as.numeric(filt[filt[,1]=="MTP ROWS",2])
    ncl <- as.numeric(filt[filt[,1]=="MTP COLUMNS",2])
    nfl <- as.numeric(filt[filt[,1]=="FILTERSETS",2])
    
    skip <- which(filt[,1]=="FILTERSET") -2

    ## parse header again to get filter information
    filters <- read.csv(files, skip=skip, nrow=nfl, sep = sep, dec = dec,
                         stringsAsFactors=FALSE, header=TRUE)
              


    dat <- read.csv(files, skip=headerline, sep = sep, dec = dec,
                    stringsAsFactors=FALSE, fill=TRUE)

    fidx <- grep("READING", colnames(dat)) # column with filter info

    ## get time
    tidx <- grep("TIME", dat[,fidx])
    ## convert time to seconds, SI and in-line with other formats
    time <- unlist(dat[tidx, (fidx+1):ncol(dat)])*3600
    
    ## get calibrated rows
    didx <- grep("Cal.", dat[,fidx])
    data <- as.data.frame(t(dat[didx,fidx:ncol(dat)]), stringsAsFactors=FALSE)
    colnames(data) <- dat[didx,1]

    data <- list()
    for ( i in 1:nrow(filters) ) {
        if ( verb ) 
            cat(paste0("\tloading data: \"", trimws(filters[i,2]), "\""))
        didx <- grep(paste0("Cal..*FS=",filters[i,1]), dat[,fidx])
        ## take raw values if no calibration exists
        if ( length(didx)==0 ) {
            didx <- which(dat[,fidx]==filters[i,1])
            if ( verb) cat(" - RAW VALUES")
        }
        dt <- t(dat[didx,(fidx+1):ncol(dat)])
        colnames(dt) <- dat[didx,1]
        data[[i]] <- list()
        data[[i]]$time <- time
        data[[i]]$data <- as.matrix(dt)
        if ( replace.zero )
            colnames(data[[i]]$data) <- sub("0","", colnames(data[[i]]$data))
                                            
        if ( verb ) cat("\n")
    }
    names(data) <- trimws(filters[,2])
    
    dataIDs <- names(data)

    ## filter only requested data
    if ( !missing(data.ids) ) {
        if ( verb )
            cat(paste("skipping",
                      paste(dataIDs[!dataIDs%in%data.ids],collapse=";"), "\n"))
        data <- data[dataIDs%in%data.ids]
        dataIDs <- dataIDs[dataIDs%in%data.ids]
    }
    data$dataIDs <- dataIDs
    
    data$Time <- time
    data$xids <- "Time"

    ## TODO: use this for all data parsers, although
    ## only required for viewPlate; but could be used for
    ## auto-grouping
    ## TODO: generalize? split well strings?
    data$wells$plate <- t(matrix(colnames(data[[data$dataIDs[1]]]$data),
                                 nrow=ncl,ncol=nrw))
    data$wells$rows <- toupper(letters[1:nrw])
    if ( replace.zero ) 
        data$wells$cols <- as.character(1:ncl)
    else data$wells$cols <- sprintf("%02d", 1:ncl)
        
    class(data) <- "platedata"
    ## TODO: get temperature and humidity 
    data
}

##if ( FALSE ) {
##    DATPATH <- "~/work/hhu_talks/uebung_201812/BIQ980-2018/data/BIQ980_groupI_20181206/biolector"
##    file <- "6_BIQ980_gI_20181206_REDUCTION-1.csv"
##     DATPATH <- "~/work/hhu_talks/uebung_201812/BIQ980-2018/data/BIQ980_groupII_20181207/biolector"
##    file <- "7_BIQ980_gI_20181206_REDUCTION-1.csv"
##   files <- file.path(DATPATH, file)
##    sep=";";dec="."
##}
readBioLectorProPlate <- function(files, data.ids,
                                  sep=";", dec=".",
                                  replace.zero=TRUE,
                                  verb=TRUE) {

    if ( verb )
      cat(paste("Parsing file:", files, "\n"))

    ## first parse header to get experiment information
    ## TODO: get and use more info
    dat <- read.csv(files, sep = sep, dec = dec,
                    stringsAsFactors=FALSE, header=FALSE, fill=TRUE)

    ## data start, for second parsing step
    dstart <- which(dat[,1] == "=====data=====")
   
    nrw <- as.numeric(sub(".*]","",grep("mtp-rows", dat[,1], value=TRUE)))
    ncl <- as.numeric(sub(".*]","",grep("mtp-columns", dat[,1], value=TRUE)))
    nfl <- as.numeric(sub(".*]","",grep("no_filterset", dat[,1], value=TRUE)))

    ## get filters
    filters <- cbind(1:nfl, rep(NA,nfl))
    for ( i in 1:nfl ) {
        fid <- paste0("\\[",sprintf("%02d",i),"_name\\]")
        filters[i,2] <- sub(fid,"",grep(fid,dat[,1],value=T))
    }

    ## parse again without header to get columns properly
    dat <- read.csv(files, sep = sep, dec = dec, skip = dstart,
                    stringsAsFactors=FALSE, header=FALSE, fill=TRUE)
    dend <- tail(which(dat[,1] == "M"), 1)
    
    ## parse again without header to get columns properly
    dat <- read.csv(files, sep = sep, dec = dec, skip = dstart, nrow = dend,
                    stringsAsFactors=FALSE, header=FALSE, fill=TRUE)

    ## skip last column, only NA, since lines end with ;
    dat <- dat[,-ncol(dat)]
    
    ## get time of main measurements
    ## NOTE: second T is time for user comments and system events
    tidx <- which(dat[,1] == "T")[1]
    fidx <- which(dat[tidx,] == "Time [h]->")
    ## convert time to seconds, SI and in line with other formats!
    time <- unlist(dat[tidx, (fidx+1):ncol(dat)])*3600
    
    ## get calibrated rows
    didx <- grep("Cali.", dat[,fidx])
    data <- as.data.frame(t(dat[didx,fidx:ncol(dat)]), stringsAsFactors=FALSE)
    colnames(data) <- dat[didx,1]

    data <- list()
    for ( i in 1:nrow(filters) ) {
        if ( verb ) 
            cat(paste0("\tloading data: \"", trimws(filters[i,2]), "\""))

        ## first search calibrated version
        fid <- paste0("Cali.",filters[i,2])
        didx <- grep(fid, dat[,fidx], fixed=TRUE)
        ## if not present, use raw
        if ( length(didx)==0 ) {
            didx <- grep(filters[i,2], dat[,fidx], fixed=TRUE)
            if ( verb) cat(" - RAW VALUES")
        }
        dt <- t(dat[didx,(fidx+1):ncol(dat)])
        colnames(dt) <- dat[didx,2]
        data[[i]] <- list()
        data[[i]]$time <- time
        data[[i]]$data <- as.matrix(dt)
        if ( replace.zero )
            colnames(data[[i]]$data) <- sub("0","", colnames(data[[i]]$data))
                                            
        if ( verb ) cat("\n")
    }
    names(data) <- trimws(filters[,2])
    
    dataIDs <- names(data)

    ## filter only requested data
    if ( !missing(data.ids) ) {
        if ( verb )
            cat(paste("skipping",
                      paste(dataIDs[!dataIDs%in%data.ids],collapse=";"), "\n"))
        data <- data[dataIDs%in%data.ids]
        dataIDs <- dataIDs[dataIDs%in%data.ids]
    }
    data$dataIDs <- dataIDs
    
    data$Time <- time
    data$xids <- "Time"

    ## TODO: use this for all data parsers, although
    ## only required for viewPlate; but could be used for
    ## auto-grouping
    ## TODO: generalize? split well strings?
    data$wells$plate <- t(matrix(colnames(data[[data$dataIDs[1]]]$data),
                                 nrow=ncl,ncol=nrw))
    data$wells$rows <- toupper(letters[1:nrw])
    if ( replace.zero ) 
        data$wells$cols <- as.character(1:ncl)
    else data$wells$cols <- sprintf("%02d", 1:ncl)
        
    class(data) <- "platedata"
    ## TODO: get temperature and humidity 
    data
}

#' Read Synergy Mx Plate Data
#'
#' Parses date exported from the Excel file that can be exported
#' from the Biotek Synergy Mx platereader. A parameter that often changes
#' is \code{skip}, the number of lines before the data starts,
#' here before the ID of a measurement in the first column, eg.
#' "600" for OD measurement at 600 nm.
#' @param skip lines to skip before parsing data
#' @param time.format format of the time, e.g., "%H:%M"%S" (see
#' \code{strptime}), or "numeric" if the time is provided in normal numbers
#' @param sep field separator used in the input tabular file
#' @param dec decimal operator used in data file (e.g., ',' if data export
#' was with german language settings)
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
      cat(paste("Parsing file:", files, "\n"))

    ## defaults
    if ( missing(skip) )
      skip <- 58


    indat <- read.csv(files, header=FALSE,stringsAsFactors=FALSE,
                      sep=sep, dec=dec, skip=skip)

    ## data IDs are in column 1, followed by data matrices starting
    ## in column 2; skip internal result calculation
    yidx <- c(which(indat[,1]!="" & indat[,1]!="Results"))

    ## filter only requested data
    dataIDs <- indat[yidx,1]
    yidx <- c(yidx, which(indat[,1]=="Results")) # add "Results" index
    names(yidx) <- c(dataIDs,"end")

    if ( !missing(data.ids) )
      dataIDs <- dataIDs[dataIDs%in%data.ids]

    data <- rep(list(NA),length(dataIDs))
    names(data) <- dataIDs
    for ( dataID in dataIDs  ) {

        if ( verb  )
            cat(paste0("\tloading data: \"", dataID, "\"\n"))

        ## get DATA:
        ## rows: the header is 2 rows after the ID, data starts 3 rows after
        ## and ends 2 rows before next
        i <- which(names(yidx)==dataID)
        hidx <- yidx[i]+2
        sidx <- yidx[i]+3
        eidx <- yidx[i+1]-2
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
#' Read BMG Optima/MARS vR3.01 R2 Plate Data
#'
#' Parses date from a CSV file that can be exported
#' from the MARS (vR3.01 R2) analysis software of BMG plate readers. Different
#' measurements (eg. OD and fluorescence) are exported separately
#' and should all be liste in parameter \code{files}. Later
#' versions of MARS (vXXX) have a simpler file format, and can be
#' parsed with \code{\link{readBMG2Plate}}
#' @param skip lines to skip before parsing data; if missing it will
#' be set to 5
#' @param sep field separator used in the input tabular file
#' @param dec decimal number symbol (e.g., ',' if data export was with
#' german language settings)
#' @inheritParams readPlateData
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{readPlateData}}, \code{\link{readBMG2Plate}}
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
          cat(paste("Parsing file:", file , "\n"))
        dat <- read.csv(file,header=FALSE,stringsAsFactors=FALSE,
                        skip=skip, sep=sep)

        ## TODO: scan to parse header
        ## for rownames "Well" or "Well Row", "Well Col"
        ## and "Content"; and then scan again with adjusted skip
        ## -> started in function readbmg

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
             cat(paste0("\tloading data: \"", dtyp, "\"\n"))
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
                        "in one data file:", file, "\n\t",
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

#' Read BMG Optima/MARS vXXX Plate Data
#'
#' Parses date from a CSV file that can be exported
#' from the MARS (vXXX) analysis software of BMG plate readers. Different
#' measurements (eg. OD and fluorescence). Data exported from an earlier
#' version of MARS (vR3.01 R2) can be parsed with \code{\link{readBMGPlate}}
#' @param skip lines to skip before parsing data
#' @param sep field separator used in the input tabular file
#' @param dec decimal number symbol (e.g., ',' if data export was with
#' german language settings)
#' @param time.format format of the time, e.g., "%H:%M"%S" (see
#' \code{strptime}),
#' @inheritParams readPlateData
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @seealso \code{\link{readPlateData}}, \code{\link{readBMGPlate}}
#' @export
readBMG2Plate <- function(files, data.ids, time.format=" %H h %M min",
                          skip=6, sep=";", dec=",", verb=TRUE) {



  data <- list()
  ## 1) PARSE ALL DATA FILES and collect the individual measurements
  for ( i in 1:length(files) ) {
    file <- files[i]
    #id <- names(files)[i]
    if ( verb )
      cat(paste("Parsing file:", file , "\n"))

    ## parse header and data separately, to obtain
    ## clean numeric data
    hdat <- read.csv(file,header=FALSE,stringsAsFactors=FALSE,
                     skip=skip,sep=sep, dec=dec)
    hdat <- hdat[1:2,]
    dat <- read.csv(file,header=FALSE,stringsAsFactors=FALSE,
                    skip=skip+2,sep=sep,dec=dec)

    ## last column in BMG2 files is usually empty, remove
    if ( sum(is.na(dat[,ncol(dat)]))==nrow(dat) ) {
      dat <- dat[,-ncol(dat)]
      hdat <- hdat[,-ncol(hdat)]
    }


    ## convert well info from column 1 to column names
    rownames(dat) <- sub("([A-Z])0","\\1",(dat[,1]))
    ## rm column 1 (well ID) and column 2 (sample ID, not used here)
    dat <- dat[,-(1:2)]
    hdat <- hdat[,-(1:2)]

    ## 20181209 : FIX MARS BUG: "Shake" instead of number in first timepoint
    dat[] <- apply(dat, 2, as.numeric)

    ## GET ALL DATA

    ## split data into  merged measurements, using header in row 1
    dtype <- trimws(sub("\\)","",sub("Raw Data \\(", "", hdat[1,])))
    types <- unique(dtype)

    ## collect all data
    types <- unique(dtype)
    dlst <- rep(list(NA),length(types))
    names(dlst) <- types
    for ( dtyp in types ) {
      if ( verb )
          cat(paste0("\tloading data: \"", dtyp, "\"\n"))
      tdat <- dat[,dtype==dtyp]
      htdat <- hdat[,dtype==dtyp]

      ## parse times in row 2
      ## TODO: use time.format and strptime!!
      hours <- as.numeric(sub(" h.*","",htdat[2,]))
      minutes <- as.numeric(sub(" min","",sub(".*h","",htdat[2,])))
      minutes[is.na(minutes)] <- 0
      times <- hours*3600 + minutes*60 # time in seconds

      dlst[[dtyp]] <- list(time=times,
                           data=t(as.matrix(tdat)))
    }
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

readbmg <- function(files, time.format=" %H h %M min", sep=";", dec=",", verb=TRUE) {

    data <- list()
    ## 1) PARSE ALL DATA FILES and collect the individual measurements
    for ( i in 1:length(files) ) {
        file <- files[i]
        ##id <- names(files)[i]
        if ( verb )
            cat(paste("Parsing file:", file , "\n"))

        ## 1) parse header
        dat <- read.csv(file,header=FALSE,stringsAsFactors=FALSE,
                        sep=sep, skip=5)

        widx <- grep("Well",dat[,1])
        cidx <- grep("Well Col",dat[,1])
        ridx <- grep("Well Row",dat[,1])
        iidx <- grep("Content", dat[,1])

    }

}


#' writes out all data from the `platexpress` object
#'
#' calls \code{\link{write.csv}} for each data type and writes
#' them to separate files names `<file>_<dataID>.csv`.
#' @param data \code{platexpress} data object
#' @param file either a character string naming a file or a connection open
#' for writing.  ‘""’ indicates output to the console.
#' @param yids IDs of the data to be written; if missing, all data will
#' written
#' @param xid ID of a data-set in the input data that can be used as x-axis
#' instead of the default Time vector
#' @param dtype type of the data to be plotted, default is the main 'data', but
#' e.g. 'orig' allows to plot the original data w/o processing (e.g.
#' interpolation, blank correction, etc.)
#' @param verb print messages if \code{>1}
#' @param ... arguments passed to write.csv
#' @export
writeData <- function(data, file, yids, xid, dtype="data", verb=1, ...) {

    ## which data to write?
    if ( missing(yids) ) yids <- data$dataIDs
    ## get x-axis data: time and temperature or another data set
    if ( missing(xid) )
        xid <- data$xids[1]
    x <- data[[xid]]

    for ( id in yids ) {
        file.name <- paste0(file, "_", make.names(id), ".csv")
        if ( verb ) cat(paste("writing", file.name, "\n"))
        dat <- data[[id]][[dtype]]
        dat <- cbind(x, dat)
        colnames(dat)[1] <- xid
        write.csv(dat, file.name, row.names=FALSE, ...)
    }
}
