## @importFrom grofit  grofit.control gcFitModel gcFitSpline gcBootspline
## @importFrom methods is

#### GROWTH PHASE SEGMENTATION

### SEGMENTED INTERFACE
#' Call \code{\link[segmented]{segmented}} for selected wells.
#'
#' The package \code{\link[segmented]{segmented}} splits curves
#' into linear segments, given a list of pre-defined breakpoints.
#' If the passed data is a biomass measure
#' (eg. OD), and option \code{log=TRUE} the slopes correspond to local
#' growth rate (per time unit). Missing (NA) or non-finite
#' y-values will be removed. \code{\link[segmented]{segmented}} has
#' a random initialization and sometimes fails. The function will
#' attempt \code{maxtry} of calls before giving up.
#' @param data \code{platexpress} data object
#' @param yid ID of the \code{platexpress} data to use
#' @param wells subset and plot order of wells
#' @param log use ln of the data
#' @param xid x-axis data ID in \code{platexpress} data
#' @param man apply a moving average \code{\link{ma}} with \code{n=man}
#' @param maxtry maximum number of attempts to run
#' \code{\link[segmented]{segmented}}
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
                            man=1, maxtry=5,
                            psis, psi, npsi=5, plot=FALSE, verb=0, ...) {

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
            psi <-  seq(min(X),max(X),length.out=npsi+2)[1:npsi+1]
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
        mxtry <- maxtry
        while( class(test)[1]=="try-error" & mxtry>0 ) {
            test <- try(o <- segmented::segmented(out.lm,psi=psis[[well]]))
            if ( mxtry<5 )
                warning(well, ": segmented failed in attempt:",
                        maxtry-mxtry+1, "of", maxtry)
            mxtry <- mxtry -1
        }
        

        if ( verb>0 ) cat("\n")
        if ( maxtry==0 ) {
            warning(well, ": segmented failed")
            next
        }
        segments[[well]] <- o
        if ( plot ) {
            segmented::plot.segmented(o,add=TRUE, col=2, lwd=3,
                                      shade=TRUE, rug=FALSE)
            segmented::lines.segmented(o,col=1, lwd=1)
            abline(v=confint(o)$x[,1])
            Sys.sleep(.4)
          }
    }
    class(segments) <- "segmented_list"
    invisible(segments)
}

#' Add results from  \code{\link{segmented_plate}}
#' to \code{platexpress} object.
#'
#' Calls the \code{\link[segmented:predict.segmented]{predict.segmented}} method
#' to reconstruct the piecewise linear growth curve, and adds it to the
#' \code{platexpress} data object. Option \code{add.slopes}
#' allows to add a time-course of growth rates (slopes
#' of the linear segments) instead.
#' @param data \code{platexpress} data object
#' @param fit list of \code{\link[segmented:segmented]{segmented}} objects
#' @param ID data ID for the new object
#' @param add.slopes add slopes instead of reconstructed data
#' @param ... arguments passed to \code{\link{addModel}}
#'@export
addModel_segmented <- function(data, fit, ID="y", add.slopes=FALSE, ...) {

    x <- data[[data$xids[1]]]

    if ( add.slopes ) {
        ## get breakpoints and slopes
        segs <- lapply(fit, function(o) c(min(x),
                                          segmented::confint.segmented(o)$x[,1],
                                          max(x)))
        slps <- lapply(fit, function(o) segmented::slope(o)$x[,1])


        ## add mus as xy data
        mus <- sapply(1:length(segs), function(i) {
            sg <- segs[[i]]
            sl <- slps[[i]]
            yo <- rep(sl, each=10)
            xo <- c(sapply(2:length(sg), function(j)
                seq(sg[j-1]+1e-3, sg[j], length.out=10)))
            approx(xo,yo, xout=x)$y
        })
        colnames(mus) <- names(fit)
        ## add to platexpress data to plot!
        invisible(addData(data, ID=paste(ID,sep="_"),
                          dat=mus, ...))
    } else {
        ## add modelled data
        OD_segs <- lapply(fit, function(o) predict(o,newdata=data.frame(x=x)))
        OD_segs <- do.call(cbind, OD_segs)
        ## add to platexpress data to plot!
        invisible(addData(data, ID=paste(ID,sep="_"),dat=exp(OD_segs), ...))
    }
}

### DPSEG INTERFACE
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
    class(segments) <- "dpseg_list"
    invisible(segments)
}

#' Add results from  \code{\link{dpseg_plate}}
#' to \code{platexpress} object.
#'
#' Calls the \code{\link[dpseg:predict.dpseg]{predict.dpseg}} method
#' to reconstruct the piecewise linear growth curve, and adds it to the
#' \code{platexpress} data object. Option \code{add.slopes}
#' allows to add a time-course of growth rates (slopes
#' of the linear segments) instead.
#' @param data \code{platexpress} data object
#' @param fit list of \code{\link[dpseg:dpseg]{dpseg}} objects
#' @param ID data ID for the new object
#' @param add.slopes add slopes instead of reconstructed data
#' @param ... arguments passed to \code{\link{addModel}}
#'@export
addModel_dpseg <- function(data, fit, ID="y", add.slopes=FALSE, ...) {

    if ( add.slopes ) {
        ## get breakpoints and slopes
        segs <- lapply(fit, function(x) x$segments[,c("x1","x2")])
        slps <- lapply(fit, function(x) x$segments[1:nrow(x$segments),"slope"])

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
        colnames(mus) <- names(fit)
        ## add to platexpress data to plot!
        invisible(addData(data, ID=paste(ID,sep="_"), dat=mus, ...))
    } else {
        ## add modelled data
        OD_segs <- lapply(fit, function(o) predict(o)$y)
        OD_segs <- do.call(cbind, OD_segs)
        ## add to platexpress data to plot!
        invisible(addData(data, ID=paste(ID,sep="_"), dat=exp(OD_segs), ...))
    }
}

#### GROWTH MODEL FITS

### GROFIT INTERFACE

### This would require to depend on grofit
## @examples
## data(ap12)
## grdat <- data2grofit(ap12data)
## fit <- gcFit.2(grdat$time, grdat$data)

### HIGH-LEVEL WRAPPER FOR GROFIT
#' high-level wrapper for \pkg{grofit}
#'
#' Calls \pkg{grofit} directly from \code{platexpress}
#' plate data and layout, using established settings. Please see
#' the documentation of package \pkg{grofit}
#' for details of the fitting procedure. This high-level wrapper uses
#' \code{platexpress} functions \code{\link{data2grofit}},
#' \code{\link{grofit.2.control}}, \code{\link{gcFit.2}} to convert
#' data, to set parameters and to call \pkg{grofit},
#' and then uses \code{\link{grofitResults}} to map results (lag phase lambda,
#' growth rate mu and capacity A) to the \code{plate} layout map, and
#' \code{\link{addModel}} to add modelled prediction of the fitted data to the
#' plate \code{data} object. All calculations
#' are executed by the original functions of \pkg{grofit}.
#'
#' The function returns a list containing the new \code{data} and \code{parameters}
#' objects, as well as the original \code{grofit} result object.
#' NOTE that the latter contains a lot more information on fit quality and
#' run statistics. All run parameters of \code{\link[grofit:gcFit]{gcFit}}
#' can be set via option \code{control}.
#' @param data a platexpress data set, see \code{\link{readPlateData}}
#' @param plate plate layout map, see \code{\link{readPlateMap}}, columns
#' of this map can be converted to \pkg{grofit} data
#' annotation, and results will be returned in the same order of wells
#' @param yid data ID of the data to be converted; default
#' is to take the first of \code{data$dataIDs}
#' @param amount name of a column in \code{plate} that provides a numeric
#' 'amount', to be used as a 'dose' annotation in grofit
#' @param fields column IDs in the plate layout map to be used for
#' \pkg{grofit} data annotation; if \code{plate} is not
#' missing at least 1 column is required
#' @param control an object of class \code{grofit.control} with an additional
#' parameter \code{plot}, as provided by the \code{platexpress} function
#' \code{\link{grofit.2.control}}; if present arguments \code{model.type},
#' \code{nboot.gc} and \code{plot} will be ignored (unless the latter
#' is not present)
#' @param model.type models that \code{grofit} will attempt to fit,
#' see \code{\link[grofit:grofit.control]{grofit.control}}; will be igonred
#' if \code{control} is passed directly
#' @param nboot.gc numbers of bootstrap samples,
#' see \code{\link[grofit:grofit.control]{grofit.control}}; will be igonred
#' if \code{control} is passed directly
#' @param plot set to TRUE for plots even without \code{interactive} use
#' @param interactive set to TRUE for interactive growth curve fitting
#' @param col color of the model data in the plate data object
#' @param verb print messages
#' @param ... further arguments to \code{\link{grofit.2.control}}
#'@export
callGrofit <- function(data, plate, yid, amount,
                       fields=c("strain","medium","substance"),
                       control, model.type=c("richards","logistic",
                                             "gompertz.exp","gompertz"),
                       nboot.gc=100, plot=TRUE, interactive=FALSE, verb=TRUE,
                       col="#0000FF", ...) {

  ## take the first data, if none was passed
  if ( missing(yid) ) {
    yid <- data$dataIDs[1]
    if ( verb )
      cat(paste("auto-selected data:", yid, "\n"))
  }

  ## generate input data, incl. plate annotation
  if ( missing(plate) )
    gdat <- data2grofit(data, yid=yid, verb=verb)
  else {
    ## if a plate map is available, use this for annotation!

    ## allow for full capabilities of grofit results: add dose
    if ( missing(amount) )
      gdat <- data2grofit(data, plate=plate, yid=yid,
                          wells=plate$well[!plate$blank],
                          eid=fields, verb=verb, ...)
    else
      gdat <- data2grofit(data, plate=plate, yid=yid,
                          wells=plate$well[!plate$blank],
                          dose=plate[[amount]][!plate$blank],
                          eid=fields, verb=verb, ...)
  }
  ## set grofit paraameters
  if ( missing(control) )
    control <- grofit.2.control(plot=plot, interactive=interactive,
                                suppress.messages=TRUE,
                                nboot.gc=nboot.gc,
                                model.type=model.type)
  if ( !"plot" %in% names(control) )
    control$plot <- plot

  ## call grofit; redirect plots to detail pdf
  odfit <- gcFit.2(gdat$time, gdat$data, control)

  ## get grofit results
  params <- grofitResults(odfit)
  #colnames(params) <- paste0(yid,"_",colnames(params))

  ## ... add to layout data
  if ( !missing(plate) )
    params <- params[as.character(plate[,"well"]),]

  ## add modelled data!
  data <- addModel(data, odfit, ID=paste0(yid,"_model"), col=col)

  res <- list(data=data, parameters=params, fit=odfit)

  res
}

### data2grofit: see AP12.R for example, TODO: fix example data and update file
#' interface to package \pkg{grofit}
#'
#' \code{\link{data2grofit}} : converts \code{\link{platexpress}} data to
#' \pkg{grofit} data format
#' @param data a platexpress data set, see \code{\link{readPlateData}}
#' @param yid data ID of the data to be converted for grofit, from \code{data$dataIDs}
#' @param min.time minimal time of the data to be used
#' @param max.time maximal time of the data to be used
#' @param wells column IDs of the data set to use, if missing all wells
#' are used
#' @param plate plate layout map, see \code{\link{readPlateMap}}, columns
#' of this map can be converted to \pkg{grofit} data
#' annotation
#' @param eid column IDs in the plate layout map to be used for
#' \pkg{grofit} data annotation; if missing but
#' \code{plate} is present, the columns 2 and 3 are used
#' @param dose vector of doses in each well, used as the third column of
#'  \pkg{grofit} data annotation, where it can be used for
#' dose-response calculations
#' @param verb print messages
#' @details Returns a simple list with two entries \code{time} and \code{data},
#' as required for \pkg{grofit}.
#' @author Rainer Machne \email{raim@tbi.univie.ac.at}
#' @export
data2grofit <- function(data, yid, min.time, max.time, wells, plate, eid, dose, verb=TRUE) {

    if ( missing(yid) ) {
        yid <- data$dataIDs[1]
        if ( verb )
          cat(paste("auto-selected data:", yid, "\n"))
    }
    dat <- data[[yid]]$data
    if ( missing(wells) )
        wells <- colnames(dat)
    else wells <- as.character(wells)

    dat <- data.frame(dat[,wells])

    ## expand time to full matrix
    ## TODO: use internal time and non-interpolated data?
    time <- data$Time
    if ( !missing(max.time) ) {
        dat <- dat[time<=max.time,]
        time <- time[time<=max.time]
    }
    if ( !missing(min.time) ) {
        dat <- dat[time>=min.time,]
        time <- time[time>=min.time]
    }
    time <- t(matrix(rep(time, ncol(dat)), c(length(time), ncol(dat))))

    ## well annotation
    ## use well as TestId
    if ( !missing(plate) ) {
        if ( missing(eid) ) ## add second column as default
            eid <- colnames(plate)[2]
        ## sub-selection
        idx <- match(wells, as.character(plate[,"well"]))
        annotation1 <- plate[idx,"well"]
        ## use both fields if two eids are given ..
        annotation2 <- plate[idx,eid[1]]
        if ( length(eid)>1 )
            for ( k in 2:length(eid) )
                annotation2 <- paste0(annotation2,"-",plate[idx,eid[k]])
        annotation <- data.frame(annotation1, annotation2)

    } else
        annotation <- data.frame(cbind(colnames(dat),
                                       rep("",ncol(dat))))
    ## dose information for grofit dose-response calculations
    ## TODO: this is ugly, do nicer!
    found.dose <- !missing(dose)
    if ( !found.dose )
        if ( !missing(plate) ) ## get dose info from plate layout: TODO
            if ( "amount" %in% colnames(plate) ) {
                idx <- match(wells,as.character(plate[,"well"]))
                dose <- as.numeric(plate[idx,"amount"])
                found.dose <- TRUE
            }
    if ( !found.dose )
        dose <- rep(0,ncol(dat))

    ## construct grofit data
    grdat <- data.frame(annotation, dose, t(dat))
    list(time=time, data=grdat)
}

#' parse \pkg{grofit} results
#'
#' parses the output of \code{\link{gcFit.2}} into a table
#' of the main model parameters, for each well, for both the
#' "best" model found by `grofit` as well the `spline` and `bootstrap`
#' fits. Parameter names are streamlined, ie. the ".model" suffix is removed.
#' @param fit grofit object, the result of a call to
#' \code{\link[grofit:gcFit]{gcFit}}
#' or the \code{platexpress} version \code{\link{gcFit.2}}
#' @param p the parameters to retrieve from \code{fit$gcTable},
#' the main results table of \code{\link[grofit:gcFit]{gcFit}}
#' @seealso \code{\link{data2grofit}}, \code{\link{growthratesResults}}
#' @export
grofitResults <- function(fit, p=c("lambda.model","mu.model","A.model","used.model",
                                   "lambda.spline","mu.spline","A.spline",
                                   "lambda.bt","mu.bt","A.bt")) {
    params <- fit$gcTable[,p]
    colnames(params) <-sub("used","model",sub(".model","",colnames(params)))
    rownames(params) <- as.character(fit$gcTable[,"TestId"])
    data.frame(well=rownames(params), params)
}



## adds predict data from grofit to plate data; called from addModel
addModel_gcFit <- function(data, fit, ID="model", ... ) {

    testid <- "TestId" # this requires wells being used as TestId in grofit

    ## copy first existing data set
    newdat <- data[[data$dataIDs[[1]]]]$data
    newdat[] <- NA
    models <- rep(NA, ncol(newdat))
    names(models) <- colnames(newdat)
    ## parse grofit results and use predict to add data
    for ( i in 1:ncol(newdat) ) {
        id <- colnames(newdat)[i]
        idx <- which(as.character(fit$gcTable[,testid])==id)
        if ( !length(idx) ) next
        if ( !fit$gcFittedModels[[idx]]$fitFlag ) next
        newdat[,i] <- stats::predict(fit$gcFittedModels[[idx]]$nls,
                                     newdata=data.frame(time=data$Time))
        models[i] <- fit$gcFittedModels[[idx]]$model
    }
    ## TODO: use models info somehow?
    ## add and return
    addData(data=data, ID=ID, dat= newdat,
            processing="grofit prediction", ...)
}

#' wrapper of grofit function \code{\link[grofit:grofit.control]{grofit.control}}
#' which adds the new \code{plot} switch used in \code{\link{gcFit.2}}
#' @param ... parameters passed on to
#' \code{\link[grofit:grofit.control]{grofit.control}}
#' @param interactive set to TRUE for interactive growth curve fitting
#' @param plot set to TRUE for plots even without \code{interactive} use
#' @seealso \code{\link{gcFit.2}}
#' @export
grofit.2.control <- function(interactive=FALSE, plot=TRUE, ...) {
    control <- grofit::grofit.control(...)
    control$interactive <- interactive
    control$plot <- plot
    control
}

#' hack of grofit function \code{\link[grofit:gcFit]{gcFit}} to allow plotting
#' even when option \code{interative==FALSE}
#' @param time a matrix of measurment times for each well as required by
#' \pkg{grofit}, provided by \code{\link{data2grofit}}
#' @param data the data as produced by \code{\link{data2grofit}}
#' @param control the \pkg{grofit} control structure, see
#' \code{\link[grofit:grofit]{grofit.control}}, but with an additional
#' argument "plot" (TRUE or FALSE), allowing to plot individual fits
#' even if interactive is set to FALSE.
#' @details Adding the field "plot" allows to de-activate
#' \code{control$interactive}, so the whole fitting procedure runs
#' through all wells, yet results are plotted. \code{\link{platexpress}}
#' also offers a hack of \code{\link[grofit:grofit.control]{grofit.control}}
#' names \code{\link{grofit.2.control}} that adds these options.
#' @seealso \code{\link{gcFit.2}}
#' @export
gcFit.2 <- function (time, data, control = grofit.2.control())  {

    requireNamespace("grofit")
    
    if (methods::is(control) != "grofit.control")
        stop("control must be of class grofit.control!")
    if ((dim(time)[1]) != (dim(data)[1]))
        stop("gcFit: Different number of datasets in data and time")
    if (!is.element(control$fit.opt, c("s", "m", "b"))) {
        warning("fit.opt must be set to 's', 'm' or 'b'. Changed to 'b'!")
        control$fit.opt = "b"
    }

     ## add gcFit.2 specific controls
    if ( !"plot"%in%names(control) )
      control$plot <- TRUE

    out.table <- NULL
    used.model <- NULL
    fitpara.all <- list()
    fitnonpara.all <- list()
    boot.all <- list()
    fitted.param <- NULL
    fitted.nonparam <- NULL
    bootstrap.param <- NULL
    for (i in 1:dim(data)[1]) {
        acttime <- as.numeric(as.matrix(time[i, ]))
        actwell <- as.numeric(as.matrix((data[i, -1:-3])))
        gcID <- as.matrix(data[i, 1:3])
        if ((control$suppress.messages == FALSE)) {
            cat("\n\n")
            cat(paste("= ", as.character(i),
                      ". growth curve =================================\n",
                sep = ""))
            cat("----------------------------------------------------\n")
        }
        if ((control$fit.opt == "m") || (control$fit.opt == "b")) {
            fitpara <- grofit::gcFitModel(acttime, actwell, gcID, control)
            fitpara.all[[i]] <- fitpara
        }
        else {
            fitpara <- list(raw.x = acttime, raw.y = actwell,
                gcID = gcID, fit.x = NA, fit.y = NA, parameters = list(A = NA,
                  mu = NA, lambda = NA, integral = NA), model = NA,
                nls = NA, reliable = NULL, fitFlag = FALSE, control = control)
            class(fitpara) <- "gcFitModel"
            fitpara.all[[i]] <- fitpara
        }
        if ((control$fit.opt == "s") || (control$fit.opt == "b")) {
            nonpara <- grofit::gcFitSpline(acttime, actwell, gcID, control)
            fitnonpara.all[[i]] <- nonpara
        }
        else {
            nonpara <- list(raw.x = acttime, raw.y = actwell,
                            gcID = gcID, fit.x = NA, fit.y = NA,
                            parameters = list(A=NA,mu=NA,lambda=NA,integral=NA),
                            parametersLowess = list(A=NA,mu=NA,lambda=NA),
                            spline = NA, reliable = NULL,
                            fitFlag = FALSE, control = control)
            class(nonpara) <- "gcFitSpline"
            fitnonpara.all[[i]] <- nonpara
        }
        wellname <- paste(as.character(data[i, 1]), as.character(data[i,
            2]), as.character(data[i, 3]), sep = "-")
        if ((control$plot == TRUE)) {
            if (fitpara$fitFlag == TRUE) {
                plot(fitpara, colData = "black", colModel = "blue", cex = 0.5)
                plot(nonpara, add = TRUE, raw = FALSE, colData = 0,
                  colSpline = "red", cex = 1.5)
                abline(h=summary(fitpara)$A,lwd=2,col="blue")
                abline(v=summary(fitpara)$lambda,lwd=2,col="blue")
            }
            else {
                plot(nonpara, add = FALSE, raw = TRUE, colData = "black",
                  colSpline = "red", cex = 0.5)
            }
            title(wellname)
            if (control$fit.opt == "m")
                legend(x = "bottomright", legend = fitpara$model,
                  col = "black", lty = 1)
            if (control$fit.opt == "s")
                legend(x = "bottomright", legend = "spline fit",
                  col = "red", lty = 1)
            if (control$fit.opt == "b")
                legend(x = "bottomright", legend = c(fitpara$model,
                  "spline fit"), col = c("blue", "red"), lty = c(1,
                  1))
        }
        reliability_tag <- NA
        if (control$interactive == TRUE) {
            answer <- readline("Are you satisfied (y/n)?")
            if (substr(answer, 1, 1) == "n") {
                cat("\n Tagged this well as unreliable !\n\n")
                reliability_tag <- FALSE
                fitpara.all[[i]]$reliable <- FALSE
                fitnonpara.all[[i]]$reliable <- FALSE
            }
            else {
                reliability_tag <- TRUE
                fitpara.all[[i]]$reliable <- TRUE
                fitnonpara.all[[i]]$reliable <- TRUE
                cat("Well was (more ore less) o.k.\n")
            }
        }
        else {
            reliability_tag <- TRUE
        }
        if (control$interactive == TRUE)
            graphics.off()
        if ((control$nboot.gc > 0) && (reliability_tag == TRUE)) {
            bt <- grofit::gcBootSpline(acttime, actwell, gcID, control)
            boot.all[[i]] <- bt
        }
        else {
            bt <- list(raw.x = acttime, raw.y = actwell, gcID = gcID,
                boot.x = NA, boot.y = NA, boot.gcSpline = NA,
                lambda = NA, mu = NA, A = NA, integral = NA,
                bootFlag = FALSE, control = control)
            class(bt) <- "gcBootSpline"
            boot.all[[i]] <- bt
        }
        description <- data.frame(TestId = data[i, 1], AddId = data[i,
            2], concentration = data[i, 3], reliability = reliability_tag,
            used.model = fitpara$model, log.x = control$log.x.gc,
            log.y = control$log.y.gc, nboot.gc = control$nboot.gc)
        fitted <- cbind(description, summary(fitpara), summary(nonpara),
            summary(bt))
        out.table <- rbind(out.table, fitted)
    }
    gcFit <- list(raw.time = time, raw.data = data, gcTable = out.table,
        gcFittedModels = fitpara.all, gcFittedSplines = fitnonpara.all,
        gcBootSplines = boot.all, control = control)
    class(gcFit) <- "gcFit"
    ## TODO: add fits to data
    ## wrapper: use platexpress data as input
    ## and add fits to data, aligned with cellGrowth::fitCellGrowth.2 wrapper
    gcFit
}



### cellGrowth INTERFACE

#' hack of the cellGrowth package function
#' \code{\link[cellGrowth:bandwidthCV]{bandwidthCV}}
#' to work with the \code{platexpress} data format, using adjusted default
#' values
#' TODO: select bandwidths from data$Time
#'
#' @param data platexpress data set, see \code{\link{readPlateData}}
#' @param yid data ID of the data to be converted for grofit, from
#' \code{data$dataIDs}
#' @param xid the x-axis to be used, defaults to the first available
#' (usually "Time")
#' @param wells column IDs of the data set to use, if missing all wells
#' are used
#' @param bandwidths see \code{\link[cellGrowth:bandwidthCV]{bandwidthCV}}
#' @param nFold see \code{\link[cellGrowth:bandwidthCV]{bandwidthCV}}
#' @param nWell see \code{\link[cellGrowth:bandwidthCV]{bandwidthCV}}
#' @param cutoff see \code{\link[cellGrowth:bandwidthCV]{bandwidthCV}}
#' @param scaleY see \code{\link[cellGrowth:bandwidthCV]{bandwidthCV}}
#' @seealso \code{\link[cellGrowth:bandwidthCV]{bandwidthCV}}
#'
#' @export
bandwidthCV.2 = function(data, yid, xid, wells,
                         bandwidths=seq(0.5,10, length.out=30), # hours - TODO: take 1/5th of data$Time
                         nFold=10,
                         nWell=50,
                         cutoff=0.95,
                         ##calibration=identity,
                         scaleY=identity # log2
) {

    if ( missing(xid) )
      xid <- data$xids[1]

    if ( missing(yid) ) # only use requested data
      yid <- data$dataIDs[1]

    ## function to split sample in n folds
    cv.folds = function (n ) {
        split(sample(1:n), rep(1:nFold, length = n))
    }

    ## function to predict max growth rate
    cvpred_mu = function(cv, h, x, z) {
        ls = lapply(
          cv,
          function(fo){
              fit = tryCatch(
                cellGrowth::fitCellGrowth(x = x[-fo], z = z[-fo],
                                          model = "locfit", locfit.h = h),
                warning=function(w){NA}, error=function(e){NA})
              if( !is.null(attr(fit,"maxGrowthRate")) ) {
                  rv = list(pred=stats::predict( fit, x[fo]),
                    mu = attr(fit, "maxGrowthRate"))
              }else{
                  rv = list(pred=rep(NA, length(fo)), mu = NA)
              }
              rv
          })
        list(pred = unlist(lapply(ls, function(x) x$pred)),
             mu = unlist(lapply(ls, function(x) x$mu)))
    }

    ## function to calculate error and std
    err2_mustd_well = function(dat, w, hs) {
        x = dat$Time
        y = dat[[yid]]$data[ , match(w, wells)]
        ## negative ODs?
        if ( !all(y>=0)) {
            warning("Negative values ",yid," found. If you are using a calibration function, this might be the problem. Values are set to NAs")
            y[y<0] = NA
        }

        z = scaleY( y )

        cv = cv.folds(length(x))
        pms = lapply(hs, function(h) cvpred_mu(cv,h,x,z))
        pred = sapply(pms, function(x) x$pred)
        mu = sapply(pms, function(x) x$mu)

        err = pred - z[unlist(cv)]
        list(err2 = colMeans(err^2),
             err2std = apply(err^2, 2, sd),
             mustd =  apply(mu, 2, sd))
    }

    ## do cv
    ws = sample(length(wells),nWell) # wells to perform cv on
    hs = bandwidths # possible bandwidths


    ## calculate mu error and std
    err2_mustd = lapply(ws,
      function(w) {
          cat(paste("Treating well",w,colnames(data[[yid]]$data)[w],"\n"))
          data <- err2_mustd_well(dat=data, w=wells[w], hs=hs)
      })


    err2 = sapply(err2_mustd, function(x) x$err2)
    err2std = sapply(err2_mustd, function(x) x$err2std)
    mustd = sapply(err2_mustd, function(x) x$mustd)

    ## which bandwidth are at one std of mini in average
    onestd_of_mini = sapply(1:ncol(err2),
      function(i) {
          m = which.min(err2[,i])
          err2[,i] <= err2[m,i] + err2std[m,i]
      })

    ## onestd_of_mini returns a list
    if( is.list(onestd_of_mini) )
      onestd_of_mini = do.call(cbind,onestd_of_mini)

    ## calculate "optimal" bandwidth
    bandwidth = hs[max( which(rowMeans(onestd_of_mini,na.rm=TRUE)>= cutoff))]


    return(list(bandwidth=bandwidth,
                wells=wells,
                bandwidths=hs,
                err2=err2,
                err2std=err2std,
                muStd=mustd,
                oneStdOfMini=onestd_of_mini))
}

#' batch wrapper for fitCellGrowth
#' TODO: align in/out with gcFit.2 or optionally cellGrowth::fitCellGrowths
#'
#' @param data platexpress data set, see \code{\link{readPlateData}}
#' @param yid data ID of the data to be converted for grofit, from
#' \code{data$dataIDs}
#' @param xid the x-axis to be used, defaults to the first available
#' (usually "Time")
#' @param wells column IDs of the data set to use, if missing all wells
#' are used
#' @param ... further parameters to the \code{cellGrowth} package function
#' \code{\link[cellGrowth:fitCellGrowth]{fitCellGrowth}}
#' @seealso \code{\link[cellGrowth:fitCellGrowth]{fitCellGrowth}}
#'
#' @export
fitCellGrowths.2 = function(data, yid, xid, wells, ...) {

    if ( missing(xid) )
      xid <- data$xids[1]

    if ( missing(yid) ) # only use requested data
      yid <- data$dataIDs[1]

    if ( missing(wells) )
      wells <- colnames(data[[yid]]$data)

    ##fits <- matrix(NA,ncol=4,nrow=length(wells))
    ##rownames(fits) <- wells
    fits <- NULL
    x <- data[[xid]]
    for ( well in wells ) {
        cat(paste("well", well, "\n"))
        y <- data[[yid]]$data[,well]
        fit <- cellGrowth::fitCellGrowth(x, z=y, ...)
        ##fits[well,] <- unlist(attributes(fit)[c(3,4,5,6)])
        fits <- append(fits, list(fit))

        ## TODO: align with grofit output and add to data
        ## calculate lag from pointOfMaxGrowthRate and MaxGrowthRate
        ## calculate X(0)
        ## TODO: add fits to data for plots
        ## data[[yid]]$fit <- fit...data
    }
    names(fits) <- wells
    fits
}

### growthrates INTERFACE

## TODO: wrapper to call all_easylinearfits etc?

#' interface to package \pkg{growthrates}
#'
#' converts \code{platexpress} data for
#' use with \pkg{growthrates}
#' @param data platexpress data set, see \code{\link{readPlateData}}
#' @param yid data ID of the data to be converted for growthrates, from
#' \code{data$dataIDs}
#' @param wells column IDs of the data set to use, if missing all wells
#' are used
#' @param plate plate layout map, see \code{\link{readPlateMap}}, to skip
#' blanks wells
#' annotation
#' @seealso \code{\link{growthratesResults}}, \code{\link{data2grofit}}
#' @export
data2growthrates <- function(data, yid, wells, plate) {

    if ( missing(wells) )
        wells <- colnames(data[[yid]]$data)
    ## filter by plate & rm blanks
    if ( !missing(plate) ) {
        wells <- wells[wells%in%plate[,"well"]]
        wells <- wells[!plate[match(wells, plate[,"well"]),"blank"]]
    }

    dat <- data[[yid]]$data[,wells,drop=FALSE]

  value <- c(dat)
  time <- rep(data$Time, ncol(dat))
  well <- rep(colnames(dat), each=nrow(dat))
  df <- data.frame(time=time, value=value,
                   well=factor(well, levels=colnames(dat)))
  df
}

#' parse fitted parameters from package \pkg{growthrates}
#'
#' parses the output of \code{\link{gcFit.2}} into a table
#' of the main model parameters, for each well
#' @param fit growthrates object, the result of a call to fit functions
#' from package \code{\link[growthrates]{growthrates}}, such as
#' \code{\link[growthrates:all_easylinear]{all_easylinear}}
#' @param scale.richards if the \code{beta} parameter from Richard's model is present,
#' multiply growthrate \code{mu} with \code{beta} (original stored as \code{mumax})
#' @seealso \code{\link{data2growthrates}}, \code{\link{grofitResults}}
#' @export
growthratesResults <- function(fit, scale.richards=TRUE) {
  res <- data.frame(growthrates::results(fit), stringsAsFactors = FALSE)
  ## replace names to match names in other packages
  ## lambda, mu; but keep K (A for some models, but K in Monod)
  nms <- colnames(res)
  nms <- sub("y0","X0",sub("lag","lambda",sub("mumax","mu",nms)))
  colnames(res) <- nms

  if ( scale.richards & "beta" %in% nms ) {
    ## NOTE/TODO: scaling mu from richards?
    res[,"mumax"] <- res[,"mu"]
    res[,"mu"] <- res[,"mu"] * res[,"beta"]
    ## TODO: also scale gompertz?
  }

  data.frame(res)
}

#' `predict' hack for  \pkg{growthrates} fits
#'
#' \pkg{growthrates} provides several
#' functions to fit growth rates, but their result objects have currently
#' slighly different structures/interfaces. To obtain relevant data from
#' \code{\link[growthrates:fit_easylinear]{fit_easylinear}} and
#' \code{\link[growthrates:fit_spline]{fit_spline}} the fitted
#' functions and parameters have to be directly used instead of just
#' calling \code{\link[growthrates:predict]{predict}} methods.
#' @param fit fit object from \pkg{growthrates} fitting functions, either
#' a single fit (`fit_' functions) or a list of fits from batch functions
#' (`all_' functions)
#' @param time the time points at which growth data is to be predicted
#' @export
grpredict <- function(fit, time) {

    ## NOTE: predict not implemented for easylinear and returning
    ## the full spline fit for smooth.spline

    ## call recursively for lists of fits
    if ( length(grep("^multiple_",class(fit)))==1 )
        return(lapply(fit@fits, grpredict, time=time))

    if ( class(fit)=="easylinear_fit" ) {

        xy <- fit@FUN(time, fit@par)[,1:2]
        ## NOTE: easylinear requires to add the lag phase
        xy[,1] <- xy[,1] + growthrates::coef(fit)["lag"]
        ## interpolate to requested time and convert to matrix
        xy <- matrix(unlist(approx(x=xy[,1], y=xy[,2], xout=time)),
                            ncol=2, byrow=FALSE)

    } else if ( class(fit)=="smooth.spline_fit" ) {

        ## NOTE: smooth.spline predict returns full spline fit as log(y)
        xy <- fit@FUN(time, fit@par)[,1:2]

    } else {
        xy <- growthrates::predict(fit, newdata=list(time=time))[,c("time","y")]
    }
    xy
}


## adds predict data from growthrates to plate data; called from addModel
addModel_growthrates <- function(data, fit, ID="model", ... ) {

    ## global time and new data matrix
    time <- data[[data$xids[1]]]

    ## newdat matrix, just copy first data set
    newdat <- data[[data$dataIDs[1]]]$data
    newdat[] <- NA

    ## get results
    lst <- grpredict(fit, time=time)

    ## convert to matrix
    lst <- lapply(lst, function(x) x[,2])
    newdat[,names(lst)] <- matrix(unlist(lst),
                                  ncol = length(lst), byrow = FALSE)


    ## add and return
    addData(data=data, ID=ID, dat= newdat,
            processing=paste("growthrates",class(fit),"prediction"), ...)

}


### COMMON INTERFACES for growthrates/grofit packages

## TODO: more common functions for all packages
## use function with() to pass specific arguments
## * convertData with arguments grofit or growthrates
## * fitResults to get results

#' Add fits from various growth model interfaces to plate data object.
#'
#' TODO: depending on class of `fit`.
#' @param data a platexpress data set, see \code{\link{readPlateData}}
#' @param fit result object from calls to batch model fitting functions
#' \code{\link{gcFit.2}}, any of \pkg{growthrates} batch functions (`all_'),
#' or \code{platexpress} interfaces \code{\link{segmented_plate}}, and
#' \code{\link{dpseg_plate}}.
#' @param ID ID for the new data set
#' @param ... arguments to \code{\link{addData}} (eg. \code{col} for
#' color selection)
#' @export
addModel <- function(data, fit, ID="model", ...) {

    if ( class(fit)=="gcFit" )
        return(addModel_gcFit(data=data, fit=fit, ID=ID, ...))
    else if ( class(fit)=="dpseg_list" )
        return(addModel_dpseg(data=data, fit=fit, ID=ID, ...))
    else if ( class(fit)=="segmented_list" )
        return(addModel_dpseg(data=data, fit=fit, ID=ID, ...))
    else if ( length(grep("^multiple_",class(fit)))==1 )
        return(addModel_growthrates(data=data, fit=fit, ID=ID, ...))

}


#' merge results to plate layout map
#'
#' a wrapper around R base function \code{\link[base:merge]{merge}}
#' for merging results per well into a \pkg{platexpress} layout map.
#' Usage is the same as for \code{\link[base:merge]{merge}} with
#' some defaults changed.
#' @param x a plate layout map or any data frame with a "well" column
#' specified by argument \code{by}
#' @param y results: any data frame with a "well" column
#' specified by argument \code{by}, see \code{\link[base:merge]{merge}}
#' @param ID prefix added to to the merged \code{y} data
#' @param by column name(s) by which to merge, 
#' @param all.x keep all rows of \code{x}, see \code{\link[base:merge]{merge}}
#' @param sort sort the results on the \code{by} column? See
#' \code{\link[base:merge]{merge}}
#' @param ... further arguments to \code{\link[base:merge]{merge}}
#' @export
mergeResults <- function(x, y, ID, by = "well", 
                         all.x=TRUE, sort=FALSE, ...) {

    if ( !missing(ID) )
        colnames(y)[2:ncol(y)] <- paste0(ID, "_",colnames(y)[2:ncol(y)])
    merge(x=x, y=y, by="well", all.x=all.x, sort=sort, ...)

}
