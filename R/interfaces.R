## @importFrom grofit  grofit.control gcFitModel gcFitSpline gcBootspline
## @importFrom methods is

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
    ## TODO: addModel, using
    ## this data to add modelled data to the plate data object
    ## by calling the appropriate grofit function obtained from used.model
    params <- fit$gcTable[,p]
    colnames(params) <-sub("used","model",sub(".model","",colnames(params)))
    rownames(params) <- as.character(fit$gcTable[,"TestId"])
    data.frame(well=rownames(params), params)
}



## adds predict data from grofit to plate data; called from addModel
addGrofitModel <- function(data, fit, ID="model", ... ) {

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
#' @seealso \code{\link{growthratesResults}}, \code{\link{data2grofit}}
#' @export
data2growthrates <- function(data, yid) {

  ## TODO: use blank and amount columns in plate layout to filter

  value <- c(data[[yid]]$data)
  time <- rep(data$Time, ncol(data[[yid]]$data))
  well <- rep(colnames(data[[yid]]$data), each=nrow(data[[yid]]$data))
  df <- data.frame(time=time, value=value,
                   well=factor(well, levels=colnames(data[[yid]]$data)))
  df
}

#' parse fitted parameters from package \pkg{growthrates}
#'
#' parses the output of \code{\link{gcFit.2}} into a table
#' of the main model parameters, for each well
#' @param fit growthrates object, the result of a call to
#' \code{\link[grofit:gcFit]{gcFit}}
#' @param scale.richards if the \code{beta} parameter from Richard's model is present,
#' multiply growthrate \code{mu} with \code{beta} (original stored as \code{mumax})
#' or the \code{platexpress} version \code{\link{gcFit.2}}
#' @seealso \code{\link{data2growthrates}}, \code{\link{grofitResults}}
#' @export
growthratesResults <- function(fit, scale.richards=TRUE) {
  res <- data.frame(growthrates::results(fit), stringsAsFactors = FALSE)
  ## replace names to match names in other packages
  ## lambda, mu and A
  nms <- colnames(res)
  nms <- sub("K","A",sub("y0","X0",sub("lag","lambda",sub("mumax","mu",nms))))
  colnames(res) <- nms

  if ( scale.richards & "beta" %in% nms ) {
    ## NOTE/TODO: scaling mu from richards?
    res[,"mumax"] <- res[,"mu"]
    res[,"mu"] <- res[,"mu"] * res[,"beta"]
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
        xy <- growthrates::predict(fit, newdata=list(time=time))[,1:2]
    }
    xy
}


## adds predict data from growthrates to plate data; called from addModel
addGrowthratesModel <- function(data, fit, ID="model", ... ) {

    ## global time and new data matrix
    time <- data[[data$xids[1]]]

    ## get results
    lst <- grpredict(fit, time=time)

    ## convert to matrix
    lst <- lapply(lst, function(x) x[,2])
    newdat <- matrix(unlist(lst), ncol = length(lst), byrow = FALSE)
    colnames(newdat) <- names(lst)

    ## add and return
    addData(data=data, ID=ID, dat= newdat,
            processing=paste("growthrates",class(fit),"prediction"), ...)

}


### COMMON INTERFACES for growthrates/grofit packages

## TODO: more common functions for all packages
## use function with() to pass specific arguments
## * convertData with arguments grofit or growthrates
## * fitResults to get results

#' Add \pkg{grofit}/\pkg{growthrates} fits to plate data object
#'
#'
#' @param data a platexpress data set, see \code{\link{readPlateData}}
#' @param fit result object from calls to batch model fitting functions
#' from \pkg{grofit} or \pkg{growthrates}, ie. result of a call to
#' \code{\link{gcFit.2}} or any of \pkg{growthrates} batch functions (`all_')
#' @param ID ID for the new data set
#' @param ... arguments to \code{\link{addData}} (eg. \code{col} for
#' color selection)
#' @export
addModel <- function(data, fit, ID="model", ...) {

    if ( class(fit)=="gcFit" )
        return(addGrofitModel(data=data, fit=fit, ID=ID, ...))
    else if ( length(grep("^multiple_",class(fit)))==1 )
        return(addGrowthratesModel(data=data, fit=fit, ID=ID, ...))

}
