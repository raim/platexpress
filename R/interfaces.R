## @importFrom grofit  grofit.control gcFitModel gcFitSpline gcBootspline
## @importFrom methods is

### GROFIT INTERFACE

#' hack of grofit function \code{\link[grofit:grofit.control]{grofit.control}}
#' which adds the new "plot" switch used in \code{\link{gcFit.2}}
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
#' without interaction
#' @param time a matrix of measurment times for each well as required by
#' \code{\link[grofit:grofit]{grofit}}, provided by \code{\link{data2grofit}}
#' @param data the data as produced by \code{\link{data2grofit}}
#' @param control the \code{\link[grofit:grofit]{grofit}} control structure, see
#' \code{\link[grofit:grofit]{grofit.control}}, but with an additional
#' entry "plot", which is set to TRUE or FALSE
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
        fit.opt = "b"
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

#' hack of cellGrowth bandwidthCV to work with platexpress data format
#' TODO: select bandwidths from data$Time
#' @export
bandwidthCV.2 = function(data, did, mid, wells,
  bandwidths=seq(0.5,10, length.out=30), # hours - TODO: take 1/5th of data$Time
  nFold=10,
  nWell=50,
  cutoff=0.95,
  ##calibration=identity,
  scaleY=identity # log2
  ) {

    if ( missing(mid) )
      mid <- data$mids[1]
    
    if ( missing(did) ) # only use requested data 
      did <- data$dataIDs[1]
    
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
        y = dat[[did]]$data[ , match(w, wells)]
        ## negative ODs?
        if ( !all(y>=0)) {
            warning("Negative values ",did," found. If you are using a calibration function, this might be the problem. Values are set to NAs")
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
          cat(paste("Treating well",w,colnames(data[[did]]$data)[w],"\n"))
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
#' @export
fitCellGrowths.2 = function(data, did, mid, wells, ...) {

    if ( missing(mid) )
      mid <- data$mids[1]
    
    if ( missing(did) ) # only use requested data 
      did <- data$dataIDs[1]

    if ( missing(wells) )
      wells <- colnames(data[[did]]$data)

    ##fits <- matrix(NA,ncol=4,nrow=length(wells))
    ##rownames(fits) <- wells
    fits <- NULL
    x <- data[[mid]]
    for ( well in wells ) {
        cat(paste("well", well, "\n"))
        y <- data[[did]]$data[,well]
        fit <- cellGrowth::fitCellGrowth(x, z=y, ...)
        ##fits[well,] <- unlist(attributes(fit)[c(3,4,5,6)])
        fits <- append(fits, list(fit))

        ## TODO: align with grofit output and add to data
        ## calculate lag from pointOfMaxGrowthRate and MaxGrowthRate
        ## calculate X(0)
        ## TODO: add fits to data for plots
        ## data[[did]]$fit <- fit...data
    }
    names(fits) <- wells
    fits
}
 
