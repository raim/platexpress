## @importFrom grofit  grofit.control gcFitModel gcFitSpline gcBootspline
## @importFrom methods is

#' hack of grofit function \code{\link[grofit:grofit]{gcFit}} to allow plotting
#' without interaction
#' @param time a matrix of measurment times for each well as required by
#' \code{\link[grofit:grofit]{grofit}}, provided by \code{\link{data2grofit}}
#' @param data the data as produced by \code{\link{data2grofit}}
#' @param control the \code{\link[grofit:grofit]{grofit}} control structure, see
#' \code{\link[grofit:grofit]{grofit.control}}, but with an additional
#' entry "plot", which is set to TRUE or FALSE
#' @details Adding the field "plot" allows to de-activate
#' \code{control$interactive}, so the whole fitting procedure runs
#' through all wells, yet results are plotted.
#' @export
gcFit.2 <- function (time, data, control = grofit::grofit.control()) 
{

    
    
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
            cat(paste("= ", as.character(i), ". growth curve =================================\n", 
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
                gcID = gcID, fit.x = NA, fit.y = NA, parameters = list(A = NA, 
                  mu = NA, lambda = NA, integral = NA), parametersLowess = list(A = NA, 
                  mu = NA, lambda = NA), spline = NA, reliable = NULL, 
                fitFlag = FALSE, control = control)
            class(nonpara) <- "gcFitSpline"
            fitnonpara.all[[i]] <- nonpara
        }
        wellname <- paste(as.character(data[i, 1]), as.character(data[i, 
            2]), as.character(data[i, 3]), sep = "-")
        if ((control$plot == TRUE)) {
            if (fitpara$fitFlag == TRUE) {
                plot(fitpara, colData = "black", colModel = "green", cex = 0.5)
                plot(nonpara, add = TRUE, raw = FALSE, colData = 0, 
                  colSpline = "red", cex = 1.5)
                ## TODO: why is it two
                abline(v=summary(fitpara)$lambda) #fitpara$parameters$lambda)
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
                  "spline fit"), col = c("green", "red"), lty = c(1, 
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
    gcFit
}
