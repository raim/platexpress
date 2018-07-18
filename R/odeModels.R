#' Monod's Growth Model
#'
#' System of two differential equations describing bacterial growth as an
#' autocatalytic process of substrate to biomass conversion
#'
#' The model is given as a system of two differential equations:
#'
#' \eqn{mu = mumax * s/(s+K)}
#' \deqn{ds/dt = phi*(sin-s) - (mu/Y)*y}
#' \deqn{dy/dt = (mu-phi)*y}
#'
#' @param time actual time (for the ode) resp. vector of simulation time steps.
#' @param init named vector with state of the system
#' @param parms parameters of Monod's growth model:
#'   \itemize{
#'      \item \code{y} initial biomass concentration (cell number),
#'      \item \code{s} initial substrate concentration
#'      \item \code{mumax} maximum growth rate (1/time),
#'      \item \code{K} saturation or half-maximal growth constant
#'      \item \code{Y} the yield factor (fraction of substrate converted to biomass
#'      \item \code{phi} culture dilution rate
#'      \item \code{sin} substrate concentration in medium influx
#'   }
#' @param \dots placeholder for additional parameters (for user-extended versions of this function)
#'
#' @return
#'
#' For \code{ode_monod}: matrix containing the simulation outputs.
#' The return value of has also class \code{deSolve}.
#'
#' For \code{grow_monod}: vector of dependent variable (\code{y}) and
#'   its log-transformed values (\code{log_y}):
#'
#' \itemize{
#' \item \code{time} time of the simulation
#' \item \code{y} cell concentration
## \item \code{log_y} natural log of total cell concentration
#' }
#'
#' @details Function \code{ode_monod} is the system of differential equations,
#' whereas \code{grow_monod} runs a numerical simulation over time.
#'
#'
#' @references
#'
#' Monod, J (1950). La technique de la culture continue, th\'eorie et applications. Annales de Institute Pasteur Paris 79(4), 390--410.
#'
#' @examples
#'
#' time <- seq(0, 30, length=200)
#' parms <- c(phi=0, sin=0, mumax=1, K=1, Y=0.5)
#' y0    <- c(s=1, y=0.01)
#' out   <- deSolve::ode(y0, time, ode_monod, parms)
#'
#' plot(out)
#'
#' o <- grow_monod(0:100, c(y0, parms))
#' plot(o)
#'
#' @family growth models
#'
#' @rdname grow_monod
#' @export ode_monod
ode_monod <- function (time, init, parms, ...) {

    with(as.list(c(parms, init)), {
        mu <- mumax * s/(s+K)
        ds <- phi*(sin-s) - (mu/Y)*y 
        dy <- (mu-phi)*y
        list(c(ds, dy)) #, log_y=unname(log(y))))
    })
    
}


#' Simple monod growth model
#' @family growth models
#' @rdname grow_monod
#' @export
grow_monod <- function(time, parms, ...) {

    ## assign parameters and solve differential equations
    init    <- parms[c("s", "y")]
    parms <- parms[c("phi", "sin", "mumax", "K", "Y")]
    out <- deSolve::ode(init, time, platexpress::ode_monod, parms = parms)
}
## attach names of parameters as attributes
attr(grow_monod, "fname") <- c("grow_monod")
attr(grow_monod, "pnames") <- c("s", "y", "phi", "sin", "mumax", "K", "Y")
class(grow_monod) <- c("growthmodel", "function")


#' Anabolism vs. Catabolism Growth Model
#' @param time actual time (for the ode) resp. vector of simulation time steps.
#' @param init named vector with the initial state of the system
#' @param parms parameters of the anabolism vs. catabolism growth model
#'   \itemize{
#'      \item \code{y} initial biomass concentration (cell number),
#'      \item \code{s} initial substrate concentration
#'      \item \code{atp} initial ATP concentration (fraction of ATP+ADP=1)
#'      \item \code{phi} culture dilution rate
#'      \item \code{sin} substrate concentration in medium influx
#'      \item \code{mumax_ab} maximum growth rate (1/time),
#'      \item \code{mumax_cd} maximum catabolic rate (1/time),
#'      \item \code{mumax_m} maintenance rate concentration/time
#'      \item \code{n_ab} anabolic ATP stoichiometry (n ATP / C-mol substrate)
#'      \item \code{n_cd} catabolic ATP stoichiometry (n ATP / C-mol substrate)
#'      \item \code{Kaab} anabolic ATP saturation or half-maximal rate constant
#'      \item \code{Ksab} anabolic substrate saturation or half-maximal rate constant
#'      \item \code{Kacd} catabolic ATP saturation or half-maximal rate constant
#'      \item \code{Kscd} catabolic substrate saturation or half-maximal rate constant
#'      \item \code{Cc} and \code{Vc} cellular carbon content and volumes
#'   }
#' @param \dots placeholder for additional parameters (for user-extended versions of this function)
#'
#' @return
#'
#' For \code{ode_anacat}: matrix containing the simulation outputs.
#' The return value of has also class \code{deSolve}.
#'
#' For \code{grow_anacat}: vector of dependent variable (\code{y}) and
#'   its log-transformed values (\code{log_y}):
#'
#' \itemize{
#' \item \code{time} time of the simulation
#' \item \code{y} cell concentration
## \item \code{log_y} natural log of total cell concentration
#' }
#'
#' @details Function \code{ode_anacat} is the system of differential equations,
#' whereas \code{grow_anacat} runs a numerical simulation over time.
#'
#' @family growth models
#' @rdname grow_anacat
#' @export 
ode_anacat <- function (time, init, parms, ...) {

    with(as.list(c(parms, init)), {
        adp <- axp - atp 
        mu_ab <- mumax_ab * s/(s+Ksab) * atp/(atp+Kaab)
        mu_cd <- mumax_cd * s/(s+Kscd) * adp/(adp+Kacd)
        mu_m <- mumax_m
        ds <- phi*(sin-s) - (mu_ab + mu_cd)*y 
        dy <- (mu_ab-phi)*y
        datp = (n_cd*mu_cd - n_ab*mu_ab - mu_m) *Cc/Vc - mu_ab*atp
        list(c(ds, dy, datp))#, log_y=unname(log(y))))
    })
    
}

#' @family growth models
#' @rdname grow_anacat
#' @export
grow_anacat <- function(time, parms, ...) {

    ## TODO: define output variables to fit
    ## OD .. from cells/OD, dry weight per cell, and C/cell
    ## cCc: 47% of dryweight is C 
    ## http://bionumbers.hms.harvard.edu/bionumber.aspx?id=100649&ver=4
    ## cpod: 5e8 cells/OD
    
    ## assign parameters and solve differential equations
    init    <- parms[c("s", "y", "atp")]
    parms <- parms[c("phi", "sin",
                     "mumax_ab", "Kaab", "Ksab", "n_ab",
                     "mumax_cd", "Kacd", "Kscd", "n_cd",
                     "mumax_m", "Cc", "Vc", "axp")]
    out <- deSolve::ode(init, time, platexpress::ode_anacat, parms = parms)
}
## attach names of parameters as attributes
attr(grow_anacat, "fname") <- c("grow_anacat")
attr(grow_anacat, "pnames") <- c("s", "y", "atp",
                                 "phi", "sin",
                                 "mumax_ab", "Kaab", "Ksab",
                                 "mumax_cd", "Kacd", "Kscd",
                                 "mumax_m", "Cc", "Vc", "axp")
class(grow_anacat) <- c("growthmodel", "function")


#' Regulated Anabolism vs. Catabolism Growth Model
#' @param time actual time (for the ode) resp. vector of simulation time steps.
#' @param init named vector with the initial state of the system
#' @param parms parameters of the anabolism vs. catabolism growth model
#'   \itemize{
#'      \item \code{y} initial biomass concentration (cell number),
#'      \item \code{s} initial substrate concentration
#'      \item \code{atp} initial ATP concentration (fraction of ATP+ADP=1)
#'      \item \code{phi} culture dilution rate
#'      \item \code{sin} substrate concentration in medium influx
#'      \item \code{mumax_ab} maximum growth rate (1/time),
#'      \item \code{mumax_cd} maximum catabolic rate (1/time),
#'      \item \code{mumax_m} maintenance rate concentration/time
#'      \item \code{n_ab} anabolic ATP stoichiometry (n ATP / C-mol substrate)
#'      \item \code{n_cd} catabolic ATP stoichiometry (n ATP / C-mol substrate)
#'      \item \code{Kaab} anabolic ATP saturation or half-maximal rate constant
#'      \item \code{Ksab} anabolic substrate saturation or half-maximal rate constant
#'      \item \code{Kacd} catabolic ATP saturation or half-maximal rate constant
#'      \item \code{Kscd} catabolic substrate saturation or half-maximal rate constant
#'      \item \code{Cc} and \code{Vc} cellular carbon content and volumes
#'   }
#' @param \dots placeholder for additional parameters (for user-extended versions of this function)
#'
#' @return
#'
#' For \code{ode_anacatr}: matrix containing the simulation outputs.
#' The return value of has also class \code{deSolve}.
#'
#' For \code{grow_anacatr}: vector of dependent variable (\code{y}) and
#'   its log-transformed values (\code{log_y}):
#'
#' \itemize{
#' \item \code{time} time of the simulation
#' \item \code{y} cell concentration
## \item \code{log_y} natural log of total cell concentration
#' }
#'
#' @details Function \code{ode_anacatr} is the system of differential equations,
#' whereas \code{grow_anacatr} runs a numerical simulation over time.
#'
#' @family growth models
#' @rdname grow_anacatr
#' @export 
ode_anacatr <- function (time, init, parms, ...) {

    with(as.list(c(parms, init)), {

        adp <- axp - atp

        mu_ab <- mumax_ab * pab * s/(s+Ksab) * atp/(atp+Kaab)
        mu_cd <- mumax_cd * pcd * s/(s+Kscd) * adp/(adp+Kacd)
        mu_m  <- mumax_m

        dpab  <- (kab * atp/(atp+Ktab) - dab * pab)*Cc/Vc -mu_ab*pab
        dpcd  <- (kcd * adp/(adp+Ktcd) - dcd * pcd)*Cc/Vc -mu_ab*pcd
        
        ds    <- phi*(sin-s) - (mu_ab + mu_cd)*y 
        dy    <- (mu_ab-phi)*y
        datp  <- (n_cd*mu_cd - n_ab*mu_ab - mu_m)*Cc/Vc - mu_ab*atp
        
        list(c(ds, dy, datp, dpab, dpcd))#, log_y=unname(log(y))))
    })
    
}

#' @family growth models
#' @rdname grow_anacatr
#' @export
grow_anacatr <- function(time, parms, ...) {

    ## TODO: define output variables to fit
    ## OD .. from cells/OD, dry weight per cell, and C/cell
    ## cCc: 47% of dryweight is C 
    ## http://bionumbers.hms.harvard.edu/bionumber.aspx?id=100649&ver=4
    ## cpod: 5e8 cells/OD
    
    ## assign parameters and solve differential equations
    init    <- parms[c("s", "y", "atp","pab","pcd")]
    parms <- parms[c("phi", "sin",
                     "kab", "Ktab", "dab",
                     "kcd", "Ktcd", "dcd",
                     "mumax_ab", "Kaab", "Ksab", "n_ab",
                     "mumax_cd", "Kacd", "Kscd", "n_cd",
                     "mumax_m", "Cc", "Vc", "axp")]
    out <- deSolve::ode(init, time, platexpress::ode_anacatr, parms = parms)
}
## attach names of parameters as attributes
attr(grow_anacatr, "fname") <- c("grow_anacatr")
attr(grow_anacatr, "pnames") <- c("s", "y", "atp","pab","pcd",
                                 "kab", "Ktab", "dab",
                                 "kcd", "Ktcd", "dcd",
                                 "phi", "sin",
                                 "mumax_ab", "Kaab", "Ksab",
                                 "mumax_cd", "Kacd", "Kscd",
                                 "mumax_m", "Cc", "Vc", "axp")
class(grow_anacatr) <- c("growthmodel", "function")
