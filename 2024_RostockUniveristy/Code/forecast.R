

#' Fit ARIMA model to univariate time series
#' 

Arima <- function(y, order=c(0, 0, 0), seasonal=c(0, 0, 0), xreg=NULL, include.mean=TRUE,
                  include.drift=FALSE, include.constant, lambda=model$lambda, biasadj=FALSE,
                  method=c("CSS-ML", "ML", "CSS"), model=NULL, x=y, ...) {
  # Remove outliers near ends
  # j <- time(x)
  # x <- na.contiguous(x)
  # if(length(j) != length(x))
  #    warning("Missing values encountered. Using longest contiguous portion of time series")
  
  series <- deparse(substitute(y))
  
  origx <- y
  if (!is.null(lambda)) {
    x <- BoxCox(x, lambda)
    lambda <- attr(x, "lambda")
    
    if (is.null(attr(lambda, "biasadj"))) {
      attr(lambda, "biasadj") <- biasadj
    }
  }
  
  if (!is.null(xreg)) {
    if(!is.numeric(xreg))
      stop("xreg should be a numeric matrix or a numeric vector")
    xreg <- as.matrix(xreg)
    if (is.null(colnames(xreg))) {
      colnames(xreg) <- if (ncol(xreg) == 1) "xreg" else paste("xreg", 1:ncol(xreg), sep = "")
    }
  }
  
  if (!is.list(seasonal)) {
    if (frequency(x) <= 1) {
      seasonal <- list(order = c(0, 0, 0), period = NA)
      if(length(x) <= order[2L])
        stop("Not enough data to fit the model")
    } else {
      seasonal <- list(order = seasonal, period = frequency(x))
      if(length(x) <= order[2L] + seasonal$order[2L] * seasonal$period)
        stop("Not enough data to fit the model")
    }
  }
  
  if (!missing(include.constant)) {
    if (include.constant) {
      include.mean <- TRUE
      if ((order[2] + seasonal$order[2]) == 1) {
        include.drift <- TRUE
      }
    }
    else {
      include.mean <- include.drift <- FALSE
    }
  }
  if ((order[2] + seasonal$order[2]) > 1 & include.drift) {
    warning("No drift term fitted as the order of difference is 2 or more.")
    include.drift <- FALSE
  }
  
  if (!is.null(model)) {
    tmp <- arima2(x, model, xreg = xreg, method = method)
    xreg <- tmp$xreg
    tmp$fitted <- NULL
    tmp$lambda <- model$lambda
  }
  else {
    if (include.drift) {
      xreg <- `colnames<-`(cbind(drift = 1:length(x), xreg),
                           make.unique(c("drift", if(is.null(colnames(xreg)) && !is.null(xreg)) rep("", NCOL(xreg)) else colnames(xreg))))
    }
    if (is.null(xreg)) {
      suppressWarnings(tmp <- stats::arima(x = x, order = order, seasonal = seasonal, include.mean = include.mean, method = method, ...))
    } else {
      suppressWarnings(tmp <- stats::arima(x = x, order = order, seasonal = seasonal, xreg = xreg, include.mean = include.mean, method = method, ...))
    }
  }
  
  # Calculate aicc & bic based on tmp$aic
  npar <- length(tmp$coef[tmp$mask]) + 1
  missing <- is.na(tmp$residuals)
  firstnonmiss <- head(which(!missing),1)
  lastnonmiss <- tail(which(!missing),1)
  n <- sum(!missing[firstnonmiss:lastnonmiss])
  nstar <- n - tmp$arma[6] - tmp$arma[7] * tmp$arma[5]
  tmp$aicc <- tmp$aic + 2 * npar * (nstar / (nstar - npar - 1) - 1)
  tmp$bic <- tmp$aic + npar * (log(nstar) - 2)
  tmp$series <- series
  tmp$xreg <- xreg
  tmp$call <- match.call()
  tmp$lambda <- lambda
  tmp$x <- origx
  
  # Adjust residual variance to be unbiased
  if (is.null(model)) {
    tmp$sigma2 <- sum(tmp$residuals ^ 2, na.rm = TRUE) / (nstar - npar + 1)
  }
  out <- structure(tmp, class = c("forecast_ARIMA", "ARIMA", "Arima"))
  out$fitted <- fitted.Arima(out)
  out$series <- series
  return(out)
}

# Return one-step fits

fitted.Arima <- function(object, h = 1, ...) {
  if (h == 1) {
    x <- getResponse(object)
    if (!is.null(object$fitted)) {
      return(object$fitted)
    }
    else if (is.null(x)) {
      # warning("Fitted values are unavailable due to missing historical data")
      return(NULL)
    }
    else if (is.null(object$lambda)) {
      return(x - object$residuals)
    }
    else {
      fits <- InvBoxCox(BoxCox(x, object$lambda) - object$residuals, object$lambda, NULL, object$sigma2)
      return(fits)
    }
  }
  else {
    return(hfitted(object = object, h = h, FUN = "Arima", ...))
  }
}


#' Forecasting using ARIMA or ARFIMA models
#'
#' Returns forecasts and other information for univariate ARIMA models.
#'
#' For \code{Arima} or \code{ar} objects, the function calls
#' \code{\link[stats]{predict.Arima}} or \code{\link[stats]{predict.ar}} and
#' constructs an object of class "\code{forecast}" from the results. For
#' \code{fracdiff} objects, the calculations are all done within
#' \code{\link{forecast.fracdiff}} using the equations given by Peiris and
#' Perera (1988).
#'
#' @param object An object of class "\code{Arima}", "\code{ar}" or
#' "\code{fracdiff}". Usually the result of a call to
#' \code{\link[stats]{arima}}, \code{\link{auto.arima}},
#' \code{\link[stats]{ar}}, \code{\link{arfima}} or
#' \code{\link[fracdiff]{fracdiff}}.
#' @param h Number of periods for forecasting. If \code{xreg} is used, \code{h}
#' is ignored and the number of forecast periods is set to the number of rows
#' of \code{xreg}.
#' @param level Confidence level for prediction intervals.
#' @param fan If \code{TRUE}, level is set to \code{seq(51,99,by=3)}. This is
#' suitable for fan plots.
#' @param xreg Future values of an regression variables (for class \code{Arima}
#' objects only). A numerical vector or matrix of external regressors; it should not be a data frame.
#' @param bootstrap If \code{TRUE}, then prediction intervals computed using
#' simulation with resampled errors.
#' @param npaths Number of sample paths used in computing simulated prediction
#' intervals when \code{bootstrap=TRUE}.
#' @param ... Other arguments.
#' @inheritParams forecast.ts
#'
#' @return An object of class "\code{forecast}".
#'
#' The function \code{summary} is used to obtain and print a summary of the
#' results, while the function \code{plot} produces a plot of the forecasts and
#' prediction intervals.
#'
#' The generic accessor functions \code{fitted.values} and \code{residuals}
#' extract useful features of the value returned by \code{forecast.Arima}.
#'
#' An object of class "\code{forecast}" is a list containing at least the
#' following elements: \item{model}{A list containing information about the
#' fitted model} \item{method}{The name of the forecasting method as a
#' character string} \item{mean}{Point forecasts as a time series}
#' \item{lower}{Lower limits for prediction intervals} \item{upper}{Upper
#' limits for prediction intervals} \item{level}{The confidence values
#' associated with the prediction intervals} \item{x}{The original time series
#' (either \code{object} itself or the time series used to create the model
#' stored as \code{object}).} \item{residuals}{Residuals from the fitted model.
#' That is x minus fitted values.} \item{fitted}{Fitted values (one-step
#' forecasts)}
#' @author Rob J Hyndman
#' @seealso \code{\link[stats]{predict.Arima}},
#' \code{\link[stats]{predict.ar}}, \code{\link{auto.arima}},
#' \code{\link{Arima}}, \code{\link[stats]{arima}}, \code{\link[stats]{ar}},
#' \code{\link{arfima}}.
#' @references Peiris, M. & Perera, B. (1988), On prediction with fractionally
#' differenced ARIMA models, \emph{Journal of Time Series Analysis},
#' \bold{9}(3), 215-220.
#' @keywords ts
#' @aliases forecast.forecast_ARIMA
#' @examples
#' fit <- Arima(WWWusage,c(3,1,0))
#' plot(forecast(fit))
#'
#' library(fracdiff)
#' x <- fracdiff.sim( 100, ma=-.4, d=.3)$series
#' fit <- arfima(x)
#' plot(forecast(fit,h=30))
#'
#' @export
forecast.Arima <- function(object, h=ifelse(object$arma[5] > 1, 2 * object$arma[5], 10),
                           level=c(80, 95), fan=FALSE, xreg=NULL, lambda=object$lambda, bootstrap=FALSE, npaths=5000, biasadj=NULL, ...) {
  # Check whether there are non-existent arguments
  all.args <- names(formals())
  user.args <- names(match.call())[-1L] # including arguments passed to 3 dots
  check <- user.args %in% all.args
  if (!all(check)) {
    error.args <- user.args[!check]
    warning(sprintf("The non-existent %s arguments will be ignored.", error.args))
  }
  
  use.drift <- is.element("drift", names(object$coef))
  x <- object$x <- getResponse(object)
  usexreg <- (use.drift | is.element("xreg", names(object))) # | use.constant)
  
  if (!is.null(xreg) && usexreg) {
    if(!is.numeric(xreg))
      stop("xreg should be a numeric matrix or a numeric vector")
    xreg <- as.matrix(xreg)
    if (is.null(colnames(xreg))) {
      colnames(xreg) <- if (ncol(xreg) == 1) "xreg" else paste("xreg", 1:ncol(xreg), sep = "")
    }
    
    origxreg <- xreg <- as.matrix(xreg)
    h <- nrow(xreg)
  }
  else {
    if(!is.null(xreg)){
      warning("xreg not required by this model, ignoring the provided regressors")
      xreg <- NULL
    }
    origxreg <- NULL
  }
  
  if (fan) {
    level <- seq(51, 99, by = 3)
  } else {
    if (min(level) > 0 & max(level) < 1) {
      level <- 100 * level
    } else if (min(level) < 0 | max(level) > 99.99) {
      stop("Confidence limit out of range")
    }
  }
  level <- sort(level)
  
  if (use.drift) {
    n <- length(x)
    #missing <- is.na(x)
    #firstnonmiss <- head(which(!missing),1)
    #n <- length(x) - firstnonmiss + 1
    if (!is.null(xreg)) {
      xreg <- `colnames<-`(cbind(drift = (1:h) + n, xreg),
                           make.unique(c("drift", if(is.null(colnames(xreg)) && !is.null(xreg)) rep("", NCOL(xreg)) else colnames(xreg))))
    } else {
      xreg <- `colnames<-`(as.matrix((1:h) + n), "drift")
    }
  }
  
  # Check if data is constant
  if (!is.null(object$constant)) {
    if (object$constant) {
      pred <- list(pred = rep(x[1], h), se = rep(0, h))
    } else {
      stop("Strange value of object$constant")
    }
  }
  else if (usexreg) {
    if (is.null(xreg)) {
      stop("No regressors provided")
    }
    object$call$xreg <- getxreg(object)
    if (NCOL(xreg) != NCOL(object$call$xreg)) {
      stop("Number of regressors does not match fitted model")
    }
    if(!identical(colnames(xreg), colnames(object$call$xreg))){
      warning("xreg contains different column names from the xreg used in training. Please check that the regressors are in the same order.")
    }
    pred <- predict(object, n.ahead = h, newxreg = xreg)
  }
  else {
    pred <- predict(object, n.ahead = h)
  }
  
  # Fix time series characteristics if there are missing values at end of series, or if tsp is missing from pred
  if (!is.null(x)) {
    tspx <- tsp(x)
    nx <- max(which(!is.na(x)))
    if (nx != length(x) | is.null(tsp(pred$pred)) | is.null(tsp(pred$se))) {
      tspx[2] <- time(x)[nx]
      start.f <- tspx[2] + 1 / tspx[3]
      pred$pred <- ts(pred$pred, frequency = tspx[3], start = start.f)
      pred$se <- ts(pred$se, frequency = tspx[3], start = start.f)
    }
  }
  
  # Compute prediction intervals
  nint <- length(level)
  if (bootstrap) # Compute prediction intervals using simulations
  {
    sim <- matrix(NA, nrow = npaths, ncol = h)
    for (i in 1:npaths)
      sim[i, ] <- simulate(object, nsim = h, bootstrap = TRUE, xreg = origxreg, lambda = lambda)
    lower <- apply(sim, 2, quantile, 0.5 - level / 200, type = 8)
    upper <- apply(sim, 2, quantile, 0.5 + level / 200, type = 8)
    if (nint > 1L) {
      lower <- t(lower)
      upper <- t(upper)
    }
    else {
      lower <- matrix(lower, ncol = 1)
      upper <- matrix(upper, ncol = 1)
    }
  }
  else { # Compute prediction intervals via the normal distribution
    lower <- matrix(NA, ncol = nint, nrow = length(pred$pred))
    upper <- lower
    for (i in 1:nint) {
      qq <- qnorm(0.5 * (1 + level[i] / 100))
      lower[, i] <- pred$pred - qq * pred$se
      upper[, i] <- pred$pred + qq * pred$se
    }
    if (!is.finite(max(upper))) {
      warning("Upper prediction intervals are not finite.")
    }
  }
  colnames(lower) <- colnames(upper) <- paste(level, "%", sep = "")
  lower <- ts(lower)
  upper <- ts(upper)
  tsp(lower) <- tsp(upper) <- tsp(pred$pred)
  method <- arima.string(object, padding = FALSE)
  seriesname <- if (!is.null(object$series)) {
    object$series
  }
  else if (!is.null(object$call$x)) {
    object$call$x
  }
  else {
    object$call$y
  }
  fits <- fitted.Arima(object)
  if (!is.null(lambda) & is.null(object$constant)) { # Back-transform point forecasts and prediction intervals
    pred$pred <- InvBoxCox(pred$pred, lambda, biasadj, pred$se^2)
    if (!bootstrap) { # Bootstrapped intervals already back-transformed
      lower <- InvBoxCox(lower, lambda)
      upper <- InvBoxCox(upper, lambda)
    }
  }
  return(structure(
    list(
      method = method, model = object, level = level,
      mean = future_msts(x, pred$pred), lower = future_msts(x, lower),
      upper = future_msts(x, upper), x = x, series = seriesname,
      fitted = copy_msts(x, fits), residuals = copy_msts(x, residuals.Arima(object))
    ),
    class = "forecast"
  ))
}

## residuals ARIMA
residuals.Arima <- function(object, type=c("innovation", "response", "regression"), h=1, ...) {
  type <- match.arg(type)
  if (type == "innovation") {
    object$residuals
  } else if (type == "response") {
    getResponse(object) - fitted(object, h = h)
  } else {
    x <- getResponse(object)
    if (!is.null(object$lambda)) {
      x <- BoxCox(x, object$lambda)
    }
    xreg <- getxreg(object)
    # Remove intercept
    if (is.element("intercept", names(object$coef))) {
      xreg <- cbind(rep(1, length(x)), xreg)
    }
    # Return errors
    if (is.null(xreg)) {
      return(x)
    } else {
      norder <- sum(object$arma[1:4])
      return(ts(
        c(x - xreg %*% as.matrix(object$coef[(norder + 1):length(object$coef)])),
        frequency = frequency(x), start = start(x)
      ))
    }
  }
}

# Simulate ARIMA model starting with observed data x
# Some of this function is borrowed from the arima.sim() function in the stats package.
# Note that myarima.sim() does simulation conditional on the values of observed x, whereas
# arima.sim() is unconditional on any observed x.

myarima.sim <- function(model, n, x, e, ...) {
  start.innov <- residuals(model)
  innov <- e
  data <- x
  # Remove initial NAs
  first.nonmiss <- which(!is.na(x))[1]
  if (first.nonmiss > 1) {
    tsp.x <- tsp(x)
    start.x <- tsp.x[1] + (first.nonmiss - 1) / tsp.x[3]
    x <- window(x, start = start.x)
    start.innov <- window(start.innov, start = start.x)
  }
  model$x <- x
  n.start <- length(x)
  x <- ts(c(start.innov, innov), start = 1 - n.start, frequency = model$seasonal.period)
  flag.noadjust <- FALSE
  if (is.null(tsp(data))) {
    data <- ts(data, frequency = 1, start = 1)
  }
  if (!is.list(model)) {
    stop("'model' must be list")
  }
  if (n <= 0L) {
    stop("'n' must be strictly positive")
  }
  p <- length(model$ar)
  q <- length(model$ma)
  d <- 0
  D <- model$seasonal.difference
  m <- model$seasonal.period
  if (!is.null(ord <- model$order)) {
    if (length(ord) != 3L) {
      stop("'model$order' must be of length 3")
    }
    if (p != ord[1L]) {
      stop("inconsistent specification of 'ar' order")
    }
    if (q != ord[3L]) {
      stop("inconsistent specification of 'ma' order")
    }
    d <- ord[2L]
    if (d != round(d) || d < 0) {
      stop("number of differences must be a positive integer")
    }
  }
  if (p) {
    minroots <- min(Mod(polyroot(c(1, -model$ar))))
    if (minroots <= 1) {
      stop("'ar' part of model is not stationary")
    }
  }
  if (length(model$ma)) {
    # MA filtering
    x <- stats::filter(x, c(1, model$ma), method = "convolution", sides = 1L)
    x[seq_along(model$ma)] <- 0
  }
  ## AR "filtering"
  len.ar <- length(model$ar)
  
  if (length(model$ar) && (len.ar <= length(data))) {
    if ((D != 0) && (d != 0)) {
      diff.data <- diff(data, lag = 1, differences = d)
      diff.data <- diff(diff.data, lag = m, differences = D)
    } else if ((D != 0) && (d == 0)) {
      diff.data <- diff(data, lag = model$seasonal.period, differences = D)
    } else if ((D == 0) && (d != 0)) {
      diff.data <- diff(data, lag = 1, differences = d)
    } else {
      diff.data <- data
    }
    
    x.new.innovations <- x[(length(start.innov) + 1):length(x)]
    x.with.data <- c(diff.data, x.new.innovations)
    
    for (i in (length(diff.data) + 1):length(x.with.data))
    {
      lagged.x.values <- x.with.data[(i - len.ar):(i - 1)]
      ar.coefficients <- model$ar[length(model$ar):1]
      sum.multiplied.x <- sum((lagged.x.values * ar.coefficients)[abs(ar.coefficients) > .Machine$double.eps])
      x.with.data[i] <- x.with.data[i] + sum.multiplied.x
    }
    
    x.end <- x.with.data[(length(diff.data) + 1):length(x.with.data)]
    x <- ts(x.end, start = 1, frequency = model$seasonal.period)
    flag.noadjust <- TRUE
  } else if (length(model$ar)) # but data too short
  {
    # AR filtering for all other cases where AR is used.
    x <- stats::filter(x, model$ar, method = "recursive")
  }
  if ((d == 0) && (D == 0) && (flag.noadjust == FALSE)) # Adjust to ensure end matches approximately
  {
    # Last 20 diffs
    if (n.start >= 20) {
      xdiff <- (model$x - x[1:n.start])[n.start - (19:0)]
    } else {
      xdiff <- model$x - x[1:n.start]
    }
    # If all same sign, choose last
    if (all(sign(xdiff) == 1) || all(sign(xdiff) == -1)) {
      xdiff <- xdiff[length(xdiff)]
    } else { # choose mean.
      xdiff <- mean(xdiff)
    }
    x <- x + xdiff
  }
  if ((n.start > 0) && (flag.noadjust == FALSE)) {
    x <- x[-(1:n.start)]
  }
  
  ######## Undo all differences
  
  if ((D > 0) && (d == 0)) {
    # Seasonal undifferencing, if there is no regular differencing
    i <- length(data) - D * m + 1
    seasonal.xi <- data[i:length(data)]
    length.s.xi <- length(seasonal.xi)
    x <- diffinv(x, lag = m, differences = D, xi = seasonal.xi)[-(1:length.s.xi)]
  } else if ((d > 0) && (D == 0)) {
    # Regular undifferencing, if there is no seasonal differencing
    x <- diffinv(x, differences = d, xi = data[length(data) - (d:1) + 1])[-(1:d)]
  } else if ((d > 0) && (D > 0)) {
    # Undifferencing for where the differencing is both Seasonal and Non-Seasonal
    # Regular first
    delta.four <- diff(data, lag = m, differences = D)
    regular.xi <- delta.four[(length(delta.four) - D):length(delta.four)]
    x <- diffinv(x, differences = d, xi = regular.xi[length(regular.xi) - (d:1) + 1])[-(1:d)]
    
    # Then seasonal
    i <- length(data) - D * m + 1
    seasonal.xi <- data[i:length(data)]
    length.s.xi <- length(seasonal.xi)
    x <- diffinv(x, lag = m, differences = D, xi = seasonal.xi)
    x <- x[-(1:length.s.xi)]
  }
  
  x <- ts(x[1:n], frequency = frequency(data), start = tsp(data)[2] + 1 / tsp(data)[3])
  return(x)
}

## simulate ARIMA
simulate.Arima <- function(object, nsim = length(object$x), seed = NULL, xreg = NULL, future = TRUE, bootstrap = FALSE, innov = NULL, lambda = object$lambda, ...) {
  # Error check:
  if (object$arma[7] < 0) {
    stop("Value for seasonal difference is < 0. Must be >= 0")
  } else if ((sum(object$arma[c(3, 4, 7)]) > 0) && (object$arma[5] < 2)) {
    stop("Invalid value for seasonal period")
  }
  
  # Check if data is included
  x <- object$x <- getResponse(object)
  if (is.null(x)) {
    n <- 0
    future <- FALSE
    if (nsim == 0L) {
      nsim <- 100
    }
  } else {
    if (is.null(tsp(x))) {
      x <- ts(x, frequency = 1, start = 1)
    }
    if (!is.null(lambda)) {
      x <- BoxCox(x, lambda)
      lambda <- attr(x, "lambda")
    }
    n <- length(x)
  }
  
  # Check xreg
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    nsim <- nrow(xreg)
  }
  use.drift <- is.element("drift", names(object$coef))
  usexreg <- (!is.null(xreg) | use.drift | !is.null(object$xreg))
  xm <- oldxm <- 0
  if (use.drift) {
    # Remove existing drift column
    if (NCOL(xreg) == 1 && all(diff(xreg) == 1)) {
      xreg <- NULL
    } else if (!is.null(colnames(xreg))) {
      xreg <- xreg[, colnames(xreg) != "drift", drop = FALSE]
    }
    # Create new drift column
    xreg <- cbind(drift = as.matrix(seq(nsim) + n * future), xreg)
  }
  # Check xreg has the correct dimensions
  if (usexreg) {
    if (is.null(xreg)) {
      stop("xreg argument missing")
    } else if (is.null(object$xreg)) {
      stop("xreg not required")
    } else if (NCOL(xreg) != NCOL(object$xreg)) {
      stop("xreg has incorrect dimension.")
    }
  }
  
  ######## Random Seed Code
  if (is.null(innov)) {
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }
    if (is.null(seed)) {
      RNGstate <- .Random.seed
    } else {
      R.seed <- .Random.seed
      set.seed(seed)
      RNGstate <- structure(seed, kind = as.list(RNGkind()))
      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
  } else {
    nsim <- length(innov)
  }
  ######## End Random seed code
  
  # Check for seasonal ARMA components and set flag accordingly. This will be used later in myarima.sim()
  flag.s.arma <- (sum(object$arma[c(3, 4)]) > 0)
  # Check for Seasonality in ARIMA model
  if (sum(object$arma[c(3, 4, 7)]) > 0) {
    # return(simulateSeasonalArima(object, nsim=nsim, seed=seed, xreg=xreg, future=future, bootstrap=bootstrap, ...))
    if (sum(object$model$phi) == 0) {
      ar <- NULL
    } else {
      ar <- as.double(object$model$phi)
    }
    if (sum(object$model$theta) == 0) {
      ma <- NULL
    } else {
      ma <- as.double(object$model$theta)
    }
    order <- c(length(ar), object$arma[6], length(ma))
    
    if (future) {
      model <- list(
        order = order, ar = ar, ma = ma, sd = sqrt(object$sigma2), residuals = residuals(object),
        seasonal.difference = object$arma[7], seasonal.period = object$arma[5], flag.seasonal.arma = flag.s.arma,
        seasonal.order = object$arma[c(3, 7, 4)]
      )
    } else {
      model <- list(order = order, ar = ar, ma = ma, sd = sqrt(object$sigma2), residuals = residuals(object))
    }
    flag.seasonal.diff <- (object$arma[7] > 0)
  } else {
    #### Non-Seasonal ARIMA specific code: Set up the model
    order <- object$arma[c(1, 6, 2)]
    if (order[1] > 0) {
      ar <- object$model$phi[1:order[1]]
    } else {
      ar <- NULL
    }
    if (order[3] > 0) {
      ma <- object$model$theta[1:order[3]]
    } else {
      ma <- NULL
    }
    if (object$arma[2] != length(ma)) {
      stop("MA length wrong")
    } else if (object$arma[1] != length(ar)) {
      stop("AR length wrong")
    }
    
    if (future) {
      model <- list(
        order = object$arma[c(1, 6, 2)], ar = ar, ma = ma, sd = sqrt(object$sigma2), residuals = residuals(object),
        seasonal.difference = 0, flag.seasonal.arma = flag.s.arma, seasonal.order = c(0, 0, 0), seasonal.period = 1
      )
    } else {
      model <- list(order = object$arma[c(1, 6, 2)], ar = ar, ma = ma, sd = sqrt(object$sigma2), residuals = residuals(object))
    }
    flag.seasonal.diff <- FALSE
    ### End non-seasonal ARIMA specific code
  }
  
  if (bootstrap) {
    res <- na.omit(c(model$residuals) - mean(model$residuals, na.rm = TRUE))
    e <- sample(res, nsim, replace = TRUE)
  } else if (is.null(innov)) {
    e <- rnorm(nsim, 0, model$sd)
  } else if (length(innov) == nsim) {
    e <- innov
  } else {
    stop("Length of innov must be equal to nsim")
  }
  
  narma <- sum(object$arma[1L:4L])
  if (length(object$coef) > narma) {
    if (names(object$coef)[narma + 1L] == "intercept") {
      xreg <- cbind(intercept = rep(1, nsim), xreg)
      if (future) {
        object$xreg <- cbind(intercept = rep(1, n), object$xreg)
      }
    }
    if (!is.null(xreg)) {
      xm <- if (narma == 0) {
        drop(as.matrix(xreg) %*% object$coef)
      } else {
        drop(as.matrix(xreg) %*% object$coef[-(1L:narma)])
      }
      if (future) {
        oldxm <- if (narma == 0) {
          drop(as.matrix(object$xreg) %*% object$coef)
        } else {
          drop(as.matrix(object$xreg) %*% object$coef[-(1L:narma)])
        }
      }
    }
  }
  if (future) {
    sim <- myarima.sim(model, nsim, x - oldxm, e = e) + xm
  } else {
    if (flag.seasonal.diff) {
      zeros <- object$arma[5] * object$arma[7]
      sim <- arima.sim(model, nsim, innov = e)
      sim <- diffinv(sim, lag = object$arma[5], differences = object$arma[7])[-(1:zeros)]
      sim <- tail(sim, nsim) + xm
    } else {
      sim <- tail(arima.sim(model, nsim, innov = e), nsim) + xm
    }
    if (!is.null(x)) {
      sim <- ts(sim, start = tsp(x)[1], frequency = tsp(x)[3])
    } else {
      sim <- ts(sim, frequency = object$frequency)
    }
    
    # If model is non-stationary, then condition simulated data on first observation
    if (!is.null(x) & (model$order[2] > 0 || flag.seasonal.diff)) {
      sim <- sim - sim[1] + x[1]
    }
  }
  if (!is.null(lambda)) {
    sim <- InvBoxCox(sim, lambda)
  }
  
  return(sim)
}


#' Forecast plot

plot.forecast <- function(x, include, PI=TRUE, showgap = TRUE, shaded=TRUE, shadebars=(length(x$mean) < 5),
                          shadecols=NULL, col=1, fcol=4, pi.col=1, pi.lty=2, ylim=NULL, main=NULL, xlab="",
                          ylab="", type="l", flty = 1, flwd = 2, ...) {
  if (is.element("x", names(x))) { # Assume stored as x
    xx <- x$x
  } else {
    xx <- NULL
  }
  if (is.null(x$lower) || is.null(x$upper) || is.null(x$level)) {
    PI <- FALSE
  }
  else if (!is.finite(max(x$upper))) {
    PI <- FALSE
  }
  
  if (!shaded) {
    shadebars <- FALSE
  }
  if (is.null(main)) {
    main <- paste("Forecasts from ", x$method, sep = "")
  }
  if (PI) {
    x$upper <- as.matrix(x$upper)
    x$lower <- as.matrix(x$lower)
  }
  
  if (is.element("lm", class(x$model)) && !is.element("ts", class(x$mean))) # Non time series linear model
  {
    plotlmforecast(
      x, PI = PI, shaded = shaded, shadecols = shadecols, col = col, fcol = fcol, pi.col = pi.col, pi.lty = pi.lty,
      ylim = ylim, main = main, xlab = xlab, ylab = ylab, ...
    )
    if (PI) {
      return(invisible(list(mean = x$mean, lower = as.matrix(x$lower), upper = as.matrix(x$upper))))
    } else {
      return(invisible(list(mean = x$mean)))
    }
  }
  
  # Otherwise assume x is from a time series forecast
  n <- length(xx)
  if (n == 0) {
    include <- 0
  } else if (missing(include)) {
    include <- length(xx)
  }
  
  # Check if all historical values are missing
  if (n > 0) {
    if (sum(is.na(xx)) == length(xx)) {
      n <- 0
    }
  }
  if (n > 0) {
    xx <- as.ts(xx)
    freq <- frequency(xx)
    strt <- start(xx)
    nx <- max(which(!is.na(xx)))
    xxx <- xx[1:nx]
    include <- min(include, nx)
    
    if (!showgap) {
      lastObs <- x$x[length(x$x)]
      lastTime <- time(x$x)[length(x$x)]
      x$mean <- ts(c(lastObs, x$mean), start = lastTime, frequency = freq)
      x$upper <- ts(rbind(lastObs, x$upper), start = lastTime, frequency = freq)
      x$lower <- ts(rbind(lastObs, x$lower), start = lastTime, frequency = freq)
    }
  }
  else {
    freq <- frequency(x$mean)
    strt <- start(x$mean)
    nx <- include <- 1
    xx <- xxx <- ts(NA, frequency = freq, end = tsp(x$mean)[1] - 1 / freq)
    
    if (!showgap) {
      warning("Removing the gap requires historical data, provide this via model$x. Defaulting showgap to TRUE.")
    }
  }
  
  pred.mean <- x$mean
  
  if (is.null(ylim)) {
    ylim <- range(c(xx[(n - include + 1):n], pred.mean), na.rm = TRUE)
    if (PI) {
      ylim <- range(ylim, x$lower, x$upper, na.rm = TRUE)
    }
  }
  npred <- length(pred.mean)
  tsx <- is.ts(pred.mean)
  if (!tsx) {
    pred.mean <- ts(pred.mean, start = nx + 1, frequency = 1)
    type <- "p"
  }
  plot(
    ts(c(xxx[(nx - include + 1):nx], rep(NA, npred)), end = tsp(xx)[2] + (nx - n) / freq + npred / freq, frequency = freq),
    xlab = xlab, ylim = ylim, ylab = ylab, main = main, col = col, type = type, ...
  )
  if (PI) {
    if (is.ts(x$upper)) {
      xxx <- time(x$upper)
    }
    else {
      xxx <- tsp(pred.mean)[1] - 1 / freq + (1:npred) / freq
    }
    idx <- rev(order(x$level))
    nint <- length(x$level)
    if (is.null(shadecols)) {
      # require(colorspace)
      if (min(x$level) < 50) { # Using very small confidence levels.
        shadecols <- rev(colorspace::sequential_hcl(100)[x$level])
      } else { # This should happen almost all the time. Colors mapped to levels.
        shadecols <- rev(colorspace::sequential_hcl(52)[x$level - 49])
      }
    }
    if (length(shadecols) == 1) {
      if (shadecols == "oldstyle") { # Default behaviour up to v3.25.
        shadecols <- heat.colors(nint + 2)[switch(1 + (nint > 1), 2, nint:1) + 1]
      }
    }
    for (i in 1:nint)
    {
      if (shadebars) {
        for (j in 1:npred)
        {
          polygon(
            xxx[j] + c(-0.5, 0.5, 0.5, -0.5) / freq, c(rep(x$lower[j, idx[i]], 2), rep(x$upper[j, idx[i]], 2)),
            col = shadecols[i], border = FALSE
          )
        }
      }
      else if (shaded) {
        polygon(
          c(xxx, rev(xxx)), c(x$lower[, idx[i]], rev(x$upper[, idx[i]])),
          col = shadecols[i], border = FALSE
        )
      }
      else if (npred == 1) {
        lines(c(xxx) + c(-0.5, 0.5) / freq, rep(x$lower[, idx[i]], 2), col = pi.col, lty = pi.lty)
        lines(c(xxx) + c(-0.5, 0.5) / freq, rep(x$upper[, idx[i]], 2), col = pi.col, lty = pi.lty)
      }
      else {
        lines(x$lower[, idx[i]], col = pi.col, lty = pi.lty)
        lines(x$upper[, idx[i]], col = pi.col, lty = pi.lty)
      }
    }
  }
  if (npred > 1 && !shadebars && tsx) {
    lines(pred.mean, lty = flty, lwd = flwd, col = fcol)
  } else {
    points(pred.mean, col = fcol, pch = 19)
  }
  if (PI) {
    invisible(list(mean = pred.mean, lower = x$lower, upper = x$upper))
  } else {
    invisible(list(mean = pred.mean))
  }
}

# Find xreg matrix in an Arima object

getxreg <- function(z) {
  # Look in the obvious place first
  if (is.element("xreg", names(z))) {
    return(z$xreg)
  }
  # Next most obvious place
  else if (is.element("xreg", names(z$coef))) {
    return(eval.parent(z$coef$xreg))
  }
  # Now check under call
  else if (is.element("xreg", names(z$call))) {
    return(eval.parent(z$call$xreg))
  }
  # Otherwise check if it exists
  else {
    armapar <- sum(z$arma[1:4]) + is.element("intercept", names(z$coef))
    npar <- length(z$coef)
    if (npar > armapar) {
      stop("It looks like you have an xreg component but I don't know what it is.\n  Please use Arima() or auto.arima() rather than arima().")
    } else { # No xreg used
      return(NULL)
    }
  }
}

arima.string <- function(object, padding=FALSE) {
  order <- object$arma[c(1, 6, 2, 3, 7, 4, 5)]
  m <- order[7]
  result <- paste("ARIMA(", order[1], ",", order[2], ",", order[3], ")", sep = "")
  if (m > 1 && sum(order[4:6]) > 0) {
    result <- paste(result, "(", order[4], ",", order[5], ",", order[6], ")[", m, "]", sep = "")
  }
  if (padding && m > 1 && sum(order[4:6]) == 0) {
    result <- paste(result, "         ", sep = "")
    if (m <= 9) {
      result <- paste(result, " ", sep = "")
    } else if (m <= 99) {
      result <- paste(result, "  ", sep = "")
    } else {
      result <- paste(result, "   ", sep = "")
    }
  }
  if (!is.null(object$xreg)) {
    if (NCOL(object$xreg) == 1 && is.element("drift", names(object$coef))) {
      result <- paste(result, "with drift        ")
    } else {
      result <- paste("Regression with", result, "errors")
    }
  }
  else {
    if (is.element("constant", names(object$coef)) || is.element("intercept", names(object$coef))) {
      result <- paste(result, "with non-zero mean")
    } else if (order[2] == 0 && order[5] == 0) {
      result <- paste(result, "with zero mean    ")
    } else {
      result <- paste(result, "                  ")
    }
  }
  if (!padding) {
    # Strip trailing spaces
    result <- gsub("[ ]*$", "", result)
  }
  return(result)
}

# Copy msts attributes from x to y
copy_msts <- function(x, y) {
  if(NROW(x) > NROW(y)) {
    # Pad y with initial NAs
    if(NCOL(y) == 1) {
      y <- c(rep(NA, NROW(x) - NROW(y)), y)
    } else {
      y <- rbind(matrix(NA, ncol=NCOL(y), nrow = NROW(x) - NROW(y)), y)
    }
  } else if(NROW(x) != NROW(y)) {
    stop("x and y should have the same number of observations")
  }
  if(NCOL(y) > 1) {
    class(y) <- c("mts", "ts", "matrix")
  } else {
    class(y) <- "ts"
  }
  if("msts" %in% class(x))
    class(y) <- c("msts", class(y))
  attr <- attributes(x)
  attributes(y)$tsp <- attr$tsp
  attributes(y)$msts <- attr$msts
  return(y)
}


# Copy msts attributes from x to y shifted to forecast period
future_msts <- function(x, y) {
  if(NCOL(y) > 1) {
    class(y) <- c("mts", "ts", "matrix")
  } else {
    class(y) <- "ts"
  }
  if("msts" %in% class(x))
    class(y) <- c("msts", class(y))
  attr <- attributes(x)
  attr$tsp[1:2] <- attr$tsp[2] + c(1,NROW(y))/attr$tsp[3]
  attributes(y)$tsp <- attr$tsp
  attributes(y)$msts <- attr$msts
  return(y)
}

## get Response functions

getResponse <- function(object, ...) UseMethod("getResponse")

#' @rdname getResponse
#' @export
getResponse.default <- function(object, ...) {
  if (is.list(object)) {
    return(object$x)
  } else {
    return(NULL)
  }
}

#' @rdname getResponse
#' @export
getResponse.lm <- function(object, ...) {
  if(!is.null(object[["x"]])){
    object[["x"]]
  }
  else{
    responsevar <- deparse(formula(object)[[2]])
    model.frame(object$model)[, responsevar]
  }
}

#' @rdname getResponse
#' @export
getResponse.Arima <- function(object, ...) {
  if (is.element("x", names(object))) {
    x <- object$x
  } else {
    series.name <- object$series
    if (is.null(series.name)) {
      return(NULL)
    } else {
      x <- try(eval.parent(parse(text = series.name)), silent = TRUE)
      if (is.element("try-error", class(x))) { # Try one level further up the chain
        x <- try(eval.parent(parse(text = series.name), 2), silent = TRUE)
      }
      if (is.element("try-error", class(x))) { # Give up
        return(NULL)
      }
    }
  }
  return(as.ts(x))
}


