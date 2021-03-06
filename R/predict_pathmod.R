#' Predict method for fitted formative-reflective models
#'
#' Once the formative-reflective model has been estimated with the \code{pathmod()} function,
#' you can produce predictions based on old or new data.
#' @param object an object of class \code{pathmod} generated by the function \code{pathmod()}.
#' If no data were supplied to the function \code{pathmod()} through the \code{data} parameter,
#' then the parameters
#' @param newdata a data frame containing the \eqn{x}-variables of the model named as they had been
#' named in the \code{xgroups} argument passed to the \code{pathmod()} function. If you do not supply
#' this argument but data have been passed to the function \code{pathmod()} through the \code{data}
#' parameter, then the function returs the fitted values.
#' @param xmean numeric vector of \eqn{x}-variable means. If no data and no \code{xmean} parameter were
#' supplied to the \code{pathmod()} call, then you have to assign this argument, otherwise this
#' parameter is not needed. If you assign it when \code{xmean} is available in \code{object},
#' then this value is used instead of the one present in \code{object}.
#' @param ymean numeric vector of \eqn{y}-variable means. If no data and no \code{ymean} parameter were
#' supplied to the \code{pathmod()} call, then you have to assign this argument, otherwise this
#' parameter is not needed. If you assign it when \code{ymean} is available in \code{object},
#' then this value is used instead of the one present in \code{object}.
#' @param xsd numeric vector of \eqn{x}-variable standard deviations. If no data and no \code{xsd} parameter were
#' supplied to the \code{pathmod()} call, then you have to assign this argument, otherwise this
#' parameter is not needed. If you assign it when \code{xsd} is available in \code{object},
#' then this value is used instead of the one present in \code{object}.
#' @param ysd numeric vector of \eqn{y}-variable means. If no data and no \code{ysd} parameter were
#' supplied to the \code{pathmod()} call, then you have to assign this argument, otherwise this
#' parameter is not needed. If you assign it when \code{ysd} is available in \code{object},
#' then this value is used instead of the one present in \code{object}.
#' @param ... further arguments passed to or from other methods
#' @return A list with slots:
#' \describe{
#'   \item{\code{xis}}{3-dimensional array with a matrix of \eqn{\xi}-variables for each value of \eqn{\alpha}}
#'   \item{\code{etas}}{3-dimensional array with a matrix of \eqn{\eta}-variables for each value of \eqn{\alpha}}
#'   \item{\code{xhats}}{3-dimensional array with a matrix of \eqn{x}-variable predictions for each value of \eqn{\alpha}}
#'   \item{\code{yhats}}{3-dimensional array with a matrix of \eqn{y}-variable predictions for each value of \eqn{\alpha}}
#' }
#' @examples
#' data(russett)
#' xg <- list(c("gini", "farm", "rent"), c("gnpr", "labo"))
#' yg <- list(c("inst", "ecks", "death"))
#' pm <- pathmod(xg, yg, data = russett[1:40, ])
#' pr <- predict(pm, russett[41:47, ])
#' @export
predict.pathmod <- function(object, newdata = NULL, xmean = NULL, ymean = NULL, xsd = NULL, ysd = NULL, ...) {
  # no newdata: returns insample xi, eta, xhat, yhat (one for each alpha)
  if (is.null(newdata)) {
    if (is.null(object$xis)) stop("no data supplied to the `pathmod()` function call or to the `newdata` parameter.")
    ncolx <- length(unlist(object$xgroups))
    ncoly <- length(unlist(object$ygroups))
    n <- nrow(object$xi[, , 1])
    na <- length(object$alphas)
    xhat <- array(0.0, c(n, ncolx, na),
                  dimnames = list(rownames(object$xis), unlist(object$xgroups), names(object$xis[1, 1, ])))
    yhat <- array(0.0, c(n, ncoly, na),
                  dimnames = list(rownames(object$xis), unlist(object$ygroups), names(object$xis[1, 1, ])))
    for (i in 1:na) {
      xh <- with(object,
                 xis[, , i] %*% solve(omegas[[i]] %*% Sxx %*% t(omegas[[i]])) %*% omegas[[i]] %*% Sxx)
      xhat[, , i] <- xh * rep(object$xsd, rep.int(nrow(xh), ncol(xh))) +
                     rep(object$xmean, rep.int(nrow(xh), ncol(xh)))
      yh <- with(object, etas[, , i] %*% t(lambdas[[i]]))
      yhat[, , i] <- yh * rep(object$ysd, rep.int(nrow(yh), ncol(yh))) +
                     rep(object$ymean, rep.int(nrow(yh), ncol(yh)))
    }
    return(list(xis = object$xis, etas = object$etas, xhats = xhat, yhats = yhat))
  }
  # manage means and scales
  if (!is.null(xmean)) object$xmean <- xmean
  else if (is.null(object$xmean)) stop("no xmean parameter available in the pathmod object, please supply it.")
  if (!is.null(ymean)) object$ymean <- ymean
  else if (is.null(object$ymean)) stop("no ymean parameter available in the pathmod object, please supply it.")
  if (!is.null(xsd)) object$xsd <- xsd
  else if (is.null(object$xsd)) stop("no xsd parameter available in the pathmod object, please supply it.")
  if (!is.null(ysd)) object$ysd <- ysd
  else if (is.null(object$ysd)) stop("no ysd parameter available in the pathmod object, please supply it.")
  # use newdata to generate xi, eta, xhat, yhat (one for each alpha)
  anames <- paste0("alpha = ", object$alphas)
  xdata <- newdata[, unlist(object$xgroups), drop = FALSE]
  xdata <- scale(xdata, center = object$xmean, scale = object$xsd)
  n <- nrow(xdata)
  na <- length(object$alphas)
  xis <- array(0.0, c(n, object$p, na),
               dimnames = list(rownames(newdata), paste0("xi", 1:object$p), anames))
  etas <- array(0.0, c(n, object$q, na),
                dimnames = list(rownames(newdata), paste0("eta", 1:object$q), anames))
  xhat <- array(0.0, c(n, ncol(xdata), na),
                dimnames = list(rownames(newdata), unlist(object$xgroups), anames))
  yhat <- array(0.0, c(n, length(unlist(object$ygroups)), na),
                dimnames = list(rownames(newdata), unlist(object$ygroups), anames))
  for (i in 1:na) {
    xis[, , i] <- xdata %*% t(object$omegas[[i]])
    etas[, , i] <- xis[, , i] %*% t(object$gammas[[i]])
    xh <- xis[, , i] %*%
      solve(object$omegas[[i]] %*% object$Sxx %*% t(object$omegas[[i]])) %*%
      object$omegas[[i]] %*% object$Sxx
    xhat[, , i] <- xh * rep(object$xsd, rep.int(nrow(xh), ncol(xh))) +
      rep(object$xmean, rep.int(nrow(xh), ncol(xh)))
    yh <- etas[, , i] %*% t(object$lambdas[[i]])
    yhat[, , i] <- yh * rep(object$ysd, rep.int(nrow(yh), ncol(yh))) +
      rep(object$ymean, rep.int(nrow(yh), ncol(yh)))
  }
  list(xis = xis, etas = etas, xhats = xhat, yhats = yhat)
}
