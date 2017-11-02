#' Plot method for the rotation of latent variables in a fitted formative-reflective model
#'
#' Plot the correlations between latent variables as alpha changes.
#' @param x an object of class \code{lvcor}
#' @param lv a string: \code{"xi"} to compute exogenous latent variables' rotation, \code{"eta"} to compute
#' endogenous latent variables rotation
#' @param ... passed to the standard plot function used inside the function to draw a scatter plot
#' @return The function produces a plot and does not return anything.
#' @export
plot.lvcor <- function(x, lv = c("xi", "eta"), ...) {
  lv <- match.arg(lv[1], c("xi", "eta"))
  if (lv == "xi") {
    xi <- matrix(t(apply(x$xicor, 3, diag)), nrow = length(x$alphas))
    plot(x$alphas, xi[, 1], type = "n", ylim = c(0, 1),
         xlab = expression(alpha),  ylab = "Correlation",
         main = "Exogenous latent variable rotation", ...)
    for (i in 1:ncol(xi)) {
      lines(x = x$alphas,  y = xi[, i], col = i)
      points(x = x$alphas, y = xi[, i], pch = i)
    }
    legend_expressions <- sapply(1:ncol(xi), function(i) {
      as.expression(substitute(A [B],
                               list(A = as.name("xi"), B = as.name(i))))
    })
    legend("bottomleft", legend = legend_expressions, lty = 1, col = 1:ncol(xi), pch = 1:ncol(xi), ncol = ncol(xi))
    grid()
  } else {
    eta <- matrix(t(apply(x$etacor, 3, diag)), nrow = length(x$alphas))
    plot(x$alphas, eta[, 1], type = "n", ylim = c(0, 1),
         xlab = expression(alpha),  ylab = "Correlation",
         main = "Endogenous latent variable rotation")
    for (i in 1:ncol(eta)) {
      lines(x = x$alphas,  y = eta[, i], col = i)
      points(x = x$alphas, y = eta[, i], pch = i)
    }
    legend_expressions <- sapply(1:ncol(eta), function(i) {
      as.expression(substitute(A [B],
                               list(A = as.name("eta"), B = as.name(i))))
    })
    legend("bottomleft", legend = legend_expressions, lty = 1, col = 1:ncol(eta), pch = 1:ncol(eta), ncol = ncol(eta))
    grid()
  }
}
