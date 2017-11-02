#' Plot method for fitted formative-reflective models
#'
#' Plot the coefficients of determination of a \code{pathmod} object for all the values in \code{alphas}.
#' @param x an object of class \code{pathmod}
#' @param ... passed to the standard plot function used inside the function to draw a scatter plot
#' @return It produces a plot, it returns nothing.
#' @export
#' @importFrom graphics grid legend lines plot points
plot.pathmod <- function(x, ...) {
  xr2 <- 1 - sapply(x$xlosses, mean)
  yr2 <- 1 - sapply(x$ylosses, mean)
  plot(x$alphas, xr2, type = "n", ylim = c(0, 1),
       xlab = expression(alpha), ylab = expression(R^2),
       main = "Coefficients of determination", ...)
  lines(x$alphas, xr2, col = "blue")
  lines(x$alphas, yr2, col = "red")
  points(x$alphas, xr2, pch = 0, col = "blue")
  points(x$alphas, yr2, pch = 1, col = "red")
  legend("bottomleft", legend = c(expression(R[x]^2), expression(R[y]^2)),
         col = c("blue", "red"), lty = 1, pch = 0:1, ncol = 2)
  grid()
}
