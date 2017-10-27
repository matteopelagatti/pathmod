#' Reduced rank regression using covariance matrices.
#'
#' Reduced rank regression solves the problem of estimating the matrices \eqn{A} and \eqn{B} in the linear model
#' \deqn{y = A B' x + \epsilon}
#' where \eqn{y} is a \eqn{p}-vector of response variables,
#' \eqn{x} is a \eqn{q}-vector of regressors and
#' \eqn{A} and \eqn{B} are coefficient matrices with dimensions \eqn{p x r} and \eqn{q x r}, respectively
#' (\eqn{1 \le r \le min(p, q)}).
#' The function is used internally by the package and so the passed arguments are not checked.
#'
#' @param Syy variance-covariance matrix of the \eqn{y} vector.
#' @param Syx covariance matrix of \eqn{y} and \eqn{x}.
#' @param Sxx variance-covariance matrix of the \eqn{x} vector.
#' @param r   rank of reduced rank regression.
#' @return A list with the following slots:
#' \describe{
#'   \item{\code{A}}{Matrix \eqn{A}}
#'   \item{\code{B}}{Matrix \eqn{B}}
#'   \item{\code{eigenvalues}}{Eigenvalues of the first \code{r} canonical correlations between \eqn{y} and \eqn{x}.}
#' }
#' @examples
#' S <- matrix(c(
#'   1.0, 0.2, 0.8, 0.0, 0.9, 0.9,
#'   0.2, 1.0, 0.0, 0.0, 0.2, 0.2,
#'   0.8, 0.0, 1.0, 0.2, 0.8, 0.8,
#'   0.0, 0.0, 0.2, 1.0, 0.1, 0.1,
#'   0.9, 0.2, 0.8, 0.1, 1.0, 0.8,
#'   0.9, 0.2, 0.8, 0.1, 0.8, 1.0),
#'   6,6,byrow=TRUE)
#' rrrbycov(S[4:6, 4:6], S[4:6, 1:3], S[1:3, 1:3], 1)
rrrbycov <- function(Syy, Syx, Sxx, r = 1) {
  vr <- seq.int(r)
  nx <- ncol(Sxx)
  invSyy <- chol2inv(chol(Syy))
  ge <- geigen::geigen(t(Syx) %*% invSyy %*% Syx, Sxx, symmetric = TRUE)
  B <- ge$vectors[, nx + 1 - vr, drop = FALSE]
  colnames(B) <- paste0("dim", vr)
  A <- Syx %*% B %*% chol2inv(chol(t(B) %*% Sxx %*% B))
  colnames(A) <- paste0("dim", vr)
  list(A = A, B = B, eigenvalues = ge$values[nx + 1 - vr])
}
