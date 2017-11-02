#' Compute rotation of latent variables in fitted formative-reflective models
#'
#' It computes the correlations (rotation) among latent variables as alpha changes.
#' @param pm an object of class \code{pathmod}
#' @return An object of class \code{lvcor} with slots:
#' \describe{
#'     \item{\code{alphas}}{a numeric vector with the values of \eqn{\alpha}}
#'     \item{\code{xicor}}{\eqn{p} x \eqn{p} x \code{length(alphas)} array with the correlation between the exogenous latent
#'     variables for \eqn{\alpha = 0} and the other values of \eqn{\alpha}}
#'     \item{\code{etacor}}{\eqn{q} x \eqn{q} x \code{length(alphas)} array with the correlation between the endogenous latent
#'     variables for \eqn{\alpha = 0} and the other values of \eqn{\alpha}}
#' }
#' @export
lvcor <- function(pm) {
  na <- length(pm$alphas)
  if (na < 2) stop("I need a pathmod object with at least two values in the argument alphas")

  OmSxx <- pm$omegas[[1]] %*% pm$Sxx
  invOmSxxOmt <- diag(sqrt(diag(OmSxx %*% t(pm$omegas[[1]])))^(-1))
  xicor <- array(0, c(pm$p, pm$p, na),
                 dimnames = list(paste0("xi", 1:pm$p, "(0)"), paste0("xi", 1:pm$p, "(alpha)"), paste("alpha =", pm$alphas)))
  for (i in seq.int(na)) {
    xicor[, , i] <- invOmSxxOmt %*% OmSxx %*% t(pm$omegas[[i]]) %*%
      diag(sqrt(diag(pm$omegas[[i]] %*% pm$Sxx %*% t(pm$omegas[[i]])))^(-1), pm$p, pm$p)
  }

  GaOmSxx <-pm$gammas[[1]] %*%  pm$omegas[[1]] %*% pm$Sxx
  invGaOmSxxOmGat <- diag(sqrt(diag(GaOmSxx %*% t(pm$gammas[[1]] %*% pm$omegas[[1]])))^(-1))
  etacor <- array(0, c(pm$q, pm$q, na),
                  dimnames = list(paste0("eta", 1:pm$q, "(0)"), paste0("eta", 1:pm$q, "(alpha)"), paste("alpha =", pm$alphas)))
  for (i in seq.int(na)) {
    etacor[, , i] <- invGaOmSxxOmGat %*% GaOmSxx %*% t(pm$gammas[[i]] %*% pm$omegas[[i]]) %*%
      diag(sqrt(diag(pm$gammas[[i]] %*% pm$omegas[[i]] %*% pm$Sxx %*% t(pm$gammas[[i]] %*% pm$omegas[[i]])))^(-1), pm$q, pm$q)
  }
  structure(list(alphas = pm$alphas, xicor = xicor, etacor = etacor), class = "lvcor")
}
