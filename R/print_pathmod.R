#' Print method for fitted formative-reflective models
#'
#' Print the most important results created by the function \code{pathmod}.
#' @param x an object of class \code{pathmod}
#' @param round number of decimal places for the numbers to be printed (mostly correlations and \eqn{R^2})
#' @param ... other parameters that are passed to the standard print function used to print the matrices
#' @return It only prints text, it returns nothing.
#' @export
print.pathmod <- function(x, round = 3, ...) {
  # matrix with coefficients of determination
  M <- matrix(0, length(x$alphas), x$p + x$q + 4,
              dimnames = list(rep("", length(x$alphas)),
                              c("alpha", "X-side", "Y-side", "R2(alpha)", paste0("x", 1:x$p), paste0("y", 1:x$q))))
  M[, 1] <- x$alphas
  for (i in seq.int(length(x$alphas))) {
    M[i, -1] <- 1 - c(mean(x$xlosses[[i]]), mean(x$ylosses[[i]]),
                      (1 - x$alphas[i])*mean(x$xlosses[[i]]) + x$alphas[i]*mean(x$ylosses[[i]]),
                      x$xlosses[[i]], x$ylosses[[i]])
  }
  # matrix with correlation between LVs and MVs and between LVs and LVs
  cor_xi_x <- array(0, c(x$p, ncol(x$Sxx), length(x$alphas)),
                    dimnames = list(paste0("xi", 1:x$p), colnames(x$Sxx), paste("alpha = ", x$alphas)))
  cor_eta_xi <- array(0, c(x$q, x$p, length(x$alphas)),
                      dimnames = list(paste0("eta", 1:x$q), paste0("xi", 1:x$p), paste("alpha = ", x$alphas)))
  cor_y_eta <- array(0, c(ncol(x$Syy), x$q, length(x$alphas)),
                     dimnames = list(colnames(x$Syy), paste0("eta", 1:x$q), paste("alpha =", x$alphas)))
  for (i in seq.int(length(x$alphas))) {
    OmSxx       <- x$omegas[[i]] %*% x$Sxx
    OmSxxOmt    <- OmSxx %*% t(x$omegas[[i]])
    invDiagOSOt <- diag(sqrt(diag(OmSxxOmt))^(-1))
    invDiagGOSOtGt <- diag(sqrt(diag(x$gammas[[i]] %*% OmSxxOmt %*% t(x$gammas[[i]])))^(-1),
                           nrow(x$gamma[[i]]), nrow(x$gamma[[i]]))
    cor_xi_x[, , i] <- invDiagOSOt %*% OmSxx %*% diag(sqrt(diag(x$Sxx))^(-1))
    cor_eta_xi[, , i] <- invDiagGOSOtGt %*% x$gammas[[i]] %*% OmSxxOmt %*% invDiagOSOt
    cor_y_eta[, , i]  <- invDiagGOSOtGt %*% x$gammas[[i]] %*% x$omegas[[i]] %*% t(x$Syx) %*% diag(sqrt(diag(x$Syy))^(-1))
  }
  cat("------- pathmod output -------\n")
  cat("\nCoefficients of determination\n")
  print(round(M, round), ...)
  cat("\nCorrelations between exogenous manifest and latent variables\n")
  print(round(cor_xi_x, round), ...)
  cat("\nCorrelations between latent variables\n")
  print(round(cor_eta_xi, round), ...)
  cat("\nCorrelations between endogenous latent and manifest variables\n")
  print(round(cor_y_eta, round), ...)
}
