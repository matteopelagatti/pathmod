#' Extract latent variables in a formative-reflective model
#'
#' Apply the method of Fattore, Pelagatti, Vittadini (2017) to extract latent variables in a formative-reflective model.
#' @param xgroups a list containing \eqn{p} char vectors with the names of the variables of each \eqn{x}-group.
#' @param ygroups a list containing \eqn{q} char vectors with the names of the variables of each \eqn{y}-group.
#' @param data the data frame with the observations (optional, but needed if you want the latent
#' variable scores).
#' @param cov the named covariance matrix of the variables (if you specify the \code{data} argument and
#' leave it \code{null}, then \code{cov} is assigned the correlation matrix of \code{data}).
#' @param alphas a numeric vector with values in the range [0, 1] for the \eqn{\alpha}-weights.
#' @param xmean optional parameter to center the \eqn{x}-variables taken from \code{data}, otherwise equal to their sample means
#' @param ymean optional parameter to center the \eqn{y}-variables taken from \code{data}, otherwise equal to their sample means
#' @param xsd optional parameter to scale the \eqn{x}-variables taken from \code{data}, otherwise equal to their sample standard deviations
#' @param ysd optional parameter to scale the \eqn{y}-variables taken from \code{data}, otherwise equal to their sample standard deviations
#' @return An object of class \code{pathmod} with slots:
#' \describe{
#'     \item{\code{omegas}}{list of matrices \eqn{\Omega}, one for each element in \code{alpha}}
#'     \item{\code{gammas}}{list of matrices \eqn{\Gamma}, one for each element in \code{alpha}}
#'     \item{\code{lambdas}}{list of matrices \eqn{\Lambda}, one for each element in \code{alpha}}
#'     \item{\code{xlosses}}{vector of \eqn{x}-side losses, one for each element in \code{alpha}}
#'     \item{\code{ylosses}}{vector of \eqn{y}-side losses, one for each element in \code{alpha}}
#'     \item{\code{alphas}}{same as input argument (values of \eqn{\alpha})}
#'     \item{\code{p}}{number of \eqn{x}-groups}
#'     \item{\code{q}}{number of \eqn{y}-groups}
#'     \item{\code{Sxx}}{variance-covariance matrix of exogenous variables}
#'     \item{\code{Syx}}{covariance matrix with covariances between endogenous and exogenous variables}
#'     \item{\code{Syy}}{variance-covariance matrix of endogenous variables}
#'     \item{\code{xis}}{three-dimensional array with scores of the exogenous latent variables, one for each element in \code{alpha} (available only if \code{data} is assigned)}
#'     \item{\code{etas}}{three-dimensional array with scores of the endogenous latent variables, one for each element in \code{alpha} (available only if \code{data} is assigned)}
#'     \item{\code{xgroups}}{the argument \code{xgroups}}
#'     \item{\code{ygroups}}{the argument \code{ygroups}}
#'     \item{\code{xmean}}{the supplied argument \code{xmean} or the means of the \eqn{x}-variables in \code{data}}
#'     \item{\code{ymean}}{the supplied argument \code{ymean} or the means of the \eqn{y}-variables in \code{data}}
#'     \item{\code{xsd}}{the supplied argument \code{xsd} or the standard deviations of the \eqn{x}-variables in \code{data}}
#'     \item{\code{ysd}}{the supplied argument \code{ysd} or the standard deviations of the \eqn{y}-variables in \code{data}}
#' }
#' @examples
#' data(russett)
#' xg <- list(c("gini", "farm", "rent"), c("gnpr", "labo"))
#' yg <- list(c("inst", "ecks", "death"))
#' pm <- pathmod(xg, yg, data = russett)
#' @export
#' @importFrom stats cor optim
pathmod <- function(xgroups, ygroups, data = NULL, cov = NULL, alphas = seq(0, 10)/10,
                    xmean = NULL, ymean = NULL, xsd = NULL, ysd = NULL) {
  # --- input checks ---
  if (is.null(cov) & is.null(data)) {
    stop("at least one of the arguments cov and data have to be assigned")
  }
  if (!is.null(cov)) {
    if (!is.matrix(cov) || !isSymmetric(cov))
      stop("the third argument, cov, must be a symmetric variance-covariance (usually correlation) matrix")
  }
  if (!is.null(data)) {
    if (!is.data.frame(data))
      stop("the 4th argument, data, must be a data.frame")
    if (is.null(cov)) cov <- cor(data[, c(unlist(xgroups), unlist(ygroups))])
  }
  if (!is.null(cov) & !is.null(data)) {
    # if (ncol(data) != ncol(cov))
    #   stop("the variables in the arguments data and cov must be the same (same number, same position)")
    if (is.null(colnames(cov)) | is.null(rownames(cov)))
      stop("The covariance matrix's rows and columns must be named as the variables in the data.frame")
  }
  if (!(is.numeric(alphas) && is.vector(alphas))) stop("argument alphas must be a numeric vector")
  if (min(alphas) < 0 || max(alphas) > 1) stop("argument alphas' elements must be in the range [0, 1]")
  # add checks for xmeans, ymeans, xsd, ysd

  # objects we will need
  alphas <- sort(alphas)
  p <- length(xgroups)
  q <- length(ygroups)
  vp <- sapply(xgroups, length)
  vp0 <- vp - 1
  vq <- sapply(ygroups, length)
  vq0 <- vq - 1
  nx <- sum(vp)
  ny <- sum(vq)
  cumvp <- c(0, cumsum(vp))
  cumvp0 <- c(0, cumsum(vp0))
  cumvq <- c(0, cumsum(vq))
  cumvq0 <- c(0, cumsum(vq0))
  xndx <- unlist(xgroups)
  yndx <- unlist(ygroups)
  Sxx <- cov[xndx, xndx]
  Syy <- cov[yndx, yndx]
  Syx <- cov[yndx, xndx]
  S   <- cov[c(xndx, yndx), c(xndx, yndx)]
  lSxx  <- lapply(xgroups, function(v) cov[v, v])
  lSyy  <- lapply(ygroups, function(v) cov[v, v])
  lSyix <- lapply(ygroups, function(v) cov[v , xndx])
  trx <- sapply(lSxx, tr)
  try <- sapply(lSyy, tr)
  mOmega  <- matrix(0, p, nx)
  mGamma  <- matrix(0, q, p)
  mLambda <- matrix(0, ny, q)
  vlx <- numeric(p)
  vly <- numeric(q)
  # --- objective function ---
  obj <- function(pars, alpha) {
    lx <- 0
    ly <- 0
    for (i in seq.int(p)) {
      omega <- c(1, pars[(cumvp0[i] + 1):cumvp0[i + 1]])
      mOmega[i, (cumvp[i] + 1):cumvp[i + 1]] <<- omega
      Sxx_omega <- lSxx[[i]] %*% omega
      vlx[i] <<- tr(lSxx[[i]] - Sxx_omega %*% chol2inv(chol(t(omega) %*% Sxx_omega)) %*% t(Sxx_omega)) / trx[i]
      lx <- lx +  vlx[i] / p
    }
    mOmegat <- t(mOmega)
    Sxixi <- mOmega %*% Sxx %*% mOmegat
    for (i in seq.int(q)) {
      Syxi <- lSyix[[i]] %*% mOmegat
      rrr <- rrrbycov(lSyy[[i]], Syxi, Sxixi, 1)
      gamma <- rrr$B
      mGamma[i, ] <<- gamma
      mLambda[(cumvq[i] + 1):cumvq[i + 1], i]  <<- rrr$A
      Syxi_gamma <- Syxi %*% gamma
      vly[i] <<- tr(lSyy[[i]] - Syxi_gamma %*% chol2inv(chol(t(gamma) %*% Sxixi %*% gamma)) %*% t(Syxi_gamma)) / try[i]
      ly <- ly + vly[i] / q
    }
    (1 - alpha)*lx + alpha*ly
  }
  # Solution when alpha = 0, also used for initial values
  lomega <- vector("list", p)
  pars   <- numeric(nx - p)
  for (i in seq.int(p)) {
    lomega[[i]] <- eigen(lSxx[[i]], TRUE)$vectors[, 1]
    pars[(cumvp0[i] + 1):cumvp0[i + 1]] <- lomega[[i]][-1] / lomega[[i]][1]
  }
  # Optimizations for various alphas
  optout  <- vector("list", length(alphas))
  names(optout) <- paste("alpha =", alphas)
  lOmega  <- optout
  lGamma  <- optout
  lLambda <- optout
  lxl <- optout
  lyl <- optout
  for (i in seq.int(length(alphas))) {
    cat("Working on alpha =", alphas[i], "\n")
    optout[[i]] <- optim(pars, obj, alpha = alphas[i], method = "BFGS")
    obj(optout[[i]]$par, alphas[i])
    lOmega[[i]]  <- mOmega
    lGamma[[i]]  <- mGamma
    lLambda[[i]] <- mLambda
    lxl[[i]] <- vlx
    lyl[[i]] <- vly
    pars <- optout[[i]]$par
  }
  if (!is.null(data)) {
    # centering and scaling the data
    xdata <- data[, unlist(xgroups), drop = FALSE]
    if (!is.null(xmean) && is.null(xsd)) {
      stdx <- scale(xdata, center = xmean, scale = TRUE)
      xsd <- attr(stdx, "scaled:scale")
    } else if (!is.null(xsd) && is.null(xmean)) {
      stdx <- scale(xdata, center = TRUE, scale = xsd)
      xmean <- attr(stdx, "scaled:center")
    } else if (!is.null(xmean) && !is.null(xsd)) {
      stdx <- scale(xdata, center = xmean, scale = xsd)
    } else {
      stdx <- scale(data[, unlist(xgroups)])
      xmean <- attr(stdx, "scaled:center")
      xsd <- attr(stdx, "scaled:scale")
    }
    ydata <- data[, unlist(ygroups), drop = FALSE]
    if (is.null(ymean)) ymean <- colMeans(ydata)
    if (is.null(ysd)) ysd <- sqrt(colMeans((ydata - rep(ymean, rep.int(nrow(ydata), ncol(ydata))))^2))
    # computation of xi and eta for each alpha
    xis  <- array(0, c(nrow(stdx), p, length(alphas)),
                  dimnames = list(rownames(stdx), paste0("xi", 1:p), paste("alpha =", alphas)))
    etas <- array(0, c(nrow(stdx), q, length(alphas)),
                  dimnames = list(rownames(stdx), paste0("eta", 1:q), paste("alpha =", alphas)))
    for (i in 1:length(alphas)) {
      xis[, , i]  <- stdx %*% t(lOmega[[i]])
      etas[, , i] <- xis[, , i] %*% t(lGamma[[i]])
    }
  } else {
    xis <- etas <- NULL
  }
  structure(list(omegas = lOmega,
                 gammas = lGamma,
                 lambdas = lLambda,
                 xlosses = lxl,
                 ylosses = lyl,
                 alphas = alphas,
                 p = p,
                 q = q,
                 Sxx = Sxx,
                 Syx = Syx,
                 Syy = Syy,
                 xis = xis,
                 etas = etas,
                 xgroups = xgroups,
                 ygroups = ygroups,
                 xmean = xmean,
                 ymean = ymean,
                 xsd = xsd,
                 ysd = ysd),
            class = "pathmod")
}
