#' Trace of a square matrix.
#'
#' @param M A square matrix.
#' @return The trace (sum of elements on the diagonal) of \code{M}.
#' @examples
#' M <- matrix(1:9, 3, 3)
#' tr(M)
#' @export
tr <- function(M) sum(diag(M))
