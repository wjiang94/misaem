#' Summarizing Fits for miss.glm
#'
#' Summary for class \code{miss.glm}.
#' @param object an object of class "\code{miss.glm}", usually, a result of a call to \code{\link{miss.glm}}.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class "\code{summary.miss.glm}", a list with following components:
#' \item{coefficients}{The matrix of coefficients and standard errors}
#' \item{loglikelihood}{Observed log-likelihood.}
#' \item{call}{the component from \code{object}.}
#' \item{formula}{the component from \code{object}.}
#' @examples
#' ## For examples see example(miss.glm)
#' @export

summary.miss.glm <- function (object, ...){

  coef.table <- matrix(, 0L, 2L)
  dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error"))
  coef.table <- cbind(Estimate = object$coefficients, `Std. Error` = object$s.err)
  ans <- list(coefficients = coef.table, loglikelihood = object$ll, call=object$call, formula=object$formula)
  class(ans) <- "summary.miss.glm"

  return(ans)
}
