#' Summarizing Fits for miss.lm
#'
#' Summary for class \code{miss.lm}.
#' @param object an object of class "\code{miss.lm}", usually, a result of a call to \code{\link{miss.lm}}.
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class "\code{summary.miss.lm}", a list with following components:
#' \item{coefficients}{The matrix of coefficients and standard errors.}
#' \item{loglikelihood}{Observed log-likelihood.}
#' \item{call}{the component from \code{object}.}
#' \item{formula}{the component from \code{object}.}
#' @examples
#' ## For examples see example(miss.lm)
#' @export

summary.miss.lm <- function (object, ...){

  coef.table <- matrix(, 0L, 2L)
  dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error"))
  coef.table <- cbind(Estimate = object$coefficients, `Std. Error` = object$s.err)
  ans <- list(coefficients = coef.table, loglikelihood = object$ll, call=object$call, formula=object$formula)
  class(ans) <- "summary.miss.lm"

  return(ans)
}
