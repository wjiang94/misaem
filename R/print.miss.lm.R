#' Print miss.lm
#'
#' Print results for class \code{miss.lm}.
#' @param x an object of class "\code{miss.lm}", usually, a result of a call to \code{\link{miss.lm}}.
#' @param digits minimal number of significant digits.
#' @param ... further arguments passed to or from other methods.
#' @examples
#' ## For examples see example(miss.lm)
#' @export


print.miss.lm <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:  ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(length(coef(x))) {
    cat("Coefficients")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits),
                  print.gap = 2, quote = FALSE)
  } else cat("No coefficients\n\n")
  cat("Standard error estimates")
  if(length(x$s.err)) {
    if(is.character(co <- x$contrasts))
      cat("  [contrasts: ",
          apply(cbind(names(co),co), 1L, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(x$s.err, digits = digits),
                  print.gap = 2, quote = FALSE)
  } else cat("No standard error estimates \n\n")
  cat("Log-likelihood:", format(signif(x$ll, digits)))
  cat("\n")
  invisible(x)
}