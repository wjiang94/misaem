#' Print Summary of miss.lm
#'
#' Print results for class \code{summary.miss.lm}.
#' @param x an object of class "\code{summary.miss.lm}", usually, a result of a call to \code{\link{summary.miss.lm}}.
#' @param digits minimal number of significant digits.
#' @param ... further arguments passed to or from other methods.
#' @examples
#' ## For examples see example(miss.lm)
#' @export
print.summary.miss.lm <-
  function (x, digits = max(3L, getOption("digits") - 3L), ...)
  {
    cat("\nCall:\n",
        paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
    if(length(coef(x))) {
      cat("Coefficients")
      cat(":\n")
      print.default(format(x$coefficients, digits = digits),
                    print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n\n")

    cat("Log-likelihood: ", format(x$loglikelihood, digits = max(4L, digits + 1L)),
        "\n", sep = "")
    cat("\n")
    invisible(x)
  }
