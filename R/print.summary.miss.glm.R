#' Print Summary of miss.glm
#'
#' Print results for class \code{summary.miss.glm}.
#' @param x an object of class "\code{summary.miss.glm}", usually, a result of a call to \code{\link{summary.miss.glm}}.
#' @param digits minimal number of significant digits.
#' @param ... further arguments passed to or from other methods.
#' @examples
#' ## For examples see example(miss.glm)
#' @export
print.summary.miss.glm <-
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
