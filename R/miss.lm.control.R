#' Auxiliary for Controlling Fitting
#'
#' Auxiliary function for \code{\link{miss.lm}} fitting. Typically only used internally by \code{\link{miss.lm.fit}}.
#' @param maxruns maximum number of iterations. The default is maxruns = 500.
#' @param tol_em the tolerance to stop EM. The default is tol_em = 1e-4.
#' @param print_iter logical indicating if output should be produced for each iteration.
#' @return A list with components named as the arguments.
#' @examples
#' ## For examples see example(miss.lm)


miss.lm.control <- function (maxruns=500,tol_em=1e-7, print_iter=TRUE) 
{
  if (!is.numeric(tol_em) || tol_em <= 0) 
    stop("value of 'tol_em' must be > 0")
  if (!is.numeric(maxruns) || maxruns <= 0) 
    stop("maximum number of iterations must be > 0")
  # if (sum(subsets %in% 1:p) < length(subsets))  {
  #   stop("Error: index of selected variables must be in the range of covariates.")
  # }
  # if (length(unique(subsets)) != length(subsets)){
  #   stop("Error: index of selected variables must not be repeated.")
  # }
  list(tol_em = tol_em, maxruns = maxruns,  print_iter= print_iter)
}