#' Auxiliary for Controlling Fitting
#'
#' Auxiliary function for \code{\link{miss.glm}} fitting. Typically only used internally by \code{\link{miss.glm.fit}}.
#' @param maxruns maximum number of iterations. The default is maxruns = 500.
#' @param tol_em the tolerance to stop SAEM. The default is tol_em = 1e-7.
#' @param nmcmc the MCMC length. The default is nmcmc = 2.
#' @param tau rate \eqn{\tau}{\tau} in the step size \eqn{(k-k_{1})^{-\tau}}{(k-k1)^(-\tau)}. The default is tau = 1.
#' @param k1 number of first iterations \eqn{k_{1}}{k1} in the step size \eqn{(k-k_{1})^{-\tau}}{(k-k1)^(-\tau)}. The default is k1=50.
#' @param seed an integer as a seed set for the random generator.
#' @param subsets Index of selected covariates if any. The default is all the covariates.
#' @param print_iter logical indicating if output should be produced for each iteration.
#' @param var_cal logical indicating if the variance of estimated parameters should be calculated.
#' @param ll_obs_cal logical indicating if the observed log-likelihood should be calculated.
#' @return A list with components named as the arguments.
#' @examples
#' ## For examples see example(miss.glm)


miss.glm.control <- function (maxruns=500,tol_em=1e-7,nmcmc=2,tau=1,k1=50, subsets =NA, seed=NA, print_iter=TRUE, var_cal=TRUE, ll_obs_cal=TRUE)
{
  if (!is.numeric(tol_em) || tol_em <= 0)
    stop("value of 'tol_em' must be > 0")
  if (!is.numeric(maxruns) || maxruns <= 0)
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(nmcmc) || nmcmc <= 0)
    stop("Number of MCMC in each iteration must be > 0")
  if (!is.numeric(tau) || tau <= 0)
    stop("Step-size rate must be > 0")
  if (!is.numeric(k1) || k1 <= 0)
    stop("Number of first k1 iterations must be > 0")
  # if (sum(subsets %in% 1:p) < length(subsets))  {
  #   stop("Error: index of selected variables must be in the range of covariates.")
  # }
  # if (length(unique(subsets)) != length(subsets)){
  #   stop("Error: index of selected variables must not be repeated.")
  # }
  list(tol_em = tol_em, maxruns = maxruns, nmcmc = nmcmc, tau = tau, k1=k1, seed = seed, subsets = subsets, print_iter= print_iter, var_cal=var_cal, ll_obs_cal=ll_obs_cal)
}
