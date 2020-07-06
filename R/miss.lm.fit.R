#' Fitting Linear Regression Model with Missing Values
#'
#' This function is used inside \code{miss.lm} to fit linear regression model with missing values, by EM algorithm.
#' @param x design matrix with missingness \eqn{N \times p}{N * p}.
#' @param y response vector \eqn{N \times 1}{N * 1}.
#' @param control a list of parameters for controlling the fitting process. For \code{miss.lm.fit} this is passed to \code{\link{miss.lm.control}}.
#' @return a list with following components:
#' \item{coefficients}{Estimated \eqn{\beta}{\beta}.}
#' \item{ll}{Observed log-likelihood.}
#' \item{s.resid}{Estimated standard error for residuals.}
#' \item{s.err}{Standard error for estimated parameters.}
#' \item{mu.X}{Estimated \eqn{\mu}{\mu}.}
#' \item{Sig.X}{Estimated \eqn{\Sigma}{\Sigma}.}
#' @import norm MASS mice
#' @importFrom methods is
#' @examples
#' ## For examples see example(miss.lm)

miss.lm.fit <- function (x, y,
                         control = list())

{
  control <- do.call("miss.lm.control", control)
  maxruns = control$maxruns
  tol_em = control$tol_em
  print_iter = control$print_iter

  xnames <- dimnames(x)[[2L]]
  x = as.matrix(x[,-1])

  #delete rows completely missing
  if(any(apply(is.na(x),1,sum)==p)){
    i_allNA=which(apply(is.na(x),1,sum)==p)
    x = x[-i_allNA,]
    y = y[-i_allNA]
  }
  if(any((is.na(y))==TRUE)){
    i_YNA=which(is.na(y)==TRUE)
    x = x[-i_YNA,]
    y = y[-i_YNA]
  }

  # Check parameter consistency
  if (is(x, "data.frame")){
    x <- as.matrix(x)
  }
  if (sum(sapply(x, is.numeric)) < ncol(x) ||
      sum(sapply(y, is.numeric)) < length(y)) {
    stop("Error: parameters 'x' and 'y' should be numeric.")
  }
  # Estimate extended mean and covariance using EM
  n <- nrow(x)
  p <- ncol(x) + 1
  s <- prelim.norm(cbind(y, x))
  thetahat <- em.norm(s, showits=print_iter, maxits=maxruns, criterion=tol_em)
  pars <- getparam.norm(s, thetahat)
  ll <- loglik.norm(s, thetahat)
  ptm <- Sys.time()
  # Calculate regression estimates
  b.est <- c(pars$mu[1] -
               pars$sigma[1,2:p] %*% solve(pars$sigma[2:p,2:p]) %*% pars$mu[2:p],
             pars$sigma[1,2:p] %*% solve(pars$sigma[2:p,2:p]))
  esigma.est <- as.vector(sqrt(
    pars$sigma[1,1] - b.est[2:p] %*% pars$sigma[2:p,2:p] %*% b.est[2:p]))
  Gram.est <- rbind(rep(0, p), cbind(rep(0, p - 1), pars$sigma[2:p,2:p])) +
    c(1, pars$mu[2:p]) %*% t(c(1, pars$mu[2:p]))
  sdb.est <- sqrt(diag(esigma.est^2 * solve(Gram.est * n)))
  time_run=Sys.time() - ptm


  names(b.est) <- xnames
  names(sdb.est) <- xnames

  list(coefficients = b.est, s.err = sdb.est, s.resid = esigma.est, ll= ll, Sig.X = pars$sigma, mu.X = pars$mu)
}
