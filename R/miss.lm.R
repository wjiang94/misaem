#' Statistical Inference for Linear Regression Models with Missing Values
#'
#' This function is used to perform statistical inference for linear regression model with missing values, by algorithm EM.
#' @param formula an object of class "\code{\link[stats]{formula}}" : a symbolic description of the linear regression model to be fitted.
#' @param data an optional data frame containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment}(\code{formula}), typically the environment from which \code{miss.lm} is called.
#' @param control a list of parameters for controlling the fitting process. For \code{miss.lm.fit} this is passed to \code{\link{miss.lm.control}}.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return An object of class "\code{miss.lm}": a list with following components:
#' \item{coefficients}{Estimated \eqn{\beta}{\beta}.}
#' \item{ll}{Observed log-likelihood.}
#' \item{s.resid}{Estimated standard error for residuals.}
#' \item{s.err}{Standard error for estimated parameters.}
#' \item{mu.X}{Estimated \eqn{\mu}{\mu}.}
#' \item{Sig.X}{Estimated \eqn{\Sigma}{\Sigma}.}
#' \item{call}{the matched call.}
#' \item{formula}{the formula supplied.}
#' @import mvtnorm stats mice
#' @examples
#' # Generate complete data
#' set.seed(1)
#' mu.X <- c(1, 1)
#' Sigma.X <- matrix(c(1, 1, 1, 4), nrow = 2)
#' n <- 50
#' p <- 2
#' X.complete <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma.X) +
#'               matrix(rep(mu.X,n), nrow=n, byrow = TRUE)
#' b <- c(2, 3, -1)
#' sigma.eps <- 0.25
#' y <- cbind(rep(1, n), X.complete) %*% b + rnorm(n, 0, sigma.eps)
#'
#' # Add missing values
#' p.miss <- 0.10
#' patterns <- runif(n*p)<p.miss #missing completely at random
#' X.obs <- X.complete
#' X.obs[patterns] <- NA
#'
#' # Estimate regression using EM
#' df.obs = data.frame(y,X.obs)
#' miss.list = miss.lm(y~., data=df.obs)
#' print(miss.list)
#' print(summary(miss.list))
#' summary(miss.list)$coef
#' @export

miss.lm <- function (formula,  data,
                      control = list(...), ...
)
{
  call <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass #otherwise X will omit NA
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  control <- do.call("miss.lm.control", control)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf)
  else matrix(,NROW(Y), 0L)

  fit <- eval(call('miss.lm.fit',
                   x = X, y = Y,
                   control = control))
  fit <- c(fit, list(call = call, formula = formula))
  class(fit) <- "miss.lm"
  return(fit)
}
