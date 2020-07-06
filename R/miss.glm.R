#' Statistical Inference for Logistic Regression Models with Missing Values
#'
#' This function is used to perform statistical inference for logistic regression model with missing values, by algorithm SAEM.
#' @param formula an object of class "\code{\link[stats]{formula}}" : a symbolic description of the logistic regression model to be fitted.
#' @param data an optional data frame containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment}(\code{formula}), typically the environment from which \code{miss.glm} is called.
#' @param control a list of parameters for controlling the fitting process. For \code{miss.glm.fit} this is passed to \code{\link{miss.glm.control}}.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return An object of class "\code{miss.glm}": a list with following components:
#' \item{coefficients}{Estimated \eqn{\beta}{\beta}.}
#' \item{ll}{Observed log-likelihood.}
#' \item{var.covar}{Variance-covariance matrix for estimated parameters.}
#' \item{s.err}{Standard error for estimated parameters.}
#' \item{mu.X}{Estimated \eqn{\mu}{\mu}.}
#' \item{Sig.X}{Estimated \eqn{\Sigma}{\Sigma}.}
#' \item{call}{the matched call.}
#' \item{formula}{the formula supplied.}
#' @import mvtnorm stats
#' @examples
#' # Generate dataset
#' N <- 100  # number of subjects
#' p <- 3     # number of explanatory variables
#' mu.star <- rep(0,p)  # mean of the explanatory variables
#' Sigma.star <- diag(rep(1,p)) # covariance
#' beta.star <- c(1, 1,  0) # coefficients
#' beta0.star <- 0 # intercept
#' beta.true = c(beta0.star,beta.star)
#' X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star) +
#'               matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
#' p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
#' y <- as.numeric(runif(N)<p1)
#'
#' # Generate missingness
#' p.miss <- 0.10
#' patterns <- runif(N*p)<p.miss #missing completely at random
#' X.obs <- X.complete
#' X.obs[patterns] <- NA
#'
#' df.obs = data.frame(y,X.obs)
#'
#' # SAEM
#' miss.list = miss.glm(y~., data=df.obs, print_iter=FALSE,seed=100)
#' print(miss.list)
#' print(summary(miss.list))
#' summary(miss.list)$coef
#' @export

miss.glm <- function (formula,  data,
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
  control <- do.call("miss.glm.control", control)
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

  fit <- eval(call('miss.glm.fit',
                   x = X, y = Y,
                   control = control))
  fit <- c(fit, list(call = call, formula = formula))
  class(fit) <- "miss.glm"
  return(fit)
}
