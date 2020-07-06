# Predict Method for miss.lm Fits
#'
#' Prediction on test with missing values for the logistic regression model.
#' @param object a fitted object of class inheriting from "miss.lm".
#' @param newdata a data frame in which to look for variables with which to predict. It can contain missing values.
#' @param seed  An integer as a seed set for the random generator.
#' @param ... Further arguments passed to or from other methods.
#' @return
#' \item{pr.y}{The prediction result for linear regression.}
#' @import norm MASS mice
#' @importFrom methods is
#' @examples
#'
#'
#'
#'
#' # Generate complete data
#' set.seed(1)
#' mu.X <- c(1, 1)
#' Sigma.X <- matrix(c(1, 1, 1, 4), nrow = 2)
#' n <- 50 # train set size
#' p <- 2 # number of covariates
#' X.complete <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma.X) +
#'               matrix(rep(mu.X,n), nrow=n, byrow = TRUE)
#' b <- c(2, 3, -1)
#' sigma.eps <- 0.25
#' y <- cbind(rep(1, n), X.complete) %*% b +
#'   rnorm(n, 0, sigma.eps)
#'
#' # Add missing values
#' p.miss <- 0.10
#' patterns <- runif(n*p)<p.miss #missing completely at random
#' X.obs <- X.complete
#' X.obs[patterns] <- NA
#' # Estimate regression using EM
#' df.obs = data.frame(y ,X.obs)
#' miss.list = miss.lm(y~., data=df.obs)
#'
#' # Generate new dataset for prediction
#' nt <- 20
#' Xt <- matrix(rnorm(nt*p), nrow=nt)%*%chol(Sigma.X)+
#'   matrix(rep(mu.X,nt), nrow=nt, byrow = TRUE)
#' # Generate missingness in new dataset
#' patterns <- runif(nt*p)<p.miss
#' Xt.obs <- Xt
#' Xt.obs[patterns] <- NA
#'
#' # Prediction with missing values
#' miss.pred = predict(miss.list, data.frame(Xt.obs))
#' print(miss.pred)
#' @export

predict.miss.lm <- function(object, newdata = NULL, seed = NA, ...)
{
  if (!is.na(seed))
    set.seed(seed)

  X.new = newdata
  mu.em = object$mu.X
  sig2.em = object$Sig.X
  beta.em = object$coef


  # Check data consistency
  if (is(X.new, "data.frame")){
    X.new <- as.matrix(X.new)
  }
  if (!is.matrix(X.new)){
    stop("Error: parameter 'X.new' should be either a matrix or a data frame.")
  }
  if (sum(sapply(X.new, is.numeric)) < ncol(X.new)) {
    stop("Error: parameter 'X.new should be numeric'.")
  }
  # Prepare structure
  X.prep <- cbind(rep(NA, nrow(X.new)), X.new)
  X.prep <- t(t(X.prep) - mu.em)
  Inv.Sigma.tmp <- solve(sig2.em)
  # Impute
  X.pred <- t(apply(X.prep, 1, imputeEllP, Inv.Sigma.tmp))
  X.pred <- t(t(X.pred) + mu.em)

  return(pr.y=X.pred[,1])
}
