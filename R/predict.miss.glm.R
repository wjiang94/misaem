# Predict Method for miss.glm Fits
#'
#' Prediction on test with missing values for the logistic regression model.
#' @param object a fitted object of class inheriting from "miss.glm".
#' @param newdata a data frame in which to look for variables with which to predict. It can contain missing values.
#' @param seed  An integer as a seed set for the random generator.
#' @param method The name of method to deal with missing values in test set. It can be 'map'(maximum a posteriori) or 'impute' (imputation by conditional expectation). Default is 'map'.
#' @param ... Further arguments passed to or from other methods.
#' @return
#' \item{pr.saem}{The prediction result for logistic regression: the probability of response y=1.}
#' @import mvtnorm stats MASS
#' @importFrom methods is
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
#'
#' # Generate new dataset for prediction
#' Nt <- 20
#' Xt <- matrix(rnorm(Nt*p), nrow=Nt)%*%chol(Sigma.star)+
#'   matrix(rep(mu.star,Nt), nrow=Nt, byrow = TRUE)

#' # Generate missingness in new dataset
#' patterns <- runif(Nt*p)<p.miss
#' Xt.obs <- Xt
#' Xt.obs[patterns] <- NA
#'
#' # Prediction with missing values
#' miss.prob = predict(miss.list, data.frame(Xt.obs), method='map')
#' print(miss.prob)
#' @export

predict.miss.glm <- function(object, newdata = NULL, seed = NA, method='map', ...)
{
  if (!is.na(seed))
    set.seed(seed)

  X.test = newdata
  mu.saem = object$mu.X
  sig2.saem = object$Sig.X
  beta.saem = object$coef


  #judge
  if (is(X.test, "data.frame")){
    X.test <- as.matrix(X.test)
  }
  if (method == "MAP" | method == "Map") {
    method <- "map"
  }
  if (method == "Impute" | method == "IMPUTE") {
    method <- "impute"
  }

  rindic = as.matrix(is.na(X.test))

  if(sum(rindic)!=0){
    if(method=='impute'){
      for(i in 1:dim(X.test)[1]){
        if(sum(rindic[i,])!=0){
          miss_col = which(rindic[i,]==TRUE)
          x2 = X.test[i,-miss_col]
          mu1 = mu.saem[miss_col]
          mu2 = mu.saem[-miss_col]
          sigma11 = sig2.saem[miss_col,miss_col]
          sigma12 = sig2.saem[miss_col,-miss_col]
          sigma22 = sig2.saem[-miss_col,-miss_col]
          sigma21 = sig2.saem[-miss_col,miss_col]
          mu_cond = mu1+sigma12 %*% solve(sigma22) %*% (x2-mu2)
          X.test[i,miss_col] =mu_cond
        }
      }
      tmp <- as.matrix(cbind.data.frame(rep(1,dim(X.test)[1]),X.test)) %*% as.matrix(beta.saem)
      pr.saem <- 1/(1+(1/exp(tmp)))

    }else if(method=='map'){
      pr2 =rep(0,dim(X.test)[1])
      mc.size = 100
      X.test = data.matrix(X.test)
      for(i in 1:dim(X.test)[1]){
        x=X.test[i,]
        if(sum(rindic[i,])==0){
          pr2[i]=log_reg(y=1,x=c(1,x),beta.saem,iflog=FALSE)
        }
        else{
          miss_col = which(rindic[i,]==TRUE)
          x2 = X.test[i,-miss_col]
          mu1 = mu.saem[miss_col]
          mu2 = mu.saem[-miss_col]
          sigma11 = sig2.saem[miss_col,miss_col]
          sigma12 = sig2.saem[miss_col,-miss_col]
          sigma22 = sig2.saem[-miss_col,-miss_col]
          sigma21 = sig2.saem[-miss_col,miss_col]
          mu_cond = mu1+sigma12 %*% solve(sigma22)%*%(x2-mu2)
          sigma_cond = sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
          x1_all=mvrnorm(n = mc.size, mu_cond, sigma_cond, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
          p1=0
          for(j in 1:mc.size){
            x[miss_col] =x1= x1_all[j,]
            p1 = p1 + log_reg(y=1,x=c(1,x),beta.saem,iflog=FALSE)
          }
          pr2[i] =p1/mc.size
        }
      }
      pr.saem = as.matrix(pr2)
    } else {
      stop("Error: There is no such method. Method should be 'map' or 'impute'. ")
    }
  } else {
    tmp <- as.matrix(cbind.data.frame(rep(1,dim(X.test)[1]),X.test)) %*% as.matrix(beta.saem)
    pr.saem <- 1/(1+(1/exp(tmp)))
  }


  return(pr.saem)
}
