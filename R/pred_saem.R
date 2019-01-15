#' pred_saem
#'
#' Prediction on test with missing values for the logistic regression model.
#' @param X.test Design matrix in test set.
#' @param beta.saem Estimated \eqn{\beta}{\beta} by SAEM.
#' @param sig2.saem Estimated \eqn{\Sigma}{\Sigma} by SAEM.
#' @param mu.saem Estimated \eqn{\mu}{\mu} by SAEM.
#' @param seed  An integer as a seed set for the radom generator. The default value is 200.
#' @param method The name of method to deal with missing values in test set. It can be 'map'(maximum a posteriori) or 'impute' (imputation by conditional expectation). Default is 'map'.
#' @return
#' \item{pr.saem}{The prediction result for logistic regression: the probability of response y=1.}
#' @import mvtnorm stats MASS
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

#' # Generate missingness
#' p.miss <- 0.10
#' patterns <- runif(N*p)<p.miss #missing completely at random
#' X.obs <- X.complete
#' X.obs[patterns] <- NA
#'
#' # SAEM
#' list.saem = miss.saem(X.obs,y)
#'
#' # Generate test set with missingness
#' Nt = 50
#' X.test <- matrix(rnorm(Nt*p), nrow=Nt)%*%chol(Sigma.star)+
#'               matrix(rep(mu.star,Nt), nrow=Nt, byrow = TRUE)
#' p1 <- 1/(1+exp(-X.test%*%beta.star-beta0.star))
#' y.test <- as.numeric(runif(Nt)<p1)
#'
#' # Prediction on test set
#' pr.saem <- pred_saem(X.test, list.saem$beta, list.saem$mu, list.saem$sig2)
#' pred.saem <- (pr.saem>0.5)*1
#' table(y.test, pred.saem)
#' @export


pred_saem = function(X.test,beta.saem,mu.saem,sig2.saem,seed=200,method='map'){

  #judge
  if (class(X.test) == "data.frame") {
    X.test <- as.matrix(X.test)
  }
  if (sum(sapply(X.test, is.numeric)) < ncol(X.test)) {
    stop("Error: the variables should be numeric.")
  }
  if (method == "MAP" | method == "Map") {
    method <- "map"
  }
  if (method == "Impute" | method == "IMPUTE") {
    method <- "impute"
  }
  set.seed(seed)
  rindic = as.matrix(is.na(X.test))

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
  return(pr.saem)
}



