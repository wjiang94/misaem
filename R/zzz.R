#' @importFrom utils packageDescription
#' @noRd
.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    pdesc <- packageDescription(pkgname)
    packageStartupMessage('')
    packageStartupMessage(pdesc$Package,
                          " "
                          , pdesc$Version,
                          " by "
                          ,pdesc$Author)
    packageStartupMessage(paste0('-> For help, please type: vignette(\''
                                 ,pkgname,
                                 '\')'))
    packageStartupMessage('')
  }}


#' louis_lr_saem
#'
#' Used in main function miss.saem. Calculate the variance of estimated parameters for logistic regression model with missing data, using Monte Carlo version of Louis formula.
#' @param beta Estimated parameter of logistic regression model.
#' @param mu Estimated parameter \eqn{\mu}.
#' @param Sigma Estimated parameter \eqn{\Sigma}.
#' @param Y Response vector \eqn{N \times 1}{N * 1}
#' @param X.obs Design matrix with missingness \eqn{N \times p}{N * p}
#' @param pos_var Index of selected covariates.
#' @param rindic Missing pattern of X.obs. If a component in X.obs is missing, the corresponding position in rindic is 1; else 0.
#' @param whichcolXmissing The column index in covariate containing at least one missing observation.
#' @param mc.size Monte Carlo sampling size.
#' @import stats
#' @return Variance of estimated \eqn{\beta}.
#' @examples
#' # Generate dataset
#' N <- 50  # number of subjects
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
#' # Louis formula to obtain variance of estimates
#' V_obs = louis_lr_saem(beta.true,mu.star,Sigma.star,y,X.obs)
#' @export
louis_lr_saem = function(beta,mu,Sigma,Y,X.obs,pos_var=1:ncol(X.obs),rindic=as.matrix(is.na(X.obs)),whichcolXmissing=(1:ncol(rindic))[apply(rindic,2,sum)>0],mc.size=2){
  n=dim(X.obs)[1]
  p = length(pos_var)
  beta = beta[c(1,pos_var+1)]
  mu = mu[pos_var]
  Sigma = Sigma[pos_var,pos_var]
  X.obs = as.matrix(X.obs[,pos_var])
  rindic = as.matrix(rindic[,pos_var])

  X.mean = as.matrix(X.obs)
  for(i in 1:ncol(X.mean)){X.mean[is.na(X.mean[,i]), i] <- mean(X.mean[,i], na.rm = TRUE)}
  X.sim <- X.mean
  G = D = I_obs = matrix(0,ncol = p+1,nrow = p+1)
  Delta = matrix (0,ncol=1,nrow=p+1)
  S.inv <- solve(Sigma)

  for (i in (1:n)) {
    jna <- which(is.na(X.obs[i,]))
    njna <- length(jna)
    if (njna==0) {
      x = matrix(c(1,X.sim[i,]))
      exp_b=exp(beta%*%x)[1]
      d2l = -x%*%t(x)*(exp_b/(1+exp_b)^2)
      I_obs = I_obs - d2l}
    if (njna>0) {
      xi <- X.sim[i,]
      Oi <- solve(S.inv[jna,jna])
      mi <- mu[jna]
      lobs <- beta[1]
      if (njna<p) {
        jobs <- setdiff(1:p,jna)
        mi <- mi - (xi[jobs] - mu[jobs])%*%S.inv[jobs,jna]%*%Oi
        lobs <- lobs + sum(xi[jobs]*beta[jobs+1])
      }
      cobs <- exp(lobs)
      xina <- xi[jna]
      betana <- beta[jna+1]
      for (m in (1:mc.size)) {
        xina.c <- mi + rnorm(njna)%*%chol(Oi)
        if (Y[i]==1)
          alpha <- (1+exp(-sum(xina*betana))/cobs)/(1+exp(-sum(xina.c*betana))/cobs)
        else
          alpha <- (1+exp(sum(xina*betana))*cobs)/(1+exp(sum(xina.c*betana))*cobs)
        if (runif(1) < alpha){
          xina <- xina.c
        }
        X.sim[i,jna] <- xina
        x = matrix(c(1,X.sim[i,]))
        exp_b=exp(beta%*%x)[1]
        dl = x*(Y[i]-exp_b/(1+exp_b))
        d2l = -x%*%t(x)*(exp_b/(1+exp_b)^2)
        D = D + 1/m * ( d2l - D)
        G = G + 1/m * ( dl %*% t(dl) - G)
        Delta = Delta + 1/m * ( dl - Delta)
      }
      I_obs = I_obs - (D+G-Delta%*%t(Delta))
    }
  }
  V_obs = solve(I_obs)
  return(V_obs)
}


#' log_reg
#'
#' Calculate the likelihood or log-likelihood for one observation of logistic regression model .
#' @param y Response value (0 or 1).
#' @param x Covariate vector of dimension \eqn{p \times 1}{p*1}.
#' @param beta Estimated parameter of logistic regression model.
#' @param iflog If TRUE, log_reg calculate the log-likelihood; else likelihood.
#' @return Likelihood or log-likelihood.
#' @examples
#' res = log_reg(1,c(1,2,3),c(1,-1,1))
#' @export
log_reg <- function(y,x,beta,iflog=TRUE){
  res <- y*(x%*%beta) - log(1+exp(x%*%beta))
  if(iflog==TRUE)
    return(res)
  else
    return(exp(res))
}


#' likelihood_saem
#'
#' Used in main function miss.saem. Calculate the observed log-likelihood for logistic regression model with missing data, using Monte Carlo version of Louis formula.
#' @param beta Estimated parameter of logistic regression model.
#' @param mu Estimated parameter \eqn{\mu}.
#' @param Sigma Estimated parameter \eqn{\Sigma}.
#' @param Y Response vector \eqn{N \times 1}{N * 1}
#' @param X.obs Design matrix with missingness \eqn{N \times p}{N * p}
#' @param rindic Missing pattern of X.obs. If a component in X.obs is missing, the corresponding position in rindic is 1; else 0.
#' @param whichcolXmissing The column index in covariate containing at least one missing observation.
#' @param mc.size Monte Carlo sampling size.
#' @import mvtnorm stats MASS
#' @return Observed log-likelihood.
#' @examples
#' # Generate dataset
#' N <- 50  # number of subjects
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
#' # Observed log-likelihood
#' ll_obs = likelihood_saem(beta.true,mu.star,Sigma.star,y,X.obs)
#' @export

likelihood_saem = function(beta,mu,Sigma,Y,X.obs,rindic=as.matrix(is.na(X.obs)),whichcolXmissing=(1:ncol(rindic))[apply(rindic,2,sum)>0],mc.size=2){
  n=dim(X.obs)[1]
  p=dim(X.obs)[2]
  lh=0
  for(i in 1:n){
    y=Y[i]
    x=X.obs[i,]
    if(sum(rindic[i,])==0){
      lr = log_reg(y=y,x=c(1,x),beta,iflog=TRUE)
      lh=lh+lr
    }
    else{
      miss_col = which(rindic[i,]==TRUE)
      x2 = x[-miss_col]
      mu1 = mu[miss_col]
      mu2 = mu[-miss_col]
      sigma11 = Sigma[miss_col,miss_col]
      sigma12 = Sigma[miss_col,-miss_col]
      sigma22 = Sigma[-miss_col,-miss_col]
      sigma21 = Sigma[-miss_col,miss_col]
      mu_cond = mu1+sigma12 %*% solve(sigma22)%*%(x2-mu2)
      sigma_cond = sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
      x1_all=mvrnorm(n = mc.size, mu_cond, sigma_cond, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
      lh_mis1=0
      for(j in 1:mc.size){
        x[miss_col] =x1= x1_all[j,]
        lh_mis1 = lh_mis1 + log_reg(y=y,x=c(1,x),beta,iflog=FALSE)
      }
      lr = log(lh_mis1/mc.size)
      lh = lh + lr
    }
  }
  return(ll=lh)
}



#' combinations
#'
#' Given all the possible patterns of missingness.
#' @param p Dimension of covariates.
#' @return A matrix containing all the possible missing patterns. Each row indicates a pattern of missingness. "1" means "observed", 0 means "missing".
#' @examples
#' comb = combinations(5)
#' @export

combinations = function(p){
  comb = NULL
  if (p<20) {
    for( i in 1:p) comb = rbind(cbind(1,comb),cbind(0,comb))
    return(comb)
  }
  else {stop("Error: the dimension of dataset is too large to possibly block your computer. Better try with number of variables smaller than 20.")}
}



#' Function for imputing single point for linear regression model
#' @param point A single observation containing missing values.
#' @param Sigma.inv Inverse of estimated \eqn{\Sigma}{\Sigma}.
#' @return Imputed observation.

imputeEllP <- function(point, Sigma.inv){
  point.new <- point
  d <- length(point)
  index <- which(is.na(point))
  if (length(index) == 1){
    point.new[index] <- -point.new[-index] %*% Sigma.inv[index,-index] /
      Sigma.inv[index,index]
  }else{
    A <- Sigma.inv[index,index]
    b <- -Sigma.inv[index,(1:d)[-index], drop = FALSE] %*% point[-index]
    point.new[index] <- solve(A) %*% b
  }
  return (point.new)
}
