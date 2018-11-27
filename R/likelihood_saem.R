#' likelihood_saem
#'
#' Used in main function miss.saem. Caculate the observed log-likelihood for logistic regression model with missing data, using Monte Carlo version of Louis formula.
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
#' @export


likelihood_saem = function(beta,mu,Sigma,Y,X.obs,rindic,whichcolXmissing,mc.size){
  n=dim(X.obs)[1]
  p=dim(X.obs)[2]
  lh=0
  lh1=0
  lh2=0
  for(i in 1:n){
    y=Y[i]
    x=X.obs[i,]
    #if(sum(rindic[i,])==0){lh=lh+log_reg(y=y,x=c(1,x),beta,iflog=FALSE)+dmvnorm(x, mean = mu, sigma = Sigma, log = FALSE)}
    if(sum(rindic[i,])==0){
      lr = log_reg(y=y,x=c(1,x),beta,iflog=TRUE)
      nm = dmvnorm(x, mean = mu, sigma = Sigma, log = TRUE)
      lh=lh+lr+nm
      lh1=lh1+lr
      lh2=lh2+nm
      #if(nm <= -100){cat(sprintf('i='),i,sprintf('logLik(lr)='),lr,sprintf('logLik(nm)='),nm,'\n')}
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
      lh1 = lh1 + lr
      if(length(x2)>1){nm = dmvnorm(x2, mean = mu2, sigma = sigma22, log = TRUE)}
      if(length(x2)==1){nm = dnorm(x2, mean=mu2, sd=sqrt(sigma22), log = TRUE)}
      lh2 = lh2 + nm
      lh = lh + lr + nm
      #if(nm <= -100){cat(sprintf('i='),i,sprintf('logLik(lr)='),lr,sprintf('logLik(nm)='),nm,'\n')}
    }
  }
  return(ll=lh)
}
