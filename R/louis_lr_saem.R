#' louis_lr_saem
#'
#' Caculate the variance of estimated parameters for logistic regression model with missing data, using Monte Carlo version of Louis formula.
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
#' @export


louis_lr_saem = function(beta,mu,Sigma,Y,X.obs,pos_var,rindic,whichcolXmissing,mc.size){
  n=dim(X.obs)[1]
  #p=dim(X.obs)[2]
  p = length(pos_var)
  beta = beta[c(1,pos_var+1)]
  mu = mu[pos_var]
  Sigma = Sigma[pos_var,pos_var]
  X.obs = X.obs[,pos_var]
  rindic = rindic[,pos_var]

  X.mean = X.obs
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
