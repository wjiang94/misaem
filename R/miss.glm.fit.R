#' Fitting Logistic Regression Models with Missing Values
#'
#' This function is used inside \code{miss.glm} to fit logistic regression model with missing values, by algorithm SAEM.
#' @param x design matrix with missingness \eqn{N \times p}{N * p}.
#' @param y response vector \eqn{N \times 1}{N * 1}.
#' @param control a list of parameters for controlling the fitting process. For \code{miss.glm.fit} this is passed to \code{\link{miss.glm.control}}.
#' @return a list with following components:
#' \item{coefficients}{Estimated \eqn{\beta}{\beta}.}
#' \item{ll}{Observed log-likelihood.}
#' \item{var.covar}{Variance-covariance matrix for estimated parameters.}
#' \item{s.err}{Standard error for estimated parameters.}
#' \item{mu.X}{Estimated \eqn{\mu}{\mu}.}
#' \item{Sig.X}{Estimated \eqn{\Sigma}{\Sigma}.}
#' @import mvtnorm stats
#' @examples
#' ## For examples see example(miss.glm)

miss.glm.fit <- function (x, y,
          control = list())

{
  control <- do.call("miss.glm.control", control)
  xnames <- dimnames(x)[[2L]]
  x = as.matrix(x[,-1])

  maxruns = control$maxruns
  tol_em = control$tol_em
  nmcmc = control$nmcmc
  tau = control$tau
  k1 = control$k1
  seed = control$seed
  print_iter = control$print_iter
  var_cal = control$var_cal
  ll_obs_cal = control$ll_obs_cal
  subsets = control$subsets

  if (!is.na(seed))
    set.seed(seed)

  p=ncol(x)

  if (is.na(subsets[1]))
    subsets = 1:p

  if (sum(subsets %in% 1:p) < length(subsets))  {
    stop("Error: index of selected variables must be in the range of covariates.")
  }

  if (length(unique(subsets)) != length(subsets)){
    stop("Error: index of selected variables must not be repeated.")
  }

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

  if (sum(sapply(x, is.numeric)) < ncol(x)) {
    stop("Error: the variables should be numeric.")
  }
  if (sum(y==1) +  sum(y==0) < nrow(x)) {
    stop("Error: y must be coded by 0 or 1, and there is no missing data in y.")
  }

  n=length(y)


  rindic = as.matrix(is.na(x))
  if(sum(rindic)>0){
    whichcolmissing = (1:ncol(rindic))[apply(rindic,2,sum)>0]
    missingcols = length(whichcolmissing)
  }
  if(sum(rindic)==0){missingcols=0}


  if(missingcols>0){
    k=0
    cstop=0.1
    seqbeta = matrix(NA,nrow=ncol(x)+1,ncol=(maxruns+1))
    seqbeta_avg = matrix(NA,nrow=ncol(x)+1,ncol=(maxruns+1))

    X.mean = x
    for(i in 1:ncol(X.mean)){
      X.mean[is.na(X.mean[,i]), i] <- mean(X.mean[,i], na.rm = TRUE)
    }
    X.sim <- X.mean

    mu = apply(X.mean,2,mean)
    Sigma = var(X.mean)*(n-1)/n
    beta= rep(0,p+1)
    beta[c(1,subsets+1)]= glm(y~ X.mean[,subsets],family=binomial(link='logit'))$coef

     if(print_iter==TRUE){
      cat(sprintf('Iteration of SAEM: \n'))
    }
    while ((cstop>tol_em)*(k<maxruns)|(k<20)){
      k = k+1
      beta.old = beta

      if(k <k1){gamma <- 1}else{gamma <- 1/(k-(k1-1))^tau}

      S.inv <- solve(Sigma)

      for (i in (1:n)) {
        jna <- which(is.na(x[i,]))
        njna <- length(jna)
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
          if(cobs==0){cobs=.Machine$double.xmin}
          if(cobs==Inf){cobs=.Machine$double.xmax}

          xina <- xi[jna]
          betana <- beta[jna+1]
          for (m in (1:nmcmc)) {
            xina.c <- mi + rnorm(njna)%*%chol(Oi)

            if (y[i]==1)
              alpha <- (1+exp(-sum(xina*betana))/cobs)/(1+exp(-sum(xina.c*betana))/cobs)
            else
              alpha <- (1+exp(sum(xina*betana))*cobs)/(1+exp(sum(xina.c*betana))*cobs)
            if (runif(1) < alpha){
              xina <- xina.c
            }
          }
          X.sim[i,jna] <- xina
        }
      }
      beta_new= rep(0,p+1)
      beta_new[c(1,subsets+1)]= glm(y~ X.sim[,subsets],family=binomial(link='logit'))$coef

      beta <- (1-gamma)*beta + gamma*beta_new
      cstop = sum((beta-beta.old)^2)

      mu <- (1-gamma)*mu + gamma*colMeans(X.sim)
      Sigma <- (1-gamma)*Sigma + gamma*cov(X.sim)

      seqbeta[,k]=beta.old

      if(k==1){
        seqbeta_avg[,k]=beta.old
      }else{
        seqbeta_avg[,k]= 1/k*rowSums(seqbeta[,1:k])
      }

      if(print_iter==TRUE & k %% 50 == 0){
        cat(sprintf('%i ', k, '...'))
      }
    }
    var_obs = ll = std_obs =NULL
    if(var_cal==TRUE){
      var_obs = louis_lr_saem(beta,mu,Sigma,y,x,subsets,rindic,whichcolmissing,mc.size=100)
      std_obs <- sqrt(diag(var_obs))
    }
    if(ll_obs_cal==TRUE){
      ll= likelihood_saem(beta,mu,Sigma,y,x,rindic,whichcolmissing,mc.size=100)
    }
  }
  if(missingcols==0){
    x = matrix(x,nrow=n)
    data.complete <- data.frame(y=y,x)
    model.complete <- glm(y ~. ,family=binomial(link='logit'),data=data.complete)
    mu = apply(x,2,mean)
    Sigma = var(x)*(n-1)/n
    beta <- model.complete$coefficients
    var_obs = ll = ll1 =ll2= std_obs =seqbeta_avg= seqbeta=NULL
    if(var_cal==TRUE){
      P <- predict(model.complete, type = "response")
      W <- diag(P*(1-P))
      X <- model.matrix(model.complete)

      var_obs <- solve(t(X)%*%W%*%X)
      std_obs <- sqrt(diag(var_obs))
    }
    if(ll_obs_cal==TRUE){
      ll = likelihood_saem(beta,mu,Sigma,y,x,rindic,whichcolmissing,mc.size=100)
    }
  }

  beta = beta[c(1,subsets+1)]
  names(beta) <- xnames[c(1,subsets+1)]
  names(std_obs) <- xnames[c(1,subsets+1)]
  list(coefficients = beta, s.err = std_obs, var.covar = var_obs, ll = ll[1,1], Sig.X = Sigma, mu.X = mu
       )
}
