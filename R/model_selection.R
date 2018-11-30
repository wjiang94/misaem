#' model_selection
#'
#' Model selection for the logistic regression model with missing data.
#' @param X.obs Design matrix with missingness \eqn{N \times p}{N * p}
#' @param y Response vector \eqn{N \times 1}{N * 1}
#' @param seed  An integer as a seed set for the radom generator. The default value is 200.
#' @return A list with components
#' \item{subset_choose}{The index of variates included in the best model selected.}
#' \item{beta}{Estimated \eqn{\beta}{\beta} for the best model.}
#' \item{sig2}{Estimated \eqn{\Sigma}{\Sigma} for the best model.}
#' \item{mu}{Estimated \eqn{\mu}{\mu} for the best model.}
#' @import mvtnorm stats
#' @examples
#' # Generate dataset
#' N <- 40  # number of subjects
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

#' # model selection for SAEM
#' list.saem.select = model_selection(X.obs,y)
#' print(list.saem.select$subset_choose)
#' print(list.saem.select$beta)
#' @export


model_selection = function(X.obs,y,seed=200){
  #judge
  if (class(X.obs) == "data.frame") {
    X.obs <- as.matrix(X.obs)
  }
  if (sum(sapply(X.obs, is.numeric)) < ncol(X.obs)) {
    stop("Error: the variables should be numeric.")
  }
  if (sum(y==1) +  sum(y==0) < nrow(X.obs)) {
    stop("Error: y must be coded by 0 or 1, and there is no missing data in y.")
  }

  set.seed(seed)
  N=dim(X.obs)[1]
  p=dim(X.obs)[2]
  subsets=combinations(p)
  ll = matrix(-Inf,nrow=p,ncol=p)
  subsets1 = subsets[rowSums(subsets)==1,]
  for (j in 1:(nrow(subsets1))){
    pos_var=which(subsets1[j,]==1)
    list.saem.subset=miss.saem(data.matrix(X.obs),y,pos_var,maxruns=500,tol_em=1e-7,nmcmc=2,print_iter=FALSE,ll_obs_cal=TRUE)
    ll[1,j] = list.saem.subset$ll
  }
  id = BIC = rep(0,p)
  subsetsi=subsets1
  SUBSETS = matrix(-Inf,nrow=p,ncol=p)
  for(i in 2:p){
    nb.x = i-1
    nb.para = (nb.x + 1) + p + p*p
    id[i-1] = d = which.max(ll[i-1,])
    pos_var=which(subsetsi[d,]==1)
    BIC[i-1] = -2*ll[i-1,d]+ nb.para * log(N)
    SUBSETS[i-1,]=subsetsi[d,]
    if(i==2){subsetsi = subsets[(rowSums(subsets)==i) & (subsets[,pos_var]==i-1),]}
    if(i>2){subsetsi = subsets[(rowSums(subsets)==i) & (rowSums(subsets[,pos_var])==i-1),]}
    if(i<p){
      for (j in 1:(nrow(subsetsi))){
        pos_var=which(subsetsi[j,]==1)
        list.saem.subset=miss.saem(data.matrix(X.obs),y,pos_var,maxruns=500,tol_em=1e-7,nmcmc=2,print_iter=FALSE,ll_obs_cal=TRUE)
        ll[i,j] = list.saem.subset$ll
      }
    }
  }
  SUBSETS[p,]=rep(1,p)
  list.saem.subset=miss.saem(data.matrix(X.obs),y,1:p,maxruns=500,tol_em=1e-7,nmcmc=2,print_iter=FALSE,ll_obs_cal=TRUE)
  ll[p,1] = list.saem.subset$ll
  nb.x = p
  nb.para = (nb.x + 1) + p + p*p
  BIC[p] = -2*ll[p,1]+ nb.para * log(N)
  subset_choose = which(SUBSETS[which.min(BIC),]==1)
  list.saem.subset=miss.saem(data.matrix(X.obs),y,subset_choose,maxruns=500,tol_em=1e-7,k1=2,nmcmc=2,print_iter=TRUE)
  beta.saem.train = list.saem.subset$beta
  mu.saem.train = list.saem.subset$mu
  sig2.saem.train = list.saem.subset$sig2
  return(list(subset_choose=subset_choose,beta=beta.saem.train,sig2=sig2.saem.train,mu=mu.saem.train))
}
