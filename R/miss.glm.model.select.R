#' miss.glm.model.select
#'
#' Model selection for the logistic regression model with missing data.
#' @param X Design matrix with missingness \eqn{N \times p}{N * p}
#' @param Y Binary response vector \eqn{N \times 1}{N * 1}
#' @param seed  An integer as a seed set for the random generator. The default value is 200.
#' @return An object of class "\code{miss.glm}".
#' @import mvtnorm stats
#' @importFrom methods is
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
#' Y <- as.numeric(runif(N)<p1)

#' # Generate missingness
#' p.miss <- 0.10
#' patterns <- runif(N*p)<p.miss #missing completely at random
#' X <- X.complete
#' X[patterns] <- NA

#' # model selection for SAEM
#' miss.model = miss.glm.model.select(Y,X,seed=100)
#' print(miss.model)
#' @export


miss.glm.model.select = function(Y, X, seed=NA){

  if (!is.na(seed))
    set.seed(seed)

  N=dim(X)[1]
  p=dim(X)[2]
  xnames <- dimnames(data.frame(X))[[2L]]

  #judge
  if (is(X, "data.frame")){
    X <- as.matrix(X)
  }
  if (sum(sapply(X, is.numeric)) < ncol(X)) {
    stop("Error: the variables should be numeric.")
  }
  if (sum(Y==1) +  sum(Y==0) < nrow(X)) {
    stop("Error: Y must be coded by 0 or 1, and there is no missing data in Y.")
  }

  df = data.frame(Y, X)
  names(df) = c('Y', xnames)

  subsets=combinations(p)
  ll = matrix(-Inf,nrow=p,ncol=p)
  subsets1 = subsets[rowSums(subsets)==1,]
  for (j in 1:(nrow(subsets1))){
    pos_var=which(subsets1[j,]==1)
    list.saem.subset=miss.glm(Y~., data=df,print_iter=FALSE,subsets=pos_var)
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
        list.saem.subset=miss.glm(Y~., data=df,print_iter=FALSE,subsets=pos_var)
        ll[i,j] = list.saem.subset$ll
      }
    }
  }
  SUBSETS[p,]=rep(1,p)
  list.saem.subset=miss.glm(Y~., data=df,print_iter=FALSE,subsets=1:p)
  ll[p,1] = list.saem.subset$ll
  nb.x = p
  nb.para = (nb.x + 1) + p + p*p
  BIC[p] = -2*ll[p,1]+ nb.para * log(N)
  subset_choose = which(SUBSETS[which.min(BIC),]==1)
  list.saem.subset=miss.glm(Y~., data=df,print_iter=FALSE,subsets=subset_choose)
  return(miss.model = list.saem.subset)
}
