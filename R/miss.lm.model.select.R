#' miss.lm.model.select
#'
#' Model selection for the linear regression model with missing data.
#' @param X Design matrix with missingness \eqn{N \times p}{N * p}
#' @param Y Response vector \eqn{N \times 1}{N * 1}
#' @return An object of class "\code{miss.lm}".
#' @import mvtnorm stats
#' @importFrom methods is
#' @examples
#' # Generate complete data
#' set.seed(1)
#' mu.X <- c(1, 1)
#' Sigma.X <- matrix(c(1, 1, 1, 4), nrow = 2)
#' n <- 50
#' p <- 2
#' X.complete <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma.X) +
#'               matrix(rep(mu.X,n), nrow=n, byrow = TRUE)
#' b <- c(2, 0, -1)
#' sigma.eps <- 0.25
#' y <- cbind(rep(1, n), X.complete) %*% b + rnorm(n, 0, sigma.eps)
#'
#' # Add missing values
#' p.miss <- 0.10
#' patterns <- runif(n*p)<p.miss #missing completely at random
#' X.obs <- X.complete
#' X.obs[patterns] <- NA
#'
#' # model selection
#' miss.model = miss.lm.model.select(y, X.obs)
#' print(miss.model)
#' @export


miss.lm.model.select = function(Y, X){
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
  subsets=combinations(p)
  ll = matrix(-Inf,nrow=p,ncol=p)
  subsets1 = subsets[rowSums(subsets)==1,]
  for (j in 1:(nrow(subsets1))){
    pos_var=which(subsets1[j,]==1)
    df = data.frame(Y,X[,pos_var])
    names(df) = c('Y', xnames[pos_var])
    list.em.subset=miss.lm(Y~., data=df, print_iter=FALSE)
    ll[1,j] = list.em.subset$ll
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
        df = data.frame(Y,X[,pos_var])
        names(df) = c('Y', xnames[pos_var])
        list.em.subset=miss.lm(Y~., data=df, print_iter=FALSE)
        ll[i,j] = list.em.subset$ll
      }
    }
  }
  SUBSETS[p,]=rep(1,p)
  df = data.frame(Y,X)
  names(df) = c('Y', xnames)
  list.em.subset=miss.lm(Y~., data=df, print_iter=FALSE)
  ll[p,1] = list.em.subset$ll
  nb.x = p
  nb.para = (nb.x + 1) + p + p*p
  BIC[p] = -2*ll[p,1]+ nb.para * log(N)
  subset_choose = which(SUBSETS[which.min(BIC),]==1)
  df = data.frame(Y,X[,subset_choose])
  names(df) = c('Y', xnames[subset_choose])
  list.em.subset=miss.lm(Y~., data=df, print_iter=FALSE)
  return(miss.model = list.em.subset)
}
