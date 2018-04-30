N <- 1000  # number of subjects
p <- 5     # number of explanatory variables
mu.star <- 1:p  #rep(0,p)  # mean of the explanatory variables
sd <- 1:p # rep(1,p) # standard deviations
C <- matrix(c(   # correlation matrix
  1,   0.8, 0,   0,   0,
  0.8, 1,   0,   0,   0,
  0,   0,   1,   0.3, 0.6,
  0,   0,   0.3, 1,   0.7,
  0,   0,   0.6, 0.7, 1
), nrow=p)
Sigma.star <- diag(sd)%*%C%*%diag(sd) # variance-covariance matrix of the explanatory variables
beta.star <- c(0.5, -0.3, 1, 0, -0.6) # coefficients
beta0.star <- -0.2  # intercept
beta.true = c(beta0.star,beta.star)

X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star) + matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
y <- as.numeric(runif(N)<p1)
data.complete <- data.frame(y=y,X.complete)
# ------- generating missing data
X.obs.MAR <- X.obs <- X.complete
p.miss =0.1
patterns = runif(N*p)<p.miss
X.obs[patterns] <- NA
for(i in c(2,4,5)){
  z <- cbind(y,X.complete[,c(1,3)])%*%matrix(sample(-5:5, 3, replace=T),ncol=1)        # linear combination
  pr <- 1/(1+exp(-z))         # pass through an inv-logit function
  r <- rbinom(N,1,pr)      # bernoulli response variable
  X.obs.MAR[r==0,i]<-NA
}
