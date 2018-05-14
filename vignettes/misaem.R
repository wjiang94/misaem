## ------------------------------------------------------------------------
library(misaem)

## ------------------------------------------------------------------------
# Generate dataset
N <- 1000  # number of subjects
p <- 5     # number of explanatory variables
mu.star <- 1:p  #rep(0,p)  # mean of the explanatory variables
sd <- 1:p # rep(1,p) # standard deviations
C <- matrix(c(   # correlation matrix
1,   0.8, 0,   0,   0,
0.8, 1,   0,   0,   0,
0,   0,   1,   0.3, 0.6,
0,   0,   0.3, 1,   0.7,
0,   0,   0.6, 0.7, 1), nrow=p)
Sigma.star <- diag(sd)%*%C%*%diag(sd) # variance-covariance matrix of the explanatory variables
beta.star <- c(0.5, -0.3, 1, 0, -0.6) # coefficients
beta0.star <- -0.2  # intercept
beta.true = c(beta0.star,beta.star)
X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star)
             + matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
y <- as.numeric(runif(N)<p1)

# Generate missingness
p.miss <- 0.10
patterns <- runif(N*p)<p.miss #missing completely at random
X.obs <- X.complete
X.obs[patterns] <- NA
# SAEM
list.saem = miss.saem(X.obs,y)
print(list.saem$beta)

