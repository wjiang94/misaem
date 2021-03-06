% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miss.glm.model.select.R
\name{miss.glm.model.select}
\alias{miss.glm.model.select}
\title{miss.glm.model.select}
\usage{
miss.glm.model.select(Y, X, seed = NA)
}
\arguments{
\item{Y}{Binary response vector \eqn{N \times 1}{N * 1}}

\item{X}{Design matrix with missingness \eqn{N \times p}{N * p}}

\item{seed}{An integer as a seed set for the random generator. The default value is 200.}
}
\value{
An object of class "\code{miss.glm}".
}
\description{
Model selection for the logistic regression model with missing data.
}
\examples{
# Generate dataset
N <- 40  # number of subjects
p <- 3     # number of explanatory variables
mu.star <- rep(0,p)  # mean of the explanatory variables
Sigma.star <- diag(rep(1,p)) # covariance
beta.star <- c(1, 1,  0) # coefficients
beta0.star <- 0 # intercept
beta.true = c(beta0.star,beta.star)
X.complete <- matrix(rnorm(N*p), nrow=N)\%*\%chol(Sigma.star) +
              matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
p1 <- 1/(1+exp(-X.complete\%*\%beta.star-beta0.star))
Y <- as.numeric(runif(N)<p1)
# Generate missingness
p.miss <- 0.10
patterns <- runif(N*p)<p.miss #missing completely at random
X <- X.complete
X[patterns] <- NA
# model selection for SAEM
miss.model = miss.glm.model.select(Y,X,seed=100)
print(miss.model)
}
