#' log_reg
#'
#' Caculate the likelihood or log-likelihood for one observation of logistic regression model .
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
