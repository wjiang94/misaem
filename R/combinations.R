#' combinations
#'
#' Given all the possible patterns of missingness.
#' @param p Dimension of covariates.
#' @return A matrix containing all the possible missing patterns. Each row indicates a pattern of missingness. "1" means "observed", 0 means "missing".
#' @export



combinations = function(p){
  comb = NULL
  if (p<15) {
    for( i in 1:p) comb = rbind(cbind(1,comb),cbind(0,comb))
    return(comb)
  }
  else {stop("Error: this value will probably block your computer, try on your own risk")}

}
