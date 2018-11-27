#' combinations
#'
#' Given all the possible patterns of missingness.
#' @param p Dimension of covariates.
#' @return A matrix containing all the possible missing patterns. Each row indicates a pattern of missingness. "1" means "observed", 0 means "missing".
#' @examples
#' comb = combinations(5)
#' @export



combinations = function(p){
  comb = NULL
  if (p<20) {
    for( i in 1:p) comb = rbind(cbind(1,comb),cbind(0,comb))
    return(comb)
  }
  else {stop("Error: the dimension of dataset is too large to possibly block your computer. Better try with number of variables smaller than 20.")}

}
