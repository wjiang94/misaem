#' Mean
#'
#' This computes the mean without na
#' @param x a numerique vector
#' @return the mean
#' @import magrittr
#' @importFrom stats na.omit
#' @examples
#' mean_no_na(c(4,5))
#' @export
mean_no_na <- function(x){
  x <- x %>% na.omit()
  res <- sum(x)/length(x)
  return(res)
}
