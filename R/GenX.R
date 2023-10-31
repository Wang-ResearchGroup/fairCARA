#' Generate covariates
#'
#' This function takes sample size argument and generates a matrix of covariates
#'
#' @param n Sample size
#' @return Covariates matrix
#' @export
#'
#'
#'
GenX <- function(n){
  X <- sample(1:2, size = n, replace = TRUE, prob = true_p)
  return(X)
}
