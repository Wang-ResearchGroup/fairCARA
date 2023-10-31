#' Generate subgroup indices
#'
#' This function takes covariates mateix as an argument and generates a matrix of subgroup indices
#'
#' @param X Covariates matrix
#' @return Subgroup membership matrix
#' @export
#'
#'
#'

GenS <- function(X){

  S1 <- (X ==1)
  S2 <- (X==2)

  S <- cbind(S1,S2)

  colnames(S) <- c("A","B")

  return(S)
}


