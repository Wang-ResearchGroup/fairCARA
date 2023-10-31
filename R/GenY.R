#' Generate outcomes
#'
#' This function takes sample size, covariates matrix, treatments, and subgroup indices as arguments and generates a vector of outcomes
#'
#' @param n Covariates matrix
#' @param X Covariates matrix
#' @param Tr Covariates matrix
#' @param S Covariates matrix
#' @return Outcome vector
#' @export
#'
#'
#'

GenY <- function(n,X,Tr,S){

  Y <- NULL

  for(i in 1:n){

    # treatment arm
    if(Tr[i]==1){
      if(S[i,1]){
        Y[i] <- rnorm(1,mu1.t,sd1.t)
      }else{
        Y[i] <- rnorm(1,mu2.t,sd2.t)
      }
    }


    # control arm
    if(Tr[i]==0){
      if(S[i,1]){
        Y[i] <- rnorm(1,mu1.c,sd1.c)
      }else{
        Y[i] <- rnorm(1,mu2.c,sd2.c)
      }
    }

  }

  return(Y)
}
