#' Compute optimal treatment assignment probabilities
#'
#' This function takes sample estimates and optimization constraints to generate a vector of treatment assignment probabilities
#'
#' @param tau Estimates subgroup treatment effects
#' @param sigma1_vec Group-level conditional variance under treatment arm
#' @param sigma0_vec Group-level conditional variance under control arm
#' @param p_vec Subgroup proportions
#' @param N_t Number of subjects enrolled up to Stage t
#' @param envy_c Envyfreeness parameter
#' @param lower_c Lower bound on treatment assignment probabilities
#' @param upper_c Upper bound on treatment assignment probabilities
#' @return Allocation vector
#' @export
#'
SubAlloc <- function(tau,sigma1_vec, sigma0_vec, p_vec, N_t,
                     envy_c,lower_c,upper_c){

  sigma11 <- sigma1_vec[1]
  sigma21 <- sigma1_vec[2]

  sigma10 <- sigma0_vec[1]
  sigma20 <- sigma0_vec[2]

  p1 <- p_vec[1]
  p2 <- p_vec[2]

  # objective function
  eval_f0 <- function(x){
    # res <- 1/2*(sigma11^2/x[1]+sigma10^2/(1-x[1])) +
    #   1/2*(sigma21^2/x[2]+sigma20^2/(1-x[2]))
    res <- p1*(sigma11^2/x[1]+sigma10^2/(1-x[1])) +
      p2*(sigma21^2/x[2]+sigma20^2/(1-x[2]))
    return( res )
  }

  delta_N <- sqrt(log(N_t)/N_t)

  # constraint function
  eval_g0 <- function(x) {
    constr <- c(x[1]-x[2]-envy_c,
                -envy_c-x[1]+x[2],
                -delta_N-log(x[1]/(1-x[1]))*tau[1],
                -delta_N-log(x[2]/(1-x[2]))*tau[2]
                #-log(x[1]/(1-x[1]))*tau[1],
                #-log(x[2]/(1-x[2]))*tau[2]
    )
    return( constr )
  }

  # Solve using NLOPT_LN_COBYLA without gradient information
  res1 <- nloptr( x0=c(0.01,0.01),
                  eval_f=eval_f0,
                  lb = c(lower_c,lower_c),
                  ub = c(upper_c,upper_c),
                  eval_g_ineq = eval_g0,
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                              "xtol_rel"=1.0e-8))

  # optimal e*
  e_star <- res1$solution[1:2]

  names(e_star) <- c("A","B")

  return(e_star)
}

