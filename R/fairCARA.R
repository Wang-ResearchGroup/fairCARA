#' Fair CARA experiments simulation
#'
#' This function simulate fair CARA experiments and return treatment effect estimates
#'
#' @param m Estimates subgroup treatment effects
#' @param n Group-level conditional variance under treatment arm
#' @param envy_c Envyfreeness parameter
#' @return
#' \item{tau_opt}{Estimated treatment effects under proposed design}
#' \item{tau_cr}{Estimated treatment effects under complete randomization}
#' \item{tau_db}{Estimated treatment effects under doubly-adaptive biased coin design}
#' \item{sub_alloc_opt}{Treatment allocation under proposed design}
#' \item{sub_alloc_cr}{Treatment allocation under complete randomization}
#' \item{sub_alloc_db}{Treatment allocation under doubly-adaptive biased coin design}
#' \item{ate_opt}{Estimated treatment effect under the proposed design}
#' \item{ate_cr}{Estimated treatment effect under complete randomization}
#' \item{ate_db}{Estimated treatment effect under doubly adaptive biased coin design}
#' \item{cover_opt}{Coverage under the proposed design}
#' \item{cover_cr}{Coverage under complete randomization}
#' \item{cover_db}{Coverage under doubly adaptive biased coin design}
#' \item{power_opt}{Power under the proposed design}
#' \item{power_cr}{Power under complete randomization}
#' \item{power_db}{Power under doubly adaptive biased coin design}
#' @export
#'

## Monte Carlo samples
fairCARA <- function(m,n,envy_c){

  # First stage
  X_1 <- GenX(m)
  S_1 <- GenS(X_1)

  # Stage 1: Assign treatment randomly
  T_1 <- rbinom(m,1,0.5)
  # Generate Y_1
  Y_1 <-  GenY(m,X_1,T_1,S_1)
  e_1 <- rep(1/2, 2) # randomly assign treatment with e=1/2
  p_1 <- colSums(S_1)/m

  # estimate subgroup ATE: use PSWeight
  tau_1 <- sd_1.t <- sd_1.c <- NULL

  for (j in 1:ncol(S_1)){
    dat1 <- as.data.frame(cbind(Y_1[S_1[,j]],T_1[S_1[,j]],X_1[S_1[,j]]))
    names(dat1) <- c("Y","Tr","X")

    tau_1[j] <-  sum(dat1$Y*dat1$Tr)/sum(dat1$Tr) - sum(dat1$Y*(1-dat1$Tr))/sum(1-dat1$Tr)
    sd_1.t[j]<- sd(dat1$Y[dat1$Tr==1])
    sd_1.c[j]<- sd(dat1$Y[dat1$Tr==0])
  }

  names(tau_1) <- c("A","B")
  names(sd_1.t) <- names(sd_1.c) <-  c("A","B")

  tau_old <- tau_1
  sd_old.t <- sd_old_db.t <- sd_1.t
  sd_old.c <- sd_old_db.c <- sd_1.c
  S_old <- S_old.db <- S_1
  n_old <- m
  T_old <- T_old.cr <- T_old.db <- T_1
  X_old <- X_1
  Y_old <- Y_old.cr <- Y_old.db <- Y_1

  p_old <- p_1

  ## fully sequential
  for(i in 1:n){ # start of fully sequential assignments

    X_i <- GenX(1)
    S_i <- GenS(X_i)
    S_new <- rbind(S_old,S_i)
    n_new <- n_old +1

    ## (1) Proposed design --------------------------------

    group_name <- colnames(S_i)[S_i]

    e_opt <- SubAlloc(tau_old,sd_old.t,sd_old.c, p_old,n_old,
                      envy_c,lower_c=0.01,upper_c=0.99)

    e_opt_i <- e_opt[group_name] # allocation for subject i

    # current allocation
    e_current_i <- sum(T_old[S_old[,group_name]])/length(T_old[S_old[,group_name]])

    T_i <- ifelse(e_current_i < e_opt_i,1,0)

    Y_i <-  GenY(1,X_i,T_i,S_i) # observe outcome of subject i

    # update old info
    T_old <- c(T_old, T_i)
    Y_old <- c(Y_old, Y_i)
    S_old <- S_new
    X_old <- c(X_old, X_i)
    n_old <- n_new

    # update tau and sd
    tau_old <- sd_old.t <-sd_old.c <- NULL
    for (j in 1:ncol(S_old)){
      dat1 <- as.data.frame(cbind(Y_old[S_old[,j]],T_old[S_old[,j]],X_old[S_old[,j]]))
      names(dat1) <- c("Y","Tr","X")
      ps <- sum(dat1$Tr)/nrow(dat1)

      tau_old[j] <-  sum(dat1$Y*dat1$Tr)/sum(dat1$Tr) - sum(dat1$Y*(1-dat1$Tr))/sum(1-dat1$Tr)
      sd_old.t[j]<- sd(dat1$Y[dat1$Tr==1])
      sd_old.c[j]<- sd(dat1$Y[dat1$Tr==0])
    }

    names(tau_old) <- c("A","B")
    names(sd_old.t) <- names(sd_old.c)<- c("A","B")


    sub_alloc_opt <- c(sum(T_old[S_old[,1]])/sum(S_old[,1]),
                       sum(T_old[S_old[,2]])/sum(S_old[,2]))

    ## (2) Complete randomization--------------------------------
    e_cr <- c(1/2,1/2)
    names(e_cr) <- c('A','B')
    T_i.cr <- rbinom(1,1,1/2)
    Y_i.cr <-  GenY(1,X_i,T_i.cr,S_i)
    # update old info
    T_old.cr <- c(T_old.cr, T_i.cr)
    Y_old.cr <- c(Y_old.cr, Y_i.cr)

    # update tau and sd
    tau_old.cr <- sd_old.cr.t <- sd_old.cr.c <- NULL
    for (j in 1:ncol(S_old)){
      dat1 <- as.data.frame(cbind(Y_old.cr[S_old[,j]],T_old.cr[S_old[,j]],X_old[S_old[,j]]))
      names(dat1) <- c("Y","Tr","X")

      tau_old.cr[j] <-  sum(dat1$Y*dat1$Tr)/sum(dat1$Tr) - sum(dat1$Y*(1-dat1$Tr))/sum(1-dat1$Tr)
      sd_old.cr.t[j]<- sd(dat1$Y[dat1$Tr==1])
      sd_old.cr.c[j]<- sd(dat1$Y[dat1$Tr==0])
    }

    names(tau_old.cr) <- c("A","B")
    names(sd_old.cr.t) <- names(sd_old.cr.c) <-  c("A","B")

    sub_alloc_cr <- c(sum(T_old.cr[S_old[,1]])/sum(S_old[,1]),
                      sum(T_old.cr[S_old[,2]])/sum(S_old[,2]))

    ## (3) DBCD --------------------------------------------
    e_db <- sd_old_db.t/(sd_old_db.t + sd_old_db.c)
    e_db_i <- e_db[group_name] # allocation for subject i

    # current allocation
    e_db_current_i <- sum(T_old.db[S_old.db[,group_name]])/length(T_old.db[S_old.db[,group_name]])

    T_i.db <- ifelse(e_db_current_i < e_db_i,1,0)

    Y_i.db <-  GenY(1,X_i,T_i.db,S_i) # observe outcome of subject i

    # update old info
    T_old.db <- c(T_old.db, T_i.db)
    Y_old.db <- c(Y_old.db, Y_i.db)
    S_old.db <- S_new

    # update tau and sd
    tau_old.db <- sd_old_db.t <-sd_old_db.c <- NULL
    for (j in 1:ncol(S_old)){
      dat1 <- as.data.frame(cbind(Y_old.db[S_old[,j]],T_old.db[S_old[,j]],X_old[S_old[,j]]))
      names(dat1) <- c("Y","Tr","X")

      tau_old.db[j] <-  sum(dat1$Y*dat1$Tr)/sum(dat1$Tr) - sum(dat1$Y*(1-dat1$Tr))/sum(1-dat1$Tr)
      sd_old_db.t[j]<- sd(dat1$Y[dat1$Tr==1])
      sd_old_db.c[j]<- sd(dat1$Y[dat1$Tr==0])
    }

    names(tau_old.db) <- c("A","B")
    names(sd_old_db.t) <- names(sd_old_db.c)<- c("A","B")

    sub_alloc_db <- c(sum(T_old.db[S_old[,1]])/sum(S_old[,1]),
                      sum(T_old.db[S_old[,2]])/sum(S_old[,2]))


  } # end of fully sequential assignments



  p_old <- colSums(S_old)/(nrow(S_old))

  ## Proposed design ------------------------------------
  tau_opt <- tau_old
  ate_opt <- tau_opt[1]*p_old[1] + tau_opt[2]*p_old[2]
  nu_opt <- 1/p_old*(sd_old.t^2/sub_alloc_opt + sd_old.c^2/(1-sub_alloc_opt))
  var_opt <-  sum(p_old^2*nu_opt) + (p_old[1]*(tau_opt[1]- ate_opt)^2 + p_old[2]*(tau_opt[2]- ate_opt)^2)
  sd_opt <- sqrt(var_opt)

  hi_opt <- ate_opt + 1.96*sd_opt/sqrt(m+n)
  lo_opt <- ate_opt - 1.96*sd_opt/sqrt(m+n)
  cover_opt <- (lo_opt < tau_true_all) & (hi_opt > tau_true_all)
  power_opt <- (lo_opt>0)


  ## Complete randomization ------------------------------------
  tau_cr <- tau_old.cr
  ate_cr <- tau_cr[1]*p_old[1] + tau_cr[2]*p_old[2]
  nu_cr <- 1/p_old*(sd_old.cr.t^2/sub_alloc_cr + sd_old.cr.c^2/(1-sub_alloc_cr))
  var_cr <-  sum(p_old^2*nu_cr) + (p_old[1]*(tau_cr[1]- ate_cr)^2 + p_old[2]*(tau_cr[2]- ate_cr)^2)
  sd_cr <- sqrt(var_cr)

  hi_cr <- ate_cr + 1.96*sd_cr/sqrt(m+n)
  lo_cr <- ate_cr - 1.96*sd_cr/sqrt(m+n)
  cover_cr <- (lo_cr < tau_true_all) & (hi_cr > tau_true_all)
  power_cr <- (lo_cr>0)

  ## DBCD ------------------------------------
  tau_db <- tau_old.db
  ate_db <- tau_db[1]*p_old[1] + tau_db[2]*p_old[2]
  nu_db <- 1/p_old*(sd_old_db.t^2/sub_alloc_db + sd_old_db.c^2/(1-sub_alloc_db))
  var_db <-  sum(p_old^2*nu_db) + (p_old[1]*(tau_db[1]- ate_db)^2 + p_old[2]*(tau_db[2]- ate_db)^2)
  sd_db <- sqrt(var_db)

  hi_db <- ate_db + 1.96*sd_db/sqrt(m+n)
  lo_db <- ate_db - 1.96*sd_db/sqrt(m+n)
  cover_db <- (lo_db < tau_true_all) & (hi_db > tau_true_all)
  power_db <- (lo_db>0)

  return(c(tau_opt,tau_cr,tau_db,
           sub_alloc_opt,sub_alloc_cr,sub_alloc_db,
           ate_opt,ate_cr,ate_db,
           cover_opt,cover_cr,cover_db,
           power_opt,power_cr,power_db))
}


