library(copula)
library(gumbel)

# Function to randomize the dataset
# Given a sorted (artificial) competing risks data.frame
# which contains the event time and status.
# Effectively, onyl the status 1 and 2 will be randomized.
#
rand.new <- function(dataset){
  n.eff <- sum(dataset$status %in% c(1,2))
  # randomly re-assign event types 1 and 2 
  dataset[which(dataset$status %in% c(1,2)),2] <- rbinom(n.eff,1,0.5)+1
  return(dataset)
}

# function to convert a paired censored dataset into a competing risks dataset
# n:     number of pairs
# times: event or censoring times
# cens:  censoring indicator: 1=uncensored, 0=censored
# tau:   Late (but not too late) final evaluation time.
#        Exceeding times are set to tau and uncensored.
#
paired2CR.new <- function(n, times1, cens1, times2, cens2, tau){
  times1 <- pmin(times1, tau)
  cens1[times1==tau] <- 1
  times2 <- pmin(times2, tau)
  cens2[times2==tau] <- 1
  times <- pmin(times1, times2)
  type <- ifelse((times1 < times2 & cens1) | (times1 == times2 & cens1 & !cens2), 1, 
                 ifelse((times2 < times1 & cens2) | (times1 == times2 & cens2 & !cens1), 2,
                        ifelse(times2 == times1 & cens1 & cens2, 3, 0)))
  return(data.frame(time=times, status=type)) 
}



# Actual function for computing the relative treatment effect estimator, 
# and the corresponding variance estimator.
#  dataf: a dataframe which is sorted according to the event times.
rte <- function(dataf){
  # unique sorted event times
  times <- unique(dataf$time[dataf$status>0])
  nt <- length(times)
  
  if(nt==0) return(c(0.5,999999999))
  
  atrisk <-  sapply(times, function(t) sum(dataf$time >= t))
  dN1 <- sapply(times, function(t) sum(dataf$status[dataf$time == t]==1))
  dN2 <- sapply(times, function(t) sum(dataf$status[dataf$time == t]==2))
  dN3 <- sapply(times, function(t) sum(dataf$status[dataf$time == t]==3))
  dN <- dN1+dN2+dN3
  
  # Kaplan-Meier estimator
  S <- cumprod(1-dN/atrisk)
  
  # Aalen-Johansen estimators
  # left-continuous version of Kaplan-Meier estimator in integral
  F1 <- cumsum(c(1,S[-nt]) * dN1/atrisk)
  F2 <- cumsum(c(1,S[-nt]) * dN2/atrisk)
  F3 <- cumsum(c(1,S[-nt]) * dN3/atrisk)
  
  num1 <- (atrisk - dN1)*dN1
  num2 <- (atrisk - dN2)*dN2
  num3 <- (atrisk - dN3)*dN3
  atriskMdN  <- atrisk - dN
  denom <- 1/(atrisk*atriskMdN^2)
  w1 <- num1*denom
  w2 <- num2*denom
  w3 <- num3*denom
  
  # Most likely, some truncation at tau, i.e., type3 event and NAE-jump of size 1 for type 3.
  # Also, assumed no atom at tau for the event times.
  if(dN[nt] == atrisk[nt]){
    var1 <- sum((F1[-nt]+S[-nt]-F1[nt-1])^2 * w1[-nt] + (F1[-nt]-F1[nt-1])^2 * (w2[-nt]+w3[-nt]))
    var2 <- sum((F2[-nt]+S[-nt]-F2[nt-1])^2 * w2[-nt] + (F2[-nt]-F2[nt-1])^2 * (w1[-nt]+w3[-nt]))
    cov12 <- sum((F2[nt-1]-F2[-nt])*(F1[nt-1]-F1[-nt]) * w3[-nt] - (S[-nt]+F2[-nt]-F2[nt-1])*(F1[nt-1]-F1[-nt])*w2[-nt] - (S[-nt]+F1[-nt]-F1[nt-1])*(F2[nt-1]-F2[-nt])*w1[-nt])
    return(c(F2[nt] + 0.5*F3[nt], 0.25*var1 + 0.25*var2 - 0.5*cov12))
  }else{
    var2 <- sum((F2+S-F2[nt])^2 * w2 + (F2-F2[nt])^2 * (w1+w3))
    var3 <- sum((F3+S-F3[nt])^2 * w3 + (F3-F3[nt])^2 * (w1+w2))
    cov23 <- sum((F2[nt]-F2)*(F3[nt]-F3) * w1 - (S+F2-F2[nt])*(F3[nt]-F3)*w2 - (S+F3-F3[nt])*(F2[nt]-F2)*w3)
    return(c(F2[nt] + 0.5*F3[nt], var2 + 0.25*var3 + cov23))
  }
}


# function to compute the relative treatment effect, variance estimate, 
#         p-values, confidence intervals
# n:          number of pairs
# times:      event or censoring times
# cens:       censoring indicator: 1=uncensored, 0=censored
# tau:        Late (but not too late) final evaluation time.
#             Exceeding times are set to tau and uncensored.
# alpha:      nominal significance level
# BS.iter:    number of bootstrap and randomization iterations.
# test.bool:  boolean. The tests will only be conducted when TRUE
# ci.bool:    boolean. The confidence intervals will only be computed when TRUE
# alt:        specifies the type of hypothesis test and confidence intervals.
#             Either of "two.sided" (default), "greater", or "less".
#             For more details, see the functions  one.test()  and  one.ci()  below.
#
rel_treat_eff.new <- function(n, times1, cens1, times2, cens2, tau, alpha=0.05,
                              rte_hyp=0.5, BS.iter=100, test.bool=FALSE, 
                              ci.bool=FALSE, alt="two.sided"){
  dataset <- paired2CR.new(n, times1, cens1, times2, cens2, tau)
  dataset <- dataset[order(dataset$time),]
  
  rte_cov <- rte(dataset)
  
  if(! alt %in% c("two.sided", "greater", "less", "g", "l")){
    print("Wrong specification of alt. Returning just the point and variance estimates.")
    return(rte_cov)
  }
  
  if(BS.iter==0) {
    return(rte_cov)
  }

  if((rte_cov[1]<0.000001) & (alt %in% c("g", "greater")))  
    return(list(RTE=rte_cov[1], RTE.var=0, 
       p.val.gauss=1, CI.gauss = c(0, 0), 
       p.val.rand=1, CI.rand = c(0, 0),
       p.val.bs=1, CI.bs = c(0, 0),
       p.val.phi.gauss=1, CI.phi.gauss = c(0, 0), 
       p.val.phi.rand=1, CI.phi.rand = c(0, 0),
       p.val.phi.bs=1, CI.phi.bs = c(0, 0)))
  
  if((rte_cov[1]<0.000001) & (alt %in% c("l", "less")))  
    return(list(RTE=rte_cov[1], RTE.var=0, 
      p.val.gauss=0, CI.gauss = c(0, 0), 
      p.val.rand=0, CI.rand = c(0, 0),
      p.val.bs=0, CI.bs = c(0, 0),
      p.val.phi.gauss=0, CI.phi.gauss = c(0, 0), 
      p.val.phi.rand=0, CI.phi.rand = c(0, 0),
      p.val.phi.bs=0, CI.phi.bs = c(0, 0)))
  
  if((rte_cov[1]<0.000001) & (alt == "two.sided"))  
    return(list(RTE=rte_cov[1], RTE.var=0, 
                p.val.gauss=0, CI.gauss = c(0, 0), 
                p.val.rand=0, CI.rand = c(0, 0),
                p.val.bs=0, CI.bs = c(0, 0),
                p.val.phi.gauss=0, CI.phi.gauss = c(0, 0), 
                p.val.phi.rand=0, CI.phi.rand = c(0, 0),
                p.val.phi.bs=0, CI.phi.bs = c(0, 0)))
  
  return(list(RTE=rte_cov[1], RTE.var=0, 
              p.val.gauss=0, CI.gauss = c(1, 1), 
              p.val.rand=0, CI.rand = c(1, 1),
              p.val.bs=0, CI.bs = c(1, 1),
              p.val.phi.gauss=0, CI.phi.gauss = c(1, 1), 
              p.val.phi.rand=0, CI.phi.rand = c(1, 1),
              p.val.phi.bs=0, CI.phi.bs = c(1, 1)))
  
  if((abs(rte_cov[1]-1)<0.000001) & (alt %in% c("g", "greater")))  
    return(list(RTE=rte_cov[1], RTE.var=0, 
      p.val.gauss=0, CI.gauss = c(1, 1), 
      p.val.rand=0, CI.rand = c(1, 1),
      p.val.bs=0, CI.bs = c(1, 1),
      p.val.phi.gauss=0, CI.phi.gauss = c(1, 1), 
      p.val.phi.rand=0, CI.phi.rand = c(1, 1),
      p.val.phi.bs=0, CI.phi.bs = c(1, 1)))
  
  if((abs(rte_cov[1]-1)<0.000001) & (alt == "two.sided"))  
    return(list(RTE=rte_cov[1], RTE.var=0, 
                p.val.gauss=0, CI.gauss = c(1, 1), 
                p.val.rand=0, CI.rand = c(1, 1),
                p.val.bs=0, CI.bs = c(1, 1),
                p.val.phi.gauss=0, CI.phi.gauss = c(1, 1), 
                p.val.phi.rand=0, CI.phi.rand = c(1, 1),
                p.val.phi.bs=0, CI.phi.bs = c(1, 1)))
  
  
  # randomization
  rand_p <- replicate(BS.iter, rand_rel_treat_eff.new(dataset))

  # bootstrap
  bs_p <- replicate(BS.iter, bs_rel_treat_eff.new(dataset, n, rte_cov[1]))
  
  p_val_gauss <- 1
  CI_gauss_lower <- 0
  CI_gauss_upper <- 1
  p_val_rand <- 1
  CI_rand_lower <- 0
  CI_rand_upper <- 1
  p_val_bs <- 1
  CI_bs_lower <- 0
  CI_bs_upper <- 1
  p_val_phi_gauss <- 1
  CI_phi_gauss_lower <- 0
  CI_phi_gauss_upper <- 1
  p_val_phi_rand <- 1
  CI_phi_rand_lower <- 0
  CI_phi_rand_upper <- 1
  p_val_phi_bs <- 1
  CI_phi_bs_lower <- 0
  CI_phi_bs_upper <- 1
  
  
  if(test.bool){
    T_0 <- (rte_cov[1]-rte_hyp)/sqrt(rte_cov[2])
    
    # the studentization with 0.5*log(0.5) is probably much better than p*log(p)
    # otherwise, for small or large values of p the test statistic could tend to have small values.
    T_phi <- (log(-log(rte_cov[1]))-log(-log(rte_hyp)))/sqrt(rte_cov[2])*rte_cov[1]*log(rte_cov[1])
    
    if(alt=="greater"  | alt=="g"){
      # 1-sided p-values (right-tailed)
      p_val_gauss <- 1 - pnorm(T_0) 
      p_val_bs <- mean(bs_p[1,] >= T_0)
      p_val_rand <- mean(rand_p[1,] >= T_0)
      
      # right because trafo phi is decreasing, but we also standardize with something negative.
      p_val_phi_gauss <- 1-pnorm(T_phi)
      p_val_phi_bs <- mean(bs_p[2,] >= T_phi)
      p_val_phi_rand <- mean(rand_p[2,] >= T_phi)
    }
    
    if(alt=="less"  | alt=="l"){
      # 1-sided p-value (left-tailed)
      p_val_gauss <- pnorm(T_0) 
      p_val_bs <- mean(bs_p[1,] <= T_0)
      p_val_rand <- mean(rand_p[1,] <= T_0)
      
      p_val_phi_gauss <- pnorm(T_phi)
      p_val_phi_bs <- mean(bs_p[2,] <= T_phi)
      p_val_phi_rand <- mean(rand_p[2,] <= T_phi)
    }
    
    if(alt=="two.sided"){
      # 2-sided p-values (two-tailed)
      p_val_gauss <- 2*(1-pnorm(abs(T_0)))
      p_val_bs <- 2*min(mean(bs_p[1,] <= T_0), mean(bs_p[1,] >= T_0))
      p_val_rand <- 2*min(mean(rand_p[1,] <= T_0), mean(rand_p[1,] >= T_0))
      
      p_val_phi_gauss <- 2*(1-pnorm(abs(T_phi)))
      p_val_phi_rand <- 2*min(mean(rand_p[2,] <= T_phi), mean(rand_p[2,] >= T_phi))
      p_val_phi_bs <- 2*min(mean(bs_p[2,] <= T_phi), mean(bs_p[2,] >= T_phi))
    }
  }
  
  if(ci.bool){
    
    if(alt=="greater"  | alt=="g"){
      # Confidence intervals ending at 1 to the right
      CI_gauss_lower <- rte_cov[1] - qnorm(1-alpha) * sqrt(rte_cov[2])
      CI_rand_lower <- rte_cov[1] - quantile(rand_p[1,], prob=1-alpha)*sqrt(rte_cov[2])
      CI_bs_lower <- rte_cov[1] - quantile(bs_p[1,], prob=1-alpha)*sqrt(rte_cov[2])
      
      CI_phi_gauss_lower <- rte_cov[1]^exp(-qnorm(1-alpha)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
      CI_phi_rand_lower <- rte_cov[1]^exp(-quantile(rand_p[2,], prob=1-alpha)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
      CI_phi_bs_lower <- rte_cov[1]^exp(-quantile(bs_p[2,], prob=1-alpha)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
    }
    
    if(alt=="less"  | alt=="l"){
      # Confidence intervals ending at 0 to the left
      CI_gauss_upper <- rte_cov[1] + qnorm(1-alpha) * sqrt(rte_cov[2])
      CI_rand_upper <- rte_cov[1] - quantile(rand_p[1,], prob=alpha)*sqrt(rte_cov[2])
      CI_bs_upper <- rte_cov[1] - quantile(bs_p[1,], prob=alpha)*sqrt(rte_cov[2])
      
      CI_phi_gauss_upper <- rte_cov[1]^exp(qnorm(1-alpha/2)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
      CI_phi_rand_upper <- rte_cov[1]^exp(-quantile(rand_p[2,], prob=alpha/2)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
      CI_phi_bs_upper <- rte_cov[1]^exp(-quantile(bs_p[2,], prob=alpha/2)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
    }
    
    if(alt=="two.sided"){
      # Confidence intervals with no restriction on any boundary
      CI_gauss_lower <- rte_cov[1] - qnorm(1-alpha/2) * sqrt(rte_cov[2])
      CI_gauss_upper <- rte_cov[1] + qnorm(1-alpha/2) * sqrt(rte_cov[2])
      CI_rand_lower <- rte_cov[1] - quantile(rand_p[1,], prob=1-alpha/2)*sqrt(rte_cov[2])
      CI_rand_upper <- rte_cov[1] - quantile(rand_p[1,], prob=alpha/2)*sqrt(rte_cov[2])
      CI_bs_lower <- rte_cov[1] - quantile(bs_p[1,], prob=1-alpha/2)*sqrt(rte_cov[2])
      CI_bs_upper <- rte_cov[1] - quantile(bs_p[1,], prob=alpha/2)*sqrt(rte_cov[2])
      
      CI_phi_gauss_lower <- rte_cov[1]^exp(-qnorm(1-alpha/2)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
      CI_phi_gauss_upper <- rte_cov[1]^exp(qnorm(1-alpha/2)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
      CI_phi_rand_lower <- rte_cov[1]^exp(-quantile(rand_p[2,], prob=1-alpha/2)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
      CI_phi_rand_upper <- rte_cov[1]^exp(-quantile(rand_p[2,], prob=alpha/2)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
      CI_phi_bs_lower <- rte_cov[1]^exp(-quantile(bs_p[2,], prob=1-alpha/2)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
      CI_phi_bs_upper <- rte_cov[1]^exp(-quantile(bs_p[2,], prob=alpha/2)*sqrt(rte_cov[2])/rte_cov[1]/log(rte_cov[1]))
    }
    
    
  }
  
  return(list(RTE=rte_cov[1], RTE.var=rte_cov[2], 
              p.val.gauss=p_val_gauss, CI.gauss = c(CI_gauss_lower, CI_gauss_upper), 
              p.val.rand=p_val_rand, CI.rand = c(CI_rand_lower, CI_rand_upper),
              p.val.bs=p_val_bs, CI.bs = c(CI_bs_lower, CI_bs_upper),
              p.val.phi.gauss=p_val_phi_gauss, CI.phi.gauss = c(CI_phi_gauss_lower, CI_phi_gauss_upper), 
              p.val.phi.rand=p_val_phi_rand, CI.phi.rand = c(CI_phi_rand_lower, CI_phi_rand_upper),
              p.val.phi.bs=p_val_phi_bs, CI.phi.bs = c(CI_phi_bs_lower, CI_phi_bs_upper)))
}


# function to randomize the dataset and then compute the normalized relative treatment effect
rand_rel_treat_eff.new <- function(dataset){
  rte_cov <- rte(rand.new(dataset))
  
  # If relative treatment effect is very close to 0 or 1.
  if(is.na((log(-log(rte_cov[1]))-log(-log(0.5)))/sqrt(rte_cov[2])*rte_cov[1]*log(rte_cov[1]))){
    if(abs(rte_cov[2])<0.000001){
      if(abs(rte_cov[1]-1) < 0.000001) return(c(Inf, Inf))
      if(abs(rte_cov[1]) < 0.000001) return(c(-Inf, -Inf))
      if(abs(rte_cov[1]-0.5) < 0.000001) return(c(0,0))
    }else{
      if(abs(rte_cov[1]-1) < 0.000001) return(c((rte_cov[1]-0.5)/sqrt(rte_cov[2]), Inf))
      if(abs(rte_cov[1]) < 0.000001) return(c((rte_cov[1]-0.5)/sqrt(rte_cov[2]), -Inf))
    }
  }
  # Return the untransformed and log-log-transformed standardized randomized relative treatment effect.
  return(c((rte_cov[1]-0.5)/sqrt(rte_cov[2]), (log(-log(rte_cov[1]))-log(-log(0.5)))/sqrt(rte_cov[2])*rte_cov[1]*log(rte_cov[1])))
}


# function to bootstrap the dataset and then compute the normalized relative treatment effect
bs_rel_treat_eff.new <- function(dataset, n, p){
  dataset.bs <- dataset[sample(1:n,n,replace=TRUE),]
  dataset.bs <- dataset.bs[order(dataset.bs$time),]
  rte_cov <- rte(dataset.bs)
  
  # If relative treatment effect is very close to 0 or 1.
  if(is.na((log(-log(rte_cov[1]))-log(-log(p)))/sqrt(rte_cov[2])*rte_cov[1]*log(rte_cov[1]))){
    if(abs(rte_cov[2])<0.000001){
      if(abs(rte_cov[1]-1) < 0.000001) return(c(Inf, Inf))
      if(abs(rte_cov[1]) < 0.000001) return(c(-Inf, -Inf))
      if(abs(rte_cov[1]-p) < 0.000001) return(c(0,0))
    }else{
      if(abs(rte_cov[1]-1) < 0.000001) return(c((rte_cov[1]-p)/sqrt(rte_cov[2]), Inf))
      if(abs(rte_cov[1]) < 0.000001) return(c((rte_cov[1]-p)/sqrt(rte_cov[2]), -Inf))
    }
  } 
  # Return the untransformed and log-log-transformed standardized bootstrapped relative treatment effect.
  return(c((rte_cov[1]-p)/sqrt(rte_cov[2]), (log(-log(rte_cov[1]))-log(-log(p)))/sqrt(rte_cov[2])*rte_cov[1]*log(rte_cov[1])))
}



dataGen <- function(n, setting){
  # setting: a vector which controls how the data is generated
  # Entries:
  # 1: copula: Gumbel-Hougaard (1), Clayton (2), indep. (3)
  # 2: parameter for copula
  # 3: S2: Exp(2) (1), mixture vs Exp (2), Go(0.1,1.4...) vs Exp(1) (3)
  # 4: censoring interval max (using U[0,max]-distribution)

  Z <- switch(setting[1],
              rgumbel(n, setting[2], dim = 2, method=1),
              rCopula(n, copula=claytonCopula(setting[2])),
              matrix(runif(2*n),nc=2)
  )
  Z <- rbind(Z, cbind(runif(0), rep(0,0)), cbind(rep(0,0), runif(0)))
  
  rootGompertz <- c(3.02949, 3.1625, 1)
  qGompertzRoot <- function(z, shape=0.6, scale=rootGompertz[setting[1]]) return(log(1-1/shape*log(1-z))/scale)  
  rootMix <- c(1.3255, 1.298, 1)
  qmix <- function(z, rate1=3, rate2=rootMix[setting[1]], p=0.5)
    return(sapply(z,function(y) uniroot(function(x){y - p*pexp(x,rate1) - (1-p)*pexp(x,rate2)}, interval=c(0,100))$root))
  
  T1 <- switch(setting[3],
               qexp(Z[,1], rate=2),
               qexp(Z[,1], rate=2),
               qGompertzRoot(Z[,1]))
  
  T2 <- switch(setting[3],
               qexp(Z[,2], rate=2),
               qmix(Z[,2]),
               qexp(Z[,2], rate=3)
  )
  C1 <- runif(n, 0, setting[4])
  C2 <- runif(n, 0, setting[4])
  
  T1[T1==0] <- 0.000000000000000001
  T2[T2==0] <- 0.000000000000000001
  C1[C1==0] <- 0.000000000000000001
  C2[C2==0] <- 0.000000000000000001
  
  return(cbind(T1,T2,C1,C2))
}



# different simulation settings:
settings <- matrix(c(
  c(1,5,1,1.1),
  c(1,5,1,1.6),
  c(1,5,1,2.7),
  c(1,5,2,1.1),
  c(1,5,2,1.6),
  c(1,5,2,2.7),
  c(2,-0.6,1,1.1),
  c(2,-0.6,1,1.6),
  c(2,-0.6,1,2.7),
  c(2,-0.6,2,1.1),
  c(2,-0.6,2,1.6),
  c(2,-0.6,2,2.7),
  c(3,0,1,1.1),
  c(3,0,1,1.6),
  c(3,0,1,2.7),
  c(3,0,2,1.1),
  c(3,0,2,1.6),
  c(3,0,2,2.7),
  # 19 till 27
  # for Gompertz vs Exp marginals
  c(1,5,3,0.7),
  c(1,5,3,1),
  c(1,5,3,1.75),
  c(2,-0.6,3,0.7),
  c(2,-0.6,3,1),
  c(2,-0.6,3,1.75),
  c(3,0,3,0.7),
  c(3,0,3,1),
  c(3,0,3,1.75)), 
  byrow=T, nc=4)


  
  # Conducts the tests. Returns the value of the relative treatment effect estimator,
  # and the (one-sided) p-values for the asymptotic, randomization, bootstrap test
  # untransformed and log-log-transformed.
  # itnr: iteration number (in simulations)
  # B: number of bootstrap and randomization iterations
  # alpha: significance level
  # n: sample size
  # setting: a vector of size 4 which specifies the data generation.
  #           Possible choices are given in the vector settings.
  # tau: end of study time
  # X1: (censored) times under treatment 1; if NA, data will be generated.
  # X2: (censored) times under treatment 2.
  # d1: uncensoring indicators for treatment 1.
  # d2: uncensoring indicators for treatment 2.
  # alt: the alternative hypothesis. Either of 
  #     "two.sided" (default), i.e.,  H_a: relative treatment effect != 0.5,
  #     "greater", i.e.,              H_a: relative treatment effect > 0.5,
  #     "less", i.e.,                 H_a: relative treatment effect < 0.5.
  #
one.test <- function(itnr, B, alpha=0.05, n, setting=c(1,1,1,1), tau, X1=NULL, X2=NULL, d1=NULL, d2=NULL, alt="two.sided"){
  
  # data generation
  if(is.null(X1)){
    data <- dataGen(n, setting)
    T1 <- data[,1]
    C1 <- data[,3]
    X1 <- pmin(T1,C1,rep(tau,n))
    d1 <- ((X1==T1) | (X1==tau)) 
    
    T2 <- data[,2]
    C2 <- data[,4]
    X2 <- pmin(T2,C2,rep(tau,n))
    d2 <- ((X2==T2) | (X2==tau)) 
  }else{
    X1 <- pmin(X1,tau)
    d1 <- (d1 | (X1==tau))
    X2 <- pmin(X2,tau)
    d2 <- (d2 | (X2==tau)) 
  }
  
  results <- rel_treat_eff.new(n, times1=X1, cens1=d1, times2=X2, cens2=d2, tau=tau, alpha=alpha, BS.iter=B, test.bool=TRUE, ci.bool=FALSE, alt=alt)
  
  return(c(results[[1]], results[[3]], results[[5]], results[[7]], results[[9]], results[[11]], results[[13]]))
}


# Computes the confidence intervals. Returns the value of the relative treatment effect estimator,
# and the (one-sided) p-values for the asymptotic, randomization, bootstrap-based
# (1-alpha)x100% confidence intervals:  untransformed and log-log-transformed.
# itnr: iteration number (in simulations)
# B: number of bootstrap and randomization iterations
# alpha: significance level
# n: sample size
# setting: a vector of size 4 which specifies the data generation.
#           Possible choices are given in the vector settings.
# tau: end of study time
# X1: (censored) times under treatment 1; if NA, data will be generated.
# X2: (censored) times under treatment 2.
# d1: uncensoring indicators for treatment 1.
# d2: uncensoring indicators for treatment 2.
# factor1: a factor multiplied to the event times of treatment group 1 to create
#   settings in which the relative treatment effect differs from 0.5.
# alt: the type of confidence interval. Either of 
#     "two.sided" (default), i.e.,  CI = [L,U],
#     "greater", i.e.,              CI = [L,1],
#     "less", i.e.,                 CI = [0,U].
#
one.CI <- function(itnr, B, alpha=0.05, n, setting=c(1,1,1,1), tau, X1=NULL, X2=NULL, d1=NULL, d2=NULL, alt="two.sided", factor1=1){

    # data generation
  if(is.null(X1)){
    data <- dataGen(n, setting)
    T1 <- data[,1] * factor1
    C1 <- data[,3]
    X1 <- pmin(T1,C1,rep(tau,n))
    d1 <- ((X1==T1) | (X1==tau)) 
    
    T2 <- data[,2]
    C2 <- data[,4]
    X2 <- pmin(T2,C2,rep(tau,n))
    d2 <- ((X2==T2) | (X2==tau)) 
  }else{
    X1 <- pmin(X1,tau)
    d1 <- (d1 | (X1==tau)) 
    X2 <- pmin(X2,tau)
    d2 <- (d2 | (X2==tau)) 
  }
  
  results <- rel_treat_eff.new(n, times1=X1, cens1=d1, times2=X2, cens2=d2, tau=tau, alpha=alpha, BS.iter=B, plotting=FALSE, test.bool=FALSE, ci.bool=TRUE, alt=alt)
  
  return(c(results[[4]], results[[6]], results[[8]], results[[10]], results[[12]], results[[14]], results[[3]], results[[5]], results[[7]], results[[9]], results[[11]], results[[13]], results[[1]]))
}


### For (parallel) simulation of the data sets and conducting the tests.
if(FALSE){
  library(parallel)
  
  # under Linux/Unix:
  cl <- makeCluster(16, type="FORK", outfile="errors.txt")
  # under Windows:
  #cl <- makeCluster(16, type="PSOCK", outfile="errors.txt")
  clusterExport(cl=cl, c('dataGen', 'one.test', 'one.ci', 'rte', 'rgumbel', 'rCopula', 'claytonCopula', 'rel_treat_eff.new', 'paired2CR.new', 'rand_rel_treat_eff.new', 'rand.new', 'bs_rel_treat_eff.new', 'rCopula', 'claytonCopula'))
  # index for the specification of simulation setting
  index <- 4
  iseed <- round(123456+index*100)
  clusterSetRNGStream(cl, iseed = iseed)
  n <- c(25,50,75,100,125,150)
  tau <- 1
  B <- 2000
  iter <- 5000
  alpha <- 0.05
  
  for(n.ind in 1:length(n)){
    
    results <- parSapply(cl,1:iter, one.test, B,alpha,n=n[n.ind] ,settings[index,],tau)
    clusterExport(cl=cl, 'results')
    
    L.res <- list(res=results, means=NA, n=n[n.ind], setting=settings[index,], setting.ind <- index, iseed=iseed, tau=tau)
    save(L.res, file=paste("~/", gsub("xx", index, "RTE_setting_xx_n"), n[n.ind],".RData",sep=""))
    
    # results matrix will be arranged in the following form:
    # columns: no trafo / trafo
    # rows: gauss / rand / bs
    # repeated three times (alpha=.01,.05,.1) and combined by extending columns
    
    results2 <- matrix(rowMeans(results[-1,] <= alpha/5), nc=2, byrow=F)
    results2 <- cbind(results2, matrix(rowMeans(results[-1,] <= alpha), nc=2, byrow=F))
    results2 <- cbind(results2, matrix(rowMeans(results[-1,] <= alpha*2), nc=2, byrow=F))
    print(n.ind)
    rownames(results2) <- c("gauss", "rand", "bs")
    colnames(results2) <- c("lin 1%", "trafo 1%", "lin 5%", "trafo 5%", "lin 10%", "trafo 10%")
    print(results2)
    L.res <- list(res=results, means=results2, RTE_mean = mean(results[1,]), n=n[n.ind], setting=settings[index,], setting.ind <- index, iseed=iseed, tau=tau)
    save(L.res, file=paste("~/", gsub("xx", index, "RTE_setting_xx_n"), n[n.ind],".RData",sep=""))
  }
  
  stopCluster(cl)
}



### APPLICATION to diabetic retinopathy data set 

if(FALSE){
  library(survival)
  data(diabetic)
  
  diabetes <- diabetic
  
  # Y for young
  # A for adult
  
  diabetesY <- subset(diabetes, age < 20)
  diabetesA <- subset(diabetes, age >= 20)
  
  # treated eyes
  treat.A <- subset(diabetesA,trt==1)
  untreat.A <- subset(diabetesA,trt==0)
  
  TimeA1 <- treat.A$time
  StatusA1 <- treat.A$status
  TimeA2 <- untreat.A$time
  StatusA2 <- untreat.A$status
  
  # untreated eyes
  treat.Y <- subset(diabetesY,trt==1)
  untreat.Y <- subset(diabetesY,trt==0)
  
  TimeY1 <- treat.Y$time
  StatusY1 <- treat.Y$status
  TimeY2 <- untreat.Y$time
  StatusY2 <- untreat.Y$status
  
  # Results
  # two-sided analysis
  set.seed(1234)
  res_A <- rel_treat_eff.new(length(TimeA1), TimeA1, StatusA1, TimeA2, StatusA2, tau=60, BS.iter=2000, test.bool=T, ci.bool=T, alt="two.sided")
  res_A
  set.seed(1234)
  res_Y <- rel_treat_eff.new(length(TimeY1), TimeY1, StatusY1, TimeY2, StatusY2, tau=60, BS.iter=2000, test.bool=T, ci.bool=T, alt="two.sided")
  res_Y
  
  # one-sided analysis (right-tailed / upper confidence interval)
  set.seed(1234)
  res_A <- rel_treat_eff.new(length(TimeA1), TimeA1, StatusA1, TimeA2, StatusA2, tau=60, BS.iter=2000, test.bool=T, ci.bool=T, alt="greater")
  res_A
  set.seed(1234)
  res_Y <- rel_treat_eff.new(length(TimeY1), TimeY1, StatusY1, TimeY2, StatusY2, tau=60, BS.iter=2000, test.bool=T, ci.bool=T, alt="greater")
  res_Y
  
  # one-sided analysis (left-tailed / lower confidence interval)
  set.seed(1234)
  res_A <- rel_treat_eff.new(length(TimeA1), TimeA1, StatusA1, TimeA2, StatusA2, tau=60, BS.iter=2000, test.bool=T, ci.bool=T, alt="less")
  res_A
  set.seed(1234)
  res_Y <- rel_treat_eff.new(length(TimeY1), TimeY1, StatusY1, TimeY2, StatusY2, tau=60, BS.iter=2000, test.bool=T, ci.bool=T, alt="less")
  res_Y
}
