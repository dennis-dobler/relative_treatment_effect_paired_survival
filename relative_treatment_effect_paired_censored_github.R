# The following two lines must be executed only once.
# This is necessary for using the Greenwood-type covariance estimators of the Aalen-Johansen estimators which are implemented in the github version.
#library(devtools)
#install_github("mclements/etm")

library(copula)
library(gumbel)
library(etm)


# function to randomize the dataset
# n:     number of pairs
# times: event or censoring times
# cens:  censoring indicator: 1=uncensored, 0=censored
# tau:   Late (but not too late) final evaluation time.
#        Exceeding times are set to tau and uncensored.
#
rand <- function(n, times1, cens1, times2, cens2, tau){
  times <- c(times1, times2)
  cens <- c(cens1, cens2)
  G <- rbinom(n, 1, 0.5)
  return(list(times1=times[(1:n)+n*G], cens1=cens[(1:n)+n*G], times2=times[(1:n)+n*(1-G)], cens2=cens[(1:n)+n*(1-G)]))
}

# function to convert a paired censored dataset into a competing risks dataset
# n:     number of pairs
# times: event or censoring times
# cens:  censoring indicator: 1=uncensored, 0=censored
# tau:   Late (but not too late) final evaluation time.
#        Exceeding times are set to tau and uncensored.
#
paired2CR <- function(n, times1, cens1, times2, cens2, tau){
  times1 <- pmin(times1, tau)
  cens1[times1==tau] <- 1
  times2 <- pmin(times2, tau)
  cens2[times2==tau] <- 1
  times <- pmin(times1, times2)
  type <- ifelse((times1 < times2 & cens1) | (times1 == times2 & cens1 & !cens2), 1, 
          ifelse((times2 < times1 & cens2) | (times1 == times2 & cens2 & !cens1), 2,
          ifelse(times2 == times1 & cens1 & cens2, 3, "cens")))
  return(data.frame(id = 1:n, from=0, to=type, time=times)) 
}


# function to compute the relative treatment effect, variance estimate, 
#     p-values, confidence intervals
# n:     number of pairs
# times: event or censoring times
# cens:  censoring indicator: 1=uncensored, 0=censored
# tau:   Late (but not too late) final evaluation time.
#        Exceeding times are set to tau and uncensored.
# alpha: nominal significance level
# BS.iter: number of bootstrap and randomization iterations.
#
rel_treat_eff <- function(n, times1, cens1, times2, cens2, tau, alpha=0.05, BS.iter=100){
  dataset <- paired2CR(n, times1, cens1, times2, cens2, tau)
  tra <- matrix(ncol=4,nrow=4,FALSE)
  tra[1,2:4] <- TRUE
  
  etm_complete <- etm(data=dataset, state.names=0:3, tra=tra, cens.name="cens", s=0, t=tau, covariance=TRUE)
  
  if(is.null(etm_complete$cov)) return(list(RTE=0.5, RTE.var=0, 
                                            p.val.gauss=1, CI.gauss = c(0, 0), 
                                            p.val.rand=1, CI.rand = c(0, 0),
                                            p.val.bs=1, CI.bs = c(0, 0),
                                            p.val.phi.gauss=1, CI.phi.gauss = c(0, 0), 
                                            p.val.phi.rand=1, CI.phi.rand = c(0, 0),
                                            p.val.phi.bs=1, CI.phi.bs = c(0, 0)))
  
  latest <- length(etm_complete$est[1,1,])
  
  if(length(etm_complete$cov[1,1,latest-1])==0) return(list(RTE=0.5, RTE.var=0, 
                                           p.val.gauss=1, CI.gauss = c(0, 0), 
                                           p.val.rand=1, CI.rand = c(0, 0),
                                           p.val.bs=1, CI.bs = c(0, 0),
                                           p.val.phi.gauss=1, CI.phi.gauss = c(0, 0), 
                                           p.val.phi.rand=1, CI.phi.rand = c(0, 0),
                                           p.val.phi.bs=1, CI.phi.bs = c(0, 0)))
  
  p <- etm_complete$est[1,3,latest] + 0.5 * etm_complete$est[1,4,latest]
  
  # if p=0, don't reject; if p=1, reject one-sided test
  if(p<0.000001)  return(list(RTE=p, RTE.var=0, 
                p.val.gauss=1, CI.gauss = c(0, 0), 
                p.val.rand=1, CI.rand = c(0, 0),
                p.val.bs=1, CI.bs = c(0, 0),
                p.val.phi.gauss=1, CI.phi.gauss = c(0, 0), 
                p.val.phi.rand=1, CI.phi.rand = c(0, 0),
                p.val.phi.bs=1, CI.phi.bs = c(0, 0)))
  
  if(abs(p-1)<0.000001)  return(list(RTE=p, RTE.var=0, 
                p.val.gauss=0, CI.gauss = c(1, 1), 
                p.val.rand=0, CI.rand = c(1, 1),
                p.val.bs=0, CI.bs = c(1, 1),
                p.val.phi.gauss=0, CI.phi.gauss = c(1, 1), 
                p.val.phi.rand=0, CI.phi.rand = c(1, 1),
                p.val.phi.bs=0, CI.phi.bs = c(1, 1)))
  
  
  
  # Indices in etm function:
  # 1 for survival, 9 for 0->2 trans., 13 for 0->3 trans.
  sigma2 <- 0.25*etm_complete$cov[1,1,latest-1] + etm_complete$cov[9,9,latest-1] + 0.25*etm_complete$cov[13,13,latest-1] + etm_complete$cov[1,9,latest-1] + etm_complete$cov[9,13,latest-1] + 0.25*etm_complete$cov[1,13,latest-1]
  # Here we used the time index latest-1 because the last observation(s) is/are uncensored at tau.
  # Hence, a jump with size S(tau-) occurs, where S is the Kaplan-Meier estimator.
  
  
  # We wish to test the right-sided alternative p > 0.5.
  # Thus reject for large values of test statistic T_0.
  
  # First, for the untransformed:
  
  T_0 <- (p-0.5)/sqrt(sigma2)
  
  if(is.na(T_0))  return(list(RTE=p, RTE.var=sigma2, 
                p.val.gauss=1, CI.gauss = c(0, 1), 
                p.val.rand=1, CI.rand = c(0, 1),
                p.val.bs=1, CI.bs = c(0, 1),
                p.val.phi.gauss=1, CI.phi.gauss = c(0, 1), 
                p.val.phi.rand=1, CI.phi.rand = c(0, 1),
                p.val.phi.bs=1, CI.phi.bs = c(0, 1)))

  # two-sided p-value
  # p_val_gauss <- 2*pnorm(abs(T_0)) - 1
  # 1-sided p-value (right-tailed)
  p_val_gauss <- 1 - pnorm(T_0) 
  CI_gauss_lower <- p - qnorm(1-alpha/2) * sqrt(sigma2)
  CI_gauss_upper <- p + qnorm(1-alpha/2) * sqrt(sigma2)
  
  # randomization
  rand_p <- replicate(BS.iter, rand_rel_treat_eff(n, times1, cens1, times2, cens2, tau, p))
  # bootstrap
  bs_p <- replicate(BS.iter, bs_rel_treat_eff(n, times1, cens1, times2, cens2, tau, p))
  
  
  # CIs: equal mass (i.e. mass alpha/2 to the left and to the right) 
  # two-sided p-value
  # p_val_rand <- 2*min(mean(rand_p[1,] <= T_0), mean(rand_p[1,] >= T_0))
  # 1-sided p-value (right-tailed)
  p_val_rand <- mean(rand_p[1,] >= T_0)
  CI_rand_lower <- p - quantile(rand_p[1,], prob=1-alpha/2)*sqrt(sigma2)
  CI_rand_upper <- p - quantile(rand_p[1,], prob=alpha/2)*sqrt(sigma2)

  # two-sided p-value
  # p_val_bs <- 2*min(mean(bs_p <= T_0), mean(bs_p >= T_0))
  # 1-sided p-value (right-tailed)
  p_val_bs <- mean(bs_p[1,] >= T_0)
  CI_bs_lower <- p - quantile(bs_p[1,], prob=1-alpha/2)*sqrt(sigma2)
  CI_bs_upper <- p - quantile(bs_p[1,], prob=alpha/2)*sqrt(sigma2)

  
  # Now, for the log-log-transformation:
  T_phi <- (log(-log(p))-log(-log(0.5)))/sqrt(sigma2)*p*log(p)
  
  
  # p-values for two-sided (non-randomized) tests: asymptotic, randomization, bootstrap.
  # p_val_phi_gauss <- 2*pnorm(sqrt(n)*abs(T_phi)) -1
  # p_val_phi_rand <- 2*min(mean(rand_p[2,] <= T_phi), mean(rand_p[2,] >= T_phi))
  # p_val_phi_bs <- 2*min(mean(bs_p[2,] <= T_phi), mean(bs_p[2,] >= T_phi))
  
  # p-values for one-sided (right-tailed, non-randomized) tests: asymptotic, randomization, bootstrap.
  p_val_phi_gauss <- 1-pnorm(T_phi)
  p_val_phi_rand <- mean(rand_p[2,] >= T_phi)
  p_val_phi_bs <- mean(bs_p[2,] >= T_phi)
  
  # two-sided confidence intervals
  CI_phi_gauss_lower <- p^exp(-qnorm(1-alpha/2)*sqrt(sigma2)/p/log(p))
  CI_phi_gauss_upper <- p^exp(qnorm(1-alpha/2)*sqrt(sigma2)/p/log(p))
  CI_phi_rand_lower <- p^exp(-quantile(rand_p[2,], prob=1-alpha/2)*sqrt(sigma2)/p/log(p))
  CI_phi_rand_upper <- p^exp(-quantile(rand_p[2,], prob=alpha/2)*sqrt(sigma2)/p/log(p))
  CI_phi_bs_lower <- p^exp(-quantile(bs_p[2,], prob=1-alpha/2)*sqrt(sigma2)/p/log(p))
  CI_phi_bs_upper <- p^exp(-quantile(bs_p[2,], prob=alpha/2)*sqrt(sigma2)/p/log(p))
  
  return(list(RTE=p, RTE.var=sigma2, 
              p.val.gauss=p_val_gauss, CI.gauss = c(CI_gauss_lower, CI_gauss_upper), 
              p.val.rand=p_val_rand, CI.rand = c(CI_rand_lower, CI_rand_upper),
              p.val.bs=p_val_bs, CI.bs = c(CI_bs_lower, CI_bs_upper),
              p.val.phi.gauss=p_val_phi_gauss, CI.phi.gauss = c(CI_phi_gauss_lower, CI_phi_gauss_upper), 
              p.val.phi.rand=p_val_phi_rand, CI.phi.rand = c(CI_phi_rand_lower, CI_phi_rand_upper),
              p.val.phi.bs=p_val_phi_bs, CI.phi.bs = c(CI_phi_bs_lower, CI_phi_bs_upper)))
}


# function to randomize the dataset and then compute the normalized relative treatment effect
rand_rel_treat_eff <- function(n, times1, cens1, times2, cens2, tau, p){
  dataset_rand <- rand(n, times1, cens1, times2, cens2, tau)
  dataset_rand <- paired2CR(n, dataset_rand$times1, dataset_rand$cens1, dataset_rand$times2, dataset_rand$cens2, tau)
  tra <- matrix(ncol=4,nrow=4,FALSE)
  tra[1,2:4] <- TRUE
  etm_complete <- etm(data=dataset_rand, state.names=0:3, tra=tra, cens.name="cens", s=0, t=tau)
  
  if(is.null(etm_complete)) print("etm (rand) is NULL!")
  if(is.null(etm_complete$cov)) return(c(0,0))
  
  latest <- length(etm_complete$est[1,1,])
  
  if(length(etm_complete$cov[1,1,latest-1])==0) return(c(0,0))
  
  p_rand <- etm_complete$est[1,3,latest] + 0.5 * etm_complete$est[1,4,latest]
  sigma2_rand <- 0.25*etm_complete$cov[1,1,latest-1] + etm_complete$cov[9,9,latest-1] + 0.25*etm_complete$cov[13,13,latest-1] + etm_complete$cov[1,9,latest-1] + etm_complete$cov[9,13,latest-1] + 0.25*etm_complete$cov[1,13,latest-1]
  
  # If p_rand is very close to 0 or 1.
  if(is.na((log(-log(p_rand))-log(-log(0.5)))/sqrt(sigma2_rand)*p_rand*log(p_rand))){
    if(abs(sigma2_rand)<0.000001){
      if(abs(p_rand-1) < 0.000001) return(c(Inf, Inf))
      if(abs(p_rand) < 0.000001) return(c(-Inf, -Inf))
      if(abs(p_rand-0.5) < 0.000001) return(c(0,0))
    }else{
      if(abs(p_rand-1) < 0.000001) return(c((p_rand-0.5)/sqrt(sigma2_rand), Inf))
      if(abs(p_rand) < 0.000001) return(c((p_rand-0.5)/sqrt(sigma2_rand), -Inf))
    }
  }
  # Return the untransformed and log-log-transformed standardized randomized relative treatment effect.
  return(c((p_rand-0.5)/sqrt(sigma2_rand), (log(-log(p_rand))-log(-log(0.5)))/sqrt(sigma2_rand)*p_rand*log(p_rand)))
}


# function to bootstrap the dataset and then compute the normalized relative treatment effect
bs_rel_treat_eff <- function(n, times1, cens1, times2, cens2, tau, p){
  bs_indices <- sample(1:n,n,replace=TRUE)
  times1_bs <- times1[bs_indices]
  cens1_bs <- cens1[bs_indices]
  times2_bs <- times2[bs_indices]
  cens2_bs <- cens2[bs_indices]
  dataset_bs <- paired2CR(n, times1_bs, cens1_bs, times2_bs, cens2_bs, tau)
  tra <- matrix(ncol=4,nrow=4,FALSE)
  tra[1,2:4] <- TRUE
  etm_complete <- etm(data=dataset_bs, state.names=0:3, tra=tra, cens.name="cens", s=0, t=tau)
  
  if(is.null(etm_complete)) print("etm (bs) is NULL!")
  if(is.null(etm_complete$cov)) return(c(0,0))
  
  latest <- length(etm_complete$est[1,1,])
  
  if(length(etm_complete$cov[1,1,latest-1])==0) return(c(0,0))
  
  p_bs <- etm_complete$est[1,3,latest] + 0.5 * etm_complete$est[1,4,latest]
  sigma2_bs <- 0.25*etm_complete$cov[1,1,latest-1] + etm_complete$cov[9,9,latest-1] + 0.25*etm_complete$cov[13,13,latest-1] + etm_complete$cov[1,9,latest-1] + etm_complete$cov[9,13,latest-1] + 0.25*etm_complete$cov[1,13,latest-1]
  
  # if p_bs is very close to 0 or 1
  if(is.na((log(-log(p_bs))-log(-log(p)))/sqrt(sigma2_bs)*p_bs*log(p_bs))){
    if(abs(sigma2_bs)<0.000001){
      if(abs(p_bs-1) < 0.000001) return(c(Inf, Inf))
      if(abs(p_bs) < 0.000001) return(c(-Inf, -Inf))
      if(abs(p_bs-p) < 0.000001) return(c(0,0))
    }else{
      if(abs(p_bs-1) < 0.000001) return(c((p_bs-p)/sqrt(sigma2_bs), Inf))
      if(abs(p_bs) < 0.000001) return(c((p_bs-p)/sqrt(sigma2_bs), -Inf))
    }
  }
  
  # Return the untransformed and log-log-transformed standardized bootstrapped relative treatment effect.
  return(c((p_bs-p)/sqrt(sigma2_bs), (log(-log(p_bs))-log(-log(p)))/sqrt(sigma2_bs)*p_bs*log(p_bs)))
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
  Z <- rbind(Z, cbind(runif(n), rep(0,n)), cbind(rep(0,n), runif(n)))
  
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


one.test <- function(itnr, B, alpha=0.05, n, setting=settings[18,], tau, X1=NULL, X2=NULL, d1=NULL, d2=NULL){
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
  
  results <- rel_treat_eff(n, times1=X1, cens1=d1, times2=X2, cens2=d2, tau=tau, alpha=alpha, BS.iter=B, plotting=FALSE)
  
  return(c(results[[1]], results[[3]], results[[5]], results[[7]], results[[9]], results[[11]], results[[13]]))
}


### For simulation of the data sets and conducting the tests.
if(FALSE){
  library(parallel)
  cl <- makeCluster(16, type="FORK", outfile="errors.txt")
  #cl <- makeCluster(16, type="PSOCK", outfile="errors.txt")
  clusterExport(cl=cl, c('dataGen', 'one.test', 'rgumbel', 'rCopula', 'claytonCopula', 'rel_treat_eff', 'paired2CR', 'etm', 'rand_rel_treat_eff', 'rand', 'bs_rel_treat_eff', 'rCopula', 'claytonCopula'))
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
  library(timereg)
  data(diabetes)
  
  # Y for young
  # A for adult
  
  diabetesY <- subset(diabetes, adult == 1)
  diabetesA <- subset(diabetes, adult == 2)
  
  treat.A <- diabetesA[2*(1:83)-1,]
  untreat.A <- diabetesA[2*(1:83),]
  
  timeA1 <- treat.A$time
  statusA1 <- treat.A$status
  timeA2 <- untreat.A$time
  statusA2 <- untreat.A$status
  
  
  treat.Y <- diabetesY[2*(1:114)-1,]
  untreat.Y <- diabetesY[2*(1:114),]
  
  timeY1 <- treat.Y$time
  statusY1 <- treat.Y$status
  timeY2 <- untreat.Y$time
  statusY2 <- untreat.Y$status
  
  set.seed(1234)
  rel_treat_eff(length(timeA1), timeA1, statusA1, timeA2, statusA2, tau=60, BS.iter=2000)
  set.seed(1234)
  rel_treat_eff(length(timeY1), timeY1, statusY1, timeY2, statusY2, tau=60, BS.iter=2000)
}