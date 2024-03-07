# Actual function for computing the relative treatment effect estimator,
# and the corresponding variance estimator.
#  dataf: a dataframe which is sorted according to the event times.
rte2 <- function(dataf){
  # unique sorted event times
  times <- unique(dataf$time[dataf$status>0])
  nt <- length(times)
  
  if(nt==0) return(c(0.5,999999999))
  
  atrisk <-  sapply(times, function(t) sum(dataf$time >= t))[-nt]
  dN1 <- sapply(times, function(t) sum(dataf$status[dataf$time == t]==1))[-nt]
  dN2 <- sapply(times, function(t) sum(dataf$status[dataf$time == t]==2))[-nt]
  dN3 <- sapply(times, function(t) sum(dataf$status[dataf$time == t]==3))[-nt]
  dN <- dN1+dN2+dN3
  
  atrisk.inv <- 1/atrisk
  
  # Kaplan-Meier estimator
  S <- cumprod(1-dN*atrisk.inv)
  
  # Aalen-Johansen estimators
  # left-continuous version of Kaplan-Meier estimator in integral
  F1 <- cumsum(c(1,S[-(nt-1)]) * dN1*atrisk.inv)
  F2 <- cumsum(c(1,S[-(nt-1)]) * dN2*atrisk.inv)
  F2MF1 <- F2-F2[nt-1]-F1+F1[nt-1]
  
  return(c(0.5 * (1 + F2[nt-1] - F1[nt-1]), 0.25*sum((-(F2MF1*dN + S*(dN2-dN1))^2*atrisk.inv + F2MF1^2*dN + 2*S*F2MF1*(dN2-dN1) + S^2*(dN1+dN2))/(atrisk-dN)^2) ))
}