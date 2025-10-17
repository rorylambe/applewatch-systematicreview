# This file contains the example code provided by Tipton & Shuster
# The meta-analyses we conducted for heart rate and blood oxygen saturation
# are based on this sample code.

#--------------------------------------------------------
#Functions for calculating and comparing Bland-altman Meta-analyses
#--------------------------------------------------------

meta <- function(Te,V_T) {
  #T #vector of estimates (either logS2 or bias)
  #V_T #vector of conditional variances
  #method #RE or RVE
  
  m <- length(Te)
  wt_FE <- 1/V_T
  T_FE <- sum(Te*wt_FE)/sum(wt_FE)
  Q <- sum(wt_FE*(Te - T_FE)^2)
  S1 <- sum(wt_FE)
  S2 <- sum(wt_FE^2)
  o2 <- (Q - (m - 1))/(S1 - S2/S1)
  
  wt_RE <- 1/(V_T + o2)
  T_RE <- sum(Te*wt_RE)/sum(wt_RE)
  V_T_RE_mod <- 1/sum(wt_RE)
  V_T_RE_rve <- (m/(m-1))*sum(wt_RE^2*(Te - T_RE)^2)/(sum(wt_RE))^2
  
  c(m,T_RE,o2,V_T_RE_mod, V_T_RE_rve)
  
}

loa_maker <- function(bias,V_bias,logs2,V_logs2) {
  
  bias_row <- meta(bias,V_bias)
  logs2_row <- meta(logs2,V_logs2)
  
  bias_mean <- bias_row[2]
  sd2_est <- exp(logs2_row[2])
  tau_est <- bias_row[3]
  
  LOA_L <- bias_mean - 2*sqrt(sd2_est + tau_est)
  LOA_U <- bias_mean + 2*sqrt(sd2_est + tau_est)
  
  m <- bias_row[1]
  tcrit <- qt(1-.05/2,m-1)
  B1 <- sd2_est^2/(sd2_est + tau_est)
  B2 <- tau_est^2/(sd2_est + tau_est)
  
  wt <- 1/V_bias
  #	S1 <- sum(wt)
  #	S2 <- sum(wt^2)
  #	S3 <- sum(wt^3)
  #	A0 <- 2*(m-1)/(S1-S2/S1)^2
  #	A1 <- 4/(S1 - S2/S1)
  #	A2 <- 2*(S2-2*S3/S1+S2^2/S1^2)/(S1-S2/S1)^2
  #	V_logT2 <- A0/tau_est^2 + A1/tau_est + A2
  
  V_logT2 <- 2/sum((V_bias + tau_est)^(-2))
  
  V_LOA_mod <- bias_row[4] + B1*logs2_row[4] + B2*V_logT2
  V_LOA_rve <- bias_row[5] + B1*logs2_row[5] + B2*V_logT2
  
  CI_L_mod <- LOA_L - tcrit*sqrt(V_LOA_mod)
  CI_U_mod <- LOA_U + tcrit*sqrt(V_LOA_mod)
  
  CI_L_rve <- LOA_L - tcrit*sqrt(V_LOA_rve)
  CI_U_rve <- LOA_U + tcrit*sqrt(V_LOA_rve)
  
  c(m, bias_mean,sqrt(sd2_est), tau_est, LOA_L, LOA_U, CI_L_mod, CI_U_mod, CI_L_rve, CI_U_rve)
  
} 

#--------------------------------------------------------
#Doing Bland-Altman Meta-analyses
#--------------------------------------------------------

#bias = vector of study specific bias estimates
#V_bias = vector of variances of study specific bias estimates (i.e., s2*/n)

#logs2 = vector of study specific estimates, where logs2 = log(s2*) + 1/(n-1) and s2* is bias adjusted
#V_logs2 = vector of variances (i.e., 2/(n-1))

out <- loa(bias,V_bias,logs2,V_logs2)
names(out) <- c("studies","bias","sd","tau2","LOA_L","LOA_U","CI_Lm","CI_Um","CI_Lr","CI_Ur")

#LOA_L, LOA_U = LOA lower and upper bounds
#CI_Lm, CI_Um = model-based random-effects meta-analysis estimate of CI for LOA
#CI_Lr, CI_Ur = robust variance estimation meta-analysis estimate of CI for LOA