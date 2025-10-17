# Bland-Altman meta-analysis for heart rate

#--------------------------------------------------------
# Heart rate (all conditions)
# Tipton & Shuster Bland-Altman LoA meta-analysis
#--------------------------------------------------------

library(ggplot2)
library(metafor)
library(dplyr)

#--------------------------------------------------------
# Load data
#--------------------------------------------------------
df <- read.csv("hr_MD.csv")
head(df)
dim(df) 

#--------------------------------------------------------
# Drop studies at high risk of bias (sensitivity analysis)
#--------------------------------------------------------
 #df <- df %>%
   #filter(RoB != "high")

 #dim(df)

#--------------------------------------------------------
# Subgroup analysis by optical HR sensor generation
#--------------------------------------------------------
# table(df$optical_hr_sensor)

# df <- df %>%
#   filter(optical_hr_sensor == "third_gen") # change per generation

# dim(df)

#--------------------------------------------------------
# Subgroup analysis by cardiac/patient vs healthy
#--------------------------------------------------------
# table(df$population)

# df <- df%>%
#   filter(population != "healthy")

# dim(df)
#--------------------------------------------------------
# Compute derived values
#--------------------------------------------------------
n  <- df$sample_size
md <- df$mean_diff
s2 <- df$sd_diff_calc^2    # reported sample variance of differences

# Bias-adjusted variance (Tipton & Shuster, eq. 2)
df$s2_star   <- s2 + (md^2)/n

# Variances needed for meta-analysis
df$V_bias    <- df$s2_star / n
df$logs2     <- log(df$s2_star) + 1/(n - 1)
df$V_logs2   <- 2/(n - 1)

#--------------------------------------------------------
# Meta-analysis helper functions
#--------------------------------------------------------

meta <- function(Te, V_T) {
  m      <- length(Te)
  wt_FE  <- 1 / V_T
  T_FE   <- sum(Te * wt_FE) / sum(wt_FE)
  Q      <- sum(wt_FE * (Te - T_FE)^2)
  S1     <- sum(wt_FE)
  S2     <- sum(wt_FE^2)
  tau2   <- max((Q - (m - 1)) / (S1 - S2 / S1), 0)
  
  wt_RE  <- 1 / (V_T + tau2)
  T_RE   <- sum(Te * wt_RE) / sum(wt_RE)
  V_mod  <- 1 / sum(wt_RE)
  V_rve  <- (m / (m - 1)) * sum(wt_RE^2 * (Te - T_RE)^2) / (sum(wt_RE))^2
  
  c(m = m, T_RE = T_RE, tau2 = tau2, V_mod = V_mod, V_rve = V_rve)
}

#--------------------------------------------------------
# LoA maker function from Tipton & Shuster paper
#--------------------------------------------------------
loa_maker <- function(bias, V_bias, logs2, V_logs2) {
  bias_row  <- meta(bias, V_bias)
  logs2_row <- meta(logs2, V_logs2)
  
  m         <- bias_row["m"]
  bias_mean <- bias_row["T_RE"]
  sd2_est   <- exp(logs2_row["T_RE"])
  tau2      <- bias_row["tau2"]
  
  LOA_L     <- bias_mean - 2 * sqrt(sd2_est + tau2)
  LOA_U     <- bias_mean + 2 * sqrt(sd2_est + tau2)
  
  if (m < 3) warning("Few studies: t-quantile and RVE unstable.")
  tcrit     <- qt(0.975, m - 1)
  
  # 95% CI for pooled mean bias 
  CI_bias_lower <- bias_mean - tcrit * sqrt(bias_row["V_mod"])
  CI_bias_upper <- bias_mean + tcrit * sqrt(bias_row["V_mod"])
  
  B1        <- sd2_est^2 / (sd2_est + tau2)
  B2        <- tau2^2    / (sd2_est + tau2)
  V_logT2   <- 2 / sum((V_bias + tau2)^(-2))
  
  V_LOA_mod <- bias_row["V_mod"] + B1 * logs2_row["V_mod"] + B2 * V_logT2
  V_LOA_rve <- bias_row["V_rve"] + B1 * logs2_row["V_rve"] + B2 * V_logT2
  
  CI_L_mod  <- LOA_L - tcrit * sqrt(V_LOA_mod)
  CI_U_mod  <- LOA_U + tcrit * sqrt(V_LOA_mod)
  CI_L_rve  <- LOA_L - tcrit * sqrt(V_LOA_rve)
  CI_U_rve  <- LOA_U + tcrit * sqrt(V_LOA_rve)
  
  c(studies = m, bias = bias_mean, sd = sqrt(sd2_est), tau2 = tau2,
    LOA_L = LOA_L, LOA_U = LOA_U,
    CI_L_mod = CI_L_mod, CI_U_mod = CI_U_mod,
    CI_L_rve = CI_L_rve, CI_U_rve = CI_U_rve,
    CI_bias_lower = CI_bias_lower, CI_bias_upper = CI_bias_upper)
}

#--------------------------------------------------------
# Apply LoA maker function to data
#--------------------------------------------------------
out <- loa_maker(df$mean_diff, df$V_bias, df$logs2, df$V_logs2)

# Assign column names to the output
names(out) <- c(
  "studies", "bias", "sd", "tau2",
  "LOA_L", "LOA_U",
  "CI_L_mod", "CI_U_mod",
  "CI_L_rve", "CI_U_rve",
  "CI_bias_lower", "CI_bias_upper"
)

print(out)

