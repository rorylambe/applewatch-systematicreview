# Meta-analysis of AFib detection
# Bivariate meta-analysis using Reitsma

#--------------------------------------------------------
# Install / load required packages
#--------------------------------------------------------
# install.packages(c("readr", "dplyr", "mada", "meta", "metafor"))
library(readr)
library(dplyr)
library(mada)
library(meta)
library(metafor)

#--------------------------------------------------------
# Load data
#--------------------------------------------------------
data <- read.csv('afib_data1.csv')
head(data)

#--------------------------------------------------------
# Drop studies at high risk of bias (sensitivity analysis)
#--------------------------------------------------------
# data <- data %>%
#   filter(RoB != "high")

# dim(data)

#--------------------------------------------------------
# Subgroup analysis by optical HR sensor generation / ECG app 1.0 vs 2.0
#--------------------------------------------------------
# table(data$optical_hr_sensor)

# keep only studies with third gen HR sensor (Series 6 onward), and ECG app 2.0
# data <- data%>%
#   filter(optical_hr_sensor == "second_gen")

# dim(data)

#--------------------------------------------------------

# descriptive stats, incl sens, spec, false-pos rate
madad(data)

# bivariate meta-analysis with reitsma
model <- (fit.reitsma <- reitsma(data))

summary(model)

# SROC curve
plot(fit.reitsma, sroclwd = 2,
     main = "AFib detection: SROC curve (bivariate model)")
points(fpr(data), sens(data), pch = 2)
legend("bottomright", c("Individual study data", "Summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "Conf. region"), lwd = c(2,1))


