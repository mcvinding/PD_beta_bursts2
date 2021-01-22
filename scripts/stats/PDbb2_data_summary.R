# PD beta burst data summary

#Load data
data <- read.csv('X:\\PD_longrest\\groupanalysis\\data_log.csv')

########################################################################
# DEMOGRAPHICS BY GROUP

# Age

# Sex

# Clinical scores



########################################################################
# DATA CLEANING SUMMARY
# N ICA comp remove
range(data$ica_remove)
mean(data$ica_remove)

# Data segments
range(data$data_length)
mean(data$data_length)
median(data$data_length)
sd(data$data_length)

range(data$drop_pct)
mean(data$drop_pct)
median(data$drop_pct)
