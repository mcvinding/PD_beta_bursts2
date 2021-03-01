# Demographics and data summaries
library(read)
library(dplyr)
library(xlsx)

########################################################################
# Get summary of all subjects
########################################################################
# Load data
subj.data <- read.xlsx2('X://PD_long//subj_data//subj_data_anonymised.xlsx',1)
subj.data$subj <- as.factor(paste("0", subj.data$Study_id, sep=""))

check.data <- read.xlsx2('X://PD_long//subj_data//subjects_and_dates.xlsx', 1)
alldat <- merge(subj.data, check.data, by.x="subj", by.y="id")
alldat$group <- as.factor(ifelse(alldat$Type=="Control" | alldat$Type=="Control ", "control", "patient"))
alldat$sex <- as.factor(ifelse(alldat$Sex.=="M" | alldat$Type=="M ", "M", "F"))
alldat$age <- as.numeric(alldat$Age.)

# Get "Old" subjects
load(file='C://Users//Mikkel//Documents//betabursts//subj_data//alldata.RData')
sdata <- data.frame(subj      = as.factor(paste("0",alldata$MEG_ID, sep="")),
                    group     = alldata$Sub_type,
                    sex       = alldata$sex,
                    age       = alldata$age)

sdata <- subset(sdata, sdata$subj %in% check.data$id)
summary(sdata$group)

# Merge
alldat <- bind_rows(alldat, sdata)

# N subj
summary(alldat$group)

# Age by group
aggregate(alldat$age, by=list(alldat$group), range)
aggregate(alldat$age, by=list(alldat$group), mean)
aggregate(alldat$age, by=list(alldat$group), sd)

# Sex by group
xtabs(~sex+group, alldat)

########################################################################
# Get summary of included subjects
########################################################################
load('X://PD_longrest//groupanalysis//alldata_subj2.Rdata')

# N subj
summary(alldata$group)

# Sex by group
xtabs(~sex+group, alldata)
chisq.test(xtabs(~sex+group, alldata))

# Age by group
aggregate(alldata$age, by=list(alldata$group), range)
aggregate(alldata$age, by=list(alldata$group), mean)
aggregate(alldata$age, by=list(alldata$group), sd)
t.test(alldata$age~alldata$group)

# Disease duration
mean(alldata$pd.dur, na.rm=TRUE)
median(alldata$pd.dur, na.rm=TRUE)
sd(alldata$pd.dur, na.rm=TRUE)

# LEDD
mean(alldata$ledd, na.rm=TRUE)
median(alldata$ledd, na.rm=TRUE)
sd(alldata$ledd, na.rm=TRUE)

# MDS UPDRS-III
mean(alldata$UPDRS, na.rm=TRUE)
median(alldata$UPDRS, na.rm=TRUE)
sd(alldata$UPDRS, na.rm=TRUE)

# MOCA
aggregate(alldata$MoCA, by=list(alldata$group), range)
aggregate(alldata$MoCA, by=list(alldata$group), mean)
aggregate(alldata$MoCA, by=list(alldata$group), sd)
t.test(alldata$MoCA~alldata$group)


########################################################################
# DATA CLEANING SUMMARY
# N ICA comp remove
range(alldata$ica_remove)
mean(alldata$ica_remove)
median(alldata$ica_remove)
wilcox.test(alldata$ica_remove[alldata$group=="patient"], alldata$ica_remove[alldata$group=="control"], alternative = "two.sided")

# Data segments
range(alldata$drop_pct)
mean(alldata$drop_pct)
median(alldata$drop_pct)
sd(alldata$drop_pct)

range(alldata$data_length)
mean(alldata$data_length)
median(alldata$data_length)
sd(alldata$data_length)
t.test(alldata$data_length~alldata$group)

#END