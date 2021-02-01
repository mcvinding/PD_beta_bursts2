###########################################################################################
# Import data to R and save for further analysis and plotting.
# Subject metadata, clinical data, N events data per subject
# Arrange into dataframe with one subject per row
###########################################################################################
library(R.matlab)
library(xlsx)
library(freesurfer)

# Define paths
wrkdir <- "X://PD_longrest//groupanalysis"
setwd(wrkdir)

# Create array to mark rejected subjects
check.data <- read.xlsx2('X://PD_long//subj_data//subjects_and_dates.xlsx', 1)
check.data$id <- as.factor(check.data$i)
check.data <- subset(check.data, rest_ec=="1")

###########################################################################################
# %%% IMPORT SUBJECT DATA %%%
###########################################################################################
# Get "New" subjects
subj_data <- read.xlsx2('X://PD_long//subj_data//subj_data_anonymised.xlsx', 1)
subj_data$age <- as.numeric(as.character(subj_data$Age.))
subj_data$sex <- as.factor(ifelse(subj_data$Sex.== 'M' | subj_data$Sex. == "M ", 'M', "F"  ))
subj_data$group <- as.factor(ifelse(subj_data$Type== 'Control' | subj_data$Type == "Control ", 'control', "patient"))
subj_data$subj <- as.factor(paste("0", subj_data$Study_id, sep=""))
subj_data$Years_Edu <- as.numeric(as.character(subj_data$Years_Edu))
subj_data$FAB <- as.numeric(as.character(subj_data$FAB))
subj_data$MoCA <- as.numeric(as.character(subj_data$MoCA))
subj_data$BDI <- as.numeric(as.character(subj_data$BDI))

# Arrange all subject data
sdata1 <- data.frame(subj      = subj_data$subj,
                     group     = subj_data$group,
                     sex       = subj_data$sex,
                     age       = subj_data$age,
                     edu_years = subj_data$Years_Edu,
                     FAB       = subj_data$FAB,
                     MoCA      = subj_data$MoCA,
                     BDI       = subj_data$BDI)

# Get "Old" subjects
load(file='C://Users//Mikkel//Documents//betabursts//subj_data//alldata.RData')
sdata2 <- data.frame(subj      = paste("0",alldata$MEG_ID, sep=""),
                     group     = alldata$Sub_type,
                     sex       = alldata$sex,
                     age       = alldata$age,
                     edu_years = alldata$edu_yrs,
                     FAB       = rep(NA, 40),
                     MoCA      = alldata$MoCA,
                     BDI       = rep(NA, 40))

# Combine 
sdata <- rbind(sdata1, sdata2)

# Reject rejected subjects
sdata <- subset(sdata, sdata$subj %in% check.data$id)

# Save
save(sdata, file='X://PD_longrest//groupanalysis//sdata.Rdata')

# (re)load
load(file='X://PD_longrest//groupanalysis//sdata.Rdata')

###########################################################################################
# %%% IMPORT CLINICAL TEST DATA %%%
###########################################################################################
utemp <- read.csv('X://PD_long//subj_data//UPDRS_PD_MEG_subscales.csv', sep=";")
udata1 <- data.frame(subj     = as.factor(paste("0",utemp$Study_id, sep="")),
                     HY.stage = as.factor(as.numeric(as.character(utemp$H_Y_STAGE))),
                     UPDRS    = as.numeric(as.character(utemp$UPDRS_TOTAL)),
                     U.F1     = as.numeric(as.character(utemp$F1)),
                     U.F2     = as.numeric(as.character(utemp$F2)),
                     U.F3     = as.numeric(as.character(utemp$F3)),
                     U.F4     = as.numeric(as.character(utemp$F4)),
                     U.F5     = as.numeric(as.character(utemp$F5)),
                     U.F6     = as.numeric(as.character(utemp$F6)),
                     U.F7     = as.numeric(as.character(utemp$F7)))

utemp2 <- read.xlsx('C://Users//Mikkel//Documents//betabursts//groupanalysis//UPDRS_raw.xlsx',1,header=T)
utemp2 <- subset(utemp2, session=="2")
udata2 <- data.frame(id     = as.factor(utemp2$id),
                     UPDRS  = utemp2$Total,
                     U.F1   = as.numeric(as.character(utemp2$F1)),
                     U.F2   = as.numeric(as.character(utemp2$F2)),
                     U.F3   = as.numeric(as.character(utemp2$F3)),
                     U.F4   = as.numeric(as.character(utemp2$F4)),
                     U.F5   = as.numeric(as.character(utemp2$F5)),
                     U.F6   = as.numeric(as.character(utemp2$F6)),
                     U.F7   = as.numeric(as.character(utemp2$F7)))

load(file='C://Users//Mikkel//Documents//betabursts//subj_data//alldata.RData')
hy.temp = data.frame(subj     = as.factor(paste("0",alldata$MEG_ID, sep="")),
                     HY.stage = as.factor(alldata$HY_stage),
                     id       = alldata$ID)

udata2 <- merge(udata2, hy.temp, by=c("id"))
udata2 <- subset(udata2, select=-c(id))
udata <- rbind(udata1, udata2)

## MERGE
clindata <- merge(sdata, udata, by.x="subj", by.y="subj")

# SAVE
save(clindata, file='X://PD_longrest//groupanalysis//clindata.Rdata')
write.csv(clindata, file='X://PD_longrest//groupanalysis//clindata.csv')

# (re)load
load(file='X://PD_longrest//groupanalysis//clindata.Rdata')

###########################################################################################
# %%% IMPORT FS ROI STATS %%%
###########################################################################################
roi.thick <- data.frame(subj=sdata$subj,
                        thick = rep(NaN, length(sdata$subj)),
                        hemi = rep('lh', length(sdata$subj)))

for (ii in 1:length(sdata$subj)){
  subj = sdata$subj[ii]
  tmp.lh <- read_fs_table(paste('X://PD_longrest//meg_data//',subj,'//',subj,'.lh.sensmotor.stats', sep=""), head=FALSE)
  roi.thick$thick[roi.thick$subj==subj] <- tmp.lh$V5
  rm(tmp.lh)
}

roi.thick$hemi <- as.factor(roi.thick$hemi)

# Save
save(roi.thick, file='X://PD_longrest//groupanalysis//thickdata.Rdata')

# Reload
load(file='X://PD_longrest//groupanalysis//thickdata.Rdata')

###########################################################################################
# %%% IMPORT PSD DATA %%%
###########################################################################################

fooof.data <- read.csv('X://PD_longrest//groupanalysis//fooof_df2.csv', sep=";")
fooof.data$subj <- as.factor(paste("0",fooof.data$subj, sep=""))

###########################################################################################
# %%% IMPORT N EVENT DATA %%%
###########################################################################################
# Read N event data
temp <- readMat("neve_u_m2_data2.mat")

subj <- as.factor(unlist(temp$subjects))

# XXX: SORT !!!!
temp2 <- readMat("neve_b_m2_data2.mat")
temp3 <- readMat("neve_b_pc_data2.mat")
temp6 <- readMat("neve_u_pc_data2.mat")

neve.data <- data.frame(nevent.b.m2=temp2$nevent.b.m2,
                        nevent.u.m2=temp$nevent.u.m2,
                        subj=subj)

# Get data log and BPM
data.log <- read.csv('X:\\PD_longrest\\groupanalysis\\data_log.csv', sep=";")
data.log$subj <- as.factor(data.log$subj)
data.log$subj <- paste("0",data.log$subj, sep="")

neve.data <- merge(neve.data, data.log, by="subj", all=FALSE)
neve.data$data_length.min <- neve.data$data_length/60
neve.data$nevent.b.m2.min <- round(neve.data$nevent.b.m2/neve.data$data_length.min)
neve.data$nevent.u.m2.min <- round(neve.data$nevent.u.m2/neve.data$data_length.min)

# Save
save(ndata, file='X://PD_longrest//groupanalysis//neve.data.Rdata')

# (re)load
load(file='X://PD_longrest//groupanalysis//neve.data.Rdata')

###########################################################################################
# %%% LEDD and disease duration
###########################################################################################
tmp <- read.xlsx('X:\\PD_long\\subj_data\\ptns_medication.xlsx', 1)
tmp$yr1 <- as.numeric(substr(tmp$DATE, start=2, stop=5))
tmp$yr2 <- as.numeric(tmp$Initial_diagnosis)
tmp$pd.dur <- tmp$yr1-tmp$yr2

tmp1.dat <- data.frame(subj=as.factor(paste('0',as.character(tmp$NATID), sep="")),
                       ledd=tmp$LEDD,
                       pd.dur=tmp$pd.dur)

load(file='C://Users//Mikkel//Documents//betabursts//subj_data//alldata.RData')

tmp2.dat <- data.frame(subj=as.factor(paste('0',as.character(alldata$MEG_ID), sep="")),
                       ledd=alldata$LEDD,
                       pd.dur=alldata$disease_dur)

med.dat <- rbind(tmp1.dat, tmp2.dat)

###########################################################################################
# %%% COLLECT ALL DATA INTO ONE %%%
###########################################################################################
# Merge data frames
ndata <- merge(neve.data, sdata, by="subj", all=TRUE)
tmpdat1 <- merge(ndata, roi.thick, by=c("subj"), all=TRUE)
tmpdat2 <- merge(tmpdat1, fooof.data, by=c("subj"), all=TRUE)
tmpdat3 <- merge(tmpdat2, med.dat, by=c("subj"), all=TRUE)
alldata <- merge(tmpdat3, udata, by=c("subj"), all.x=TRUE)
alldata <- subset(alldata, alldata$subj %in% check.data$id)

save(alldata, file='X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
save(alldata, file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj2.Rdata')
write.csv(alldata, file='X://PD_longrest//groupanalysis//alldata_subj2.csv')

#END