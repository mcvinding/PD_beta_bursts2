# Import data to R and save for further analysis and plotting.
library(R.matlab)
library(xlsx)
library(freesurfer)

# Define paths
wrkdir <- "X://PD_longrest//groupanalysis"
setwd(wrkdir)

###########################################################################################
# %%% IMPORT SUBJECT DATA %%%
###########################################################################################
# Get "New" subjects
# subj_data <- read.xlsx2('C://Users//Mikkel//Documents//PDbb2//groupanalysis//subj_data_anonymised.xlsx', 1)
subj_data <- read.xlsx2('X://PD_long//subj_data//subj_data_anonymised.xlsx', 1)
subj_data$age <- as.numeric(as.character(subj_data$Age.))
subj_data$sex <- as.factor(ifelse(subj_data$Sex.== 'M' | subj_data$Sex. == "M ", 'M', "F"  ))
subj_data$group <- as.factor(ifelse(subj_data$Type== 'Control' | subj_data$Type == "Control ", 'control', "patient"))
subj_data$subj <- as.factor(paste("0", subj_data$Study_id, sep=""))

# Arrange all subject data
sdata1 <- data.frame(subj=subj_data$subj,
                     group=subj_data$group,
                     sex=subj_data$sex,
                     age=subj_data$age)

# Get "Old" subjects
load(file='C://Users//Mikkel//Documents//betabursts//subj_data//alldata.RData')
sdata2 <- data.frame(subj=paste("0",alldata$MEG_ID, sep=""),
                         group=alldata$Sub_type,
                         sex=alldata$sex,
                         age=alldata$age)

# Combine 
sdata <- rbind(sdata1, sdata2)

# Reject rejected subjects
check.data <- read.xlsx2('X://PD_long//subj_data//subjects_and_dates.xlsx', 1)
check.data <- subset(check.data, rest_ec==1)

sdata <- subset(sdata, sdata$subj %in% check.data$id)

# Save
save(sdata, file='X://PD_longrest//groupanalysis//sdata.Rdata')

# Reload
load(file='X://PD_longrest//groupanalysis//sdata.Rdata')

###########################################################################################
# %%% IMPORT FS ROI STATS %%%
###########################################################################################

roi.lh.thick <- data.frame(subj=sdata$subj,
                           thick = rep(NaN, length(sdata$subj)),
                           hemi = rep('lh', length(sdata$subj)))
                                      
roi.rh.thick <- data.frame(subj=sdata$subj,
                           thick = rep(NaN, length(sdata$subj)),
                           hemi = rep('rh', length(sdata$subj)))

for (ii in 1:length(sdata$subj)){
  subj = sdata$subj[ii]
  tmp.lh <- read_fs_table(paste('X://PD_longrest//meg_data//',subj,'//',subj,'.lh.sensmotor.stats', sep=""), head=FALSE)
  roi.lh.thick$thick[roi.lh.thick$subj==subj] <- tmp.lh$V5
  
  tmp.rh <- read_fs_table(paste('X://PD_longrest//meg_data//',subj,'//',subj,'.rh.sensmotor.stats', sep=""), head=FALSE)
  roi.rh.thick$thick[roi.rh.thick$subj==subj] <- tmp.rh$V5
  
  rm(tmp.lh, tmp.rh)
}

roi.thick <- rbind(roi.lh.thick, roi.rh.thick)
roi.thick$hemi <- as.factor(roi.thick$hemi)

# Save
save(roi.thick, file='X://PD_longrest//groupanalysis//thickdata.Rdata')

# Reload
load(file='X://PD_longrest//groupanalysis//thickdata.Rdata')

# Selct only left hemi for now!
lh.roi.thick <- subset(roi.thick, hemi=='lh')

###########################################################################################
# %%% IMPORT N EVENT DATA %%%
###########################################################################################
# Read N event data
temp <- readMat("neve_b_m1_data.mat")
# hemi <- as.factor(unlist(temp$hemiN))
# nevent.b.m1 <- temp$nevent.b.m1
subj <- as.factor(unlist(temp$subjects))

temp2 <- readMat("neve_b_m2_data.mat")
temp3 <- readMat("neve_b_pc_data.mat")
temp4 <- readMat("neve_u_m1_data.mat")
temp5 <- readMat("neve_u_m2_data.mat")
temp6 <- readMat("neve_u_pc_data.mat")

neve.data <- data.frame(nevent.b.m1=temp$nevent.b.m1,
                        nevent.b.m2=temp2$nevent.b.m2,
                        nevent.b.pc=temp3$nevent.b.pc,
                        nevent.u.m1=temp4$nevent.u.m1,
                        nevent.u.m2=temp5$nevent.u.m2,
                        nevent.u.pc=temp6$nevent.u.pc,
                        subj=subj)
neve.data$nevent.min <- round(neve.data$nevent/3)

# Merge data frames
ndata <- merge(neve.data, sdata, by="subj", all=FALSE)
ndata <- merge(ndata, lh.roi.thick, by=c("subj"))  # Add hemi later

# Save
save(ndata, file='X://PD_longrest//groupanalysis//ndata_all.Rdata')
# write.csv(ndata, file='X://PD_longrest//groupanalysis//ndata.csv')

###########################################################################################
# %%% IMPORT EVENT LENGTH DATA %%%
###########################################################################################

# Read event length data
temp <- readMat("leneve_b_m1.mat")
# hemi <- as.factor(unlist(temp$hemi))
leneve <- temp$len.b.m1
subj <- as.factor(unlist(temp$sub.b.m1))

len.data.b.m1 <- data.frame(leneve=leneve,
                            subj=subj)
len.data.b.m1$leneve.ms <- len.data$lenev*1000

# Merge data frames
ldata <- merge(len.data, sdata, by="subj", all=FALSE)
ldata <- merge(ldata, roi.thick, by=c("subj", "hemi"))

# Save
save(ldata, file='X://PD_longrest//groupanalysis//ldata.Rdata')
write.csv(ldata, file='X://PD_longrest//groupanalysis//ldata.csv')

## IMPORT EVENT POW DATA
# Read event power data
temp <- readMat("maxeve_data.mat")
hemi <- as.factor(unlist(temp$hemi))
maxeve <- temp$maxeve
subj <- as.factor(unlist(temp$subj))

max.data <- data.frame(maxeve=maxeve,
                       subj=subj,
                       hemi=hemi)

# Merge data frames
mdata <- merge(max.data, sdata, by="subj", all=FALSE)
mdata <- merge(mdata, roi.thick, by=c("subj", "hemi"))

# Save
save(mdata, file='X://PD_longrest//groupanalysis//mdata.Rdata')
write.csv(mdata, file='X://PD_longrest//groupanalysis//mdata.csv')

## IMPORT EVENT INTERVAL DATA
# Read event interval data
temp <- readMat("toeeve_data.mat")
hemi <- as.factor(unlist(temp$hemi))
toeeve <- temp$toeeve
subj <- as.factor(unlist(temp$subjs))

toe.data <- data.frame(toeeve=toeeve,
                       subj=subj,
                       hemi=hemi)

toe.data$toeeve.ms <- toe.data$toeeve*1000

# Merge data frames
tdata <- merge(toe.data, sdata, by="subj", all=FALSE)
tdata <- merge(tdata, roi.thick, by=c("subj", "hemi"))

# Save
save(tdata, file='X://PD_longrest//groupanalysis//tdata.Rdata')
write.csv(tdata, file='X://PD_longrest//groupanalysis//tdata.csv')


# MERGE ALL

bbdata <- data.frame(leneve=leneve,
                     toeeve=toeeve,
                     maxeve=maxeve,
                     subj=subj,
                     hemi=hemi)

bbdata <- merge(bbdata, sdata, by="subj", all=FALSE)
bbdata <- merge(bbdata, roi.thick, by=c("subj", "hemi"))
bbdata$leneve.ms <- bbdata$leneve*1000
bbdata$toeeve.ms <- bbdata$toeeve*1000

# Save
save(bbdata, file='X://PD_longrest//groupanalysis//bbdata.Rdata')
write.csv(bbdata, file='X://PD_longrest//groupanalysis//bbdata.csv')

#END