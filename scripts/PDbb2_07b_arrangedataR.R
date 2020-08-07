# Import data to R and save for further analysis and plotting.
library(R.matlab)
library(xlsx)
library(freesurfer)

# Define paths
wrkdir <- "X://PD_longrest//groupanalysis"
setwd(wrkdir)

# IMPORT SUBJECT DATA
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
save(ndata, file='X://PD_longrest//groupanalysis//sdata.Rdata')

## IMPORT FS ROI STATS
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


## IMPORT N EVENT DATA
# Read N event data
temp <- readMat("neve_data.mat")
hemi <- unlist(temp$hemi)
neve <- temp$neve
subj <- unlist(temp$subjs)

neve.data <- data.frame(nevent=neve,
                        subj=subj,
                        hemi=hemi)
neve.data$nevent.min <- round(neve.data$nevent/3)

# Merge data frames
ndata <- merge(neve.data, sdata, by="subj", all=FALSE)
ndata <- merge(ndata, roi.thick, by=c("subj", "hemi"))

# Save
# save(ndata, file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//ndata.Rdata')
save(ndata, file='X://PD_longrest//groupanalysis//ndata.Rdata')
write.csv(ndata, file='X://PD_longrest//groupanalysis//ndata.csv')

## IMPORT EVENT LENGTH DATA
# Read event length data
temp <- readMat("leneve_data.mat")
hemi <- unlist(temp$hemi)
leneve <- temp$leneve
subj <- unlist(temp$subjs)

len.data <- data.frame(leneve=leneve,
                       subj=subj,
                       hemi=hemi)

len.data$leneve.ms <- len.data$lenev*1000

# Merge data frames
ldata <- merge(len.data, sdata, by="subj", all=FALSE)
ldata <- merge(ldata, roi.thick, by=c("subj", "hemi"))

# Save
save(ldata, file='X://PD_longrest//groupanalysis//ldata.Rdata')
write.csv(ldata, file='X://PD_longrest//groupanalysis//ldata.csv')

## IMPORT EVENT POW DATA
# Read event power data
temp <- readMat("maxeve_data.mat")
hemi <- unlist(temp$hemi)
maxeve <- temp$maxeve
subj <- unlist(temp$subj)

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
hemi <- unlist(temp$hemi)
toeeve <- temp$toeeve
subj <- unlist(temp$subjs)

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

#END