###########################################################################################
# Import data to R and save for further analysis and plotting.
# Subject metadata, clinical data, event length, peak amplitude and time until next event (tue)
# data per individual event.
# Arrange into dataframe with one event per row
###########################################################################################
library(R.matlab)
library(xlsx)

# Define paths
wrkdir <- "X://PD_longrest//groupanalysis"
setwd(wrkdir)

###########################################################################################
# %%% Load SUBJECT METADATA %%%
###########################################################################################
# load(file='X://PD_longrest//groupanalysis//alldata_subj.Rdata')
# load(file='X://PD_longrest//groupanalysis//sdata.Rdata')
load(file='X://PD_longrest//groupanalysis//clindata.Rdata')
load(file='X://PD_longrest//groupanalysis//thickdata.Rdata')
lh.roi.thick <- subset(roi.thick, hemi=='lh') # Selct only left hemi for now!

###########################################################################################
# %%% IMPORT EVENT LENGTH DATA %%%
###########################################################################################
temp <- readMat("leneve_u_m2.mat")
leneve <- temp$len.u.m2
subj <- as.factor(unlist(temp$sub.u.m2))
# idxr <- as.factor(1:length(leneve))
# 
# len.data <- data.frame(leneve=leneve,
#                        subj=subj,
#                        idxr=idxr)
# len.data$leneve.ms <- len.data$lenev*1000

###########################################################################################
# %%%  IMPORT EVENT POW DATA %%%
###########################################################################################
temp <- readMat("maxeve_u_m2.mat")
maxeve <- temp$max.u.m2
subj <- as.factor(unlist(temp$sub.u.m2))
# idxr <- as.factor(1:length(maxeve))
# 
# max.data <- data.frame(maxeve=maxeve,
#                        subj=subj,
#                        idxr=idxr)

###########################################################################################
# %%% IMPORT EVENT INTERVAL DATA %%%
###########################################################################################
temp <- readMat("tueeve_u_m2.mat")
tueeve <- temp$tue.u.m2
subj <- as.factor(unlist(temp$sub.u.m2))
# idxr <- as.factor(1:length(maxeve))
# 
# tue.data <- data.frame(tueeve=tueeve,
#                        subj=subj,
#                        idxr=idxr)
# tue.data$toeeve.ms <- toe.data$toeeve*1000

###########################################################################################
# %%% COLLECT ALL DATA INTO ONE %%%
###########################################################################################
bbdata <- data.frame(leneve=leneve,
                     tueeve=tueeve,
                     maxeve=maxeve,
                     subj=subj)
bbdata$leneve.ms <- bbdata$leneve*1000
bbdata$tueeve.ms <- bbdata$tueeve*1000
bbdata <- merge(bbdata, clindata, by="subj", all=FALSE)

# Save
save(bbdata, file='X://PD_longrest//groupanalysis//bbdata.Rdata')
write.csv(bbdata, file='X://PD_longrest//groupanalysis//bbdata.csv')

#END