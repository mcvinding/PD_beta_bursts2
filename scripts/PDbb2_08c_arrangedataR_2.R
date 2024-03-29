###########################################################################################
# Import data to R and save for further analysis and plotting. Subject metadata, clinical 
# data, event length, peak amplitude and time until next event (tue) data per individual 
# event. Arrange into data frame with one event per row
#
# Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease. medRxiv.org. https://doi.org/10.1101/2021.06.27.21259592
#
# @mcvinding
###########################################################################################
library(R.matlab)
library(xlsx)

# Define paths
wrkdir <- "/home/mikkel/PD_longrest/groupanalysis"
setwd(wrkdir)

###########################################################################################
# %%% Load SUBJECT METADATA %%%
###########################################################################################
# load(file='alldata_subj2.Rdata')
load(file='sdata.Rdata')
load(file='clindata.Rdata')
load(file='thickdata.Rdata')

###########################################################################################
# %%% IMPORT EVENT LENGTH DATA %%%
###########################################################################################
temp <- readMat("leneve_u_m22.mat")
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
temp <- readMat("maxeve_u_m22.mat")
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
temp <- readMat("tueeve_u_m22.mat")
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
tmp1 <- data.frame(leneve=leneve,
                     tueeve=tueeve,
                     maxeve=maxeve,
                     subj=subj)
tmp1$leneve.ms <- tmp1$leneve*1000
tmp1$tueeve.ms <- tmp1$tueeve*1000
tmp2 <- merge(tmp1, sdata, by="subj", all=FALSE)
tmp3 <- merge(tmp2, roi.thick, by="subj", all=FALSE)
clintmp <- clindata[,c(1,9:17)]
bbdata <- merge(tmp3, clintmp, by="subj", all=TRUE)

# Save
save(bbdata, file='bbdata2.Rdata')
write.csv(bbdata, file='bbdata2.csv')

#END