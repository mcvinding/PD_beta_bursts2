# Import data to R and save for further analysis and plotting.
library(R.matlab)
library(xlsx)

# Define paths
wrkdir <- "C://Users//Mikkel//Documents//PDbb2//groupanalysis"
setwd(wrkdir)

# IMPORT SUBJECT DATA
# Get "New" subjects
subj_data <- read.xlsx2('C://Users//Mikkel//Documents//PDbb2//groupanalysis//subj_data_anonymised.xlsx', 1)
subj_data$age <- as.numeric(as.character(subj_data$Age.))
subj_data$sex <- as.factor(ifelse(subj_data$Sex.== 'M' | subj_data$Sex. == "M ", 'M', "F"  ))
subj_data$group <- as.factor(ifelse(subj_data$Type== 'Control' | subj_data$Type == "Control ", 'Control', "Patient"  ))
subj_data$subj <- as.factor(paste("0", subj_data$Study_id, sep=""))

# Get "Old" subjects
# TO DO

# Arrange all subject data
sdata <- data.frame(subj=subj_data$subj,
                    group=subj_data$group,
                    sex=subj_data$sex,
                    age=subj_data$age)


# Reject rejected subjects
# TO DO

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

save(ndata, file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//ndata.Rdata')


