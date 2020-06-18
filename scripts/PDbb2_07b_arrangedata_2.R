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
subj_data$id <- as.factor(paste("0", subj_data$Study_id, sep=""))

# Get "Old" subjects
# TO DO

# Arrange all subject data
sdata <- data.frame(id=subj_data$id,
                    group=subj_data$group,
                    sex=subj_data$sex,
                    age=subj_data$age)


# Reject rejected subjects
# TO DO

# Read N event data
temp <- readMat("neve_table.mat")
PDn1 <- temp$PDn1
PDn2 <- temp$PDn2
ctrln1 <- temp$ctrln1
ctrln2 <- temp$ctrln2
nevent <- c(temp$PDn1,temp$PDn2,temp$ctrln1,temp$ctrln2)
subs <- unlist(c(temp$PD.subs,temp$PD.subs,temp$ctrl.subs,temp$ctrl.subs))
group <- c(rep("ptns",2*length(temp$PD.subs)),rep("ctrl",2*length(temp$ctrl.subs)))
session <- c(rep(c(rep("1",length(temp$PD.subs)), rep("2",length(temp$PD.subs))),2))

neve.data <- data.frame(nevent=nevent,
                        subs=subs,
                        group=group,
                        session=session)

neve.data$nevent.min <- round(neve.data$nevent/3)
