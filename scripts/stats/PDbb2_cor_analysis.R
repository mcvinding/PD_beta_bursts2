## Correlation analysis
library(dplyr)
library(corrplot)
library(Hmisc)

# Define paths
wrkdir <- "X://PD_longrest//groupanalysis"
# wrkdir <- "C://Users//Mikkel//Documents//PDbb2//groupanalysis"
setwd(wrkdir)


# Load data
load('alldata_subj2.Rdata')
load('bbdata2.Rdata')

# Get subject level bb measures
bbsum <- aggregate(cbind(leneve, tueeve, maxeve)~subj, data=bbdata, FUN=median)
alldata <- merge(alldata, bbsum, by="subj")


# Re-arrange data
cordat <- select(alldata, nevent.u.m2.min,nevent.u.m1.min,nevent.b.m2.min,nevent.b.m1.min,
                 leneve,tueeve,maxeve,
                 a_intercept,a_slope,beta_pw,beta_cf,alpha_pw,alpha_cf,
                 age,MoCA,thick,U.F1,U.F2,U.F3,U.F4,U.F5,U.F6,U.F7)

# Correlation
cormat <- rcorr(as.matrix(cordat))

cormat <- cor(cordat, use="pairwise.complete.obs")

col<- colorRampPalette(c("blue", "white", "red"))(20)
corrplot(cormat$r, method="number",p.mat=cormat$P, type="full",sig.level=0.05, insig="blank",diag=F,col=col)
corrplot(cormat$r)

# Export
write.csv(cormat, file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//corrmat2.csv')
write.csv(cordat, file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//corrdat2.csv')





