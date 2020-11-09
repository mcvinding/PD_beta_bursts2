## Correlation analysis
library(dplyr)
library(corrplot)
library(Hmisc)

# Load data
# load(file='X://PD_longrest//groupanalysis//alldata_subj.Rdata')
load('C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj.Rdata')
load('C://Users//Mikkel//Documents//PDbb2//groupanalysis//bbdata.Rdata')


# Get subject level bb measures
bbsum <- aggregate(cbind(leneve, tueeve, maxeve)~subj, data=bbdata, FUN=median)
alldata <- merge(alldata, bbsum, by="subj")


# Re-arrange data
cordat <- select(alldata, nevent.u.m2,leneve,tueeve,maxeve,
                  age,MoCA,thick,UPDRS,
                  a_intercept,a_slope,beta_pw,beta_cf,alpha_pw,alpha_cf,
                  skw.u,krt.u)

# Correlation
cormat <- rcorr(as.matrix(cordat))

cormat <- cor(cordat, use="pairwise.complete.obs")

col<- colorRampPalette(c("blue", "white", "red"))(20)
corrplot(cormat$r, method="number",p.mat=cormat$P, type="full",sig.level=0.05, insig="blank",diag=F,col=col)
corrplot(cormat$r)

# Export
write.csv(cormat, file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//corrmat.csv')
write.csv(cordat, file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//corrdat.csv')





