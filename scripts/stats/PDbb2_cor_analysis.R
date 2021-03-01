## Correlation analysis
library(dplyr)
library(corrplot)
library(Hmisc)

# Define paths
wrkdir <- "X://PD_longrest//groupanalysis"
# wrkdir <- "C://Users//Mikkel//Documents//PDbb2//groupanalysis"
setwd(wrkdir)

col<- colorRampPalette(c("blue", "white", "red"))(20)

# Load data
load('alldata_subj2.Rdata')
load('bbdata2.Rdata')

# Get subject level bb measures
bbsum <- aggregate(cbind(leneve, tueeve, maxeve)~subj, data=bbdata, FUN=median)
alldata <- merge(alldata, bbsum, by="subj")
alldata$U.F45 <- alldata$U.F4+alldata$U.F5

# Re-arrange data
cordat <- select(alldata, nevent.u.m2.min, leneve,tueeve,maxeve,
                 a_intercept,a_slope,beta_pw,beta_cf,alpha_pw,alpha_cf,
                 age,thick, MoCA,,U.F1,U.F2,U.F3,U.F45,U.F6,U.F7,ledd, pd.dur)

# Correlation
# cormat <- cor(cordat, use="pairwise.complete.obs")
cormat <- rcorr(as.matrix(cordat), type="pearson")

corrplot(cormat$r, method="number",p.mat=cormat$P, type="full", sig.level=0.05, insig="blank",diag=F,col=col, cl.pos="n")
corrplot(cormat$r, sig.level=0.05, insig="label_sig",diag=F,col=col)
corrplot(cormat$r, method="shade", shade.col=NA, tl.col="black", tl.srt=45, addCoef.col="black", 
         p.mat=cormat$P, type="full", sig.level=0.05, insig="n",diag=F,col=col, cl.pos="n")

pairs(cordat, pch=".", upper.panel=NULL)

# Export
write.csv(cormat, file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//corrmat2.csv')
write.csv(cordat, file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//corrdat2.csv')





