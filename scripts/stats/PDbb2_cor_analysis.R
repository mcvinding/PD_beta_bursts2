## Correlation analysis: all-to-all correlation between variables
#
# Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., 
#  Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor 
#  rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease [Preprint]. 
#  medRxiv.org https://doi.org/10.1101/2021.06.27.21259592
#

library(dplyr)
library(corrplot)
library(Hmisc)
library(BayesFactor)

# Define paths
wrkdir <- "X://PD_longrest//groupanalysis"
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

write.csv(cordat, file='X://PD_longrest//groupanalysis//cordat.csv')

# Correlation
correlationBF(cordat$nevent.u.m2.min, cordat$U.F45, rscale="medium") #, nullInterval=c(-0.05,0.05))
tt <- correlationBF(cordat$nevent.u.m2.min, cordat$U.F45, rscale="medium", posterior=T, iter=10000)
plot((tt[,"rho"]))
mean((tt[,'rho'])[5001:10000]>0)*2

correlationBF(cordat$nevent.u.m2.min, cordat$U.F45, rscale="medium.narrow") #, nullInterval=c(-0.05,0.05))
tt <- correlationBF(cordat$nevent.u.m2.min, cordat$U.F45, rscale="medium", posterior=T, iter=10000)
plot((tt[,"rho"]))
mean((tt[,'rho'])[5001:10000]>0)

bayes.cor.test(cordat$nevent.u.m2.min, cordat$U.F45)
bayes.cor.test(cordat$beta_pw, cordat$U.F45)

cormat <- rcorr(as.matrix(cordat), type="pearson")

# Export
write.csv(cormat$r, file='X://PD_longrest//groupanalysis//corrmat2r.csv')
write.csv(cormat$P, file='X://PD_longrest//groupanalysis//corrmat2p.csv')

#END