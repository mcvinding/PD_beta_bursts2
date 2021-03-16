# PD beta burst statistics: analysis of N events
### CLEAN UP!!!! ###
library(lme4)
library(arm)
library(ggplot2)
library(car)
source('X://PD_longrest//scripts//functions//zscore.R')

## Load data
setwd('X://PD_longrest//groupanalysis//')
load('X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
# load('C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj2.Rdata')

# ## Center and standardize variables
# alldata$age.centerd <- alldata$age-mean(alldata$age)
# alldata$ageZ <- zscore(alldata$age)
# alldata$thickZ <- zscore(alldata$thick)

# Inspect hist
ggplot( aes(x=nevent.u.m2.min, fill=group), data=alldata) +
  geom_histogram(color="black", alpha=0.6, position = 'identity', bins=25)

# Inspect ~age
ggplot(aes(x=age, y=nevent.u.m2.min, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_smooth(method=lm)

# ######################################################################################
# # LMER regression model
mod.neve.Full3 <- glm(nevent.u.m2.min ~ (group+age.centerd+sex+thickz)^3, data=alldata, family=poisson)
anova(mod.neve.Full3, test="Chisq")

# Model summary
mod.sim <- sim(mod.neve.Full3, n.sims=1000)
cf <- coef(mod.sim)
x1 <- summary(mod.neve.Full3)$coefficients
x2 <- t(apply(cf, 2, quantile, c(0.025, 0.975)))
cbind(x1[,1], x2)

## Difference
# Group
c(exp(x1[2])*100-100,
  quantile(exp(cf[,2])*100-100, c(0.025, 0.975)))

# Ctrl x age : female.
c(exp(x1[3])*100-100,
  quantile(exp(cf[,3])*100-100, c(0.025, 0.975)))

# Ctrl x age : male.
c(exp(x1[3]+x1[9])*100-100,
  quantile(exp(cf[,3]+cf[,9])*100-100, c(0.025, 0.975)))

# Ptns x age (female).
c(exp(x1[3]+x1[6])*100-100,
  quantile(exp(cf[,3]+cf[,6])*100-100, c(0.025, 0.975)))

# Ptns x age (male).
c(exp(x1[3]+x1[6]+x1[9]+x1[12])*100-100,
  quantile(exp(cf[,3]+cf[,6]+cf[,9]+cf[,12])*100-100, c(0.025, 0.975)))

# Thick (female)
c(exp(x1[5])*100-100,
  quantile(exp(cf[,5])*100-100, c(0.025, 0.975)))

# Thick (male)
c(exp(+x1[5]+x1[11])*100-100,
  quantile(exp(cf[,5]+cf[,11])*100-100, c(0.025, 0.975)))

# Age x Thick (female)
c(exp(x1[3]+x1[5]+x1[10])*100-100,
  quantile(exp(cf[,3]+cf[,5]+cf[,10])*100-100, c(0.025, 0.975)))


## Plot scatterplot amd model
agespan <- seq(min(alldata$age.centerd),max(alldata$age.centerd), 0.1)
agespan2 <- seq(min(alldata$age),max(alldata$age), 0.1)

new.dat <- data.frame(age.centerd=rep(agespan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thickZ=0,
                      age=rep(agespan2, 4))

new.dat$pred  <- exp(predict(mod.neve.Full3, new.dat, re.form=NA))
ggplot(aes(x=age, y=nevent.u.m2.min, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)



thickspan <- seq(min(alldata$thickZ),max(alldata$thickZ), 0.01)
thickspan2 <- seq(min(alldata$thick),max(alldata$thick), length.out = length(thickspan))

new.dat <- data.frame(thickZ=rep(thickspan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(thickspan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(thickspan)), 2)),
                      age.centerd=0,
                      thick=rep(thickspan2, 4))

new.dat$pred <- exp(predict(mod.neve.Full3, new.dat, re.form=NA))

ggplot(aes(x=thickZ, y=nevent.u.m2.min, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1, data=new.dat)

#END