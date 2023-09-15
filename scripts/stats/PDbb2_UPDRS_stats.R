# UPDRS stats: analysis of sensorimotor signal feature on MDS-UPDRS-III score divided into factors.
#
# Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, 
#   M., Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical 
#   sensorimotor rhythms are uniquely linked to the severity of specific symptoms in 
#   Parkinson's disease [Preprint]. medRxiv.org https://doi.org/10.1101/2021.06.27.21259592
#

library(arm)
library(lmtest)
library(bayestestR)
library(car)
library(brms)
source('/home/mikkel/PD_longrest/scripts/functions/zscore.R')

# Load data
wrkdir <- '/home/mikkel/PD_longrest/groupanalysis'
outdir <- '/home/mikkel/PD_longrest/output'
setwd(wrkdir)
load(file='alldata_subj2.Rdata')
load(file='bbdata2.Rdata')

# Prepare data
bbsum <- aggregate(cbind(leneve, tueeve, maxeve)~subj, data=bbdata, FUN=median)
alldata <- merge(alldata, bbsum, by="subj")

pd.data <- subset(alldata, group == 'patient')
pd.data <- subset(pd.data, !is.na(U.F7))
pd.data <- subset(pd.data, !is.na(U.F6))
pd.data <- subset(pd.data, !is.na(U.F3))
pd.data <- subset(pd.data, !is.na(U.F1))
pd.data$U.F45 <- pd.data$U.F4+pd.data$U.F5

pd.data$log.leneve <- log(pd.data$leneve)
pd.data$log.tueeve <- log(pd.data$tueeve)
pd.data$log.maxeve <- log(pd.data$maxeve)

# Transform data
pd.data$beta_pwz <- zscore(pd.data$beta_pw)
pd.data$alpha_pwz <- zscore(pd.data$alpha_pw)
pd.data$nevent.u.m2.minz <- zscore(pd.data$nevent.u.m2.min)
pd.data$lenevez <- zscore(pd.data$leneve)
pd.data$tueevez <- zscore(pd.data$tueeve)
pd.data$maxevez <- zscore(pd.data$maxeve)
pd.data$a_interceptz <-zscore(pd.data$a_intercept)
pd.data$a_slopez <- zscore(pd.data$a_slope)
pd.data$alpha_cfz <- zscore(pd.data$alpha_cf)
pd.data$beta_cfz <- zscore(pd.data$beta_cf)
pd.data$pd.durz  <- zscore(pd.data$pd.dur)
pd.data$leddz <- zscore(pd.data$ledd)
pd.data$U.F1z  <- zscore(pd.data$U.F1)
pd.data$U.F2z  <- zscore(pd.data$U.F2)
pd.data$U.F3z  <- zscore(pd.data$U.F3)
pd.data$U.F45z <- zscore(pd.data$U.F45)
pd.data$U.F6z  <- zscore(pd.data$U.F6)
pd.data$U.F7z  <- zscore(pd.data$U.F7)
pd.data$pd.durz  <- zscore(pd.data$pd.dur)

######################################################################################
# total MDS-UPDRS III 

Fmod.x <- lm(UPDRSz ~ nevent.u.m2.minz +lenevez + tueevez + maxevez 
             +a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+ 
             age.centerd + sex + thickz + leddz + pd.durz, 
              data=pd.data)

summary(Fmod.x)
as.data.frame(vif(Fmod.x))
qqnorm(resid(F1mod.x))
qqline(resid(F1mod.x))

#Fmod.x <- brm(UPDRSz ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
#                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
#                 age.centerd + sex + thickz + leddz + pd.durz, 
#               data=pd.data, save_pars = save_pars(all = TRUE))

# Sig tests
Fmod.neve <- update(Fmod.x, ~. -nevent.u.m2.minz)
bf.neve <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.neve))
#lrtest(Fmod.neve, Fmod.x)
Fmod.leneve <- update(Fmod.x, ~. -lenevez)
bf.leneve <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.leneve))
#lrtest(Fmod.leneve, Fmod.x)
Fmod.tueeve <- update(Fmod.x, ~. -tueevez)
bf.tueeve <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.tueeve))
#lrtest(Fmod.tueeve, Fmod.x)
Fmod.maxeve <- update(Fmod.x, ~. -maxevez)
bf.maxeve <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.maxeve))
#lrtest(Fmod.maxeve, Fmod.x)
Fmod.intcpt <- update(Fmod.x, ~. -a_interceptz)
bf.intcpt <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.intcpt))
#lrtest(Fmod.intcpt, Fmod.x)
Fmod.slope <- update(Fmod.x, ~. -a_slopez)
bf.slope <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.slope))
#lrtest(Fmod.slope, Fmod.x)
Fmod.alpapw <- update(Fmod.x, ~. -alpha_pwz)
bf.alpapw <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.alpapw))
#lrtest(Fmod.x, Fmod.alpapw)
Fmod.alpacf <- update(Fmod.x, ~. -alpha_cfz)
bf.alpacf<- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.alpacf))
#lrtest(Fmod.x, Fmod.alpacf)
Fmod.betapw <- update(Fmod.x, ~. -beta_pwz)
bf.betapw <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.betapw))
#lrtest(Fmod.x, Fmod.betapw)
Fmod.betacf <- update(Fmod.x, ~. -beta_cfz)
bf.betacf <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.betacf))
#lrtest(Fmod.x, Fmod.betacf)
Fmod.age <- update(Fmod.x, ~. -age.centerd)
bf.age <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.age))
#lrtest(Fmod.x, Fmod.age)
Fmod.sex <- update(Fmod.x, ~. -sex)
bf.sex <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.sex))
#lrtest(Fmod.x, Fmod.sex)
Fmod.thick <- update(Fmod.x, ~. -thickz)
bf.thick <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.thick))
#lrtest(Fmod.x, Fmod.thick)
Fmod.dur <- update(Fmod.x, ~. -pd.durz)
bf.dur<- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.dur))
#lrtest(Fmod.x, Fmod.dur)
Fmod.ledd <- update(Fmod.x, ~. -leddz)
bf.ledd <- as.numeric(bayesfactor_models(Fmod.x, denominator = Fmod.ledd))
#lrtest(Fmod.x, Fmod.ledd)

# SUMMARY
#mod.sim <- sim(Fmod.x, n.sims=1000)
#x1 <- coef(Fmod.x)
#x2 <- t(apply(coef(mod.sim), 2, quantile, c(0.025, 0.975)))
#sums <- cbind(x1, x2)
#sums
######################################################################################
# MDS-UPDRS-III Factor 1
F1mod.x <- lm(U.F1z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                age.centerd + sex + thickz+ leddz + pd.durz, 
              data=pd.data)

#F1mod.x <- brm(U.F1z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
#                a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
#                age.centerd + sex + thickz + leddz + pd.durz, 
#              data=pd.data, save_pars = save_pars(all = TRUE))

qqnorm(resid(F1mod.x))
qqline(resid(F1mod.x))

summary(F1mod.x)

# Sig tests
F1mod.neve <- update(F1mod.x, ~. -nevent.u.m2.minz)
lrtest(F1mod.neve, F1mod.x)
bayesfactor_models(F1mod.x, denominator = F1mod.neve)

F1mod.leneve <- update(F1mod.x, ~. -lenevez)
bayesfactor_models(F1mod.x, denominator = F1mod.leneve)
lrtest(F1mod.leneve, F1mod.x)

F1mod.tueeve <- update(F1mod.x, ~. -tueevez)
bayesfactor_models(F1mod.x, denominator = F1mod.tueeve)
lrtest(F1mod.tueeve, F1mod.x)

F1mod.maxeve <- update(F1mod.x, ~. -maxevez)
bayesfactor_models(F1mod.x, denominator = F1mod.maxeve)
lrtest(F1mod.maxeve, F1mod.x)

F1mod.intcpt <- update(F1mod.x, ~. -a_interceptz)
bayesfactor_models(F1mod.x, denominator = F1mod.intcpt)
lrtest(F1mod.intcpt, F1mod.x)

F1mod.slope <- update(F1mod.x, ~. -a_slopez)
bayesfactor_models(F1mod.x, denominator = F1mod.slope)
lrtest(F1mod.slope, F1mod.x)

F1mod.alpapw <- update(F1mod.x, ~. -alpha_pwz)
bayesfactor_models(F1mod.x, denominator = F1mod.alpapw)
lrtest(F1mod.x, F1mod.alpapw)

F1mod.alpacf <- update(F1mod.x, ~. -alpha_cfz)
bayesfactor_models(F1mod.x, denominator = F1mod.alpacf)
lrtest(F1mod.x, F1mod.alpacf)

F1mod.betapw <- update(F1mod.x, ~. -beta_pwz)
bayesfactor_models(F1mod.x, denominator = F1mod.betapw)
lrtest(F1mod.x, F1mod.betapw)

F1mod.betacf <- update(F1mod.x, ~. -beta_cfz)
bayesfactor_models(F1mod.x, denominator = F1mod.betacf)
lrtest(F1mod.x, F1mod.betacf)

F1mod.age <- update(F1mod.x, ~. -age.centerd)
bayesfactor_models(F1mod.x, denominator = F1mod.age)
lrtest(F1mod.x, F1mod.age)

F1mod.sex <- update(F1mod.x, ~. -sex)
bayesfactor_models(F1mod.x, denominator = F1mod.sex)
lrtest(F1mod.x, F1mod.sex)

F1mod.thick <- update(F1mod.x, ~. -thickz)
bayesfactor_models(F1mod.x, denominator = F1mod.thick)
lrtest(F1mod.x, F1mod.thick)

F1mod.dur <- update(F1mod.x, ~. -pd.durz)
bf.dur<- as.numeric(bayesfactor_models(F1mod.x, denominator = F1mod.dur))
#lrtest(Fmod.x, Fmod.dur)
F1mod.ledd <- update(F1mod.x, ~. -leddz)
bf.ledd <- as.numeric(bayesfactor_models(F1mod.x, denominator = F1mod.ledd))

# SUMMARY
mod.sim <- sim(F1mod.x, n.sims=1000)
x1 <- coef(F1mod.x)
x2 <- t(apply(coef(mod.sim), 2, quantile, c(0.025, 0.975)))
sums1 <- cbind(x1, x2)
sums1

######################################################################################
# MDS-UPDRS-III Factor 2
F2mod.x <- lm(U.F2z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                age.centerd + sex + thickz+ leddz + pd.durz, 
              data=pd.data)

qqnorm(resid(F2mod.x))
qqline(resid(F2mod.x))

F2mod.neve <- update(F2mod.x, ~. -nevent.u.m2.minz)
bayesfactor_models(F2mod.x, denominator = F2mod.neve)
lrtest(F2mod.x, F2mod.neve)

F2mod.leneve <- update(F2mod.x, ~. -lenevez)
bayesfactor_models(F2mod.x, denominator = F2mod.leneve)
lrtest(F2mod.x, F2mod.leneve)

F2mod.tueeve <- update(F2mod.x, ~. -tueevez)
bayesfactor_models(F2mod.x, denominator = F2mod.tueeve)
lrtest(F2mod.x, F2mod.tueeve)

F2mod.maxeve <- update(F2mod.x, ~. -maxevez)
bayesfactor_models(F2mod.x, denominator = F2mod.maxeve)
lrtest(F2mod.x, F2mod.maxeve)

F2mod.intcpt <- update(F2mod.x, ~. -a_interceptz)
bayesfactor_models(F2mod.x, denominator = F2mod.intcpt)
lrtest(F2mod.x, F2mod.intcpt)

F2mod.slope <- update(F2mod.x, ~. -a_slopez)
bayesfactor_models(F2mod.x, denominator = F2mod.slope)
lrtest(F2mod.x, F2mod.slope)

F2mod.alpapw <- update(F2mod.x, ~. -alpha_pwz)
bayesfactor_models(F2mod.x, denominator = F2mod.alpapw)
lrtest(F2mod.x, F2mod.alpapw)

F2mod.alpacf <- update(F2mod.x, ~. -alpha_cfz)
bayesfactor_models(F2mod.x, denominator = F2mod.alpacf)
lrtest(F2mod.x, F2mod.alpacf)

F2mod.betapw <- update(F2mod.x, ~. -beta_pwz)
bayesfactor_models(F2mod.x, denominator = F2mod.betapw)
lrtest(F2mod.x, F2mod.betapw)

F2mod.betacf <- update(F2mod.x, ~. -beta_cfz)
bayesfactor_models(F2mod.x, denominator = F2mod.betacf)
lrtest(F2mod.x, F2mod.betacf)

F2mod.age <- update(F2mod.x, ~. -age.centerd)
bayesfactor_models(F2mod.x, denominator = F2mod.age)
lrtest(F2mod.x, F2mod.age)

F2mod.sex <- update(F2mod.x, ~. -sex)
bayesfactor_models(F2mod.x, denominator = F2mod.sex)
lrtest(F2mod.x, F2mod.sex)

F2mod.thick <- update(F2mod.x, ~. -thickz)
bayesfactor_models(F2mod.x, denominator = F2mod.thick)
lrtest(F2mod.x, F2mod.thick)

F2mod.dur <- update(F2mod.x, ~. -pd.durz)
bf.dur<- as.numeric(bayesfactor_models(F2mod.x, denominator = F2mod.dur))
#lrtest(Fmod.x, Fmod.dur)
F1mod.ledd <- update(F1mod.x, ~. -leddz)
bf.ledd <- as.numeric(bayesfactor_models(F2mod.x, denominator = F2mod.ledd))

# SUMMARY
mod2.sim <- sim(F2mod.x, n.sims=1000)
x1 <- coef(F2mod.x)
x2 <- t(apply(coef(mod2.sim), 2, quantile, c(0.025, 0.975)))
sums2 <- cbind(x1, x2)
sums2

######################################################################################
# MDS-UPDRS-III Factor 3
F3mod.x <- lm(U.F3z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                age.centerd + sex + thickz+ leddz + pd.durz, 
              data=pd.data)

qqnorm(resid(F3mod.x))
qqline(resid(F3mod.x))

F3mod.neve <- update(F3mod.x, .~. -nevent.u.m2.minz)
bayesfactor_models(F3mod.x, denominator = F3mod.neve)
lrtest(F3mod.x, F3mod.neve)

F3mod.leneve <- update(F3mod.x, ~. -lenevez)
bayesfactor_models(F3mod.x, denominator = F3mod.leneve)
lrtest(F3mod.x, F3mod.leneve)

F3mod.tueeve <- update(F3mod.x, ~. -tueevez)
bayesfactor_models(F3mod.x, denominator = F3mod.tueeve)
lrtest(F3mod.x, F3mod.tueeve)

F3mod.maxeve <- update(F3mod.x, ~. -maxevez)
bayesfactor_models(F3mod.x, denominator = F3mod.maxeve)
lrtest(F3mod.x, F3mod.maxeve)

F3mod.intcpt <- update(F3mod.x, ~. -a_interceptz)
bayesfactor_models(F3mod.x, denominator = F3mod.intcpt)
lrtest(F3mod.x, F3mod.intcpt)

F3mod.slope <- update(F3mod.x, ~. -a_slopez)
bayesfactor_models(F3mod.x, denominator = F3mod.slope)
lrtest(F3mod.x, F3mod.slope)

F3mod.alpapw <- update(F3mod.x, ~. -alpha_pwz)
bayesfactor_models(F3mod.x, denominator = F3mod.alpapw)
lrtest(F3mod.x, F3mod.alpapw)

F3mod.alpacf <- update(F3mod.x, ~. -alpha_cfz)
bayesfactor_models(F3mod.x, denominator = F3mod.alpacf)
lrtest(F3mod.x, F3mod.alpacf)

F3mod.betapw <- update(F3mod.x, ~. -beta_pwz)
bayesfactor_models(F3mod.x, denominator = F3mod.betapw)
lrtest(F3mod.x, F3mod.betapw)

F3mod.betacf <- update(F3mod.x, ~. -beta_cfz)
bayesfactor_models(F3mod.x, denominator = F3mod.betacf)
lrtest(F3mod.x, F3mod.betacf)

F3mod.age <- update(F3mod.x, ~. -age.centerd)
bayesfactor_models(F3mod.x, denominator = F3mod.age)
lrtest(F3mod.x, F3mod.age)

F3mod.sex <- update(F3mod.x, ~. -sex)
bayesfactor_models(F3mod.x, denominator = F3mod.sex)
lrtest(F3mod.x, F3mod.sex)

F3mod.thick <- update(F3mod.x, ~. -thickz)
bayesfactor_models(F3mod.x, denominator = F3mod.thick)
lrtest(F3mod.x, F3mod.thick)

# SUMMARY
mod3.sim <- sim(F3mod.x, n.sims=1000)
x1 <- summary(F3mod.x)$coefficients
x2 <- t(apply(coef(mod3.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1[,1], x2)
sums

######################################################################################
## MDS-UPDRS-III Factor 4+5 combined
# ML model
F45mod.x <- lm(U.F45z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                 age.centerd + sex + thickz  + leddz + pd.durz, 
              data=pd.data)

qqnorm(resid(F45mod.x))
qqline(resid(F45mod.x))

# brms model
F45mod.x <- brm(U.F45z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                 age.centerd + sex + thickz + pd.durz, 
               data=pd.data, save_pars = save_pars(all = TRUE), iter=20000, cores=4)

summary(F45mod.x)

F45mod.neve <- update(F45mod.x, .~. -nevent.u.m2.minz, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.neve)
lrtest(F45mod.x, F45mod.neve)

F45mod.leneve <- update(F45mod.x, ~. -lenevez, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.leneve)
lrtest(F45mod.x, F45mod.leneve)

F45mod.tueeve <- update(F45mod.x, ~. -tueevez, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.tueeve)
lrtest(F45mod.x, F45mod.tueeve)

F45mod.maxeve <- update(F45mod.x, ~. -maxevez, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.maxeve)
lrtest(F45mod.x, F45mod.maxeve)

F45mod.intcpt <- update(F45mod.x, ~. -a_interceptz, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.intcpt)
lrtest(F45mod.x, F45mod.intcpt)

F45mod.slope <- update(F45mod.x, ~. -a_slopez, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.slope)
lrtest(F45mod.x, F45mod.slope)

F45mod.alpapw <- update(F45mod.x, ~. -alpha_pwz, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.alpapw)
lrtest(F45mod.x, F45mod.alpapw)

F45mod.alpacf <- update(F45mod.x, ~. -alpha_cfz, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.alpacf)
lrtest(F45mod.x, F45mod.alpacf)

F45mod.betapw <- update(F45mod.x, ~. -beta_pwz, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.betapw)
lrtest(F45mod.x, F45mod.betapw)

F45mod.betacf <- update(F45mod.x, ~. -beta_cfz, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.betacf)
lrtest(F45mod.x, F45mod.betacf)

F45mod.age <- update(F45mod.x, ~. -age.centerd, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.age)
lrtest(F45mod.x, F45mod.age)

F45mod.sex <- update(F45mod.x, ~. -sex, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.sex)
lrtest(F45mod.x, F45mod.sex)

F45mod.thick <- update(F45mod.x, ~. -thickz, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.thick)
lrtest(F45mod.x, F45mod.thick)

F45mod.dur <- update(F45mod.x, ~. -pd.dur, cores=4)
bayesfactor_models(F45mod.x, denominator = F45mod.dur)
lrtest(F45mod.x, F45mod.thick)

# SUMMARY
mod45.sim <- sim(F45mod.x, n.sims=1000)
x1 <- coef(F45mod.x)
x2 <- t(apply(coef(mod45.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1, x2)
round(sums, digits=3)

######################################################################################
## MDS-UPDRS-III Factor 4 (RIGHT bradykinesia)

F4mod.x <- lm(U.F4z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                 age.centerd + sex + thickz  + pd.durz + leddz, 
               data=pd.data)
summary(F4mod.x)
qqnorm(resid(F4mod.x))
qqline(resid(F4mod.x))

#F4mod.x <- brm(U.F4z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
#                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
#                 age.centerd + sex + thickz+ leddz + pd.durz, 
#               data=pd.data, save_pars = save_pars(all = TRUE))

F4mod.neve <- update(F4mod.x, .~. -nevent.u.m2.minz)
bayesfactor_models(F4mod.x, denominator = F4mod.neve)
lrtest(F4mod.x, F4mod.neve)

F4mod.leneve <- update(F4mod.x, ~. -lenevez)
bayesfactor_models(F4mod.x, denominator = F4mod.leneve)
lrtest(F4mod.x, F4mod.leneve)

F4mod.tueeve <- update(F4mod.x, ~. -tueevez, cores=4)
bayesfactor_models(F4mod.x, denominator = F4mod.tueeve)
lrtest(F4mod.x, F4mod.tueeve)

F4mod.maxeve <- update(F4mod.x, ~. -maxevez)
bayesfactor_models(F4mod.x, denominator = F4mod.maxeve)
lrtest(F4mod.x, F4mod.maxeve)

F4mod.intcpt <- update(F4mod.x, ~. -a_interceptz)
bayesfactor_models(F4mod.x, denominator = F4mod.intcpt)
lrtest(F4mod.x, F4mod.intcpt)

F4mod.slope <- update(F4mod.x, ~. -a_slopez)
bayesfactor_models(F4mod.x, denominator = F4mod.slope)
lrtest(F4mod.x, F4mod.slope)

F4mod.alpapw <- update(F4mod.x, ~. -alpha_pwz)
bayesfactor_models(F4mod.x, denominator = F4mod.alpapw)
lrtest(F4mod.x, F4mod.alpapw)

F4mod.alpacf <- update(F4mod.x, ~. -alpha_cfz)
bayesfactor_models(F4mod.x, denominator = F4mod.alpacf)
lrtest(F4mod.x, F4mod.alpacf)

F4mod.betapw <- update(F4mod.x, ~. -beta_pwz)
bayesfactor_models(F4mod.x, denominator = F4mod.betapw)
lrtest(F4mod.x, F4mod.betapw)

F4mod.betacf <- update(F4mod.x, ~. -beta_cfz)
bayesfactor_models(F4mod.x, denominator = F4mod.betacf)
lrtest(F4mod.x, F4mod.betacf)

F4mod.age <- update(F4mod.x, ~. -age.centerd)
bayesfactor_models(F4mod.x, denominator = F4mod.age)
lrtest(F4mod.x, F4mod.age)

F4mod.sex <- update(F4mod.x, ~. -sex, cores=4)
bayesfactor_models(F4mod.x, denominator = F4mod.sex)
lrtest(F4mod.x, F4mod.sex)

F4mod.thick <- update(F4mod.x, ~. -thickz, cores=4)
bayesfactor_models(F4mod.x, denominator = F4mod.thick)
lrtest(F4mod.x, F4mod.thick)

F4mod.dur <- update(F4mod.x, ~. -pd.durz)
bayesfactor_models(F4mod.x, denominator = F4mod.dur)
lrtest(F4mod.x, F4mod.thick)

F4mod.ledd <- update(F4mod.x, ~. -leddz)
bayesfactor_models(F4mod.x, denominator = F4mod.ledd)
lrtest(F4mod.x, F4mod.ledd)

# SUMMARY
mod4.sim <- sim(F4mod.x, n.sims=1000)
x1 <- coef(F4mod.x)
x2 <- t(apply(coef(mod4.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1, x2)
round(sums, digits=3)

######################################################################################
## MDS-UPDRS-III Factor 5 (LEFT bradykinesia)

F5mod.x <- lm(U.F5z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                 age.centerd + sex + thickz + leddz + pd.durz , 
               data=pd.data)

#summary(F5mod.x)
#as.data.frame(vif(F5mod.x))
#qqnorm(resid(F5mod.x))
#qqline(resid(F5mod.x))

#F5mod.x <- brm(U.F5z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
#                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
#                 age.centerd + sex + thickz+ leddz + pd.durz, 
#               data=pd.data, iter=20000, save_pars = save_pars(all = TRUE))

F5mod.neve <- update(F5mod.x, .~. -nevent.u.m2.minz)
bf.neve5 <- as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.neve))
#lrtest(F5mod.x, F5mod.neve)

F5mod.leneve <- update(F5mod.x, ~. -lenevez)
bf.leneve5 <- as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.leneve))
#lrtest(F5mod.x, F5mod.leneve)

F5mod.tueeve <- update(F5mod.x, ~. -tueevez)
bf.tueve5 <- as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.tueeve))
#lrtest(F5mod.x, F5mod.tueeve)

F5mod.maxeve <- update(F5mod.x, ~. -maxevez)
bf.maxeve5 <-as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.maxeve))
#lrtest(F5mod.x, F5mod.maxeve)

F5mod.intcpt <- update(F5mod.x, ~. -a_interceptz)
bf.intcpt5 <-as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.intcpt))
#lrtest(F5mod.x, F5mod.intcpt)

F5mod.slope <- update(F5mod.x, ~. -a_slopez)
bf.slope5 <-as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.slope))
#lrtest(F5mod.x, F5mod.slope)

F5mod.alpapw <- update(F5mod.x, ~. -alpha_pwz)
bf.alpapw5 <-as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.alpapw))
#lrtest(F5mod.x, F5mod.alpapw)

F5mod.alpacf <- update(F5mod.x, ~. -alpha_cfz)
bf.alpacf5 <-as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.alpacf))
#lrtest(F5mod.x, F5mod.alpacf)

F5mod.betapw <- update(F5mod.x, ~. -beta_pwz)
bf.betapw5 <-as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.betapw))
#lrtest(F5mod.x, F5mod.betapw)

F5mod.betacf <- update(F5mod.x, ~. -beta_cfz)
bf.betacf5 <-as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.betacf))
#lrtest(F5mod.x, F5mod.betacf)

F5mod.age <- update(F5mod.x, ~. -age.centerd)
bf.age5 <-as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.age))
#lrtest(F5mod.x, F5mod.age)

F5mod.sex <- update(F5mod.x, ~. -sex)
bf.sex5 <- as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.sex))
#lrtest(F5mod.x, F5mod.sex)

F5mod.thick <- update(F5mod.x, ~. -thickz)
bf.thick5 <- as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.thick))
#lrtest(F5mod.x, F5mod.thick)

F5mod.dur <- update(F5mod.x, ~. -pd.durz)
bf.dur5 <- as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.dur))
#lrtest(F5mod.x, F5mod.dur)

F5mod.ledd <- update(F5mod.x, ~. -leddz)
bf.ledd5 <- as.numeric(bayesfactor_models(F5mod.x, denominator = F5mod.ledd))
#lrtest(F5mod.x, F5mod.ledd)

# SUMMARY
#mod5.sim <- sim(F5mod.x, n.sims=1000)
#x1 <- coef(F5mod.x)
#x2 <- t(apply(coef(mod5.sim), 2, quantile, c(0.025, 0.975)))
#sums <- cbind(x1, x2)
#round(sums, digits=3)
######################################################################################
# MDS-UPDRS-III Factor 6
F6mod.x <- lm(U.F6z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                age.centerd + sex + thickz + leddz + pd.durz, 
              data=pd.data)

qqnorm(resid(F6mod.x))
qqline(resid(F6mod.x))

F6mod.neve <- update(F6mod.x, .~. -nevent.u.m2.minz)
bayesfactor_models(F6mod.x, denominator = F6mod.neve)
lrtest(F6mod.x, F6mod.neve)

F6mod.leneve <- update(F6mod.x, ~. -lenevez)
bayesfactor_models(F6mod.x, denominator = F6mod.leneve)
lrtest(F6mod.x, F6mod.leneve)

F6mod.tueeve <- update(F6mod.x, ~. -tueevez)
bayesfactor_models(F6mod.x, denominator = F6mod.tueeve)
lrtest(F6mod.x, F6mod.tueeve)

F6mod.maxeve <- update(F6mod.x, ~. -maxevez)
bayesfactor_models(F6mod.x, denominator = F6mod.maxeve)
lrtest(F6mod.x, F6mod.maxeve)

F6mod.intcpt <- update(F6mod.x, ~. -a_interceptz)
bayesfactor_models(F6mod.x, denominator = F6mod.intcpt)
lrtest(F6mod.x, F6mod.intcpt)

F6mod.slope <- update(F6mod.x, ~. -a_slopez)
bayesfactor_models(F6mod.x, denominator = F6mod.slope)
lrtest(F6mod.x, F6mod.slope)

F6mod.alpapw <- update(F6mod.x, ~. -alpha_pwz)
bayesfactor_models(F6mod.x, denominator = F6mod.alpapw)
lrtest(F6mod.x, F6mod.alpapw)

F6mod.alpacf <- update(F6mod.x, ~. -alpha_cfz)
bayesfactor_models(F6mod.x, denominator = F6mod.alpacf)
lrtest(F6mod.x, F6mod.alpacf)

F6mod.betapw <- update(F6mod.x, ~. -beta_pwz)
bayesfactor_models(F6mod.x, denominator = F6mod.betapw)
lrtest(F6mod.x, F6mod.betapw)

F6mod.betacf <- update(F6mod.x, ~. -beta_cfz)
bayesfactor_models(F6mod.x, denominator = F6mod.betacf)
lrtest(F6mod.x, F6mod.betacf)

F6mod.age <- update(F6mod.x, ~. -age.centerd)
bayesfactor_models(F6mod.x, denominator = F6mod.age)
lrtest(F6mod.x, F6mod.age)

F6mod.sex <- update(F6mod.x, ~. -sex)
bayesfactor_models(F6mod.x, denominator = F6mod.sex)
lrtest(F6mod.x, F6mod.sex)

F6mod.thick <- update(F6mod.x, ~. -thickz)
bayesfactor_models(F6mod.x, denominator = F6mod.thick)
lrtest(F6mod.x, F6mod.thick)

# SUMMARY
mod6.sim <- sim(F6mod.x, n.sims=1000)
x1 <- coef(F6mod.x)
x2 <- t(apply(coef(mod6.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1, x2)
sums

######################################################################################
# MDS-UPDRS-III Factor 7
F7mod.x <- lm(U.F7z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                 age.centerd + sex + thickz+ leddz + pd.durz, 
               data=pd.data)

qqnorm(resid(F7mod.x))
qqline(resid(F7mod.x))

F7mod.neve <- update(F7mod.x, .~. -nevent.u.m2.minz)
bayesfactor_models(F7mod.x, denominator = F7mod.neve)
lrtest(F7mod.x, F7mod.neve)

F7mod.leneve <- update(F7mod.x, ~. -lenevez)
bayesfactor_models(F7mod.x, denominator = F7mod.leneve)
lrtest(F7mod.x, F7mod.leneve)

F7mod.tueeve <- update(F7mod.x, ~. -tueevez)
bayesfactor_models(F7mod.x, denominator = F7mod.tueeve)
lrtest(F7mod.x, F7mod.tueeve)

F7mod.maxeve <- update(F7mod.x, ~. -maxevez)
bayesfactor_models(F7mod.x, denominator = F7mod.maxeve)
lrtest(F7mod.x, F7mod.maxeve)

F7mod.intcpt <- update(F7mod.x, ~. -a_interceptz)
bayesfactor_models(F7mod.x, denominator = F7mod.intcpt)
lrtest(F7mod.x, F7mod.intcpt)

F7mod.slope <- update(F7mod.x, ~. -a_slopez)
bayesfactor_models(F7mod.x, denominator = F7mod.slope)
lrtest(F7mod.x, F7mod.slope)

F7mod.alpapw <- update(F7mod.x, ~. -alpha_pwz)
bayesfactor_models(F7mod.x, denominator = F7mod.alpapw)
lrtest(F7mod.x, F7mod.alpapw)

F7mod.alpacf <- update(F7mod.x, ~. -alpha_cfz)
bayesfactor_models(F7mod.x, denominator = F7mod.alpacf)
lrtest(F7mod.x, F7mod.alpacf)

F7mod.betapw <- update(F7mod.x, ~. -beta_pwz)
bayesfactor_models(F7mod.x, denominator = F7mod.betapw)
lrtest(F7mod.x, F7mod.betapw)

F7mod.betacf <- update(F7mod.x, ~. -beta_cfz)
bayesfactor_models(F7mod.x, denominator = F7mod.betacf)
lrtest(F7mod.x, F7mod.betacf)

F7mod.age <- update(F7mod.x, ~. -age.centerd)
bayesfactor_models(F7mod.x, denominator = F7mod.age)
lrtest(F7mod.x, F7mod.age)

F7mod.sex <- update(F7mod.x, ~. -sex)
bayesfactor_models(F7mod.x, denominator = F7mod.sex)
lrtest(F7mod.x, F7mod.sex)

F7mod.thick <- update(F7mod.x, ~. -thickz)
bayesfactor_models(F7mod.x, denominator = F7mod.thick)
lrtest(F7mod.x, F7mod.thick)

# SUMMARY
mod7.sim <- sim(F6mod.x, n.sims=1000)
x1 <- coef(F7mod.x)
x2 <- t(apply(coef(mod7.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1, x2)
sums

# SAVE
setwd(outdir)
save("F1mod.x","F2mod.x","F3mod.x","F45mod.x","F6mod.x","F7mod.x", file="updrsmods.RData")

#END
