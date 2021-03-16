# UPDRS stats
library(arm)
library(lmtest)
source('X://PD_longrest//scripts//functions//zscore.R')

# Load data
wrkdir <- "X://PD_longrest//groupanalysis"
setwd(wrkdir)
load(file='X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
load(file='X://PD_longrest//groupanalysis//bbdata2.Rdata')

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
# pd.data$leneve.ms <- pd.data$leneve*1000
# pd.data$tueeve.ms <- pd.data$tueeve*1000
pd.data$nevent.u.m2.minz <- zscore(pd.data$nevent.u.m2.min)
pd.data$lenevez <- zscore(pd.data$leneve)
pd.data$tueevez <- zscore(pd.data$tueeve)
pd.data$maxevez <- zscore(pd.data$maxeve)
pd.data$a_interceptz <-zscore(pd.data$a_intercept)
pd.data$a_slopez <- zscore(pd.data$a_slope)
pd.data$alpha_cfz <- zscore(pd.data$alpha_cf)
pd.data$beta_cfz <- zscore(pd.data$beta_cf)
pd.data$U.F1z  <- zscore(pd.data$U.F1)
pd.data$U.F2z  <- zscore(pd.data$U.F2)
pd.data$U.F3z  <- zscore(pd.data$U.F3)
pd.data$U.F45z <- zscore(pd.data$U.F45)
pd.data$U.F6z  <- zscore(pd.data$U.F6)
pd.data$U.F7z  <- zscore(pd.data$U.F7)

######################################################################################
# MDS-UPDRS-III Factor 1
F1mod.x <- lm(U.F1z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                 age.centerd + sex + thickz, 
               data=pd.data)

qqnorm(resid(F1mod.x))
qqline(resid(F1mod.x))

summary(F1mod.x)

# Sig tests
F1mod.neve <- update(F1mod.x, ~. -nevent.u.m2.minz)
lrtest(F1mod.neve, F1mod.x)

F1mod.leneve <- update(F1mod.x, ~. -lenevez)
lrtest(F1mod.leneve, F1mod.x)

F1mod.tueeve <- update(F1mod.x, ~. -tueevez)
lrtest(F1mod.tueeve, F1mod.x)

F1mod.maxeve <- update(F1mod.x, ~. -maxevez)
lrtest(F1mod.maxeve, F1mod.x)

F1mod.intcpt <- update(F1mod.x, ~. -a_interceptz)
lrtest(F1mod.intcpt, F1mod.x)

F1mod.slope <- update(F1mod.x, ~. -a_slopez)
lrtest(F1mod.slope, F1mod.x)

F1mod.alpapw <- update(F1mod.x, ~. -alpha_pwz)
lrtest(F1mod.x, F1mod.alpapw)

F1mod.alpacf <- update(F1mod.x, ~. -alpha_cfz)
lrtest(F1mod.x, F1mod.alpacf)

F1mod.betapw <- update(F1mod.x, ~. -beta_pwz)
lrtest(F1mod.x, F1mod.betapw)

F1mod.betacf <- update(F1mod.x, ~. -beta_cfz)
lrtest(F1mod.x, F1mod.betacf)

F1mod.age <- update(F1mod.x, ~. -age.centerd)
lrtest(F1mod.x, F1mod.age)

F1mod.sex <- update(F1mod.x, ~. -sex)
lrtest(F1mod.x, F1mod.sex)

F1mod.thick <- update(F1mod.x, ~. -thickz)
lrtest(F1mod.x, F1mod.thick)

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
                age.centerd + sex + thickz, 
              data=pd.data)

qqnorm(resid(F2mod.x))
qqline(resid(F2mod.x))

F2mod.neve <- update(F2mod.x, ~. -nevent.u.m2.minz)
lrtest(F2mod.x, F2mod.neve)

F2mod.leneve <- update(F2mod.x, ~. -lenevez)
lrtest(F2mod.x, F2mod.leneve)

F2mod.tueeve <- update(F2mod.x, ~. -tueevez)
lrtest(F2mod.x, F2mod.tueeve)

F2mod.maxeve <- update(F2mod.x, ~. -maxevez)
lrtest(F2mod.x, F2mod.maxeve)

F2mod.intcpt <- update(F2mod.x, ~. -a_interceptz)
lrtest(F2mod.x, F2mod.intcpt)

F2mod.slope <- update(F2mod.x, ~. -a_slopez)
lrtest(F2mod.x, F2mod.slope)

F2mod.alpapw <- update(F2mod.x, ~. -alpha_pwz)
lrtest(F2mod.x, F2mod.alpapw)

F2mod.alpacf <- update(F2mod.x, ~. -alpha_cfz)
lrtest(F2mod.x, F2mod.alpacf)

F2mod.betapw <- update(F2mod.x, ~. -beta_pwz)
lrtest(F2mod.x, F2mod.betapw)

F2mod.betacf <- update(F2mod.x, ~. -beta_cfz)
lrtest(F2mod.x, F2mod.betacf)

F2mod.age <- update(F2mod.x, ~. -age.centerd)
lrtest(F2mod.x, F2mod.age)

F2mod.sex <- update(F2mod.x, ~. -sex)
lrtest(F2mod.x, F2mod.sex)

F2mod.thick <- update(F2mod.x, ~. -thickz)
lrtest(F2mod.x, F2mod.thick)

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
                age.centerd + sex + thickz, 
              data=pd.data)

qqnorm(resid(F3mod.x))
qqline(resid(F3mod.x))

F3mod.neve <- update(F3mod.x, .~. -nevent.u.m2.minz)
lrtest(F3mod.x, F3mod.neve)

F3mod.leneve <- update(F3mod.x, ~. -lenevez)
lrtest(F3mod.x, F3mod.leneve)

F3mod.tueeve <- update(F3mod.x, ~. -tueevez)
lrtest(F3mod.x, F3mod.tueeve)

F3mod.maxeve <- update(F3mod.x, ~. -maxevez)
lrtest(F3mod.x, F3mod.maxeve)

F3mod.intcpt <- update(F3mod.x, ~. -a_interceptz)
lrtest(F3mod.x, F3mod.intcpt)

F3mod.slope <- update(F3mod.x, ~. -a_slopez)
lrtest(F3mod.x, F3mod.slope)

F3mod.alpapw <- update(F3mod.x, ~. -alpha_pwz)
lrtest(F3mod.x, F3mod.alpapw)

F3mod.alpacf <- update(F3mod.x, ~. -alpha_cfz)
lrtest(F3mod.x, F3mod.alpacf)

F3mod.betapw <- update(F3mod.x, ~. -beta_pwz)
lrtest(F3mod.x, F3mod.betapw)

F3mod.betacf <- update(F3mod.x, ~. -beta_cfz)
lrtest(F3mod.x, F3mod.betacf)

F3mod.age <- update(F3mod.x, ~. -age.centerd)
lrtest(F3mod.x, F3mod.age)

F3mod.sex <- update(F3mod.x, ~. -sex)
lrtest(F3mod.x, F3mod.sex)

F3mod.thick <- update(F3mod.x, ~. -thickz)
lrtest(F3mod.x, F3mod.thick)

# SUMMARY
mod3.sim <- sim(F3mod.x, n.sims=1000)
x1 <- summary(F3mod.x)$coefficients
x2 <- t(apply(coef(mod3.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1[,1], x2)
sums

######################################################################################
# MDS-UPDRS-III Factor 4+5 combined
F45mod.x <- lm(U.F45z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                 age.centerd + sex + thickz, 
              data=pd.data)

qqnorm(resid(F45mod.x))
qqline(resid(F45mod.x))

F45mod.neve <- update(F45mod.x, .~. -nevent.u.m2.minz)
lrtest(F45mod.x, F45mod.neve)

F45mod.leneve <- update(F45mod.x, ~. -lenevez)
lrtest(F45mod.x, F45mod.leneve)

F45mod.tueeve <- update(F45mod.x, ~. -tueevez)
lrtest(F45mod.x, F45mod.tueeve)

F45mod.maxeve <- update(F45mod.x, ~. -maxevez)
lrtest(F45mod.x, F45mod.maxeve)

F45mod.intcpt <- update(F45mod.x, ~. -a_interceptz)
lrtest(F45mod.x, F45mod.intcpt)

F45mod.slope <- update(F45mod.x, ~. -a_slopez)
lrtest(F45mod.x, F45mod.slope)

F45mod.alpapw <- update(F45mod.x, ~. -alpha_pwz)
lrtest(F45mod.x, F45mod.alpapw)

F45mod.alpacf <- update(F45mod.x, ~. -alpha_cfz)
lrtest(F45mod.x, F45mod.alpacf)

F45mod.betapw <- update(F45mod.x, ~. -beta_pwz)
lrtest(F45mod.x, F45mod.betapw)

F45mod.betacf <- update(F45mod.x, ~. -beta_cfz)
lrtest(F45mod.x, F45mod.betacf)

F45mod.age <- update(F45mod.x, ~. -age.centerd)
lrtest(F45mod.x, F45mod.age)

F45mod.sex <- update(F45mod.x, ~. -sex)
lrtest(F45mod.x, F45mod.sex)

F45mod.thick <- update(F45mod.x, ~. -thickz)
lrtest(F45mod.x, F45mod.thick)

# SUMMARY
mod45.sim <- sim(F45mod.x, n.sims=1000)
x1 <- coef(F45mod.x)
x2 <- t(apply(coef(mod45.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1, x2)
sums

######################################################################################
# MDS-UPDRS-III Factor 6
F6mod.x <- lm(U.F6z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                age.centerd + sex + thickz, 
              data=pd.data)

qqnorm(resid(F6mod.x))
qqline(resid(F6mod.x))

F6mod.neve <- update(F6mod.x, .~. -nevent.u.m2.minz)
lrtest(F6mod.x, F6mod.neve)

F6mod.leneve <- update(F6mod.x, ~. -lenevez)
lrtest(F6mod.x, F6mod.leneve)

F6mod.tueeve <- update(F6mod.x, ~. -tueevez)
lrtest(F6mod.x, F6mod.tueeve)

F6mod.maxeve <- update(F6mod.x, ~. -maxevez)
lrtest(F6mod.x, F6mod.maxeve)

F6mod.intcpt <- update(F6mod.x, ~. -a_interceptz)
lrtest(F6mod.x, F6mod.intcpt)

F6mod.slope <- update(F6mod.x, ~. -a_slopez)
lrtest(F6mod.x, F6mod.slope)

F6mod.alpapw <- update(F6mod.x, ~. -alpha_pwz)
lrtest(F6mod.x, F6mod.alpapw)

F6mod.alpacf <- update(F6mod.x, ~. -alpha_cfz)
lrtest(F6mod.x, F6mod.alpacf)

F6mod.betapw <- update(F6mod.x, ~. -beta_pwz)
lrtest(F6mod.x, F6mod.betapw)

F6mod.betacf <- update(F6mod.x, ~. -beta_cfz)
lrtest(F6mod.x, F6mod.betacf)

F6mod.age <- update(F6mod.x, ~. -age.centerd)
lrtest(F6mod.x, F6mod.age)

F6mod.sex <- update(F6mod.x, ~. -sex)
lrtest(F6mod.x, F6mod.sex)

F6mod.thick <- update(F6mod.x, ~. -thickz)
lrtest(F6mod.x, F6mod.thick)

# SUMMARY
mod6.sim <- sim(F6mod.x, n.sims=1000)
x1 <- coef(F6mod.x)
x2 <- t(apply(coef(mod6.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1, x2)
sums

######################################################################################
# MDS-UPDRS-III Factor 7
F7mod.x <- glm(U.F7z ~ nevent.u.m2.minz + lenevez + tueevez + maxevez +
                 a_interceptz + a_slopez + alpha_pwz + alpha_cfz + beta_pwz + beta_cfz+
                 age.centerd + sex + thickz, 
               data=pd.data)

qqnorm(resid(F7mod.x))
qqline(resid(F7mod.x))

F7mod.neve <- update(F7mod.x, .~. -nevent.u.m2.minz)
lrtest(F7mod.x, F7mod.neve)

F7mod.leneve <- update(F7mod.x, ~. -lenevez)
lrtest(F7mod.x, F7mod.leneve)

F7mod.tueeve <- update(F7mod.x, ~. -tueevez)
lrtest(F7mod.x, F7mod.tueeve)

F7mod.maxeve <- update(F7mod.x, ~. -maxevez)
lrtest(F7mod.x, F7mod.maxeve)

F7mod.intcpt <- update(F7mod.x, ~. -a_interceptz)
lrtest(F7mod.x, F7mod.intcpt)

F7mod.slope <- update(F7mod.x, ~. -a_slopez)
lrtest(F7mod.x, F7mod.slope)

F7mod.alpapw <- update(F7mod.x, ~. -alpha_pwz)
lrtest(F7mod.x, F7mod.alpapw)

F7mod.alpacf <- update(F7mod.x, ~. -alpha_cfz)
lrtest(F7mod.x, F7mod.alpacf)

F7mod.betapw <- update(F7mod.x, ~. -beta_pwz)
lrtest(F7mod.x, F7mod.betapw)

F7mod.betacf <- update(F7mod.x, ~. -beta_cfz)
lrtest(F7mod.x, F7mod.betacf)

F7mod.age <- update(F7mod.x, ~. -age.centerd)
lrtest(F7mod.x, F7mod.age)

F7mod.sex <- update(F7mod.x, ~. -sex)
lrtest(F7mod.x, F7mod.sex)

F7mod.thick <- update(F7mod.x, ~. -thickz)
lrtest(F7mod.x, F7mod.thick)

# SUMMARY
mod7.sim <- sim(F6mod.x, n.sims=1000)
x1 <- coef(F7mod.x)
x2 <- t(apply(coef(mod7.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1, x2)
sums

#END