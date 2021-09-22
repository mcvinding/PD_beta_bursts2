# PD beta burst statistics: analysis of N events by Poisson regression
#
# <ref>
#

library(lme4)
library(arm)
library(ggplot2)
library(car)
library(brms)
# source('X://PD_longrest//scripts//functions//zscore.R')

## Load data
setwd('X://PD_longrest//groupanalysis//')
load('X://PD_longrest//groupanalysis//alldata_subj2.Rdata')

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

p1 <- get_prior(bf(nevent.u.m2.min ~ (group+age.centerd+sex+thickz)^3, family=poisson), data=alldata)

tst <- brm(nevent.u.m2.min ~ (group+age.centerd+sex+thickz)^3,
                              prior = p1,
                              family=poisson,
                              data   = alldata, 
                              warmup = 1000, iter   = 2000, 
                              chains = 2, inits  = "random",
                              cores  = 3)
                              
# save_pars('all')) 




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

## Save
save(mod.neve.Full3, file='mod_neveBF.RData')

#END