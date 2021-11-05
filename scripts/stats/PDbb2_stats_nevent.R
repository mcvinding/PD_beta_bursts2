# PD beta burst statistics: analysis of N events by Poisson regression
#
# Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., 
#  Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor 
#  rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease [Preprint]. 
#  medRxiv.org https://doi.org/10.1101/2021.06.27.21259592
#

library(lme4)
library(arm)
library(ggplot2)
library(car)
library(bayestestR)

## Load data
setwd('X://PD_longrest//groupanalysis//')
load('X://PD_longrest//groupanalysis//alldata_subj2.Rdata')

# ######################################################################################
# LMER regression model (ANOVA method)
mod.neve.Full3 <- glm(nevent.u.m2.min ~ (group+age.centerd+sex+thickz)^3, data=alldata, family=poisson)
anova(mod.neve.Full3, test="Chisq")

# LMER regression model (BIC method)
mod.neve.0   <- glm(nevent.u.m2.min ~ 1, data=alldata, family=poisson)
mod.neve.G   <- update(mod.neve.0, ~. + group)
mod.neve.A   <- update(mod.neve.G, ~. + age.centerd)
mod.neve.S   <- update(mod.neve.A, ~. + sex)
mod.neve.T   <- update(mod.neve.S, ~. + thickz)
mod.neve.GA  <- update(mod.neve.T, ~. + group:age.centerd)
mod.neve.GS  <- update(mod.neve.GA, ~. + group:sex)
mod.neve.GT  <- update(mod.neve.GS, ~. + group:thickz)
mod.neve.AS  <- update(mod.neve.GT, ~. + age.centerd:sex)
mod.neve.AT  <- update(mod.neve.AS, ~. + age.centerd:thickz)
mod.neve.ST  <- update(mod.neve.AT, ~. + sex:thickz)
mod.neve.GAS <- update(mod.neve.ST, ~. + group:age.centerd:sex)
mod.neve.GAT <- update(mod.neve.GAS, ~. + group:age.centerd:thickz)
mod.neve.GST <- update(mod.neve.GAT, ~. + group:sex:thickz)
mod.neve.AST <- update(mod.neve.GST, ~. + age.centerd:sex:thickz)


## BF test
bayesfactor_models(mod.neve.G, denominator = mod.neve.0)
bayesfactor_models(mod.neve.A, denominator = mod.neve.G)
bayesfactor_models(mod.neve.S, denominator = mod.neve.A)
bayesfactor_models(mod.neve.T, denominator = mod.neve.S)
bayesfactor_models(mod.neve.GA, denominator = mod.neve.T)
bayesfactor_models(mod.neve.GS, denominator = mod.neve.GA)
bayesfactor_models(mod.neve.GT, denominator = mod.neve.GS)
bayesfactor_models(mod.neve.AS, denominator = mod.neve.GT)
bayesfactor_models(mod.neve.AT, denominator = mod.neve.AS)
bayesfactor_models(mod.neve.ST, denominator = mod.neve.AT)
bayesfactor_models(mod.neve.GAS, denominator = mod.neve.ST)
bayesfactor_models(mod.neve.GAT, denominator = mod.neve.GAS)
bayesfactor_models(mod.neve.GST, denominator = mod.neve.GAT)
bayesfactor_models(mod.neve.AST, denominator = mod.neve.GST)


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
setwd('X://PD_longrest//output')
save(mod.neve.Full3, file='mod_neveBF.RData')

#END