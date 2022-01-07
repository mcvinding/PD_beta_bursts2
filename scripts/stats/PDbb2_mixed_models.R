# Mixed-model analysis of burst level data (event length, time between events, max amplitude)
#
# Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., 
#  Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor 
#  rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease [Preprint]. 
#  medRxiv.org https://doi.org/10.1101/2021.06.27.21259592
#

library(lme4)
library(ggplot2)
library(arm)
library(bayestestR)
source('X://PD_longrest//scripts//functions//zscore.R')

## Load data
load(file='X://PD_longrest//groupanalysis//bbdata2.Rdata')
load(file='X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
outdir <- 'X://PD_longrest//output'

## Transform variables
bbdata$age.centerd <- bbdata$age-mean(bbdata$age)
bbdata$log.len <- log(bbdata$leneve)
bbdata$log.tue <- log(bbdata$tueeve)
bbdata$log.max <- log(bbdata$maxeve)
bbdata$thickz <- zscore(bbdata$thick)

######################################################################################
## EVENT LENGTH
lenmod <- lmer(log.len ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2) + (1|subj), 
               data=bbdata, REML=FALSE)
qqnorm(resid(lenmod))
qqline(resid(lenmod))

# LMER regression model (BIC method)
lenmod.0 <- lmer(log.len ~ 1 + (1|subj), data=bbdata, REML=FALSE)
lenmod.G   <- update(lenmod.0, ~. +group)
lenmod.A   <- update(lenmod.G, ~. +age.centerd)
lenmod.S   <- update(lenmod.A, ~. +sex)
lenmod.T   <- update(lenmod.S, ~. +thickz)
lenmod.A2  <- update(lenmod.T, ~. +I(age.centerd^2))
lenmod.GA  <- update(lenmod.A2, ~. +group:age.centerd)
lenmod.GS  <- update(lenmod.GA, ~. +group:sex)
lenmod.GT  <- update(lenmod.GS, ~. +group:thickz)
lenmod.SA  <- update(lenmod.GT, ~. +sex:age.centerd)
lenmod.AT  <- update(lenmod.SA, ~. +age.centerd:thickz)
lenmod.ST  <- update(lenmod.AT, ~. +sex:thickz)
lenmod.GAS <- update(lenmod.ST, ~. +group:age.centerd:sex)
lenmod.GAT <- update(lenmod.GAS, ~. +group:age.centerd:thickz)
lenmod.GST <- update(lenmod.GAT, ~. +group:sex:thickz)
lenmod.AST <- update(lenmod.GST, ~. +age.centerd:sex:thickz)

bayesfactor_models(lenmod.G, denominator   = lenmod.0)
bayesfactor_models(lenmod.A, denominator   = lenmod.G)
bayesfactor_models(lenmod.S, denominator   = lenmod.A)
bayesfactor_models(lenmod.T, denominator   = lenmod.S)
bayesfactor_models(lenmod.A2, denominator  = lenmod.T)
bayesfactor_models(lenmod.GA, denominator  = lenmod.A2)
bayesfactor_models(lenmod.GS, denominator  = lenmod.GA)
bayesfactor_models(lenmod.GT, denominator  = lenmod.GS)
bayesfactor_models(lenmod.SA, denominator  = lenmod.GT)
bayesfactor_models(lenmod.AT, denominator  = lenmod.SA)
bayesfactor_models(lenmod.ST, denominator  = lenmod.AT)
bayesfactor_models(lenmod.GAS, denominator = lenmod.ST)
bayesfactor_models(lenmod.GAT, denominator = lenmod.GAS)
bayesfactor_models(lenmod.GST, denominator = lenmod.GAT)
bayesfactor_models(lenmod.AST, denominator = lenmod.GAT)

# ANOVA method
anova(lenmod.G,
      lenmod.A, 
      lenmod.S,
      lenmod.T,
      lenmod.A2,
      lenmod.GA,
      lenmod.GS,
      lenmod.GT,
      lenmod.SA,
      lenmod.AT,
      lenmod.ST,
      lenmod.GAS,
      lenmod.GAT, 
      lenmod.GST,
      lenmod.AST,
      lenmod.0,
      test="Chisq")

# Model summary
mod.sim <- sim(lenmod, n.sims=1000)
cf <- fixef(mod.sim)
x1 <- fixef(lenmod)
x2 <- t(apply(cf, 2, quantile, c(0.025, 0.975)))
cbind(x1, x2)

# Age
c(exp(x1[3])*100-100,
  quantile(exp(cf[,3])*100-100, c(0.025, 0.975)))

# Sex
c(exp(x1[4])*100-100,
  quantile(exp(cf[,4])*100-100, c(0.025, 0.975)))

# Age x CT
c(exp(x1[3]+x1[5]+x1[10])*100-100,
  quantile(exp(cf[,3]+cf[,5]+cf[,10])*100-100, c(0.025, 0.975)))


######################################################################################
## TIME UNTIL EVENT
tuemod <- lmer(log.tue ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2) + (1|subj), data=bbdata, REML=FALSE)

qqnorm(resid(tuemod))
qqline(resid(tuemod))

# LMER regression model (BIC method)
tuemod.0 <- lmer(log.tue ~ 1 + (1|subj), data=bbdata, REML=FALSE)
tuemod.G   <- update(tuemod.0, ~. +group)
tuemod.A   <- update(tuemod.G, ~. +age.centerd)
tuemod.S   <- update(tuemod.A, ~. +sex)
tuemod.T   <- update(tuemod.S, ~. +thickz)
tuemod.A2   <- update(tuemod.T, ~. +I(age.centerd^2))
tuemod.GA  <- update(tuemod.A2, ~. +group:age.centerd)
tuemod.GS  <- update(tuemod.GA, ~. +group:sex)
tuemod.GT  <- update(tuemod.GS, ~. +group:thickz)
tuemod.SA  <- update(tuemod.GT, ~. +sex:age.centerd)
tuemod.AT  <- update(tuemod.SA, ~. +age.centerd:thickz)
tuemod.ST  <- update(tuemod.AT, ~. +sex:thickz)
tuemod.GAS <- update(tuemod.ST, ~. +group:age.centerd:sex)
tuemod.GAT <- update(tuemod.GAS, ~. +group:age.centerd:thickz)
tuemod.GST <- update(tuemod.GAT, ~. +group:sex:thickz)
tuemod.AST <- update(tuemod.GST, ~. +age.centerd:sex:thickz)

bayesfactor_models(tuemod.G, denominator   = tuemod.0)
bayesfactor_models(tuemod.A, denominator   = tuemod.G)
bayesfactor_models(tuemod.S, denominator   = tuemod.A)
bayesfactor_models(tuemod.T, denominator   = tuemod.S)
bayesfactor_models(tuemod.A2, denominator  = tuemod.T)
bayesfactor_models(tuemod.GA, denominator  = tuemod.A2)
bayesfactor_models(tuemod.GS, denominator  = tuemod.GA)
bayesfactor_models(tuemod.GT, denominator  = tuemod.GS)
bayesfactor_models(tuemod.SA, denominator  = tuemod.GT)
bayesfactor_models(tuemod.AT, denominator  = tuemod.SA)
bayesfactor_models(tuemod.ST, denominator  = tuemod.AT)
bayesfactor_models(tuemod.GAS, denominator = tuemod.ST)
bayesfactor_models(tuemod.GAT, denominator = tuemod.GAS)
bayesfactor_models(tuemod.GST, denominator = tuemod.GAT)
bayesfactor_models(tuemod.AST, denominator = tuemod.GAT)

# ANOVA method
anova(tuemod.G, 
      tuemod.A, 
      tuemod.S, 
      tuemod.T, 
      tuemod.A2,
      tuemod.GA,
      tuemod.GS,
      tuemod.GT,
      tuemod.SA,
      tuemod.AT,
      tuemod.ST,
      tuemod.GAS,
      tuemod.GAT,
      tuemod.GST,
      tuemod.AST,
      tuemod.0,
      test="Chisq")

# Model summary
mod.sim <- sim(tuemod, n.sims=1000)
cf <- fixef(mod.sim)
x1 <- fixef(tuemod)
x2 <- t(apply(cf, 2, quantile, c(0.025, 0.975)))
cbind(x1, x2)
(exp(cbind(x1, x2)))*100-100

# Sex x age: female
c(exp(x1[3])*100-100,
  quantile(exp(cf[,3])*100-100, c(0.025, 0.975)))

# Sex x age: male
c(exp(x1[3]+x1[9])*100-100,
  quantile(exp(cf[,3]+cf[,9])*100-100, c(0.025, 0.975)))

######################################################################################
## MAX
maxmod <- lmer(log.max ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2) + (1|subj), data=bbdata, REML=FALSE)

qqnorm(resid(maxmod))
qqline(resid(maxmod))

# LMER regression model (BIC method)
maxmod.0 <- lmer(log.max ~ 1 + (1|subj), data=bbdata, REML=FALSE)
maxmod.G   <- update(maxmod.0, ~. +group)
maxmod.A   <- update(maxmod.G, ~. +age.centerd)
maxmod.S   <- update(maxmod.A, ~. +sex)
maxmod.T   <- update(maxmod.S, ~. +thickz)
maxmod.A2  <- update(maxmod.T, ~. +I(age.centerd^2))
maxmod.GA  <- update(maxmod.A2, ~. +group:age.centerd)
maxmod.GS  <- update(maxmod.GA, ~. +group:sex)
maxmod.GT  <- update(maxmod.GS, ~. +group:thickz)
maxmod.SA  <- update(maxmod.GT, ~. +sex:age.centerd)
maxmod.AT  <- update(maxmod.SA, ~. +age.centerd:thickz)
maxmod.ST  <- update(maxmod.AT, ~. +sex:thickz)
maxmod.GAS <- update(maxmod.ST, ~. +group:age.centerd:sex)
maxmod.GAT <- update(maxmod.GAS, ~. +group:age.centerd:thickz)
maxmod.GST <- update(maxmod.GAT, ~. +group:sex:thickz)
maxmod.AST <- update(maxmod.GST, ~. +age.centerd:sex:thickz)

bayesfactor_models(maxmod.G, denominator   = maxmod.0)
bayesfactor_models(maxmod.A, denominator   = maxmod.G)
bayesfactor_models(maxmod.S, denominator   = maxmod.A)
bayesfactor_models(maxmod.T, denominator   = maxmod.S)
bayesfactor_models(maxmod.A2, denominator  = maxmod.T)
bayesfactor_models(maxmod.GA, denominator  = maxmod.A2)
bayesfactor_models(maxmod.GS, denominator  = maxmod.GA)
bayesfactor_models(maxmod.GT, denominator  = maxmod.GS)
bayesfactor_models(maxmod.SA, denominator  = maxmod.GT)
bayesfactor_models(maxmod.AT, denominator  = maxmod.SA)
bayesfactor_models(maxmod.ST, denominator  = maxmod.AT)
bayesfactor_models(maxmod.GAS, denominator = maxmod.ST)
bayesfactor_models(maxmod.GAT, denominator = maxmod.GAS)
bayesfactor_models(maxmod.GST, denominator = maxmod.GAT)
bayesfactor_models(maxmod.AST, denominator = maxmod.GAT)

# ANOVA method
anova(maxmod.G,
      maxmod.A,
      maxmod.S,
      maxmod.T,
      maxmod.A2,
      maxmod.GA,
      maxmod.GS,
      maxmod.GT,
      maxmod.SA,
      maxmod.AT,
      maxmod.ST,
      maxmod.GAS,
      maxmod.GAT,
      maxmod.GST,
      maxmod.AST,
      maxmod.0,
      test="Chisq")

# Model summary
mod.sim <- sim(maxmod, n.sims=1000)
cf <- fixef(mod.sim)
x1 <- fixef(maxmod)
x2 <- t(apply(cf, 2, quantile, c(0.025, 0.975)))
cbind(x1, x2)
(exp(cbind(x1, x2)))*100-100

# Male Ptns / male Ctrl
c(exp(x1[2]+x1[4]+x1[7])*100-100,
  quantile(exp(cf[,2]+cf[,4]+cf[,7])*100-100, c(0.025, 0.975)))

# Male ctrl / Female Ctrl
c(exp(x1[4])*100-100,
  quantile(exp(cf[,4])*100-100, c(0.025, 0.975)))

# Male Ptns / Female pten
c(exp(x1[4]+x1[7])*100-100,
  quantile(exp(cf[,4]+cf[,7])*100-100, c(0.025, 0.975)))

### Save models
setwd(outdir)
save(lenmod, file='lenmod.RData')
save(tuemod, file='tuemod.RData')
save(maxmod, file='maxmod.RData')

#END