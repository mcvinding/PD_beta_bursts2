# Fit models for the main analysis. Analysis of signal features as as a function 
# of age, sex, group, and cortical thickness.
#
# Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., 
#  Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor 
#  rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease [Preprint]. 
#  medRxiv.org https://doi.org/10.1101/2021.06.27.21259592
#

## Overwrite output files?
overwrite = FALSE

## Libraries
library(lme4)
library(arm)
library(ggplot2)
library(car)
library(brms)
source('/home/mikkel/PD_longrest/scripts/functions/zscore.R')

## Load data
wrkdir <- '/home/mikkel/PD_longrest/groupanalysis'
outdir <- '/home/mikkel/PD_longrest/output'
setwd(outdir)

load(file='/home/mikkel/PD_longrest/groupanalysis/alldata_subj2.Rdata')
load(file='/home/mikkel/PD_longrest/groupanalysis/bbdata2.Rdata')
peaks <- read.csv("peaks_no_fooof.csv", sep = ",")
peaks$subj <- paste(rep(0,1), peaks$subj, sep = "")
alldata <- merge(alldata, peaks, by="subj")

## Transform variables
bbdata$age.centerd <- bbdata$age-mean(bbdata$age)
bbdata$log.len <- log(bbdata$leneve)
bbdata$log.tue <- log(bbdata$tueeve)
bbdata$log.max <- log(bbdata$maxeve)
bbdata$thickz <- zscore(bbdata$thick)

################################################################################
### PSD MODELS 
################################################################################

######################################################################################
# 1/f intercept

# finter.Full3 <- lm(a_intercept ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
# finter.Full2 <- lm(a_intercept ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata)
# 
# qqnorm(resid(finter.Full2))
# qqline(resid(finter.Full2))

# Bayes with BRMS
if (!file.exists('mod_finter_2.RData') || overwrite){
  print('1/f intercept model')
  finter.mod <- brm(a_intercept ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2),
                    data=alldata,
                    iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(finter.mod, file='mod_finter_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_finter_mod.RData') || overwrite){
  print('1/f intercept step-models')
  finter.0   <- brm(a_intercept ~ 1, data=alldata, iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  finter.G   <- update(finter.0, ~. +group, newdata=alldata, cores=4)
  finter.A   <- update(finter.G, ~. +age.centerd, newdata=alldata, cores=4)
  finter.S   <- update(finter.A, ~. +sex, newdata=alldata, cores=4)
  finter.T   <- update(finter.S, ~. +thickz, newdata=alldata, cores=4)
  finter.A2  <- update(finter.T, ~. +I(age.centerd^2), newdata=alldata, cores=4)
  finter.GA  <- update(finter.A2, ~. +group:age.centerd, newdata=alldata, cores=4)
  finter.GS  <- update(finter.GA, ~. +group:sex, newdata=alldata, cores=4)
  finter.GT  <- update(finter.GS, ~. +group:thickz, newdata=alldata, cores=4)
  finter.SA  <- update(finter.GT, ~. +sex:age.centerd, newdata=alldata, cores=4)
  finter.AT  <- update(finter.SA, ~. +age.centerd:thickz, newdata=alldata, cores=4)
  finter.ST  <- update(finter.AT, ~. +sex:thickz, newdata=alldata, cores=4)
  # finter.GAS <- update(finter.ST, ~. +group:age.centerd:sex)
  # finter.GAT <- update(finter.GAS, ~. +group:age.centerd:thickz)
  # finter.GST <- update(finter.GAT, ~. +group:sex:thickz)
  # finter.AST <- update(finter.GST, ~. +age.centerd:sex:thickz)

  # Save and clear (to save memory)
  save(list = ls(pattern = "^finter."), file = 'all_finter_mod.RData')
  rm(list = ls(pattern = "^finter"))
}

######################################################################################
# 1/f slope
# fslope.Full3 <- lm(a_slope ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
# fslope.Full2 <- lm(a_slope ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata)
# 
# qqnorm(resid(fslope.Full2))
# qqline(resid(fslope.Full2))

# Bayes with BRMS
if (!file.exists('mod_fslope_2.RData') || overwrite){
  print('1/f slope model')
  fslope.mod <- brm(a_slope ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), 
                    data=alldata,
                    iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(fslope.mod, file='mod_fslope_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_fslope_mod.RData') || overwrite){
  print('1/f slope step-model')
  fslope.0   <- brm(a_slope ~ 1, data=alldata,  iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  fslope.G   <- update(fslope.0, ~. +group, newdata=alldata, cores=4)
  fslope.A   <- update(fslope.G, ~. +age.centerd, newdata=alldata, cores=4)
  fslope.S   <- update(fslope.A, ~. +sex, newdata=alldata, cores=4)
  fslope.T   <- update(fslope.S, ~. +thickz, newdata=alldata, cores=4)
  fslope.A2  <- update(fslope.T, ~. +I(age.centerd^2), newdata=alldata, cores=4)
  fslope.GA  <- update(fslope.A2, ~. +group:age.centerd, newdata=alldata, cores=4)
  fslope.GS  <- update(fslope.GA, ~. +group:sex, newdata=alldata, cores=4)
  fslope.GT  <- update(fslope.GS, ~. +group:thickz, newdata=alldata, cores=4)
  fslope.SA  <- update(fslope.GT, ~. +sex:age.centerd, newdata=alldata, cores=4)
  fslope.AT  <- update(fslope.SA, ~. +age.centerd:thickz, newdata=alldata, cores=4)
  fslope.ST  <- update(fslope.AT, ~. +sex:thickz, newdata=alldata, cores=4)
  # fslope.GAS <- update(fslope.ST, ~. +group:age.centerd:sex)
  # fslope.GAT <- update(fslope.GAS, ~. +group:age.centerd:thickz)
  # fslope.GST <- update(fslope.GAT, ~. +group:sex:thickz)
  # fslope.AST <- update(fslope.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^fslope."), file = 'all_fslope_mod.RData')
  rm(list = ls(pattern = "^fslope"))
}

######################################################################################
# Beta power
# beta_pw.Full3 <- lm(beta_pw ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
# beta_pw.Full2 <- lm(beta_pw ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata)
# 
# qqnorm(resid(beta_pw.Full2))
# qqline(resid(beta_pw.Full2))

# Bayes with BRMS
if (!file.exists('mod_betapw_2.RData') || overwrite){
  print('beta power model')
  beta_pw.mod <- brm(beta_pw ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), 
                    data=alldata,
                    iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(beta_pw.mod, file='mod_betapw_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_betapw_mod.RData') || overwrite){
  print('beta power step-models')
  beta_pw.0   <- brm(beta_pw ~ 1, data=alldata,  iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  beta_pw.G   <- update(beta_pw.0, ~. +group, newdata=alldata, cores=4)
  beta_pw.A   <- update(beta_pw.G, ~. +age.centerd, newdata=alldata, cores=4)
  beta_pw.S   <- update(beta_pw.A, ~. +sex, newdata=alldata, cores=4)
  beta_pw.T   <- update(beta_pw.S, ~. +thickz, newdata=alldata, cores=4)
  beta_pw.A2  <- update(beta_pw.T, ~. +I(age.centerd^2), newdata=alldata, cores=4)
  beta_pw.GA  <- update(beta_pw.A2, ~. +group:age.centerd, newdata=alldata, cores=4)
  beta_pw.GS  <- update(beta_pw.GA, ~. +group:sex, newdata=alldata, cores=4)
  beta_pw.GT  <- update(beta_pw.GS, ~. +group:thickz, newdata=alldata, cores=4)
  beta_pw.SA  <- update(beta_pw.GT, ~. +sex:age.centerd, newdata=alldata, cores=4)
  beta_pw.AT  <- update(beta_pw.SA, ~. +age.centerd:thickz, newdata=alldata, cores=4)
  beta_pw.ST  <- update(beta_pw.AT, ~. +sex:thickz, newdata=alldata, cores=4)
  # beta_pw.GAS <- update(beta_pw.ST, ~. +group:age.centerd:sex)
  # beta_pw.GAT <- update(beta_pw.GAS, ~. +group:age.centerd:thickz)
  # beta_pw.GST <- update(beta_pw.GAT, ~. +group:sex:thickz)
  # beta_pw.AST <- update(beta_pw.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^beta_pw."), file = 'all_betapw_mod.RData')
  rm(list = ls(pattern = "^beta_pw"))
}
######################################################################################
# Beta peak freq
# beta_cf.Full3 <- lm(beta_cf ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
# beta_cf.Full2 <- lm(beta_cf ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata)
# 
# qqnorm(resid(beta_cf.Full2))
# qqline(resid(beta_cf.Full2))

# Bayes with BRMS
if (!file.exists('mod_betacf_2.RData') || overwrite){
  print('BETA CENTER FREQUENCY model')
  beta_cf.mod <- brm(beta_cf ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), 
                     data=alldata,
                     iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(beta_cf.mod, file='mod_betacf_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_betacf_mod.RData') || overwrite){
  print('BETA CENTER FREQUENCY step-models')
  beta_cf.0   <- brm(beta_cf ~ 1, data=alldata,  iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  beta_cf.G   <- update(beta_cf.0, ~. +group, newdata=alldata, cores=4)
  beta_cf.A   <- update(beta_cf.G, ~. +age.centerd, newdata=alldata, cores=4)
  beta_cf.S   <- update(beta_cf.A, ~. +sex, newdata=alldata, cores=4)
  beta_cf.T   <- update(beta_cf.S, ~. +thickz, newdata=alldata, cores=4)
  beta_cf.A2  <- update(beta_cf.T, ~. +I(age.centerd^2), newdata=alldata, cores=4)
  beta_cf.GA  <- update(beta_cf.A2, ~. +group:age.centerd, newdata=alldata, cores=4)
  beta_cf.GS  <- update(beta_cf.GA, ~. +group:sex, newdata=alldata, cores=4)
  beta_cf.GT  <- update(beta_cf.GS, ~. +group:thickz, newdata=alldata, cores=4)
  beta_cf.SA  <- update(beta_cf.GT, ~. +sex:age.centerd, newdata=alldata, cores=4)
  beta_cf.AT  <- update(beta_cf.SA, ~. +age.centerd:thickz, newdata=alldata, cores=4)
  beta_cf.ST  <- update(beta_cf.AT, ~. +sex:thickz, newdata=alldata, cores=4)
  # beta_cf.GAS <- update(beta_cf.ST, ~. +group:age.centerd:sex)
  # beta_cf.GAT <- update(beta_cf.GAS, ~. +group:age.centerd:thickz)
  # beta_cf.GST <- update(beta_cf.GAT, ~. +group:sex:thickz)
  # beta_cf.AST <- update(beta_cf.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^beta_cf."), file = 'all_betacf_mod.RData')
  rm(list = ls(pattern = "^beta_cf"))
}

######################################################################################
# Alpha power
# alpha_pw.Full3 <- lm(alpha_pw ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
# alpha_pw.Full2 <- lm(alpha_pw ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata)
# 
# qqnorm(resid(alpha_pw.Full2))
# qqline(resid(alpha_pw.Full2))

# Bayes with BRMS
if (!file.exists('mod_alphapw_2.RData') || overwrite){
  print('ALPHA power model')
  alpha_pw.mod <- brm(alpha_pw ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), 
                     data=alldata,
                     iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(alpha_pw.mod, file='mod_alphapw_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_alphapw_mod.RData') || overwrite){
  print('ALPHA power step-models')
  alpha_pw.0   <- brm(alpha_pw ~ 1, data=alldata,  iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  alpha_pw.G   <- update(alpha_pw.0, ~. +group, newdata=alldata, cores=4)
  alpha_pw.A   <- update(alpha_pw.G, ~. +age.centerd, newdata=alldata, cores=4)
  alpha_pw.S   <- update(alpha_pw.A, ~. +sex, newdata=alldata, cores=4)
  alpha_pw.T   <- update(alpha_pw.S, ~. +thickz, newdata=alldata, cores=4)
  alpha_pw.A2  <- update(alpha_pw.T, ~. +I(age.centerd^2), newdata=alldata, cores=4)
  alpha_pw.GA  <- update(alpha_pw.A2, ~. +group:age.centerd, newdata=alldata, cores=4)
  alpha_pw.GS  <- update(alpha_pw.GA, ~. +group:sex, newdata=alldata, cores=4)
  alpha_pw.GT  <- update(alpha_pw.GS, ~. +group:thickz, newdata=alldata, cores=4)
  alpha_pw.SA  <- update(alpha_pw.GT, ~. +sex:age.centerd, newdata=alldata, cores=4)
  alpha_pw.AT  <- update(alpha_pw.SA, ~. +age.centerd:thickz, newdata=alldata, cores=4)
  alpha_pw.ST  <- update(alpha_pw.AT, ~. +sex:thickz, newdata=alldata, cores=4)
  # alpha_pw.GAS <- update(alpha_pw.ST, ~. +group:age.centerd:sex)
  # alpha_pw.GAT <- update(alpha_pw.GAS, ~. +group:age.centerd:thickz)
  # alpha_pw.GST <- update(alpha_pw.GAT, ~. +group:sex:thickz)
  # alpha_pw.AST <- update(alpha_pw.GST, ~. +age.centerd:sex:thickz)

  # Save and clear (to save memory)
  save(list = ls(pattern = "^alpha_pw."), file = 'all_alphapw_mod.RData')
  rm(list = ls(pattern = "^alpha_pw"))
}

######################################################################################
# Alpha peak freq
# alpha_cf.Full3 <- lm(alpha_cf ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
# alpha_cf.Full2 <- lm(alpha_cf ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata)
# 
# qqnorm(resid(alpha_cf.Full2))
# qqline(resid(alpha_cf.Full2))

# Bayes with BRMS
if (!file.exists('mod_alphacf_2.RData') || overwrite){
  print('ALPHA CENTER FREQUENCY  model')
  alpha_cf.mod <- brm(alpha_cf ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), 
                      data=alldata,
                      iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(alpha_cf.mod, file='mod_alphacf_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_alphacf_mod.RData') || overwrite){
  print('ALPHA CENTER FREQUENCY step-models')
  alpha_cf.0   <- brm(alpha_cf ~ 1, data=alldata,  iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  alpha_cf.G   <- update(alpha_cf.0, ~. +group, newdata=alldata, cores=4)
  alpha_cf.A   <- update(alpha_cf.G, ~. +age.centerd, newdata=alldata, cores=4)
  alpha_cf.S   <- update(alpha_cf.A, ~. +sex, newdata=alldata, cores=4)
  alpha_cf.T   <- update(alpha_cf.S, ~. +thickz, newdata=alldata, cores=4)
  alpha_cf.A2  <- update(alpha_cf.T, ~. +I(age.centerd^2), newdata=alldata, cores=4)
  alpha_cf.GA  <- update(alpha_cf.A2, ~. +group:age.centerd, newdata=alldata, cores=4)
  alpha_cf.GS  <- update(alpha_cf.GA, ~. +group:sex, newdata=alldata, cores=4)
  alpha_cf.GT  <- update(alpha_cf.GS, ~. +group:thickz, newdata=alldata, cores=4)
  alpha_cf.SA  <- update(alpha_cf.GT, ~. +sex:age.centerd, newdata=alldata, cores=4)
  alpha_cf.AT  <- update(alpha_cf.SA, ~. +age.centerd:thickz, newdata=alldata, cores=4)
  alpha_cf.ST  <- update(alpha_cf.AT, ~. +sex:thickz, newdata=alldata, cores=4)
  # alpha_cf.GAS <- update(alpha_cf.ST, ~. +group:age.centerd:sex)
  # alpha_cf.GAT <- update(alpha_cf.GAS, ~. +group:age.centerd:thickz)
  # alpha_cf.GST <- update(alpha_cf.GAT, ~. +group:sex:thickz)
  # alpha_cf.AST <- update(alpha_cf.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^alpha_cf."), file = 'all_alphacf_mod.RData')
  rm(list = ls(pattern = "^alpha_cf"))
}

################################################################################
### BURST MODELS
################################################################################

# ##############################################################################
## N events by Poisson regression

# LMER regression model (ANOVA method)
# mod.neve.Full2_lm <- glm(nevent.u.m2.min ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata, family=poisson)
# mod.neve.Full3_lm <- glm(nevent.u.m2.min ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata, family=poisson)
# 
# anova(mod.neve.Full2_lm, test="Chisq")
# anova(mod.neve.Full3_lm, test="Chisq")

# Bayes with BRMS
if (!file.exists('mod_neveBF_2.RData') || overwrite){
  print('N events model')
  mod_neve <- brm(nevent.u.m2.min ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), 
                  data=alldata, family=poisson, 
                  iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(mod_neve, file='mod_neveBF_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_nevmod.RData') || overwrite){
  print('N events step-models')
  mod.neve.0   <- brm(nevent.u.m2.min ~ 1, data=alldata, family=poisson, iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  mod.neve.G   <- update(mod.neve.0, ~. + group, newdata=alldata, cores=4)
  mod.neve.A   <- update(mod.neve.G, ~. + age.centerd, newdata=alldata, cores=4)
  mod.neve.S   <- update(mod.neve.A, ~. + sex, newdata=alldata,  cores=4)
  mod.neve.T   <- update(mod.neve.S, ~. + thickz, newdata=alldata, cores=4)
  mod.neve.A2  <- update(mod.neve.T, ~. + I(age.centerd^2), newdata=alldata, cores=4)
  mod.neve.GA  <- update(mod.neve.A2, ~. + group:age.centerd, newdata=alldata, cores=4)
  mod.neve.GS  <- update(mod.neve.GA, ~. + group:sex, newdata=alldata, cores=4)
  mod.neve.GT  <- update(mod.neve.GS, ~. + group:thickz, newdata=alldata, cores=4)
  mod.neve.AS  <- update(mod.neve.GT, ~. + age.centerd:sex, newdata=alldata, cores=4)
  mod.neve.AT  <- update(mod.neve.AS, ~. + age.centerd:thickz, newdata=alldata, cores=4)
  mod.neve.ST  <- update(mod.neve.AT, ~. + sex:thickz, newdata=alldata, cores=4)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^mod.neve."), file = 'all_nevmod.RData')
  rm(list = ls(pattern = "^mod.neve"))
}

# ##############################################################################
## EVENT LENGTH
# lenmod <- lmer(log.len ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2) + (1|subj), 
#                data=bbdata, REML=FALSE)
# qqnorm(resid(lenmod))
# qqline(resid(lenmod))

# Bayes with BRMS
if (!file.exists('lenmod_2.RData') || overwrite){
  print(' EVENT LENGTH model')
  lenmod <- brm(log.len ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2) + (1|subj), 
                data=bbdata,
                iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(lenmod, file='lenmod_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_lenmod.RData') || overwrite){
  print(' EVENT LENGTH step-model')
  lenmod.0 <- brm(log.len ~ 1 + (1|subj), data=bbdata, iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  lenmod.G   <- update(lenmod.0, ~. +group, newdata=bbdata, cores=4)
  lenmod.A   <- update(lenmod.G, ~. +age.centerd, newdata=bbdata, cores=4)
  lenmod.S   <- update(lenmod.A, ~. +sex, newdata=bbdata, cores=4)
  lenmod.T   <- update(lenmod.S, ~. +thickz, newdata=bbdata, cores=4)
  lenmod.A2  <- update(lenmod.T, ~. +I(age.centerd^2), newdata=bbdata, cores=4)
  lenmod.GA  <- update(lenmod.A2, ~. +group:age.centerd, newdata=bbdata, cores=4)
  lenmod.GS  <- update(lenmod.GA, ~. +group:sex, newdata=bbdata, cores=4)
  lenmod.GT  <- update(lenmod.GS, ~. +group:thickz, newdata=bbdata, cores=4)
  lenmod.SA  <- update(lenmod.GT, ~. +sex:age.centerd, newdata=bbdata, cores=4)
  lenmod.AT  <- update(lenmod.SA, ~. +age.centerd:thickz, newdata=bbdata, cores=4)
  lenmod.ST  <- update(lenmod.AT, ~. +sex:thickz, newdata=bbdata, cores=4)
  # lenmod.GAS <- update(lenmod.ST, ~. +group:age.centerd:sex)
  # lenmod.GAT <- update(lenmod.GAS, ~. +group:age.centerd:thickz)
  # lenmod.GST <- update(lenmod.GAT, ~. +group:sex:thickz)
  # lenmod.AST <- update(lenmod.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^lenmod."), file = 'all_lenmod.RData')
  rm(list = ls(pattern = "^lenmod"))
}

######################################################################################
## TIME UNTIL EVENT mixed model
# tuemod_ln <- lmer(log.tue ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2) + (1|subj), data=bbdata, REML=FALSE)
# 
# qqnorm(resid(tuemod))
# qqline(resid(tuemod))

# Bayes with BRMS
if (!file.exists('tuemod_2.RData') || overwrite){
  print('TIME UNTIL EVENT model')
  tuemod <- brm(log.tue ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2) + (1|subj), 
                data=bbdata,
                iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(tuemod, file='tuemod_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_tuemod.RData') || overwrite){
  print('TIME UNTIL EVENT step-models')
  tuemod.0   <- brm(log.tue ~ 1 + (1|subj), data=bbdata, iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  tuemod.G   <- update(tuemod.0, ~. +group, newdata=bbdata, cores=4)
  tuemod.A   <- update(tuemod.G, ~. +age.centerd, newdata=bbdata, cores=4)
  tuemod.S   <- update(tuemod.A, ~. +sex, newdata=bbdata, cores=4)
  tuemod.T   <- update(tuemod.S, ~. +thickz, newdata=bbdata, cores=4)
  tuemod.A2   <- update(tuemod.T, ~. +I(age.centerd^2), newdata=bbdata, cores=4)
  tuemod.GA  <- update(tuemod.A2, ~. +group:age.centerd, newdata=bbdata, cores=4)
  tuemod.GS  <- update(tuemod.GA, ~. +group:sex, newdata=bbdata, cores=4)
  tuemod.GT  <- update(tuemod.GS, ~. +group:thickz, newdata=bbdata, cores=4)
  tuemod.SA  <- update(tuemod.GT, ~. +sex:age.centerd, newdata=bbdata, cores=4)
  tuemod.AT  <- update(tuemod.SA, ~. +age.centerd:thickz, newdata=bbdata, cores=4)
  tuemod.ST  <- update(tuemod.AT, ~. +sex:thickz, newdata=bbdata, cores=4)
  # tuemod.GAS <- update(tuemod.ST, ~. +group:age.centerd:sex)
  # tuemod.GAT <- update(tuemod.GAS, ~. +group:age.centerd:thickz)
  # tuemod.GST <- update(tuemod.GAT, ~. +group:sex:thickz)
  # tuemod.AST <- update(tuemod.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^tuemod."), file = 'all_tuemod.RData')
  rm(list = ls(pattern = "^tuemod"))
}

######################################################################################
## Burst max amplitude mixed model
# maxmod_ln <- lmer(log.max ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2) + (1|subj), data=bbdata, REML=FALSE)
# qqnorm(resid(maxmod))
# qqline(resid(maxmod))

# Bayes with BRMS
if (!file.exists('maxmod_2.RData') || overwrite){
  print('Burst max amplitude mixed model')
  print('Max amplitude model')
  maxmod <- brm(log.max ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2) + (1|subj), 
                data=bbdata,
                iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(maxmod, file='maxmod_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_maxmod.RData') || overwrite){
  print('Max amplitude step-models')
  maxmod.0   <- brm(log.max ~ 1 + (1|subj), data=bbdata, iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  maxmod.G   <- update(maxmod.0, ~. +group, newdata=bbdata, cores=4)
  maxmod.A   <- update(maxmod.G, ~. +age.centerd, newdata=bbdata, cores=4)
  maxmod.S   <- update(maxmod.A, ~. +sex, newdata=bbdata, cores=4)
  maxmod.T   <- update(maxmod.S, ~. +thickz, newdata=bbdata, cores=4)
  maxmod.A2  <- update(maxmod.T, ~. +I(age.centerd^2), newdata=bbdata, cores=4)
  maxmod.GA  <- update(maxmod.A2, ~. +group:age.centerd, newdata=bbdata, cores=4)
  maxmod.GS  <- update(maxmod.GA, ~. +group:sex, newdata=bbdata, cores=4)
  maxmod.GT  <- update(maxmod.GS, ~. +group:thickz, newdata=bbdata, cores=4)
  maxmod.SA  <- update(maxmod.GT, ~. +sex:age.centerd, newdata=bbdata, cores=4)
  maxmod.AT  <- update(maxmod.SA, ~. +age.centerd:thickz, newdata=bbdata, cores=4)
  maxmod.ST  <- update(maxmod.AT, ~. +sex:thickz, newdata=bbdata, cores=4)
  # maxmod.GAS <- update(maxmod.ST, ~. +group:age.centerd:sex)
  # maxmod.GAT <- update(maxmod.GAS, ~. +group:age.centerd:thickz)
  # maxmod.GST <- update(maxmod.GAT, ~. +group:sex:thickz)
  # maxmod.AST <- update(maxmod.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^maxmod."), file = 'all_maxmod.RData')
  rm(list = ls(pattern = "^maxmod"))
}


######################################################################################
# "Raw" beta power
# beta_pw.Full3 <- lm(beta_pw ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
raw_beta_pw.Full2 <- lm(raw_beta_pw ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata)
summary(raw_beta_pw.Full2)
# qqnorm(resid(beta_pw.Full2))
# qqline(resid(beta_pw.Full2))

# Bayes with BRMS
if (!file.exists('mod_rawbetapw_2.RData') || overwrite){
  print('raw beta power model')
  raw_beta_pw.mod <- brm(raw_beta_pw ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), 
                     data=alldata,
                     iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(raw_beta_pw.mod, file='mod_rawbetapw_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_rawbetapw_mod.RData') || overwrite){
  print('raw beta power step-models')
  raw_beta_pw.0   <- brm(raw_beta_pw ~ 1, data=alldata,  iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  raw_beta_pw.G   <- update(raw_beta_pw.0, ~. +group, newdata=alldata, cores=4)
  raw_beta_pw.A   <- update(raw_beta_pw.G, ~. +age.centerd, newdata=alldata, cores=4)
  raw_beta_pw.S   <- update(raw_beta_pw.A, ~. +sex, newdata=alldata, cores=4)
  raw_beta_pw.T   <- update(raw_beta_pw.S, ~. +thickz, newdata=alldata, cores=4)
  raw_beta_pw.A2  <- update(raw_beta_pw.T, ~. +I(age.centerd^2), newdata=alldata, cores=4)
  raw_beta_pw.GA  <- update(raw_beta_pw.A2, ~. +group:age.centerd, newdata=alldata, cores=4)
  raw_beta_pw.GS  <- update(raw_beta_pw.GA, ~. +group:sex, newdata=alldata, cores=4)
  raw_beta_pw.GT  <- update(raw_beta_pw.GS, ~. +group:thickz, newdata=alldata, cores=4)
  raw_beta_pw.SA  <- update(raw_beta_pw.GT, ~. +sex:age.centerd, newdata=alldata, cores=4)
  raw_beta_pw.AT  <- update(raw_beta_pw.SA, ~. +age.centerd:thickz, newdata=alldata, cores=4)
  raw_beta_pw.ST  <- update(raw_beta_pw.AT, ~. +sex:thickz, newdata=alldata, cores=4)
  # beta_pw.GAS <- update(beta_pw.ST, ~. +group:age.centerd:sex)
  # beta_pw.GAT <- update(beta_pw.GAS, ~. +group:age.centerd:thickz)
  # beta_pw.GST <- update(beta_pw.GAT, ~. +group:sex:thickz)
  # beta_pw.AST <- update(beta_pw.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^raw_beta_pw."), file = 'all_rawbetapw_mod.RData')
  rm(list = ls(pattern = "^raw_beta_pw"))
}
######################################################################################
# "Raw" beta peak freq
#beta_cf.Full3 <- lm(beta_cf ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
raw_beta_cf.Full2 <- lm(raw_beta_cf ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata)
summary(raw_beta_cf.Full2)
# qqnorm(resid(beta_cf.Full2))
# qqline(resid(beta_cf.Full2))

# Bayes with BRMS
if (!file.exists('mod_rawbetacf_2.RData') || overwrite){
  print('UNFILTERED BETA CENTER FREQUENCY model')
  raw_beta_cf.mod <- brm(raw_beta_cf ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), 
                     data=alldata,
                     iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(raw_beta_cf.mod, file='mod_rawbetacf_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_rawbetacf_mod.RData') || overwrite){
  print('UNFILTERED BETA CENTER FREQUENCY step-models')
  raw_beta_cf.0   <- brm(raw_beta_cf ~ 1, data=alldata,  iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  raw_beta_cf.G   <- update(raw_beta_cf.0, ~. +group, newdata=alldata, cores=4)
  raw_beta_cf.A   <- update(raw_beta_cf.G, ~. +age.centerd, newdata=alldata, cores=4)
  raw_beta_cf.S   <- update(raw_beta_cf.A, ~. +sex, newdata=alldata, cores=4)
  raw_beta_cf.T   <- update(raw_beta_cf.S, ~. +thickz, newdata=alldata, cores=4)
  raw_beta_cf.A2  <- update(raw_beta_cf.T, ~. +I(age.centerd^2), newdata=alldata, cores=4)
  raw_beta_cf.GA  <- update(raw_beta_cf.A2, ~. +group:age.centerd, newdata=alldata, cores=4)
  raw_beta_cf.GS  <- update(raw_beta_cf.GA, ~. +group:sex, newdata=alldata, cores=4)
  raw_beta_cf.GT  <- update(raw_beta_cf.GS, ~. +group:thickz, newdata=alldata, cores=4)
  raw_beta_cf.SA  <- update(raw_beta_cf.GT, ~. +sex:age.centerd, newdata=alldata, cores=4)
  raw_beta_cf.AT  <- update(raw_beta_cf.SA, ~. +age.centerd:thickz, newdata=alldata, cores=4)
  raw_beta_cf.ST  <- update(raw_beta_cf.AT, ~. +sex:thickz, newdata=alldata, cores=4)
  # beta_cf.GAS <- update(beta_cf.ST, ~. +group:age.centerd:sex)
  # beta_cf.GAT <- update(beta_cf.GAS, ~. +group:age.centerd:thickz)
  # beta_cf.GST <- update(beta_cf.GAT, ~. +group:sex:thickz)
  # beta_cf.AST <- update(beta_cf.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^raw_beta_cf."), file = 'all_rawbetacf_mod.RData')
  rm(list = ls(pattern = "^raw_beta_cf"))
}

######################################################################################
# "Raw" Alpha power
#raw_alpha_pw.Full3 <- lm(raw_al_pw ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
raw_alpha_pw.Full2 <- lm(raw_alpha_pw ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata)
summary(raw_alpha_pw.Full2)
# qqnorm(resid(alpha_pw.Full2))
# qqline(resid(alpha_pw.Full2))

# Bayes with BRMS
if (!file.exists('mod_rawalphapw_2.RData') || overwrite){
  print('UNFILTERED ALPHA power model')
  raw_alpha_pw.mod <- brm( raw_alpha_pw ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), 
                      data=alldata,
                      iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(raw_alpha_pw.mod, file='mod_rawalphapw_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_rawalphapw_mod.RData') || overwrite){
  print('raw ALPHA power step-models')
  raw_alpha_pw.0   <- brm( raw_alpha_pw ~ 1, data=alldata,  iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  raw_alpha_pw.G   <- update(raw_alpha_pw.0, ~. +group, newdata=alldata, cores=4)
  raw_alpha_pw.A   <- update(raw_alpha_pw.G, ~. +age.centerd, newdata=alldata, cores=4)
  raw_alpha_pw.S   <- update(raw_alpha_pw.A, ~. +sex, newdata=alldata, cores=4)
  raw_alpha_pw.T   <- update(raw_alpha_pw.S, ~. +thickz, newdata=alldata, cores=4)
  raw_alpha_pw.A2  <- update(raw_alpha_pw.T, ~. +I(age.centerd^2), newdata=alldata, cores=4)
  raw_alpha_pw.GA  <- update(raw_alpha_pw.A2, ~. +group:age.centerd, newdata=alldata, cores=4)
  raw_alpha_pw.GS  <- update(raw_alpha_pw.GA, ~. +group:sex, newdata=alldata, cores=4)
  raw_alpha_pw.GT  <- update(raw_alpha_pw.GS, ~. +group:thickz, newdata=alldata, cores=4)
  raw_alpha_pw.SA  <- update(raw_alpha_pw.GT, ~. +sex:age.centerd, newdata=alldata, cores=4)
  raw_alpha_pw.AT  <- update(raw_alpha_pw.SA, ~. +age.centerd:thickz, newdata=alldata, cores=4)
  raw_alpha_pw.ST  <- update(raw_alpha_pw.AT, ~. +sex:thickz, newdata=alldata, cores=4)
  # alpha_pw.GAS <- update(alpha_pw.ST, ~. +group:age.centerd:sex)
  # alpha_pw.GAT <- update(alpha_pw.GAS, ~. +group:age.centerd:thickz)
  # alpha_pw.GST <- update(alpha_pw.GAT, ~. +group:sex:thickz)
  # alpha_pw.AST <- update(alpha_pw.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^raw_alpha_pw."), file = 'all_rawalphapw_mod.RData')
  rm(list = ls(pattern = "^raw_alpha_pw"))
}

######################################################################################
# Alpha peak freq
#alpha_cf.Full3 <- lm(alpha_cf ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
raw_alpha_cf.Full2 <- lm(raw_alpha_cf ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), data=alldata)
summary(raw_alpha_cf.Full2 ) 
# qqnorm(resid(alpha_cf.Full2))
# qqline(resid(alpha_cf.Full2))

# Bayes with BRMS
if (!file.exists('mod_rawalphacf_2.RData') || overwrite){
  print('UNFILETRED ALPHA CENTER FREQUENCY  model')
  raw_alpha_cf.mod <- brm(raw_alpha_cf ~ (group+age.centerd+sex+thickz)^2+I(age.centerd^2), 
                      data=alldata,
                      iter=20000, cores=4, save_pars = save_pars(all = TRUE), sample_prior="yes")
  ## Save
  setwd(outdir)
  save(raw_alpha_cf.mod, file='mod_rawalphacf_2.RData')
}

# Model comparison with bridge sampling
if (!file.exists('all_alphacf_mod.RData') || overwrite){
  print('ALPHA CENTER FREQUENCY step-models')
  raw_alpha_cf.0   <- brm(alpha_cf ~ 1, data=alldata,  iter=20000, cores=4, save_pars = save_pars(all = TRUE))
  raw_alpha_cf.G   <- update(raw_alpha_cf.0, ~. +group, newdata=alldata, cores=4)
  raw_alpha_cf.A   <- update(raw_alpha_cf.G, ~. +age.centerd, newdata=alldata, cores=4)
  raw_alpha_cf.S   <- update(raw_alpha_cf.A, ~. +sex, newdata=alldata, cores=4)
  raw_alpha_cf.T   <- update(raw_alpha_cf.S, ~. +thickz, newdata=alldata, cores=4)
  raw_alpha_cf.A2  <- update(raw_alpha_cf.T, ~. +I(age.centerd^2), newdata=alldata, cores=4)
  raw_alpha_cf.GA  <- update(raw_alpha_cf.A2, ~. +group:age.centerd, newdata=alldata, cores=4)
  raw_alpha_cf.GS  <- update(raw_alpha_cf.GA, ~. +group:sex, newdata=alldata, cores=4)
  raw_alpha_cf.GT  <- update(raw_alpha_cf.GS, ~. +group:thickz, newdata=alldata, cores=4)
  raw_alpha_cf.SA  <- update(raw_alpha_cf.GT, ~. +sex:age.centerd, newdata=alldata, cores=4)
  raw_alpha_cf.AT  <- update(raw_alpha_cf.SA, ~. +age.centerd:thickz, newdata=alldata, cores=4)
  raw_alpha_cf.ST  <- update(raw_alpha_cf.AT, ~. +sex:thickz, newdata=alldata, cores=4)
  # alpha_cf.GAS <- update(alpha_cf.ST, ~. +group:age.centerd:sex)
  # alpha_cf.GAT <- update(alpha_cf.GAS, ~. +group:age.centerd:thickz)
  # alpha_cf.GST <- update(alpha_cf.GAT, ~. +group:sex:thickz)
  # alpha_cf.AST <- update(alpha_cf.GST, ~. +age.centerd:sex:thickz)
  
  # Save and clear (to save memory)
  save(list = ls(pattern = "^raw_alpha_cf."), file = 'all_rawalphacf_mod.RData')
  rm(list = ls(pattern = "^raw_alpha_cf"))

  #END
