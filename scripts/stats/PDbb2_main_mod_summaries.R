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
library(brms)

outdir <- '/home/mikkel/PD_longrest/output'
setwd(outdir)

################################################################################
### PSD MODELS 
################################################################################

######################################################################################
# 1/f intercept
load(file='all_finter_mod.RData')
load(file='mod_finter_2.RData')

# Model summary
fixef(finter.mod)
summary(finter.mod)

# ## BF test
bayesfactor_models(finter.G, denominator   = finter.0)
bayesfactor_models(finter.A, denominator   = finter.G)
bayesfactor_models(finter.S, denominator   = finter.A)
bayesfactor_models(finter.T, denominator   = finter.S)
bayesfactor_models(finter.A2, denominator  = finter.T)
bayesfactor_models(finter.GA, denominator  = finter.A2)
bayesfactor_models(finter.GS, denominator  = finter.GA)
bayesfactor_models(finter.GT, denominator  = finter.GS)
bayesfactor_models(finter.SA, denominator  = finter.GT)
bayesfactor_models(finter.AT, denominator  = finter.SA)
bayesfactor_models(finter.ST, denominator  = finter.AT)

# Hypothesis testing
finter.tst <- hypothesis(finter.mod  , c("grouppatient>0",
                                         "age.centerd>0",
                                         "sexM<0",
                                         "thickz<0",
                                         "Iage.centerdE2>0",
                                         "grouppatient:age.centerd<0",
                                         "grouppatient:sexM>0",
                                         "grouppatient:thickz>0",
                                         "age.centerd:sexM<0",
                                         "age.centerd:thickz>0",
                                         "sexM:thickz>0"), alpha=0.05)
finter.tst$hypothesis$Post.Prob2 <- ifelse(finter.tst$hypothesis$Post.Prob>0.5, (1-finter.tst$hypothesis$Post.Prob)*2, finter.tst$hypothesis$Post.Prob*2)
finter.tst$hypothesis$Post.Prob2

rm(list = ls(pattern = "^finter"))

######################################################################################
# 1/f slope
load(file='all_fslope_mod.RData')
load(file='mod_fslope_2.RData')

# Model summary
fixef(fslope.mod)
summary(fslope.mod)

# ## BF test
bayesfactor_models(fslope.G, denominator   = fslope.0)
bayesfactor_models(fslope.A, denominator   = fslope.G)
bayesfactor_models(fslope.S, denominator   = fslope.A)
bayesfactor_models(fslope.T, denominator   = fslope.S)
bayesfactor_models(fslope.A2, denominator  = fslope.T)
bayesfactor_models(fslope.GA, denominator  = fslope.A2)
bayesfactor_models(fslope.GS, denominator  = fslope.GA)
bayesfactor_models(fslope.GT, denominator  = fslope.GS)
bayesfactor_models(fslope.SA, denominator  = fslope.GT)
bayesfactor_models(fslope.AT, denominator  = fslope.SA)
bayesfactor_models(fslope.ST, denominator  = fslope.AT)
# bayesfactor_models(fslope.GAS, denominator = fslope.ST)
# bayesfactor_models(fslope.GAT, denominator = fslope.GAS)
# bayesfactor_models(fslope.GST, denominator = fslope.GAT)
# bayesfactor_models(fslope.AST, denominator = fslope.GAT)

# Hypothesis testing
fslope.tst <- hypothesis(fslope.mod  , c("grouppatient>0",
                                         "age.centerd>0",
                                         "sexM<0",
                                         "thickz<0",
                                         "Iage.centerdE2>0",
                                         "grouppatient:age.centerd<0",
                                         "grouppatient:sexM>0",
                                         "grouppatient:thickz>0",
                                         "age.centerd:sexM<0",
                                         "age.centerd:thickz>0",
                                         "sexM:thickz>0"), alpha=0.05)
fslope.tst$hypothesis$Post.Prob2 <- ifelse(fslope.tst$hypothesis$Post.Prob>0.5, (1-fslope.tst$hypothesis$Post.Prob)*2, fslope.tst$hypothesis$Post.Prob*2)
fslope.tst$hypothesis$Post.Prob2

rm(list = ls(pattern = "^fslope"))

######################################################################################
# Beta power
load(file='all_betapw_mod.RData')
load(file='mod_betapw_2.RData')

# Model summary
fixef(beta_pw.mod)
summary(beta_pw.mod)

# ## BF test
bayesfactor_models(beta_pw.G, denominator   = beta_pw.0)
bayesfactor_models(beta_pw.A, denominator   = beta_pw.G)
bayesfactor_models(beta_pw.S, denominator   = beta_pw.A)
bayesfactor_models(beta_pw.T, denominator   = beta_pw.S)
bayesfactor_models(beta_pw.A2, denominator  = beta_pw.T)
bayesfactor_models(beta_pw.GA, denominator  = beta_pw.A2)
bayesfactor_models(beta_pw.GS, denominator  = beta_pw.GA)
bayesfactor_models(beta_pw.GT, denominator  = beta_pw.GS)
bayesfactor_models(beta_pw.SA, denominator  = beta_pw.GT)
bayesfactor_models(beta_pw.AT, denominator  = beta_pw.SA)
bayesfactor_models(beta_pw.ST, denominator  = beta_pw.AT)
# bayesfactor_models(beta_pw.GAS, denominator = beta_pw.ST)
# bayesfactor_models(beta_pw.GAT, denominator = beta_pw.GAS)
# bayesfactor_models(beta_pw.GST, denominator = beta_pw.GAT)
# bayesfactor_models(beta_pw.AST, denominator = beta_pw.GAT)

# Hypothesis testing
beta_pw.tst <- hypothesis(beta_pw.mod  , c("grouppatient>0",
                                           "age.centerd>0",
                                           "sexM<0",
                                           "thickz<0",
                                           "Iage.centerdE2>0",
                                           "grouppatient:age.centerd<0",
                                           "grouppatient:sexM>0",
                                           "grouppatient:thickz>0",
                                           "age.centerd:sexM<0",
                                           "age.centerd:thickz>0",
                                           "sexM:thickz>0"), alpha=0.05)
beta_pw.tst$hypothesis$Post.Prob2 <- ifelse(beta_pw.tst$hypothesis$Post.Prob>0.5, (1-beta_pw.tst$hypothesis$Post.Prob)*2, beta_pw.tst$hypothesis$Post.Prob*2)
beta_pw.tst$hypothesis$Post.Prob2

rm(list = ls(pattern = "^beta_pw"))

######################################################################################
# Beta peak freq
load(file='all_betacf_mod.RData')
load(file='mod_betacf_2.RData')

# Model summary
fixef(beta_cf.mod)
summary(beta_cf.mod)

# ## BF test
bayesfactor_models(beta_cf.G, denominator   = beta_cf.0)
bayesfactor_models(beta_cf.A, denominator   = beta_cf.G)
bayesfactor_models(beta_cf.S, denominator   = beta_cf.A)
bayesfactor_models(beta_cf.T, denominator   = beta_cf.S)
bayesfactor_models(beta_cf.A2, denominator  = beta_cf.T)
bayesfactor_models(beta_cf.GA, denominator  = beta_cf.A2)
bayesfactor_models(beta_cf.GS, denominator  = beta_cf.GA)
bayesfactor_models(beta_cf.GT, denominator  = beta_cf.GS)
bayesfactor_models(beta_cf.SA, denominator  = beta_cf.GT)
bayesfactor_models(beta_cf.AT, denominator  = beta_cf.SA)
bayesfactor_models(beta_cf.ST, denominator  = beta_cf.AT)
# bayesfactor_models(beta_cf.GAS, denominator = beta_cf.ST)
# bayesfactor_models(beta_cf.GAT, denominator = beta_cf.GAS)
# bayesfactor_models(beta_cf.GST, denominator = beta_cf.GAT)
# bayesfactor_models(beta_cf.AST, denominator = beta_cf.GAT)

# Hypothesis testing
beta_cf.tst <- hypothesis(beta_cf.mod  , c("grouppatient>0",
                                           "age.centerd>0",
                                           "sexM<0",
                                           "thickz<0",
                                           "Iage.centerdE2>0",
                                           "grouppatient:age.centerd<0",
                                           "grouppatient:sexM>0",
                                           "grouppatient:thickz>0",
                                           "age.centerd:sexM<0",
                                           "age.centerd:thickz>0",
                                           "sexM:thickz>0"), alpha=0.05)
beta_cf.tst$hypothesis$Post.Prob2 <- ifelse(beta_cf.tst$hypothesis$Post.Prob>0.5, (1-beta_cf.tst$hypothesis$Post.Prob)*2, beta_cf.tst$hypothesis$Post.Prob*2)
beta_cf.tst$hypothesis$Post.Prob2

rm(list = ls(pattern = "^beta_cf"))

######################################################################################
# Alpha power
load(file='all_alphapw_mod.RData')
load(file='mod_alphapw_2.RData')

# Model summary
fixef(alpha_pw.mod)
summary(alpha_pw.mod)

# ## BF test
bayesfactor_models(alpha_pw.G, denominator   = alpha_pw.0)
bayesfactor_models(alpha_pw.A, denominator   = alpha_pw.G)
bayesfactor_models(alpha_pw.S, denominator   = alpha_pw.A)
bayesfactor_models(alpha_pw.T, denominator   = alpha_pw.S)
bayesfactor_models(alpha_pw.A2, denominator  = alpha_pw.T)
bayesfactor_models(alpha_pw.GA, denominator  = alpha_pw.A2)
bayesfactor_models(alpha_pw.GS, denominator  = alpha_pw.GA)
bayesfactor_models(alpha_pw.GT, denominator  = alpha_pw.GS)
bayesfactor_models(alpha_pw.SA, denominator  = alpha_pw.GT)
bayesfactor_models(alpha_pw.AT, denominator  = alpha_pw.SA)
bayesfactor_models(alpha_pw.ST, denominator  = alpha_pw.AT)
# bayesfactor_models(alpha_pw.GAS, denominator = alpha_pw.ST)
# bayesfactor_models(alpha_pw.GAT, denominator = alpha_pw.GAS)
# bayesfactor_models(alpha_pw.GST, denominator = alpha_pw.GAT)
# bayesfactor_models(alpha_pw.AST, denominator = alpha_pw.GAT)

# Hypothesis testing
alpha_pw.tst <- hypothesis(alpha_pw.mod  , c("grouppatient>0",
                                             "age.centerd>0",
                                             "sexM<0",
                                             "thickz<0",
                                             "Iage.centerdE2>0",
                                             "grouppatient:age.centerd<0",
                                             "grouppatient:sexM>0",
                                             "grouppatient:thickz>0",
                                             "age.centerd:sexM<0",
                                             "age.centerd:thickz>0",
                                             "sexM:thickz>0"), alpha=0.05)
alpha_pw.tst$hypothesis$Post.Prob2 <- ifelse(alpha_pw.tst$hypothesis$Post.Prob>0.5, (1-alpha_pw.tst$hypothesis$Post.Prob)*2, alpha_pw.tst$hypothesis$Post.Prob*2)
alpha_pw.tst$hypothesis$Post.Prob2

rm(list = ls(pattern = "^alpha_pw"))

######################################################################################
# Alpha peak freq
load(file='all_alphacf_mod.RData')
load(file='mod_alphacf_2.RData')

# Model summary
fixef(alpha_cf.mod)
summary(alpha_cf.mod)

# ## BF test
bayesfactor_models(alpha_cf.G, denominator   = alpha_cf.0)
bayesfactor_models(alpha_cf.A, denominator   = alpha_cf.G)
bayesfactor_models(alpha_cf.S, denominator   = alpha_cf.A)
bayesfactor_models(alpha_cf.T, denominator   = alpha_cf.S)
bayesfactor_models(alpha_cf.A2, denominator  = alpha_cf.T)
bayesfactor_models(alpha_cf.GA, denominator  = alpha_cf.A2)
bayesfactor_models(alpha_cf.GS, denominator  = alpha_cf.GA)
bayesfactor_models(alpha_cf.GT, denominator  = alpha_cf.GS)
bayesfactor_models(alpha_cf.SA, denominator  = alpha_cf.GT)
bayesfactor_models(alpha_cf.AT, denominator  = alpha_cf.SA)
bayesfactor_models(alpha_cf.ST, denominator  = alpha_cf.AT)
# bayesfactor_models(alpha_cf.GAS, denominator = alpha_cf.ST)
# bayesfactor_models(alpha_cf.GAT, denominator = alpha_cf.GAS)
# bayesfactor_models(alpha_cf.GST, denominator = alpha_cf.GAT)
# bayesfactor_models(alpha_cf.AST, denominator = alpha_cf.GAT)

# Hypothesis testing
alpha_cf.tst <- hypothesis(alpha_cf.mod, c("grouppatient>0",
                                           "age.centerd>0",
                                           "sexM<0",
                                           "thickz<0",
                                           "Iage.centerdE2>0",
                                           "grouppatient:age.centerd<0",
                                           "grouppatient:sexM>0",
                                           "grouppatient:thickz>0",
                                           "age.centerd:sexM<0",
                                           "age.centerd:thickz>0",
                                           "sexM:thickz>0"), alpha=0.05)
alpha_cf.tst$hypothesis$Post.Prob2 <- ifelse(alpha_cf.tst$hypothesis$Post.Prob>0.5, (1-alpha_cf.tst$hypothesis$Post.Prob)*2, alpha_cf.tst$hypothesis$Post.Prob*2)
alpha_cf.tst$hypothesis$Post.Prob2

rm(list = ls(pattern = "^alpha_cf"))

################################################################################
### BURST MODELS
################################################################################

# ##############################################################################
## N events
load(file='all_nevmod.RData')
load(file='mod_neveBF_2.RData')

# Model summary
fixef(mod_neve)
summary(mod_neve)

## BF test
bayesfactor_models(mod.neve.G, denominator = mod.neve.0)
bayesfactor_models(mod.neve.A, denominator = mod.neve.G)
bayesfactor_models(mod.neve.S, denominator = mod.neve.A)
bayesfactor_models(mod.neve.T, denominator = mod.neve.S)
bayesfactor_models(mod.neve.A2, denominator = mod.neve.T)
bayesfactor_models(mod.neve.GA, denominator = mod.neve.A2)
bayesfactor_models(mod.neve.GS, denominator = mod.neve.GA)
bayesfactor_models(mod.neve.GT, denominator = mod.neve.GS)
bayesfactor_models(mod.neve.AS, denominator = mod.neve.GT)
bayesfactor_models(mod.neve.AT, denominator = mod.neve.AS)
bayesfactor_models(mod.neve.ST, denominator = mod.neve.AT)

# Hypothesis testing
neve.tst <- hypothesis(mod_neve, c("grouppatient>0",
                                   "age.centerd>0",
                                   "sexM<0",
                                   "thickz<0",
                                   "Iage.centerdE2>0",
                                   "grouppatient:age.centerd<0",
                                   "grouppatient:sexM>0",
                                   "grouppatient:thickz>0",
                                   "age.centerd:sexM<0",
                                   "age.centerd:thickz>0",
                                   "sexM:thickz>0"), alpha=0.05)
neve.tst$hypothesis$Post.Prob2 <- ifelse(neve.tst$hypothesis$Post.Prob>0.5, (1-neve.tst$hypothesis$Post.Prob)*2, neve.tst$hypothesis$Post.Prob*2)
neve.tst$hypothesis$Post.Prob2

rm(list = ls(pattern = "^mod.neve"))

# ##############################################################################
## EVENT LENGTH
load(file='all_lenmod.RData')
load(file='lenmod_2.RData')

# Model summary
fixef(lenmod)
summary(lenmod)

## BF test
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

# Hypothesis testing
len.tst <- hypothesis(lenmod, c("grouppatient>0",
                                "age.centerd>0",
                                "sexM<0",
                                "thickz<0",
                                "Iage.centerdE2>0",
                                "grouppatient:age.centerd<0",
                                "grouppatient:sexM>0",
                                "grouppatient:thickz>0",
                                "age.centerd:sexM<0",
                                "age.centerd:thickz>0",
                                "sexM:thickz>0"), alpha=0.05)
len.tst$hypothesis$Post.Prob2 <- ifelse(len.tst$hypothesis$Post.Prob>0.5, (1-len.tst$hypothesis$Post.Prob)*2, len.tst$hypothesis$Post.Prob*2)
len.tst$hypothesis$Post.Prob2

rm(list = ls(pattern = "^lenmod"))

######################################################################################
## TIME UNTIL EVENT mixed model
load(file='all_tuemod.RData')
load(file='tuemod_2.RData')

# Model summary
fixef(tuemod)
summary(tuemod)

## BF test
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

# Hypothesis testing
tue.tst <- hypothesis(tuemod, c("grouppatient>0",
                                "age.centerd>0",
                                "sexM<0",
                                "thickz<0",
                                "Iage.centerdE2>0",
                                "grouppatient:age.centerd<0",
                                "grouppatient:sexM>0",
                                "grouppatient:thickz>0",
                                "age.centerd:sexM<0",
                                "age.centerd:thickz>0",
                                "sexM:thickz>0"), alpha=0.05)
tue.tst$hypothesis$Post.Prob2 <- ifelse(tue.tst$hypothesis$Post.Prob>0.5, (1-tue.tst$hypothesis$Post.Prob)*2, tue.tst$hypothesis$Post.Prob*2)
tue.tst$hypothesis$Post.Prob2

rm(list = ls(pattern = "^tuemod"))

######################################################################################
## Burst max amplitude mixed model
load(file='all_maxmod.RData')
load(file='maxmod_2.RData')

# Model summary
fixef(maxmod)
summary(maxmod)

## BF test
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

# Hypothesis testing
max.tst <- hypothesis(maxmod, c("grouppatient>0",
                                "age.centerd>0",
                                "sexM<0",
                                "thickz<0",
                                "Iage.centerdE2>0",
                                "grouppatient:age.centerd<0",
                                "grouppatient:sexM>0",
                                "grouppatient:thickz>0",
                                "age.centerd:sexM<0",
                                "age.centerd:thickz>0",
                                "sexM:thickz>0"), alpha=0.05)
max.tst$hypothesis$Post.Prob2 <- ifelse(max.tst$hypothesis$Post.Prob>0.5, (1-max.tst$hypothesis$Post.Prob)*2, max.tst$hypothesis$Post.Prob*2)
max.tst$hypothesis$Post.Prob2

rm(list = ls(pattern = "^maxmod"))

#END