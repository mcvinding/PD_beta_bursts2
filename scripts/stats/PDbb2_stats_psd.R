# PD beta burst statistics: analsyis of PSD output from FOOOF analysis
#
# Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., 
#  Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor 
#  rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease [Preprint]. 
#  medRxiv.org https://doi.org/10.1101/2021.06.27.21259592
#

library(ggplot2)
library(arm)
library(lmtest)
library(bayestestR)

# Load data
load(file='X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
outdir <- 'X://PD_longrest//output'

######################################################################################
# 1/f intercept
finter.Full3 <- lm(a_intercept ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
qqnorm(resid(finter.Full3))
qqline(resid(finter.Full3))

anova(finter.Full3, test="Chisq")

# LMER regression model (BIC method)
finter.0 <- lm(a_intercept ~ 1, data=alldata)
finter.G   <- update(finter.0, ~. +group)
finter.A   <- update(finter.G, ~. +age.centerd)
finter.S   <- update(finter.A, ~. +sex)
finter.T   <- update(finter.S, ~. +thickz)
finter.A2  <- update(finter.T, ~. +I(age.centerd^2))
finter.GA  <- update(finter.A2, ~. +group:age.centerd)
finter.GS  <- update(finter.GA, ~. +group:sex)
finter.GT  <- update(finter.GS, ~. +group:thickz)
finter.SA  <- update(finter.GT, ~. +sex:age.centerd)
finter.AT  <- update(finter.SA, ~. +age.centerd:thickz)
finter.ST  <- update(finter.AT, ~. +sex:thickz)
finter.GAS <- update(finter.ST, ~. +group:age.centerd:sex)
finter.GAT <- update(finter.GAS, ~. +group:age.centerd:thickz)
finter.GST <- update(finter.GAT, ~. +group:sex:thickz)
finter.AST <- update(finter.GST, ~. +age.centerd:sex:thickz)

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
bayesfactor_models(finter.GAS, denominator = finter.ST)
bayesfactor_models(finter.GAT, denominator = finter.GAS)
bayesfactor_models(finter.GST, denominator = finter.GAT)
bayesfactor_models(finter.AST, denominator = finter.GAT)

# LMER regression model (ANOVA method)
lrtest(finter.0,
       finter.G,
       finter.A,
       finter.S,
       finter.T,
       finter.A2,
       finter.GA,
       finter.GS,
       finter.GT,
       finter.SA,
       finter.AT,
       finter.ST,
       finter.GA,
       finter.GAT,
       finter.GST,
       finter.AST)

# Model summary
set.seed(1)
inter.sim <- sim(finter.Full3, n.sims=1000)
x1 <- coef(finter.Full3)
x2 <- t(apply(coef(inter.sim), 2, quantile, c(0.025, 0.975)))
cbind(x1, x2)

# Group % change
c((x1[2]/x1[1])*100,
  quantile((coef(inter.sim)[,2]/coef(inter.sim)[,1])*100, c(0.025, 0.975)))

# # Group: female
# # x1+x2 - x1 / x1
# c((x1[2]/abs(x1[1]))*100,
#   quantile(((coef(inter.sim)[,2]/abs(coef(inter.sim)[,1])))*100, c(0.025, 0.975)))
# 
# # Group: male
# # x1+x2+x4+x7 - x1+x4 / 
# c(((x1[2]+x1[7])/abs((x1[1]+x1[4])))*100,
#   quantile((coef(inter.sim)[,2]+coef(inter.sim)[,7])/abs(coef(inter.sim)[,2]+coef(inter.sim)[,1])*100, c(0.025, 0.975)))
# 
# # male-female ptns
# # x1+x2+x4+x7 - x1+x2
# c(((x1[4]+x1[7])/abs((x1[1]+x1[2])))*100,
#   quantile((coef(inter.sim)[,4]+coef(inter.sim)[,7])/abs(coef(inter.sim)[,2]+coef(inter.sim)[,1])*100, c(0.025, 0.975)))
# 
# 
# # male-female ctrl
# # x1+x4 - x1 /
# c((x1[4]/abs(x1[1]))*100,
#   quantile(((coef(inter.sim)[,4]/abs(coef(inter.sim)[,1])))*100, c(0.025, 0.975)))

######################################################################################
# 1/f slope
fslope.Full3 <- lm(a_slope ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)
 
qqnorm(resid(fslope.Full3))
qqline(resid(fslope.Full3))

# LMER regression model (BIC method)
fslope.0 <- lm(a_slope ~ 1, data=alldata)
fslope.G   <- update(fslope.0, ~. +group)
fslope.A   <- update(fslope.G, ~. +age.centerd)
fslope.S   <- update(fslope.A, ~. +sex)
fslope.T   <- update(fslope.S, ~. +thickz)
fslope.A2  <- update(fslope.T, ~. +I(age.centerd^2))
fslope.GA  <- update(fslope.A2, ~. +group:age.centerd)
fslope.GS  <- update(fslope.GA, ~. +group:sex)
fslope.GT  <- update(fslope.GS, ~. +group:thickz)
fslope.SA  <- update(fslope.GT, ~. +sex:age.centerd)
fslope.AT  <- update(fslope.SA, ~. +age.centerd:thickz)
fslope.ST  <- update(fslope.AT, ~. +sex:thickz)
fslope.GAS <- update(fslope.ST, ~. +group:age.centerd:sex)
fslope.GAT <- update(fslope.GAS, ~. +group:age.centerd:thickz)
fslope.GST <- update(fslope.GAT, ~. +group:sex:thickz)
fslope.AST <- update(fslope.GST, ~. +age.centerd:sex:thickz)

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
bayesfactor_models(fslope.GAS, denominator = fslope.ST)
bayesfactor_models(fslope.GAT, denominator = fslope.GAS)
bayesfactor_models(fslope.GST, denominator = fslope.GAT)
bayesfactor_models(fslope.AST, denominator = fslope.GAT)

# LMER regression model (ANOVA method)
lrtest(fslope.0,
       fslope.G,
       fslope.A,
       fslope.S,
       fslope.T,
       fslope.A2,
       fslope.GA,
       fslope.GS,
       fslope.GT,
       fslope.SA,
       fslope.AT,
       fslope.ST,
       fslope.GAS,
       fslope.GAT,
       fslope.GST,
       fslope.AST)

# Model summary
slope.sim <- sim(fslope.Full3, n.sims=1000)
x1 <- coef(fslope.Full3)
x2 <- t(apply(coef(mod.sim), 2, quantile, c(0.025, 0.975)))
round(cbind(x1, x2), digits=3)

# Group
c((x1[2]/x1[1])*100,
  quantile((coef(slope.sim)[,2]/coef(slope.sim)[,1])*100, c(0.025, 0.975)))
# 
# # Thickness
# c(((x1[5]/x1[1]))*100,
#   quantile(((coef(mod.sim)[,5]/coef(mod.sim)[,1]))*100, c(0.025, 0.975)))

######################################################################################
# Beta power
beta_pw.Full3 <- lm(beta_pw ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)

qqnorm(resid(beta_pw.Full3))
qqline(resid(beta_pw.Full3))

# LMER regression model (BIC method)
beta_pw.0 <- lm(beta_pw ~ 1, data=alldata)
beta_pw.G   <- update(beta_pw.0, ~. +group)
beta_pw.A   <- update(beta_pw.G, ~. +age.centerd)
beta_pw.S   <- update(beta_pw.A, ~. +sex)
beta_pw.T   <- update(beta_pw.S, ~. +thickz)
beta_pw.A2  <- update(beta_pw.T, ~. +I(age.centerd^2))
beta_pw.GA  <- update(beta_pw.A2, ~. +group:age.centerd)
beta_pw.GS  <- update(beta_pw.GA, ~. +group:sex)
beta_pw.GT  <- update(beta_pw.GS, ~. +group:thickz)
beta_pw.SA  <- update(beta_pw.GT, ~. +sex:age.centerd)
beta_pw.AT  <- update(beta_pw.SA, ~. +age.centerd:thickz)
beta_pw.ST  <- update(beta_pw.AT, ~. +sex:thickz)
beta_pw.GAS <- update(beta_pw.ST, ~. +group:age.centerd:sex)
beta_pw.GAT <- update(beta_pw.GAS, ~. +group:age.centerd:thickz)
beta_pw.GST <- update(beta_pw.GAT, ~. +group:sex:thickz)
beta_pw.AST <- update(beta_pw.GST, ~. +age.centerd:sex:thickz)

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
bayesfactor_models(beta_pw.GAS, denominator = beta_pw.ST)
bayesfactor_models(beta_pw.GAT, denominator = beta_pw.GAS)
bayesfactor_models(beta_pw.GST, denominator = beta_pw.GAT)
bayesfactor_models(beta_pw.AST, denominator = beta_pw.GAT)

# LMER regression model (ANOVA method)
lrtest(beta_pw.0,
       beta_pw.G,
       beta_pw.A,
       beta_pw.S,
       beta_pw.T,
       beta_pw.A2,
       beta_pw.GA,
       beta_pw.GS,
       beta_pw.GT,
       beta_pw.SA,
       beta_pw.AT,
       beta_pw.ST,
       beta_pw.GAS,
       beta_pw.GAT,
       beta_pw.GST,
       beta_pw.AST)

# Model summary
beta_pw.sim <- sim(beta_pw.Full3, n.sims=1000)
x1 <- coef(beta_pw.Full3)
x2 <- t(apply(coef(beta_pw.sim), 2, quantile, c(0.025, 0.975)))
round(cbind(x1, x2), digits=3)

c(((x1[2]/x1[1]))*100,
  quantile(((coef(beta_pw.sim)[,2]/coef(beta_pw.sim)[,1]))*100, c(0.025, 0.975)))

((x1[3]/x1[1]))*100
quantile(((coef(beta_pw.sim)[,3]/coef(beta_pw.sim)[,1]))*100, c(0.025, 0.975))


######################################################################################
# Beta peak freq
beta_cf.Full3 <- lm(beta_cf ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)

qqnorm(resid(beta_cf.Full3))
qqline(resid(beta_cf.Full3))

# LMER regression model (BIC method)
beta_cf.0 <- lm(beta_cf ~ 1, data=alldata)
beta_cf.G   <- update(beta_cf.0, ~. +group)
beta_cf.A   <- update(beta_cf.G, ~. +age.centerd)
beta_cf.S   <- update(beta_cf.A, ~. +sex)
beta_cf.T   <- update(beta_cf.S, ~. +thickz)
beta_cf.A2  <- update(beta_cf.T, ~. +I(age.centerd^2))
beta_cf.GA  <- update(beta_cf.A2, ~. +group:age.centerd)
beta_cf.GS  <- update(beta_cf.GA, ~. +group:sex)
beta_cf.GT  <- update(beta_cf.GS, ~. +group:thickz)
beta_cf.SA  <- update(beta_cf.GT, ~. +sex:age.centerd)
beta_cf.AT  <- update(beta_cf.SA, ~. +age.centerd:thickz)
beta_cf.ST  <- update(beta_cf.AT, ~. +sex:thickz)
beta_cf.GAS <- update(beta_cf.ST, ~. +group:age.centerd:sex)
beta_cf.GAT <- update(beta_cf.GAS, ~. +group:age.centerd:thickz)
beta_cf.GST <- update(beta_cf.GAT, ~. +group:sex:thickz)
beta_cf.AST <- update(beta_cf.GST, ~. +age.centerd:sex:thickz)

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
bayesfactor_models(beta_cf.GAS, denominator = beta_cf.ST)
bayesfactor_models(beta_cf.GAT, denominator = beta_cf.GAS)
bayesfactor_models(beta_cf.GST, denominator = beta_cf.GAT)
bayesfactor_models(beta_cf.AST, denominator = beta_cf.GAT)

lrtest(beta_cf.0,
       beta_cf.G,
       beta_cf.A,
       beta_cf.S,
       beta_cf.T,
       beta_cf.A2,
       beta_cf.GA,
       beta_cf.GS,
       beta_cf.GT,
       beta_cf.SA,
       beta_cf.AT,
       beta_cf.ST,
       beta_cf.GAS,
       beta_cf.GAT,
       beta_cf.GST,
       beta_cf.AST)

# Model summary
beta_cf.sim <- sim(beta_cf.Full3, n.sims=1000)
x1 <- coef(beta_cf.Full3)
x2 <- t(apply(coef(beta_cf.sim), 2, quantile, c(0.025, 0.975)))
round(cbind(x1, x2), digits=3)

######################################################################################
# Alpha power
alpha_pw.Full3 <- lm(alpha_pw ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)

qqnorm(resid(alpha_pw.Full3))
qqline(resid(alpha_pw.Full3))

# LMER regression model (BIC method)
alpha_pw.0 <- lm(alpha_pw ~ 1, data=alldata)
alpha_pw.G   <- update(alpha_pw.0, ~. +group)
alpha_pw.A   <- update(alpha_pw.G, ~. +age.centerd)
alpha_pw.S   <- update(alpha_pw.A, ~. +sex)
alpha_pw.T   <- update(alpha_pw.S, ~. +thickz)
alpha_pw.A2  <- update(alpha_pw.T, ~. +I(age.centerd^2))
alpha_pw.GA  <- update(alpha_pw.A2, ~. +group:age.centerd)
alpha_pw.GS  <- update(alpha_pw.GA, ~. +group:sex)
alpha_pw.GT  <- update(alpha_pw.GS, ~. +group:thickz)
alpha_pw.SA  <- update(alpha_pw.GT, ~. +sex:age.centerd)
alpha_pw.AT  <- update(alpha_pw.SA, ~. +age.centerd:thickz)
alpha_pw.ST  <- update(alpha_pw.AT, ~. +sex:thickz)
alpha_pw.GAS <- update(alpha_pw.ST, ~. +group:age.centerd:sex)
alpha_pw.GAT <- update(alpha_pw.GAS, ~. +group:age.centerd:thickz)
alpha_pw.GST <- update(alpha_pw.GAT, ~. +group:sex:thickz)
alpha_pw.AST <- update(alpha_pw.GST, ~. +age.centerd:sex:thickz)

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
bayesfactor_models(alpha_pw.GAS, denominator = alpha_pw.ST)
bayesfactor_models(alpha_pw.GAT, denominator = alpha_pw.GAS)
bayesfactor_models(alpha_pw.GST, denominator = alpha_pw.GAT)
bayesfactor_models(alpha_pw.AST, denominator = alpha_pw.GAT)

lrtest(alpha_pw.0,
       alpha_pw.G,
       alpha_pw.A,
       alpha_pw.S,
       alpha_pw.T,
       alpha_pw.A2,
       alpha_pw.GA,
       alpha_pw.GS,
       alpha_pw.GT,
       alpha_pw.SA,
       alpha_pw.AT,
       alpha_pw.ST,
       alpha_pw.GAS,
       alpha_pw.GAT,
       alpha_pw.GST,
       alpha_pw.AST)

# Model summary
alpha_pw.sim <- sim(alpha_pw.Full3, n.sims=1000)
x1 <- coef(alpha_pw.Full3)
x2 <- t(apply(coef(alpha_pw.sim), 2, quantile, c(0.025, 0.975)))
round(cbind(x1, x2), digits=3)

######################################################################################
# Alpha peak freq
alpha_cf.Full3 <- lm(alpha_cf ~ (group+age.centerd+sex+thickz)^3+I(age.centerd^2), data=alldata)

qqnorm(resid(alpha_cf.Full3))
qqline(resid(alpha_cf.Full3))

# LMER regression model (BIC method)
alpha_cf.0 <- lm(alpha_cf ~ 1, data=alldata)
alpha_cf.G   <- update(alpha_cf.0, ~. +group)
alpha_cf.A   <- update(alpha_cf.G, ~. +age.centerd)
alpha_cf.S   <- update(alpha_cf.A, ~. +sex)
alpha_cf.T   <- update(alpha_cf.S, ~. +thickz)
alpha_cf.A2  <- update(alpha_cf.T, ~. +I(age.centerd^2))
alpha_cf.GA  <- update(alpha_cf.A2, ~. +group:age.centerd)
alpha_cf.GS  <- update(alpha_cf.GA, ~. +group:sex)
alpha_cf.GT  <- update(alpha_cf.GS, ~. +group:thickz)
alpha_cf.SA  <- update(alpha_cf.GT, ~. +sex:age.centerd)
alpha_cf.AT  <- update(alpha_cf.SA, ~. +age.centerd:thickz)
alpha_cf.ST  <- update(alpha_cf.AT, ~. +sex:thickz)
alpha_cf.GAS <- update(alpha_cf.ST, ~. +group:age.centerd:sex)
alpha_cf.GAT <- update(alpha_cf.GAS, ~. +group:age.centerd:thickz)
alpha_cf.GST <- update(alpha_cf.GAT, ~. +group:sex:thickz)
alpha_cf.AST <- update(alpha_cf.GST, ~. +age.centerd:sex:thickz)

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
bayesfactor_models(alpha_cf.GAS, denominator = alpha_cf.ST)
bayesfactor_models(alpha_cf.GAT, denominator = alpha_cf.GAS)
bayesfactor_models(alpha_cf.GST, denominator = alpha_cf.GAT)
bayesfactor_models(alpha_cf.AST, denominator = alpha_cf.GAT)

lrtest(alpha_cf.0,
       alpha_cf.G,
       alpha_cf.A,
       alpha_cf.S,
       alpha_cf.T,
       alpha_cf.A2,
       alpha_cf.GA,
       alpha_cf.GS,
       alpha_cf.GT,
       alpha_cf.SA,
       alpha_cf.AT,
       alpha_cf.ST,
       alpha_cf.GAS,
       alpha_cf.GAT,
       alpha_cf.GST,
       alpha_cf.AST)

# Model summary
alpha_cf.sim <- sim(alpha_cf.Full3, n.sims=1000)
x1 <- coef(alpha_cf.Full3)
x2 <- t(apply(coef(alpha_cf.sim), 2, quantile, c(0.025, 0.975)))
round(cbind(x1, x2), digits=3)

# age: ctrl
c(x1[3],
  quantile(coef(alpha_cf.sim)[,3], c(0.025, 0.975)))

# age: Ptns
c(x1[3]+x1[6],
  quantile(coef(alpha_cf.sim)[,3]+coef(alpha_cf.sim)[,6], c(0.025, 0.975)))

# thick: ctrl
c(x1[5],
  quantile(coef(alpha_cf.sim)[,5], c(0.025, 0.975)))

# thcik: Ptns
c(x1[5]+x1[8],
  quantile(coef(alpha_cf.sim)[,5]+coef(alpha_cf.sim)[,8], c(0.025, 0.975)))

## Save models
setwd(outdir)
save(finter.Full3, file='mod_finter.RData')
save(fslope.Full3, file='mod_fslope.RData')
save(beta_pw.Full3, file='mod_betapw.RData')
save(beta_cf.Full3, file='mod_betacf.RData')
save(alpha_pw.Full3, file='mod_alphapw.RData')
save(alpha_cf.Full3, file='mod_alphacf.RData')

#END