# PD beta burst statistics: analsyis of PSD output from FOOOF analysis
library(ggplot2)
library(arm)
library(lmtest)

# Load data
load(file='X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
outdir <- 'X://PD_longrest//output'

# Inspect
ggplot(aes(x=age, y=a_intercept, color=group, shape=sex), data=alldata)+
  geom_point()+geom_smooth(method=lm)+
  theme_bw()+ggtitle('1/f intercept ~ age')

ggplot(aes(x=age, y=a_slope, color=group, shape=sex), data=alldata)+
  geom_point()+geom_smooth(method=lm)+
  theme_bw()+ggtitle('1/f slope ~ age')

ggplot(aes(x=age, y=beta_pw, color=group, shape=sex), data=alldata)+
  geom_point()+geom_smooth(method=lm)+
  theme_bw()+ggtitle('Beta power ~ age')

ggplot(aes(x=age, y=beta_cf, color=group, shape=sex), data=alldata)+
  geom_point()+geom_smooth(method=lm)+
  theme_bw()+ggtitle('Beta peak freq ~ age')

ggplot(aes(x=age, y=log(alpha_pw), color=group, shape=sex), data=alldata)+
  geom_point()+geom_smooth(method=lm)+
  theme_bw()+ggtitle('Alpha power ~ age')

ggplot(aes(x=age, y=alpha_cf, color=group, shape=sex), data=alldata)+
  geom_point()+geom_smooth(method=lm)+
  theme_bw()+ggtitle('Alpha peak freq ~ age')

######################################################################################
# 1/f intercept
finter.Full3 <- lm(a_intercept ~ (group+age.centerd+sex+thickz)^3, data=alldata)

qqnorm(resid(finter.Full3))
qqline(resid(finter.Full3))

# anova(finter.Full3, test="Chisq")

finter.AST <- update(finter.Full3, ~. -age.centerd:sex:thickz)
finter.GST <- update(finter.AST,   ~. -group:sex:thickz)
finter.GAT <- update(finter.GST,   ~. -group:age.centerd:thickz)
finter.GSA <- update(finter.GAT,   ~. -group:sex:age.centerd)
finter.ST  <- update(finter.GSA,   ~. -sex:thickz)
finter.AT  <- update(finter.ST,    ~. -age.centerd:thickz)
finter.SA  <- update(finter.AT,    ~. -sex:age.centerd)
finter.GT  <- update(finter.SA,    ~. -group:thickz)
finter.GS  <- update(finter.GT,    ~. -group:sex)
finter.GA  <- update(finter.GS,    ~. -group:age.centerd)
finter.T   <- update(finter.GA,    ~. -thickz)
finter.S   <- update(finter.T,     ~. -sex)
finter.A   <- update(finter.S,     ~. -age.centerd)
finter.G   <- update(finter.A,     ~. -group)
    
lrtest(finter.G,
       finter.A,
       finter.S,
       finter.T,
       finter.GA,
       finter.GS,
       finter.GT,
       finter.SA,
       finter.AT,
       finter.ST,
       finter.GSA,
       finter.GAT,
       finter.GST,
       finter.AST, 
       finter.Full3)

# Model summary
set.seed(1)
inter.sim <- sim(finter.Full3, n.sims=1000)
x1 <- coef(finter.Full3)
x2 <- t(apply(coef(inter.sim), 2, quantile, c(0.025, 0.975)))
cbind(x1, x2)

# Group: female
# x1+x2 - x1 / x1
c((x1[2]/abs(x1[1]))*100,
  quantile(((coef(inter.sim)[,2]/abs(coef(inter.sim)[,1])))*100, c(0.025, 0.975)))

# Group: male
# x1+x2+x4+x7 - x1+x4 / 
c(((x1[2]+x1[7])/abs((x1[1]+x1[4])))*100,
  quantile((coef(inter.sim)[,2]+coef(inter.sim)[,7])/abs(coef(inter.sim)[,2]+coef(inter.sim)[,1])*100, c(0.025, 0.975)))

# male-female ptns
# x1+x2+x4+x7 - x1+x2
c(((x1[4]+x1[7])/abs((x1[1]+x1[2])))*100,
  quantile((coef(inter.sim)[,4]+coef(inter.sim)[,7])/abs(coef(inter.sim)[,2]+coef(inter.sim)[,1])*100, c(0.025, 0.975)))


# male-female ctrl
# x1+x4 - x1 /
c((x1[4]/abs(x1[1]))*100,
  quantile(((coef(inter.sim)[,4]/abs(coef(inter.sim)[,1])))*100, c(0.025, 0.975)))

######################################################################################
# 1/f slope
fslope.Full3 <- lm(a_slope ~ (group+age.centerd+sex+thickz)^3, data=alldata)

qqnorm(resid(fslope.Full3))
qqline(resid(fslope.Full3))

fslope.AST <- update(fslope.Full3, ~. -age.centerd:sex:thickz)
fslope.GST <- update(fslope.AST,   ~. -group:sex:thickz)
fslope.GAT <- update(fslope.GST,   ~. -group:age.centerd:thickz)
fslope.GSA <- update(fslope.GAT,   ~. -group:sex:age.centerd)
fslope.ST  <- update(fslope.GSA,   ~. -sex:thickz)
fslope.AT  <- update(fslope.ST,    ~. -age.centerd:thickz)
fslope.SA  <- update(fslope.AT,    ~. -sex:age.centerd)
fslope.GT  <- update(fslope.SA,    ~. -group:thickz)
fslope.GS  <- update(fslope.GT,    ~. -group:sex)
fslope.GA  <- update(fslope.GS,    ~. -group:age.centerd)
fslope.T   <- update(fslope.GA,    ~. -thickz)
fslope.S   <- update(fslope.T,     ~. -sex)
fslope.A   <- update(fslope.S,     ~. -age.centerd)
fslope.G   <- update(fslope.A,     ~. -group)

lrtest(fslope.G,
       fslope.A,
       fslope.S,
       fslope.T,
       fslope.GA,
       fslope.GS,
       fslope.GT,
       fslope.SA,
       fslope.AT,
       fslope.ST,
       fslope.GSA,
       fslope.GAT,
       fslope.GST,
       fslope.AST, 
       fslope.Full3)

# Model summary
mod.sim <- sim(fslope.Full3, n.sims=1000)
x1 <- coef(fslope.Full3)
x2 <- t(apply(coef(mod.sim), 2, quantile, c(0.025, 0.975)))
round(cbind(x1, x2), digits=3)

# Group
c((x1[2]/x1[1])*100,
  quantile((coef(mod.sim)[,2]/coef(mod.sim)[,1])*100, c(0.025, 0.975)))

# Thickness
c(((x1[5]/x1[1]))*100,
  quantile(((coef(mod.sim)[,5]/coef(mod.sim)[,1]))*100, c(0.025, 0.975)))

######################################################################################
# Beta power
beta_pw.Full3 <- lm(beta_pw ~ (group+age.centerd+sex+thickz)^3, data=alldata)

qqnorm(resid(beta_pw.Full3))
qqline(resid(beta_pw.Full3))

beta_pw.AST <- update(beta_pw.Full3, ~. -age.centerd:sex:thickz)
beta_pw.GST <- update(beta_pw.AST,   ~. -group:sex:thickz)
beta_pw.GAT <- update(beta_pw.GST,   ~. -group:age.centerd:thickz)
beta_pw.GSA <- update(beta_pw.GAT,   ~. -group:sex:age.centerd)
beta_pw.ST  <- update(beta_pw.GSA,   ~. -sex:thickz)
beta_pw.AT  <- update(beta_pw.ST,    ~. -age.centerd:thickz)
beta_pw.SA  <- update(beta_pw.AT,    ~. -sex:age.centerd)
beta_pw.GT  <- update(beta_pw.SA,    ~. -group:thickz)
beta_pw.GS  <- update(beta_pw.GT,    ~. -group:sex)
beta_pw.GA  <- update(beta_pw.GS,    ~. -group:age.centerd)
beta_pw.T   <- update(beta_pw.GA,    ~. -thickz)
beta_pw.S   <- update(beta_pw.T,     ~. -sex)
beta_pw.A   <- update(beta_pw.S,     ~. -age.centerd)
beta_pw.G   <- update(beta_pw.A,     ~. -group)

lrtest(beta_pw.G,
       beta_pw.A,
       beta_pw.S,
       beta_pw.T,
       beta_pw.GA,
       beta_pw.GS,
       beta_pw.GT,
       beta_pw.SA,
       beta_pw.AT,
       beta_pw.ST,
       beta_pw.GSA,
       beta_pw.GAT,
       beta_pw.GST,
       beta_pw.AST, 
       beta_pw.Full3)

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
beta_cf.Full3 <- lm(beta_cf ~ (group+age.centerd+sex+thickz)^3, data=alldata)

qqnorm(resid(beta_cf.Full3))
qqline(resid(beta_cf.Full3))

beta_cf.AST <- update(beta_cf.Full3, ~. -age.centerd:sex:thickz)
beta_cf.GST <- update(beta_cf.AST,   ~. -group:sex:thickz)
beta_cf.GAT <- update(beta_cf.GST,   ~. -group:age.centerd:thickz)
beta_cf.GSA <- update(beta_cf.GAT,   ~. -group:sex:age.centerd)
beta_cf.ST  <- update(beta_cf.GSA,   ~. -sex:thickz)
beta_cf.AT  <- update(beta_cf.ST,    ~. -age.centerd:thickz)
beta_cf.SA  <- update(beta_cf.AT,    ~. -sex:age.centerd)
beta_cf.GT  <- update(beta_cf.SA,    ~. -group:thickz)
beta_cf.GS  <- update(beta_cf.GT,    ~. -group:sex)
beta_cf.GA  <- update(beta_cf.GS,    ~. -group:age.centerd)
beta_cf.T   <- update(beta_cf.GA,    ~. -thickz)
beta_cf.S   <- update(beta_cf.T,     ~. -sex)
beta_cf.A   <- update(beta_cf.S,     ~. -age.centerd)
beta_cf.G   <- update(beta_cf.A,     ~. -group)

lrtest(beta_cf.G,
       beta_cf.A,
       beta_cf.S,
       beta_cf.T,
       beta_cf.GA,
       beta_cf.GS,
       beta_cf.GT,
       beta_cf.SA,
       beta_cf.AT,
       beta_cf.ST,
       beta_cf.GSA,
       beta_cf.GAT,
       beta_cf.GST,
       beta_cf.AST, 
       beta_cf.Full3)

# Model summary
beta_cf.sim <- sim(beta_cf.Full3, n.sims=1000)
x1 <- coef(beta_cf.Full3)
x2 <- t(apply(coef(beta_cf.sim), 2, quantile, c(0.025, 0.975)))
round(cbind(x1, x2), digits=3)

######################################################################################
# Alpha power
alpha_pw.Full3 <- lm(alpha_pw ~ (group+age.centerd+sex+thickz)^3, data=alldata)

qqnorm(resid(alpha_pw.Full3))
qqline(resid(alpha_pw.Full3))

alpha_pw.AST <- update(alpha_pw.Full3, ~. -age.centerd:sex:thickz)
alpha_pw.GST <- update(alpha_pw.AST,   ~. -group:sex:thickz)
alpha_pw.GAT <- update(alpha_pw.GST,   ~. -group:age.centerd:thickz)
alpha_pw.GSA <- update(alpha_pw.GAT,   ~. -group:sex:age.centerd)
alpha_pw.ST  <- update(alpha_pw.GSA,   ~. -sex:thickz)
alpha_pw.AT  <- update(alpha_pw.ST,    ~. -age.centerd:thickz)
alpha_pw.SA  <- update(alpha_pw.AT,    ~. -sex:age.centerd)
alpha_pw.GT  <- update(alpha_pw.SA,    ~. -group:thickz)
alpha_pw.GS  <- update(alpha_pw.GT,    ~. -group:sex)
alpha_pw.GA  <- update(alpha_pw.GS,    ~. -group:age.centerd)
alpha_pw.T   <- update(alpha_pw.GA,    ~. -thickz)
alpha_pw.S   <- update(alpha_pw.T,     ~. -sex)
alpha_pw.A   <- update(alpha_pw.S,     ~. -age.centerd)
alpha_pw.G   <- update(alpha_pw.A,     ~. -group)

lrtest(alpha_pw.G,
       alpha_pw.A,
       alpha_pw.S,
       alpha_pw.T,
       alpha_pw.GA,
       alpha_pw.GS,
       alpha_pw.GT,
       alpha_pw.SA,
       alpha_pw.AT,
       alpha_pw.ST,
       alpha_pw.GSA,
       alpha_pw.GAT,
       alpha_pw.GST,
       alpha_pw.AST, 
       alpha_pw.Full3)

# Model summary
alpha_pw.sim <- sim(alpha_pw.Full3, n.sims=1000)
x1 <- coef(alpha_pw.Full3)
x2 <- t(apply(coef(alpha_pw.sim), 2, quantile, c(0.025, 0.975)))
round(cbind(x1, x2), digits=3)

######################################################################################
# Alpha peak freq
alpha_cf.Full3 <- lm(alpha_cf ~ (group+age.centerd+sex+thickz)^3, data=alldata)

qqnorm(resid(alpha_cf.Full3))
qqline(resid(alpha_cf.Full3))

alpha_cf.AST <- update(alpha_cf.Full3, ~. -age.centerd:sex:thickz)
alpha_cf.GST <- update(alpha_cf.AST,   ~. -group:sex:thickz)
alpha_cf.GAT <- update(alpha_cf.GST,   ~. -group:age.centerd:thickz)
alpha_cf.GSA <- update(alpha_cf.GAT,   ~. -group:sex:age.centerd)
alpha_cf.ST  <- update(alpha_cf.GSA,   ~. -sex:thickz)
alpha_cf.AT  <- update(alpha_cf.ST,    ~. -age.centerd:thickz)
alpha_cf.SA  <- update(alpha_cf.AT,    ~. -sex:age.centerd)
alpha_cf.GT  <- update(alpha_cf.SA,    ~. -group:thickz)
alpha_cf.GS  <- update(alpha_cf.GT,    ~. -group:sex)
alpha_cf.GA  <- update(alpha_cf.GS,    ~. -group:age.centerd)
alpha_cf.T   <- update(alpha_cf.GA,    ~. -thickz)
alpha_cf.S   <- update(alpha_cf.T,     ~. -sex)
alpha_cf.A   <- update(alpha_cf.S,     ~. -age.centerd)
alpha_cf.G   <- update(alpha_cf.A,     ~. -group)

lrtest(alpha_cf.G,
       alpha_cf.A,
       alpha_cf.S,
       alpha_cf.T,
       alpha_cf.GA,
       alpha_cf.GS,
       alpha_cf.GT,
       alpha_cf.SA,
       alpha_cf.AT,
       alpha_cf.ST,
       alpha_cf.GSA,
       alpha_cf.GAT,
       alpha_cf.GST,
       alpha_cf.AST, 
       alpha_cf.Full3)

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