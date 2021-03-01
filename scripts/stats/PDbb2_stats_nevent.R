# PD beta burst statistics: analysis of N events
### CLEAN UP!!!! ###
library(lme4)
library(arm)
library(ggplot2)
library(car)


## Load data
setwd('X://PD_longrest//groupanalysis//')
load('X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
# load('C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj2.Rdata')

## Center variables
alldata$age.centerd <- alldata$age-mean(alldata$age)
alldata$thick.centerd <- alldata$thick-mean(alldata$thick)

# Inspect hist
ggplot( aes(x=nevent.u.m2.min, fill=group), data=alldata) +
  geom_histogram(color="black", alpha=0.6, position = 'identity', bins=25)

# Inspect ~age
ggplot(aes(x=age, y=nevent.u.m2.min, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_smooth(method=lm)

# ######################################################################################
# # LMER regression model
tstmod.Full3 <- glm(nevent.u.m2.min ~ (group+age.centerd+sex+thick.centerd)^3, data=alldata, family=poisson)

tstmod.AST <- update(tstmod.Full3, ~. -age.centerd:sex:thick.centerd)
tstmod.GST <- update(tstmod.AST, ~. -group:sex:thick.centerd)
tstmod.GAT <- update(tstmod.GST, ~. -group:age.centerd:thick.centerd)
tstmod.GSA <- update(tstmod.GAT, ~. -group:sex:age.centerd)
tstmod.ST <- update(tstmod.GSA, ~. -sex:thick.centerd)
tstmod.AT <- update(tstmod.ST, ~. -age.centerd:thick.centerd)
tstmod.SA <- update(tstmod.AT, ~. -sex:age.centerd)
tstmod.GT <- update(tstmod.SA, ~. -group:thick.centerd)
tstmod.GS <- update(tstmod.GT, ~. -group:sex)
tstmod.GA <- update(tstmod.GS, ~. -group:age.centerd)
tstmod.T <- update(tstmod.GA, ~. -thick.centerd)
tstmod.S <- update(tstmod.T, ~. -sex)
tstmod.A <- update(tstmod.S, ~. -age.centerd)
tstmod.G <- update(tstmod.A, ~. -group)

anova(tstmod.G, tstmod.A, tstmod.S, tstmod.T,
      tstmod.GA, tstmod.GS, tstmod.GT, tstmod.SA, tstmod.AT, tstmod.ST,
      tstmod.GSA, tstmod.GAT, tstmod.GST, tstmod.AST,
      tstmod.Full3,
      test="Chisq")

anova(tstmod.Full3, test="Chisq")

# Model summary
mod.sim <- sim(tstmod.Full3, n.sims=1000)
x1 <- summary(tstmod.Full3)$coefficients
x2 <- t(apply(coef(mod.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1[,1], x2)

(exp(sums)-1)*100

(((x1[2]+x1[9])/x1[1]))*100
quantile((((coef(mod.sim)[,2]+coef(mod.sim)[,2])/coef(mod.sim)[,1]))*100, c(0.025, 0.975))


# Plot top model
agespan <- seq(min(alldata$age.centerd),max(alldata$age.centerd), 0.1)
agespan2 <- seq(min(alldata$age),max(alldata$age), 0.1)

new.dat <- data.frame(age.centerd=rep(agespan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thick.centerd=mean(alldata$thick.centerd),
                      age=rep(agespan2, 4))

new.dat$pred  <- exp(predict(tstmod.Full3, new.dat, re.form=NA))
ggplot(aes(x=age, y=nevent.u.m2.min, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)



thickspan <- seq(min(alldata$thick.centerd),max(alldata$thick.centerd), 0.01)
thickspan2 <- seq(min(alldata$thick),max(alldata$thick), 0.01)

new.dat <- data.frame(thick.centerd=rep(thickspan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(thickspan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(thickspan)), 2)),
                      age.centerd=mean(alldata$age.centerd),
                      thick=rep(thickspan2, 4))

new.dat$pred <- exp(predict(tstmod.Full3, new.dat, re.form=NA))

ggplot(aes(x=thick, y=nevent.u.m2.min, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1, data=new.dat)

#END