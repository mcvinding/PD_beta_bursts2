# Mixed-model analysis of burst level data
library(lme4)
library(ggplot2)
library(arm)
source('X://PD_longrest//scripts//functions//zscore.R')

## Load data
load(file='X://PD_longrest//groupanalysis//bbdata2.Rdata')
load(file='X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
# load(file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//bbdata2.Rdata')

bbdata$age.centerd <- bbdata$age-mean(bbdata$age)
bbdata$log.len <- log(bbdata$leneve)
bbdata$log.tue <- log(bbdata$tueeve)
bbdata$log.max <- log(bbdata$maxeve)
bbdata$thick.centerd <- bbdata$thick-mean(bbdata$thick)
bbdata$thickz <- zscore(bbdata$thick)
bbdata$lenz <- zscore(bbdata$leneve)


######################################################################################
## EVENT LENGTH
lenmod <- lmer(log.len ~ (group+age.centerd+sex+thickz)^3 + (1|subj), 
               data=bbdata, REML=FALSE)

qqnorm(resid(lenmod))
qqline(resid(lenmod))

lenmod.AST <- update(lenmod,     ~. -age.centerd:sex:thickz)
lenmod.GST <- update(lenmod.AST, ~. -group:sex:thickz)
lenmod.GAT <- update(lenmod.GST, ~. -group:age.centerd:thickz)
lenmod.GSA <- update(lenmod.GAT, ~. -group:sex:age.centerd)
lenmod.ST  <- update(lenmod.GSA, ~. -sex:thickz)
lenmod.AT  <- update(lenmod.ST,  ~. -age.centerd:thickz)
lenmod.SA  <- update(lenmod.AT,  ~. -sex:age.centerd)
lenmod.GT  <- update(lenmod.SA,  ~. -group:thickz)
lenmod.GS  <- update(lenmod.GT,  ~. -group:sex)
lenmod.GA  <- update(lenmod.GS,  ~. -group:age.centerd)
lenmod.T   <- update(lenmod.GA,  ~. -thickz)
lenmod.S   <- update(lenmod.T,   ~. -sex)
lenmod.A   <- update(lenmod.S,   ~. -age.centerd)
lenmod.G   <- update(lenmod.A,   ~. -group)

anova(lenmod.G,
      lenmod.A, 
      lenmod.S,
      lenmod.T,
      lenmod.GA,
      lenmod.GS,
      lenmod.GT,
      lenmod.SA,
      lenmod.AT,
      lenmod.ST,
      lenmod.GSA,
      lenmod.GAT, 
      lenmod.GST,
      lenmod.AST,
      lenmod,
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

# Plot top model
bbsum <- aggregate(cbind(leneve, tueeve, maxeve)~subj, data=bbdata, FUN=median)
alldata <- merge(alldata, bbsum, by="subj")
agespan <- seq(min(bbdata$age.centerd),max(bbdata$age.centerd), 0.1)
agespan2 <- seq(min(bbdata$age),max(bbdata$age), 0.1)

new.dat <- data.frame(age.centerd=rep(agespan,4),
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thickz=0,
                      age=rep(agespan2, 4))

new.dat$pred <- exp(predict(lenmod, new.dat, re.form=NA))
ggplot(aes(x=age, y=leneve, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)

######################################################################################
## TIME UNTIL EVENT
tuemod <- lmer(log.tue ~ (group+age.centerd+sex+thickz)^3 + (1|subj), data=bbdata, REML=FALSE)

qqnorm(resid(tuemod))
qqline(resid(tuemod))

tuemod.AST <- update(tuemod,     ~. -age.centerd:sex:thickz)
tuemod.GST <- update(tuemod.AST, ~. -group:sex:thickz)
tuemod.GAT <- update(tuemod.GST, ~. -group:age.centerd:thickz)
tuemod.GSA <- update(tuemod.GAT, ~. -group:sex:age.centerd)
tuemod.ST  <- update(tuemod.GSA, ~. -sex:thickz)
tuemod.AT  <- update(tuemod.ST,  ~. -age.centerd:thickz)
tuemod.SA  <- update(tuemod.AT,  ~. -sex:age.centerd)
tuemod.GT  <- update(tuemod.SA,  ~. -group:thickz)
tuemod.GS  <- update(tuemod.GT,  ~. -group:sex)
tuemod.GA  <- update(tuemod.GS,  ~. -group:age.centerd)
tuemod.T   <- update(tuemod.GA,  ~. -thickz)
tuemod.S   <- update(tuemod.T,   ~. -sex)
tuemod.A   <- update(tuemod.S,   ~. -age.centerd)
tuemod.G   <- update(tuemod.A,   ~. -group)

anova(tuemod.G, tuemod.A, tuemod.S, tuemod.T,
      tuemod.GA, tuemod.GS, tuemod.GT, tuemod.SA, tuemod.AT, tuemod.ST,
      tuemod.GSA, tuemod.GAT, tuemod.GST, tuemod.AST,
      tuemod,
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

# Plot
bbsum <- aggregate(cbind(leneve, tueeve, maxeve)~subj, data=bbdata, FUN=median)
alldata <- merge(alldata, bbsum, by="subj")
agespan <- seq(min(bbdata$age.centerd),max(bbdata$age.centerd), 0.1)
agespan2 <- seq(min(bbdata$age),max(bbdata$age), 0.1)

new.dat <- data.frame(age.centerd=rep(agespan,4),
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thickz=0,
                      age=rep(agespan2, 4))

tmpdat <- bbdata
tmpdat$tuepred.fx <- exp(predict(tuemod, re.form=NA))
tmpdat$tuepred.rx <- exp(predict(tuemod))

new.dat$pred <- exp(predict(tuemod, new.dat, re.form=NA))
ggplot(aes(x=age, y=tueeve, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)


ggplot(aes(x=age, y=tuepred.rx, color=group, shape=sex), data=tmpdat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1, data=new.dat)


######################################################################################
## MAX
maxmod <- lmer(log.max ~ (group+age.centerd+sex+thickz)^3 + (1|subj), data=bbdata, REML=FALSE)

qqnorm(resid(maxmod))
qqline(resid(maxmod))

maxmod.AST <- update(maxmod,     ~. -age.centerd:sex:thickz)
maxmod.GST <- update(maxmod.AST, ~. -group:sex:thickz)
maxmod.GAT <- update(maxmod.GST, ~. -group:age.centerd:thickz)
maxmod.GSA <- update(maxmod.GAT, ~. -group:sex:age.centerd)
maxmod.ST  <- update(maxmod.GSA, ~. -sex:thickz)
maxmod.AT  <- update(maxmod.ST,  ~. -age.centerd:thickz)
maxmod.SA  <- update(maxmod.AT,  ~. -sex:age.centerd)
maxmod.GT  <- update(maxmod.SA,  ~. -group:thickz)
maxmod.GS  <- update(maxmod.GT,  ~. -group:sex)
maxmod.GA  <- update(maxmod.GS,  ~. -group:age.centerd)
maxmod.T   <- update(maxmod.GA,  ~. -thickz)
maxmod.S   <- update(maxmod.T,   ~. -sex)
maxmod.A   <- update(maxmod.S,   ~. -age.centerd)
maxmod.G   <- update(maxmod.A,   ~. -group)

anova(maxmod.G, maxmod.A, maxmod.S, maxmod.T,
      maxmod.GA, maxmod.GS, maxmod.GT, maxmod.SA, maxmod.AT, maxmod.ST,
      maxmod.GSA, maxmod.GAT, maxmod.GST, maxmod.AST,
      maxmod,
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

#Plots
tmpdat$maxpred.fx <- exp(predict(maxmod, re.form=NA))
tmpdat$maxpred.rx <- exp(predict(maxmod))
ggplot(aes(x=age, y=maxpred.rx, color=group, shape=sex), data=tmpdat)+
  geom_point()+
  geom_line(aes(y = maxpred.fx, linetype=sex), size = 1)


new.dat$pred <- exp(predict(maxmod, new.dat, re.form=NA))
ggplot(aes(x=age, y=maxeve, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)

#END