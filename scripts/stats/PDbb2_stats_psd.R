# PD beta burst statistics: analsyis of PSD output from FOOOF analysis
library(ggplot2)
library(arm)

# Load data
load(file='X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
# load('C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj.Rdata')
# alldata$age.centerd <- alldata$age-mean(alldata$age)
# alldata$thick.centerd <- alldata$thick-mean(alldata$thick)

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
finter.Full3 <- glm(a_intercept ~ (group+age.centerd+sex+thick.centerd)^3, data=alldata)
anova(finter.AST, test="Chisq")


d.AD <- data.frame(counts=c(18,17,15,20,10,20,25,13,12),
                   outcome=gl(3,1,9),
                   treatment=gl(3,3))
glm1 <- glm(counts ~ outcome + treatment, family = gaussian(), data=d.AD)
glm0 <- update(glm1, . ~ 1)

anova(glm0,glm1, test="F")


finter.AST <- update(finter.Full3, ~. -age.centerd:sex:thick.centerd)
finter.GST <- update(finter.AST, ~. -group:sex:thick.centerd)
finter.GAT <- update(finter.GST, ~. -group:age.centerd:thick.centerd)
finter.GSA <- update(finter.GAT, ~. -group:sex:age.centerd)
finter.ST  <- update(finter.GSA, ~. -sex:thick.centerd)
finter.AT  <- update(finter.ST, ~. -age.centerd:thick.centerd)
finter.SA  <- update(finter.AT, ~. -sex:age.centerd)
finter.GT  <- update(finter.SA, ~. -group:thick.centerd)
finter.GS  <- update(finter.GT, ~. -group:sex)
finter.GA  <- update(finter.GS, ~. -group:age.centerd)
finter.T   <- update(finter.GA, ~. -thick.centerd)
finter.S   <- update(finter.T, ~. -sex)
finter.A   <- update(finter.S, ~. -age.centerd)
finter.G   <- update(finter.A, ~. -group)

anova(finter.G,
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
      finter.Full3, test="F")

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
inter.sim <- sim(inter.Full3, n.sims=1000)
x1 <- summary(inter.Full3)$coefficients
x2 <- t(apply(coef(inter.sim), 2, quantile, c(0.025, 0.975)))
cbind(x1[,1], x2)

# Group: female
(((x1[1]-x1[2])/x1[1]))*100
quantile((((coef(inter.sim)[,1]-coef(inter.sim)[,2])/coef(inter.sim)[,1])-1)*100, c(0.025, 0.975))

# Group: male
(((x1[1]-x1[2]-x1[3]-x1[7])/x1[1])-1)*100
quantile((((coef(inter.sim)[,1]-coef(inter.sim)[,2]-coef(inter.sim)[,3]-coef(inter.sim)[,7])/coef(inter.sim)[,1])-1)*100, c(0.025, 0.975))

# Plot top model
agespan <- seq(min(alldata$age.centerd),max(alldata$age.centerd), 0.1)
agespan2 <- seq(min(alldata$age),max(alldata$age), 0.1)

new.dat <- data.frame(age.centerd=rep(agespan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thick.centerd=mean(alldata$thick.centerd),
                      age=rep(agespan2, 4))

new.dat$pred <- predict(tstmod.Full3, new.dat, re.form=NA)
ggplot(aes(x=age, y=a_intercept, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)+
  theme_classic()


# NOT WORKING!!!
plot_ly(x=alldata$age, z=alldata$a_intercept, y=alldata$thick, type="scatter3d", mode="markers", color=alldata$group)

######################################################################################
# 1/f slope
fslope.Full3 <- glm(a_slope ~ (group+age.centerd+sex+thick.centerd)^3, data=alldata, family=gaussian)
anova(tstmod.Full3, test="Chisq")

# Model summary
mod.sim <- sim(tstmod.Full3, n.sims=1000)
x1 <- summary(tstmod.Full3)$coefficients
x2 <- t(apply(coef(mod.sim), 2, quantile, c(0.025, 0.975)))
cbind(x1[,1], x2)

# Group
((x1[2]/x1[1]))*100
quantile(((coef(mod.sim)[,2]/coef(mod.sim)[,1]))*100, c(0.025, 0.975))

# Thickness
((x1[5]/x1[1]))*100
quantile(((coef(mod.sim)[,5]/coef(mod.sim)[,1]))*100, c(0.025, 0.975))

# Plot top model
agespan <- seq(min(alldata$age.centerd),max(alldata$age.centerd), 0.1)
agespan2 <- seq(min(alldata$age),max(alldata$age), 0.1)

new.dat <- data.frame(age.centerd=rep(agespan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thick.centerd=mean(alldata$thick.centerd),
                      age=rep(agespan2, 4))

new.dat$pred <- predict(tstmod.Full3, new.dat, re.form=NA)
ggplot(aes(x=age, y=a_slope, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)

plot_ly(x=alldata$age, z=alldata$a_slope, y=alldata$thick, type="scatter3d", mode="markers", color=alldata$group, symbol=alldata$sex)


######################################################################################
# Beta power
beta_pw.Full3 <- glm(beta_pw ~ (group+age.centerd+sex+thick.centerd)^3, data=alldata, family=gaussian)
anova(tstmod.Full3, test="Chisq")

# Model summary
beta_pw.sim <- sim(beta_pw.Full3, n.sims=1000)
x1 <- summary(beta_pw.Full3)$coefficients
x2 <- t(apply(coef(beta_pw.sim), 2, quantile, c(0.025, 0.975)))
cbind(x1[,1], x2)

((x1[2]/x1[1]))*100
quantile(((coef(beta_pw.sim)[,2]/coef(beta_pw.sim)[,1]))*100, c(0.025, 0.975))

((x1[3]/x1[1]))*100
quantile(((coef(beta_pw.sim)[,3]/coef(beta_pw.sim)[,1]))*100, c(0.025, 0.975))

# Plot top model
agespan <- seq(min(alldata$age.centerd),max(alldata$age.centerd), 0.1)
agespan2 <- seq(min(alldata$age),max(alldata$age), 0.1)

new.dat <- data.frame(age.centerd=rep(agespan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thick.centerd=mean(alldata$thick.centerd),
                      age=rep(agespan2, 4))

new.dat$pred <- predict(beta_pw.Full3, new.dat, re.form=NA)
ggplot(aes(x=age, y=beta_pw, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)

######################################################################################
# Beta peak freq
beta_cf.Full3 <- glm(beta_cf ~ (group+age.centerd+sex+thick.centerd)^3, data=alldata, family=gaussian)
anova(beta_cf.Full3, test="Chisq")

# Model summary
beta_cf.sim <- sim(beta_cf.Full3, n.sims=1000)
x1 <- summary(beta_cf.Full3)$coefficients
x2 <- t(apply(coef(beta_cf.sim), 2, quantile, c(0.025, 0.975)))
cbind(x1[,1], x2)

# Plot top model
agespan <- seq(min(alldata$age.centerd),max(alldata$age.centerd), 0.1)
agespan2 <- seq(min(alldata$age),max(alldata$age), 0.1)

new.dat <- data.frame(age.centerd=rep(agespan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thick.centerd=mean(alldata$thick.centerd),
                      age=rep(agespan2, 4))

new.dat$pred <- predict(beta_cf.Full3, new.dat, re.form=NA)
ggplot(aes(x=age, y=beta_cf, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)

######################################################################################
# Alpha power
alpha_pw.Full3 <- glm(alpha_pw ~ (group+age.centerd+sex+thick.centerd)^3, data=alldata, family=gaussian)
anova(alpha_pw.Full3, test="Chisq")

# Model summary
alpha_pw.sim <- sim(alpha_pw.Full3, n.sims=1000)
x1 <- summary(alpha_pw.Full3)$coefficients
x2 <- t(apply(coef(alpha_pw.sim), 2, quantile, c(0.025, 0.975)))
cbind(x1[,1], x2)

# Plot top model
agespan <- seq(min(alldata$age.centerd),max(alldata$age.centerd), 0.1)
agespan2 <- seq(min(alldata$age),max(alldata$age), 0.1)

new.dat <- data.frame(age.centerd=rep(agespan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thick.centerd=mean(alldata$thick.centerd),
                      age=rep(agespan2, 4))

new.dat$pred <- predict(alpha_pw.Full3, new.dat, re.form=NA)
ggplot(aes(x=age, y=alpha_pw, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)

######################################################################################
# Alpha peak freq
alpha_cf.Full3 <- glm(alpha_cf ~ (group+age.centerd+sex+thick.centerd)^3, data=alldata, family=gaussian)
anova(alpha_cf.Full3, test="Chisq")

# Model summary
alpha_cf.sim <- sim(alpha_cf.Full3, n.sims=1000)
x1 <- summary(alpha_cf.Full3)$coefficients
x2 <- t(apply(coef(alpha_cf.sim), 2, quantile, c(0.025, 0.975)))
cbind(x1[,1], x2)

# Plot top model
agespan <- seq(min(alldata$age.centerd),max(alldata$age.centerd), 0.1)
agespan2 <- seq(min(alldata$age),max(alldata$age), 0.1)

new.dat <- data.frame(age.centerd=rep(agespan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thick.centerd=mean(alldata$thick.centerd),
                      age=rep(agespan2, 4))

new.dat$pred <- predict(alpha_cf.Full3, new.dat, re.form=NA)
ggplot(aes(x=age, y=alpha_cf, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_line(aes(y=pred, linetype=sex), size = 1, data=new.dat)

#END