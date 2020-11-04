# 
library(lme4)
library(ggplot2)

## Load data
# load(file='X://PD_longrest//groupanalysis//bbdata.Rdata')
load(file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//bbdata.Rdata')

bbdata$age.centered <- bbdata$age-mean(bbdata$age)
bbdata$log.len <- log(bbdata$leneve)
bbdata$log.tue <- log(bbdata$tueeve)
bbdata$log.max <- log(bbdata$maxeve)

## EVENT LENGTH
lenmod <- lmer(log.len ~ I(age.centered^2)*group*sex + age.centered*group*sex + (1|subj), data=bbdata, REML=FALSE)
lenmod.1 <- update(lenmod, ~. -I(age.centered^2):group:sex)
lenmod.2 <- update(lenmod.1, ~. -age.centered:group:sex)
lenmod.3 <- update(lenmod.2, ~. -I(age.centered^2):group)
lenmod.4 <- update(lenmod.3, ~. -age.centered:group)
lenmod.5 <- update(lenmod.4, ~. -I(age.centered^2):sex)
lenmod.6 <- update(lenmod.5, ~. -age.centered:sex)
lenmod.7 <- update(lenmod.6, ~. -group:sex) 
lenmod.8 <- update(lenmod.7, ~. -group)
lenmod.9 <- update(lenmod.8, ~. -sex)
lenmod.10 <- update(lenmod.9, ~. -I(age.centered^2))
lenmod.11 <- update(lenmod.10, ~. -age.centered)

anova(lenmod,lenmod.1,lenmod.2,lenmod.3,lenmod.4,lenmod.5,lenmod.6,lenmod.7,lenmod.8,lenmod.9,lenmod.10,lenmod.11)

tmpdat <- subset(bbdata, hemi=="lh")
tmpdat$lenpred.fx <- exp(predict(lenmod, re.form=NA))
tmpdat$lenpred.rx <- exp(predict(lenmod))

ggplot(aes(x=age, y=lenpred.rx, color=group, shape=sex), data=tmpdat)+
  geom_point()+
  geom_line(aes(y = lenpred.fx, linetype=sex), size = 1)



## TIME UNTIL EVENT
tuemod <- lmer(log.tue ~ I(age.centered^2)*group*sex + age.centered*group*sex + (1|subj), data=bbdata, REML=FALSE)
tuemod.1 <- update(tuemod, ~. -I(age.centered^2):group:sex)
tuemod.2 <- update(tuemod.1, ~. -age.centered:group:sex)
tuemod.3 <- update(tuemod.2, ~. -I(age.centered^2):group)
tuemod.4 <- update(tuemod.3, ~. -age.centered:group)
tuemod.5 <- update(tuemod.4, ~. -I(age.centered^2):sex)
tuemod.6 <- update(tuemod.5, ~. -age.centered:sex)
tuemod.7 <- update(tuemod.6, ~. -group:sex) 
tuemod.8 <- update(tuemod.7, ~. -group)
tuemod.9 <- update(tuemod.8, ~. -sex)
tuemod.10 <- update(tuemod.9, ~. -I(age.centered^2))
tuemod.11 <- update(tuemod.10, ~. -age.centered)

anova(tuemod,tuemod.1,tuemod.2,tuemod.3,tuemod.4,tuemod.5,tuemod.6,tuemod.s,tuemod.7,tuemod.8,tuemod.9,tuemod.10)

tmpdat <- bbdata
tmpdat$tuepred.fx <- exp(predict(tuemod, re.form=NA))
tmpdat$tuepred.rx <- exp(predict(tuemod))

ggplot(aes(x=age, y=tuepred.rx, color=group, shape=sex), data=tmpdat)+
  geom_point()+
  geom_line(aes(y = tuepred.fx, linetype=sex), size = 1)


## MAX
maxmod <- lmer(log.max ~ I(age.centered^2)*group*sex + age.centered*group*sex + (1|subj), data=bbdata, REML=FALSE)
maxmod.1 <- update(maxmod, ~. -I(age.centered^2):group:sex)
maxmod.2 <- update(maxmod.1, ~. -age.centered:group:sex)
maxmod.3 <- update(maxmod.2, ~. -I(age.centered^2):group)
maxmod.4 <- update(maxmod.3, ~. -age.centered:group)
maxmod.5 <- update(maxmod.4, ~. -I(age.centered^2):sex)
maxmod.6 <- update(maxmod.5, ~. -age.centered:sex)
maxmod.7 <- update(maxmod.6, ~. -group:sex) 
maxmod.8 <- update(maxmod.7, ~. -group)
maxmod.9 <- update(maxmod.8, ~. -sex)
maxmod.10 <- update(maxmod.9, ~. -I(age.centered^2))
maxmod.11 <- update(maxmod.10, ~. -age.centered)

anova(maxmod,maxmod.1,maxmod.2,maxmod.3,maxmod.4,maxmod.5,maxmod.6,maxmod.7,maxmod.8,maxmod.9,maxmod.10,maxmod.11)

tmpdat$maxpred.fx <- exp(predict(maxmod.6, re.form=NA))
tmpdat$maxpred.rx <- exp(predict(maxmod.6))

ggplot(aes(x=age, y=maxpred.rx, color=group, shape=sex), data=tmpdat)+
  geom_point()+
  geom_line(aes(y = maxpred.fx, linetype=sex), size = 1)








