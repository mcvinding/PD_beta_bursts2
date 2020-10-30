# 
library(lme4)
library(ggplot2)

## Load data
load(file='X://PD_longrest//groupanalysis//bbdata.Rdata')

bbdata$age.centered <- bbdata$age-mean(bbdata$age)
bbdata$log.len <- log(bbdata$leneve)
bbdata$log.toe <- log(bbdata$toeeve)
bbdata$log.max <- log(bbdata$maxeve)

## EVENT LENGTH
lenmod <- lmer(log.len ~ I(age.centered^2)*group*sex + age.centered*group*sex + (1|subj), data=bbdata, subset=hemi=='lh', REML=FALSE)
lenmod.1 <- update(lenmod, ~. -I(age.centered^2):group:sex)
lenmod.2 <- update(lenmod.1, ~. -age.centered:group:sex)
lenmod.3 <- update(lenmod.2, ~. -I(age.centered^2):group)
lenmod.4 <- update(lenmod.3, ~. -age.centered:group)
lenmod.5 <- update(lenmod.4, ~. -I(age.centered^2):sex)
lenmod.6 <- update(lenmod.5, ~. -age.centered:sex)
lenmod.s <- update(lenmod.6, ~. -group:sex) 
lenmod.7 <- update(lenmod.s, ~. -group)
lenmod.8 <- update(lenmod.7, ~. -sex)
lenmod.9 <- update(lenmod.8, ~. -I(age.centered^2))
lenmod.10 <- update(lenmod.9, ~. -age.centered)

anova(lenmod,lenmod.1,lenmod.2,lenmod.3,lenmod.4,lenmod.5,lenmod.6,lenmod.s,lenmod.7,lenmod.8,lenmod.9,lenmod.10)

tmpdat <- subset(bbdata, hemi=="lh")
tmpdat$lenpred.fx <- exp(predict(lenmod, re.form=NA))
tmpdat$lenpred.rx <- exp(predict(lenmod))

ggplot(aes(x=age, y=lenpred.rx, color=group, shape=sex), data=tmpdat)+
  geom_point()+
  geom_line(aes(y = lenpred.fx, linetype=sex), size = 1)



## TIME UNTIL EVENT
toemod <- lmer(log.toe ~ I(age.centered^2)*group*sex + age.centered*group*sex + (1|subj), data=bbdata, subset=hemi=='lh', REML=FALSE)
toemod.1 <- update(toemod, ~. -I(age.centered^2):group:sex)
toemod.2 <- update(toemod.1, ~. -age.centered:group:sex)
toemod.3 <- update(toemod.2, ~. -I(age.centered^2):group)
toemod.4 <- update(toemod.3, ~. -age.centered:group)
toemod.5 <- update(toemod.4, ~. -I(age.centered^2):sex)
toemod.6 <- update(toemod.5, ~. -age.centered:sex)
toemod.s <- update(toemod.6, ~. -group:sex) 
toemod.7 <- update(toemod.s, ~. -group)
toemod.8 <- update(toemod.7, ~. -sex)
toemod.9 <- update(toemod.8, ~. -I(age.centered^2))
toemod.10 <- update(toemod.9, ~. -age.centered)

anova(toemod,toemod.1,toemod.2,toemod.3,toemod.4,toemod.5,toemod.6,toemod.s,toemod.7,toemod.8,toemod.9,toemod.10)

tmpdat$toepred.fx <- exp(predict(toemod, re.form=NA))
tmpdat$toepred.rx <- exp(predict(toemod))

ggplot(aes(x=age, y=toepred.rx, color=group, shape=sex), data=tmpdat)+
  geom_point()+
  geom_line(aes(y = toepred.fx, linetype=sex), size = 1)


## MAX
maxmod <- lmer(log.max ~ I(age.centered^2)*group*sex + age.centered*group*sex + (1|subj), data=bbdata, subset=hemi=='lh', REML=FALSE)
maxmod.1 <- update(maxmod, ~. -I(age.centered^2):group:sex)
maxmod.2 <- update(maxmod.1, ~. -age.centered:group:sex)
maxmod.3 <- update(maxmod.2, ~. -I(age.centered^2):group)
maxmod.4 <- update(maxmod.3, ~. -age.centered:group)
maxmod.5 <- update(maxmod.4, ~. -I(age.centered^2):sex)
maxmod.6 <- update(maxmod.5, ~. -age.centered:sex)
maxmod.s <- update(maxmod.6, ~. -group:sex) 
maxmod.7 <- update(maxmod.s, ~. -group)
maxmod.8 <- update(maxmod.7, ~. -sex)
maxmod.9 <- update(maxmod.8, ~. -I(age.centered^2))
maxmod.10 <- update(maxmod.9, ~. -age.centered)

anova(maxmod,maxmod.1,maxmod.2,maxmod.3,maxmod.4,maxmod.5,maxmod.6,maxmod.s,maxmod.7,maxmod.8,maxmod.9,maxmod.10)

tmpdat$maxpred.fx <- exp(predict(maxmod, re.form=NA))
tmpdat$maxpred.rx <- exp(predict(maxmod))

ggplot(aes(x=age, y=maxpred.rx, color=group, shape=sex), data=tmpdat)+
  geom_point()+
  geom_line(aes(y = maxpred.fx, linetype=sex), size = 1)

anova(maxmod,maxmod.1,maxmod.2,maxmod.3,maxmod.4,maxmod.5)








