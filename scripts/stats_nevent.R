# Analysis of N events
library(lme4)
library(arm)
library(ggplot2)
# Bayes?

## Load data
load('X://PD_longrest//groupanalysis//ndata.Rdata')

## Center variables
ndata$age.centerd <- ndata$age-mean(ndata$age)
ndata$thick.centerd <- ndata$thick-mean(ndata$thick)

# * age
# * thickness

## Polynomial variables

## Quick analysis
### 3rd degree 
tstmod.a <- lmer(nevent.min ~ age.centerd + I(age.centerd^2) + I(age.centerd^3) + 
                   group +  sex + 
                   group:sex +
                   group:age.centerd +
                   sex:age.centerd +
                   group:sex:age.centerd +
                   group:I(age.centerd^2) +
                   sex:I(age.centerd^2) +
                   group:sex:I(age.centerd^2) +
                   group:I(age.centerd^3) +
                   sex:I(age.centerd^3) +
                   group:sex:I(age.centerd^3) +
                   (1|subj) + (1|hemi),  data=ndata, REML=F)

tstmod.b <- lmer(nevent.min ~ age.centerd + I(age.centerd^2) + 
                   group +  sex + 
                   group:sex +
                   group:age.centerd +
                   sex:age.centerd +
                   group:sex:age.centerd +
                   group:I(age.centerd^2) +
                   sex:I(age.centerd^2) +
                   group:sex:I(age.centerd^2) +
                   (1|subj) + (1|hemi),  data=ndata, REML=F)
tstmod.c <- lmer(nevent.min ~ age.centerd + 
                   group +  sex + 
                   group:sex +
                   group:age.centerd +
                   sex:age.centerd +
                   group:sex:age.centerd +
                   (1|subj) + (1|hemi),  data=ndata, REML=F)

anova(tstmod.a, tstmod.b, tstmod.c)

summary(tstmod)

# 2nd degree
formula <- nevent.min ~ thick.centerd + I(thick.centerd^2) + group + sex + age.centerd +
  thick.centerd:group + thick.centerd:sex + thick.centerd:age.centerd +
  thick.centerd:group:sex + thick.centerd:group:age.centerd + thick.centerd:sex:age.centerd +
  thick.centerd:group:sex:age.centerd +
  I(thick.centerd^2):group + I(thick.centerd^2):sex + I(thick.centerd^2):age.centerd +
  I(thick.centerd^2):group:sex + I(thick.centerd^2):group:age.centerd + I(thick.centerd^2):sex:age.centerd +
  I(thick.centerd^2):group:sex:age.centerd +
  (1|subj) + (1|hemi)


tstmod2 <- glmer(formula,  data=ndata, REML=FALSE, family = poisson)
tstmod2.1 <- update(tstmod2, ~. -I(thick.centerd^2):group:sex:age.centerd -thick.centerd:group:sex:age.centerd)
tstmod2.2 <- update(tstmod2.1, ~. -I(thick.centerd^2):sex:age.centerd -thick.centerd:sex:age.centerd)
tstmod2.3 <- update(tstmod2.2, ~. -I(thick.centerd^2):group:age.centerd -thick.centerd:group:age.centerd)
tstmod2.4 <- update(tstmod2.3, ~. -I(thick.centerd^2):group:sex -thick.centerd:group:sex)
tstmod2.5 <- update(tstmod2.4, ~. -I(thick.centerd^2):age.centerd -thick.centerd:age.centerd)
tstmod2.6 <- update(tstmod2.5, ~. -I(thick.centerd^2):sex -thick.centerd:sex)
tstmod2.7 <- update(tstmod2.6, ~. -I(thick.centerd^2):group -thick.centerd:group)
tstmod2.8 <- update(tstmod2.7, ~. -age.centerd)
tstmod2.9 <- update(tstmod2.8, ~. -group)
tstmod2.10 <- update(tstmod2.9, ~. -sex)
tstmod2.11 <- update(tstmod2.10, ~. -I(thick.centerd^2))
tstmod2.12 <- update(tstmod2.11, ~. -thick.centerd)

anova(tstmod2, tstmod2.1, tstmod2.2, tstmod2.3, tstmod2.4, tstmod2.5, tstmod2.6, tstmod2.7, tstmod2.8, tstmod2.9,tstmod2.10, tstmod2.11, tstmod2.12)

summary(tstmod2.5)

## Model based on age as primary predictor
formula <- nevent.min ~ age.centerd + I(age.centerd^2) + group + sex + thick.centerd +
  age.centerd:group + age.centerd:sex + age.centerd:thick.centerd +
  age.centerd:group:sex + age.centerd:group:thick.centerd + age.centerd:sex:thick.centerd +
  age.centerd:group:sex:thick.centerd +
  I(age.centerd^2):group + I(age.centerd^2):sex + I(age.centerd^2):thick.centerd +
  I(age.centerd^2):group:sex + I(age.centerd^2):group:thick.centerd + I(age.centerd^2):sex:thick.centerd +
  I(age.centerd^2):group:sex:thick.centerd +
  (1|subj) + (1|hemi)

agemod <- glmer(formula,  data=ndata, family=poisson)
agemod.1 <- update(agemod, ~. -I(age.centerd^2):group:sex:thick.centerd -age.centerd:group:sex:thick.centerd)
agemod.2 <- update(agemod.1, ~. -I(age.centerd^2):sex:thick.centerd -age.centerd:sex:thick.centerd)
agemod.3 <- update(agemod.2, ~. -I(age.centerd^2):group:thick.centerd -age.centerd:group:thick.centerd)
agemod.4 <- update(agemod.3, ~. -I(age.centerd^2):group:sex -age.centerd:group:sex)
agemod.5 <- update(agemod.4, ~. -I(age.centerd^2):thick.centerd -age.centerd:thick.centerd)
agemod.6 <- update(agemod.5, ~. -I(age.centerd^2):sex -age.centerd:sex)
agemod.7 <- update(agemod.6, ~. -I(age.centerd^2):group -age.centerd:group)
agemod.8 <- update(agemod.7, ~. -thick.centerd)
agemod.9 <- update(agemod.8, ~. -group)
agemod.10 <- update(agemod.9, ~. -sex)
agemod.11 <- update(agemod.10, ~. -I(age.centerd^2))
agemod.12 <- update(agemod.11, ~. -age.centerd)

anova(agemod, agemod.1, agemod.2, agemod.3, agemod.4, agemod.5, agemod.6, agemod.7, agemod.8, agemod.9, agemod.10, agemod.11, agemod.12)


## MISC

nmod <- lmer(nevent.min ~  age.centerd * sex * group + thick.centerd*age.centerd*group + (1|hemi) + (1|subj), data=ndata)
nmod.x3 <- lmer(nevent.min ~ I(age-mean(age)) * sex * group + (1|hemi) + (1|subj), data=ndata)

nplot <- ggplot(aes(x=age, y=nevent, color=group, shape=sex), data=ndata)+
  geom_point()+
  geom_smooth(method=lm)
nplot

t.test(ndata$nevent.min[ndata$group=="Patient"], ndata$nevent.min[ndata$group=="Control"])
