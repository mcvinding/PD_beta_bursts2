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

summary(tstmod.a)

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
formula <- nevent.min ~ age.centerd + I(age^2) + group + sex + thick.centerd +
  age.centerd:group + age.centerd:sex + age.centerd:thick.centerd +
  age.centerd:group:sex + age.centerd:group:thick.centerd + age.centerd:sex:thick.centerd +
  age.centerd:group:sex:thick.centerd +
  I(age^2):group + I(age^2):sex + I(age^2):thick.centerd +
  I(age^2):group:sex + I(age^2):group:thick.centerd + I(age^2):sex:thick.centerd +
  I(age^2):group:sex:thick.centerd +
  (1|subj) + (1|hemi)

agemod <- glmer(formula,  data=ndata, family=poisson)
agemod <- lmer(formula,  data=ndata, REML=FALSE)

agemod.1 <- update(agemod, ~. -I(age^2):group:sex:thick.centerd -age.centerd:group:sex:thick.centerd)
agemod.2 <- update(agemod.1, ~. -I(age^2):sex:thick.centerd -age.centerd:sex:thick.centerd)
agemod.3 <- update(agemod.2, ~. -I(age^2):group:thick.centerd -age.centerd:group:thick.centerd)
agemod.4 <- update(agemod.3, ~. -I(age^2):group:sex -age.centerd:group:sex)
agemod.5 <- update(agemod.4, ~. -I(age^2):thick.centerd -age.centerd:thick.centerd)
agemod.6 <- update(agemod.5, ~. -I(age^2):sex -age.centerd:sex)
agemod.7 <- update(agemod.6, ~. -I(age^2):group -age.centerd:group)
agemod.8 <- update(agemod.7, ~. -thick.centerd)
agemod.9 <- update(agemod.8, ~. -group)
agemod.10 <- update(agemod.9, ~. -sex)
agemod.11 <- update(agemod.10, ~. -I(aged^2))
agemod.12 <- update(agemod.11, ~. -age.centerd)

anova(agemod, agemod.1, agemod.2, agemod.3, agemod.4, agemod.5, agemod.6, agemod.7, agemod.8, agemod.9, agemod.10, agemod.11, agemod.12)

#############################################
# No square
formula <- nevent.min ~ age.centerd * group * sex * thick.centerd * hemi + (1|subj)

agemod <- glmer(formula,  data=ndata, family=poisson)

agemod <- lmer(formula,  data=ndata, REML=FALSE)

# 5-way
agemod.1 <- update(agemod, ~. -age.centerd:group:sex:thick.centerd:hemi)
# 4-way
agemod.2 <- update(agemod.1, ~. -age.centerd:sex:thick.centerd:hemi)
agemod.3 <- update(agemod.2, ~. -age.centerd:group:thick.centerd:hemi)
agemod.4 <- update(agemod.3, ~. -age.centerd:group:sex:hemi)
agemod.5 <- update(agemod.4, ~. -thick.centerd:group:sex:hemi)
agemod.6 <- update(agemod.5, ~. -thick.centerd:group:sex:age.centerd)
# 3-way
agemod.7 <- update(agemod.6, ~. -age.centerd:thick.centerd:hemi)
agemod.8 <- update(agemod.7, ~. -age.centerd:thick.centerd:group)
agemod.9 <- update(agemod.8, ~. -age.centerd:thick.centerd:sex)
agemod.10 <- update(agemod.9, ~. -age.centerd:sex:hemi)
agemod.11 <- update(agemod.10, ~. -age.centerd:group:hemi)
agemod.12 <- update(agemod.11, ~. -age.centerd:group:sex)
agemod.13 <- update(agemod.12, ~. -group:thick.centerd:sex)
agemod.14 <- update(agemod.13, ~. -group:thick.centerd:hemi)
agemod.15 <- update(agemod.14, ~. -group:hemi:sex)
agemod.16 <- update(agemod.15, ~. -thick.centerd:hemi:sex)
# 2-way
agemod.17 <- update(agemod.16, ~. -age.centerd:sex)
agemod.18 <- update(agemod.17, ~. -age.centerd:group)
agemod.19 <- update(agemod.18, ~. -age.centerd:thick.centerd)
agemod.20 <- update(agemod.19, ~. -age.centerd:hemi)
agemod.21 <- update(agemod.20, ~. -thick.centerd:sex)
agemod.22 <- update(agemod.21, ~. -thick.centerd:hemi)
agemod.23 <- update(agemod.22, ~. -thick.centerd:group)
agemod.24 <- update(agemod.23, ~. -sex:hemi)
agemod.25 <- update(agemod.24, ~. -sex:group)
agemod.26 <- update(agemod.25, ~. -group:hemi)
# single
agemod.27 <- update(agemod.26, ~. -thick.centerd)
agemod.28 <- update(agemod.27, ~. -sex)
agemod.29 <- update(agemod.28, ~. -group)
agemod.30 <- update(agemod.29, ~. -hemi)
agemod.31 <- update(agemod.30, ~. -age.centerd)

anova(agemod.1,
      agemod.2,
      agemod.3,
      agemod.4,
      agemod.5,
      agemod.6,
      agemod.7,
      agemod.8,
      agemod.9,
      agemod.10,
      agemod.11,
      agemod.12,
      agemod.13,
      agemod.14,
      agemod.15,
      agemod.16,
      agemod.17,
      agemod.18,
      agemod.19,
      agemod.20,
      agemod.21,
      agemod.22,
      agemod.23,
      agemod.24,
      agemod.25,
      agemod.26,
      agemod.27,
      agemod.28,
      agemod.29,
      agemod.30,
      agemod.31)

# ONLY ONE HEMI
formula <- nevent.min ~ age.centerd * group * sex * thick.centerd

agemod <- glm(formula,  data=ndata, family=poisson, subset = hemi=="lh")
agemod <- lm(formula,  data=ndata, subset = hemi=="lh")

# 4-way
agemod.1 <- update(agemod, ~. -thick.centerd:group:sex:age.centerd)
# 3-way
agemod.2 <- update(agemod.1, ~. -age.centerd:thick.centerd:group)
agemod.3 <- update(agemod.2, ~. -age.centerd:thick.centerd:sex)
agemod.4 <- update(agemod.3, ~. -age.centerd:group:sex)
agemod.5 <- update(agemod.4, ~. -group:thick.centerd:sex)
# 2-way
agemod.6 <- update(agemod.5, ~. -age.centerd:sex)
agemod.7 <- update(agemod.6, ~. -age.centerd:group)
agemod.8 <- update(agemod.7, ~. -age.centerd:thick.centerd)
agemod.9 <- update(agemod.8, ~. -thick.centerd:sex)
agemod.10 <- update(agemod.9, ~. -thick.centerd:group)
agemod.11 <- update(agemod.10, ~. -sex:group)
# single
agemod.12 <- update(agemod.11, ~. -thick.centerd)
agemod.13 <- update(agemod.12, ~. -sex)
agemod.14 <- update(agemod.13, ~. -group)
agemod.15 <- update(agemod.14, ~. -age.centerd)

anova(agemod.1, agemod)


# BAYES
library(BayesFactor)
library(brms)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.

formula <- nevent.min ~ age.centerd * group * sex * thick.centerd

bmagemod <- brm(formula, data = ndata, family = gaussian, 
                     save_all_pars = TRUE, iter = 5000)
# 4-way
bmagemod.1 <- update(bmagemod, ~. -thick.centerd:group:sex:age.centerd)
# 3-way
bmagemod.2 <- update(bmagemod.1, ~. -age.centerd:thick.centerd:group)
bmagemod.3 <- update(bmagemod.2, ~. -age.centerd:thick.centerd:sex)
bmagemod.4 <- update(bmagemod.3, ~. -age.centerd:group:sex)
bmagemod.5 <- update(bmagemod.4, ~. -group:thick.centerd:sex)
# 2-way
bmagemod.6 <- update(bmagemod.5, ~. -age.centerd:sex)
bmagemod.7 <- update(bmagemod.6, ~. -age.centerd:group)
bmagemod.8 <- update(bmagemod.7, ~. -age.centerd:thick.centerd)
bmagemod.9 <- update(bmagemod.8, ~. -thick.centerd:sex)
bmagemod.10 <- update(bmagemod.9, ~. -thick.centerd:group)
bmagemod.11 <- update(bmagemod.10, ~. -sex:group)
# single
bmagemod.12 <- update(bmagemod.11, ~. -thick.centerd)
bmagemod.13 <- update(bmagemod.12, ~. -sex)
bmagemod.14 <- update(bmagemod.13, ~. -group)
bmagemod.15 <- update(bmagemod.14, ~. -age.centerd)

bayes_factor(bmagemod,bmagemod.1)
bayes_factor(bmagemod.1,bmagemod.2)
bayes_factor(bmagemod.2,bmagemod.3)
bayes_factor(bmagemod.3,bmagemod.4)
bayes_factor(bmagemod.4,bmagemod.5)
bayes_factor(bmagemod.5,bmagemod.6)
bayes_factor(bmagemod.6,bmagemod.7)
bayes_factor(bmagemod.7,bmagemod.8)
bayes_factor(bmagemod.8,bmagemod.9)
bayes_factor(bmagemod.9,bmagemod.10)
bayes_factor(bmagemod.10,bmagemod.11)
bayes_factor(bmagemod.11,bmagemod.12)
bayes_factor(bmagemod.12,bmagemod.13)
bayes_factor(bmagemod.13,bmagemod.14)
bayes_factor(bmagemod.14,bmagemod.15)

# Hypothesis testing
hypothesis(bmagemod.9, "age.centerd>0")    
hypothesis(bmagemod.9, "grouppatient>0")
hypothesis(bmagemod.9, "age.centerd:grouppatient<0")
hypothesis(bmagemod.5, "age.centerd:sexM<0")
hypothesis(bmagemod.5, "age.centerd:grouppatient:sexM>0")
hypothesis(bmagemod, "sexM<0")
hypothesis(bmagemod, "sexM<0")
hypothesis(bmagemod, "sexM<0")
hypothesis(bmagemod, "sexM<0")
hypothesis(bmagemod, "sexM<0")



aa <- lmBF(nevent.min ~ age.centerd*group*sex*thick.centerd, data = ndata)
bb <- lmBF(nevent.min ~ age.centerd*group*sex*thick.centerd-age.centerd:thick.centerd:group, data = ndata)
bb2 <- lmBF(nevent.min ~ age.centerd+group+sex+thick.centerd+
              age.centerd:group + age.centerd:sex + age.centerd:thick.centerd + group:sex + group:thick.centerd + sex:thick.centerd+
              age.centerd:thick.centerd:group + age.centerd:thick.centerd:sex + thick.centerd:age:sex+
              age.centerd:thick.centerd:group:sex, data = ndata)


## MISC

nmod <- lmer(nevent.min ~  age.centerd * sex * group + thick.centerd*age.centerd*group + (1|hemi) + (1|subj), data=ndata)
nmod.x3 <- lmer(nevent.min ~ I(age-mean(age)) * sex * group + (1|hemi) + (1|subj), data=ndata)

nplot <- ggplot(aes(x=age, y=nevent, color=group, shape=sex), data=ndata)+
  geom_point()+
  geom_smooth(method=lm)
nplot

t.test(ndata$nevent.min[ndata$group=="patient"], ndata$nevent.min[ndata$group=="control"])
''