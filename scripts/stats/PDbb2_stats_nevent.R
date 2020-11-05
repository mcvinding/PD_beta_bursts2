# PD beta burst statistics: analysis of N events
# CLEAN UP!!!!
library(lme4)
library(arm)
library(ggplot2)
library(BayesFactor)
library(brms)
# Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
# Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
# Bayes?

## Load data
# load('X://PD_longrest//groupanalysis//alldata_all.Rdata')
load('C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj.Rdata')

## Center variables
alldata$age.centerd <- alldata$age-mean(alldata$age)
alldata$age.min <- alldata$age-min(alldata$age)
alldata$thick.centerd <- alldata$thick-mean(alldata$thick)

# Inspect hist
# ggplot( aes(x=nevent.b.m2, fill=group), data=alldata) +
#   geom_histogram(color="black", alpha=0.6, position = 'identity', bins=50)
ggplot( aes(x=nevent.u.m2, fill=group), data=alldata) +
  geom_histogram(color="black", alpha=0.6, position = 'identity', bins=50)

# Group diff
# t.test(nevent.b.m2~group, data=alldata)
t.test(nevent.u.m2~group, data=alldata)

# Inspect ~age
ggplot(aes(x=age, y=nevent.b.m2, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_smooth(method=lm)
ggplot(aes(x=age, y=nevent.u.m2, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_smooth(method=lm)

## Polynomial variables
# nmod.b.m2 <- glm(nevent.b.m2 ~ I(age.centerd^2)*age.centerd*sex*group, data=alldata, family=poisson)
nmod.u.m2 <- glm(nevent.u.m2 ~ I(age.centerd^2)*age.centerd*sex*group, data=alldata, family=poisson)
# summary(nmod.b.m2)
summary(nmod.u.m2)

anova(nmod.u.m2, test="Chisq")

new.dat <- alldata
new.dat$pred <- exp(predict(nmod.b.m2, re.form=NA))

ggplot(aes(x=age, y=nevent.b.m2, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)

######################################################################################
# TEST
tstmod.a <- glm(nevent.u.m2 ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=poisson)
tstmod2.1 <- update(tstmod.a, ~. -group:sex:I(age.centerd^2))
tstmod2.2 <- update(tstmod2.1, ~. -sex:I(age.centerd^2))
tstmod2.3 <- update(tstmod2.2, ~. -group:I(age.centerd^2))
tstmod2.4 <- update(tstmod2.3, ~. -group:sex:age.centerd)
tstmod2.5 <- update(tstmod2.4, ~. -sex:age.centerd)
tstmod2.6 <- update(tstmod2.5, ~. -group:age.centerd)
tstmod2.7 <- update(tstmod2.6, ~. -group:sex)
tstmod2.8 <- update(tstmod2.7, ~. -sex)
tstmod2.9 <- update(tstmod2.8, ~. -group)
tstmod2.10 <- update(tstmod2.9, ~. -I(age.centerd^2))
tstmod2.11 <- update(tstmod2.10, ~. -age.centerd)

anova(tstmod2.11,tstmod2.10,tstmod2.9,tstmod2.8,tstmod2.7,tstmod2.6,tstmod2.5,tstmod2.4,tstmod2.3,tstmod2.2,tstmod2.1,tstmod.a,
      test="Chisq")

# anova(tstmod.a, test="Chisq")

new.dat <- alldata
new.dat$pred <- exp(predict(tstmod2.4, re.form=NA))
ggplot(aes(x=age, y=nevent.u.m2, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)

######################################################################################
## THICK
tckmod <- lm(nevent.u.m2 ~ thick.centerd + I(thick.centerd^2) +
              group + sex + 
              group:sex +
              group:thick.centerd +
              sex:thick.centerd +
              group:sex:thick.centerd +
              group:I(thick.centerd^2) +
              sex:I(thick.centerd^2) +
              group:sex:I(thick.centerd^2),
            data=alldata)
tckmod2.1 <- update(tckmod, ~. -group:sex:I(thick.centerd^2))
tckmod2.2 <- update(tckmod2.1, ~. -sex:I(thick.centerd^2))
tckmod2.3 <- update(tckmod2.2, ~. -group:I(thick.centerd^2))
tckmod2.4 <- update(tckmod2.3, ~. -group:sex:thick.centerd)
tckmod2.5 <- update(tckmod2.4, ~. -sex:thick.centerd)
tckmod2.6 <- update(tckmod2.5, ~. -group:thick.centerd)
tckmod2.7 <- update(tckmod2.6, ~. -group:sex)
tckmod2.8 <- update(tckmod2.7, ~. -sex)
tckmod2.9 <- update(tckmod2.8, ~. -group)
tckmod2.10 <- update(tckmod2.9, ~. -I(thick.centerd^2))
tckmod2.11 <- update(tckmod2.10, ~. -thick.centerd)

anova(tckmod2.11,tckmod2.10,
      tckmod2.9,
      tckmod2.8, 
      tckmod2.7, 
      tckmod2.6, 
      tckmod2.5, 
      tckmod2.4, 
      tckmod2.3, 
      tckmod2.2,
      tckmod2.1,
      tckmod, test="F")

new.dat <- alldata
new.dat$pred <- predict(tckmod2.2, re.form=NA)
ggplot(aes(x=thick, y=nevent.u.m2, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)

###################################################################################
# Not properly implemented
# BAYES
# BF
a.a <- lmBF(nevent.b.m2 ~ age + group +  sex + 
              group:sex +
              group:age +
              sex:age +
              group:sex:age, 
            data = alldata)
a.b <- lmBF(nevent.b.m2 ~ age + group +  sex + 
              group:sex +
              group:age +
              sex:age, 
            data = alldata)
a.c <- lmBF(nevent.b.m2 ~ age + group +  sex + 
              group:sex +
              group:age, 
            data = alldata)
a.d <- lmBF(nevent.b.m2 ~ age + group +  sex + 
              group:sex, 
            data = alldata)
a.e <- lmBF(nevent.b.m2 ~ age + group +  sex, 
            data = alldata)
a.f <- lmBF(nevent.b.m2 ~ age + group, 
            data = alldata)
a.g <- lmBF(nevent.b.m2 ~ age, 
            data = alldata)


t.a <- lmBF(nevent ~ thick + group +  sex + 
              group:sex +
              group:thick +
              sex:thick +
              group:sex:thick, 
            data = nlh)
t.b <- lmBF(nevent ~ thick + group +  sex + 
              group:sex +
              group:thick +
              sex:thick, 
            data = nlh)
t.c <- lmBF(nevent ~ thick + group +  sex + 
              group:sex +
              group:thick, 
            data = nlh)
t.d <- lmBF(nevent ~ thick + group +  sex + 
              group:sex, 
            data = nlh)
t.e <- lmBF(nevent ~ thick + group +  sex, 
            data = nlh)
t.f <- lmBF(nevent ~ thick + group, 
            data = nlh)
t.g <- lmBF(nevent ~ thick, 
            data = nlh)

## brms
set.seed(666)
brmod.a <- brm(nevent ~ age + I(age^2) + group +  sex + 
                  group:sex +
                  group:age +
                  sex:age +
                  group:sex:age +
                  group:I(age^2) +
                  sex:I(age^2) +
                  group:sex:I(age^2),
                data=nlh, family=poisson,
                save_all_pars = TRUE, iter = 1000)

brmod.b <- brm(nevent ~ age + I(age^2) + group +  sex + 
                  group:sex +
                  group:age +
                  sex:age +
                  group:sex:age +
                  group:I(age^2) +
                  sex:I(age^2),
                data=nlh, family=poisson,
                save_all_pars = TRUE, iter = 1000)

br.nev2 <- brm(nevent.min ~ group+session+(1|subs), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)
br.nev1 <- brm(nevent.min ~ session+(1|subs), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)
br.nev0 <- brm(nevent.min ~ 1+(1|subs), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)




# BAYES
library(BayesFactor)
library(brms)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.

formula <- nevent.min ~ age.centerd * group * sex * thick.centerd

bmagemod <- brm(formula, data = alldata, family = gaussian, 
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



aa <- lmBF(nevent.min ~ age.centerd*group*sex*thick.centerd, data = alldata)
bb <- lmBF(nevent.min ~ age.centerd*group*sex*thick.centerd-age.centerd:thick.centerd:group, data = alldata)
bb2 <- lmBF(nevent.min ~ age.centerd+group+sex+thick.centerd+
              age.centerd:group + age.centerd:sex + age.centerd:thick.centerd + group:sex + group:thick.centerd + sex:thick.centerd+
              age.centerd:thick.centerd:group + age.centerd:thick.centerd:sex + thick.centerd:age:sex+
              age.centerd:thick.centerd:group:sex, data = alldata)
