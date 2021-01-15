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
load('X://PD_longrest//groupanalysis//alldata_subj.Rdata')
# load('C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj.Rdata')

## Center variables
alldata$age.centerd <- alldata$age-mean(alldata$age)
alldata$age.min <- alldata$age-min(alldata$age)
alldata$thick.centerd <- alldata$thick-mean(alldata$thick)

# Inspect hist
# ggplot( aes(x=nevent.b.m2, fill=group), data=alldata) +
#   geom_histogram(color="black", alpha=0.6, position = 'identity', bins=50)
ggplot( aes(x=nevent.u.m2, fill=group), data=alldata) +
  geom_histogram(color="black", alpha=0.6, position = 'identity', bins=50)

# CLEAN UP EVERYTHING IN THIS SECTION
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

new.dat <- alldata
new.dat$pred <- exp(predict(tstmod2.4, re.form=NA))
ggplot(aes(x=age, y=nevent.u.m2, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)


######################################################################################
# No square term
tstmod.8 <- glm(nevent.u.m2 ~ age.centerd*sex*group, data=alldata, family=poisson)
tstmod.7 <- update(tstmod.8, ~. -group:sex:age.centerd)
tstmod.6 <- update(tstmod.7, ~. -sex:age.centerd)
tstmod.5 <- update(tstmod.6, ~. -group:sex)
tstmod.4 <- update(tstmod.5, ~. -group:age.centerd)
tstmod.3 <- update(tstmod.4, ~. -sex)
tstmod.2 <- update(tstmod.3, ~. -age.centerd)
tstmod.1 <- update(tstmod.2, ~. -group)

anova(tstmod.1,tstmod.2,tstmod.3,tstmod.4,tstmod.5,tstmod.6,tstmod.7,tstmod.8, test="Chisq")

# Plot top model
new.dat <- alldata
new.dat$pred <- exp(predict(tstmod.8, re.form=NA))
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
a.7 <- lmBF(nevent.b.m2 ~ age.centerd + group +  sex + 
              group:sex +
              group:age.centerd +
              sex:age.centerd +
              group:sex:age.centerd, 
            data = alldata)
a.6 <- lmBF(nevent.b.m2 ~ age.centerd + group +  sex + 
              group:sex +
              group:age.centerd +
              sex:age.centerd, 
            data = alldata)
a.5 <- lmBF(nevent.b.m2 ~ age.centerd + group +  sex + 
              group:sex +
              group:age.centerd, 
            data = alldata)
a.4 <- lmBF(nevent.b.m2 ~ age.centerd + group +  sex + 
              group:sex, 
            data = alldata)
a.3 <- lmBF(nevent.b.m2 ~ age.centerd + group +  sex, 
            data = alldata)
a.2 <- lmBF(nevent.b.m2 ~ age.centerd + group, 
            data = alldata)
a.1 <- lmBF(nevent.b.m2 ~ group, 
            data = alldata)



aa <- lmBF(nevent.min ~ age.centerd*group*sex*thick.centerd, data = alldata)
bb <- lmBF(nevent.min ~ age.centerd*group*sex*thick.centerd-age.centerd:thick.centerd:group, data = alldata)
bb2 <- lmBF(nevent.min ~ age.centerd+group+sex+thick.centerd+
              age.centerd:group + age.centerd:sex + age.centerd:thick.centerd + group:sex + group:thick.centerd + sex:thick.centerd+
              age.centerd:thick.centerd:group + age.centerd:thick.centerd:sex + thick.centerd:age:sex+
              age.centerd:thick.centerd:group:sex, data = alldata)

##############################################################################
# BRMS
library(brms)
br.nev3 <- brm(nevent.u.m2 ~ age.centerd*sex*group, data=alldata, family=poisson, 
               save_all_pars=TRUE, iter = 2000)

hypothesis(br.nev3, "grouppatient>0")                         # Ptns vs Ctrl
hypothesis(br.nev3, "age.centerd>0")            
hypothesis(br.nev3, "sexM>0")                                 #*
hypothesis(br.nev3, "age.centerd:grouppatient>0")
hypothesis(br.nev3, "age.centerd:sexM<0")                     #*
hypothesis(br.nev3, "sexM:grouppatient>0")                    #*
hypothesis(br.nev3, "age.centerd:sexM:grouppatient>0")



new.dat <- alldata
new.dat$pred <- predict(br.nev3)
new.dat$preds <- new.dat$pred[,1]
ggplot(aes(x=age, y=nevent.u.m2, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y=preds, linetype=sex), size = 1)




tstmod.7 <- update(tstmod.8, ~. -group:sex:age.centerd)
tstmod.6 <- update(tstmod.7, ~. -sex:age.centerd)
tstmod.5 <- update(tstmod.6, ~. -group:sex)
tstmod.4 <- update(tstmod.5, ~. -group:age.centerd)
tstmod.3 <- update(tstmod.4, ~. -sex)
tstmod.2 <- update(tstmod.3, ~. -age.centerd)
tstmod.1 <- update(tstmod.2, ~. -group)