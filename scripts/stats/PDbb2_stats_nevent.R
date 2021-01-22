# PD beta burst statistics: analysis of N events
# CLEAN UP!!!!
library(lme4)
library(arm)
library(ggplot2)
library(BayesFactor)
library(brms)
library(car)
# Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
# Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
# Bayes?

## Load data
load('X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
<<<<<<< HEAD
# load('C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj2.Rdata')
=======
# load('C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj.Rdata')
>>>>>>> master

## Center variables
alldata$age.centerd <- alldata$age-mean(alldata$age)
# alldata$age.min <- alldata$age-min(alldata$age)
alldata$thick.centerd <- alldata$thick-mean(alldata$thick)

# Inspect hist
ggplot( aes(x=nevent.u.m2.min, fill=group), data=alldata) +
  geom_histogram(color="black", alpha=0.6, position = 'identity', bins=25)

# Inspect ~age
ggplot(aes(x=age, y=nevent.u.m2.min, color=group, shape=sex), data=alldata)+
  geom_point()+
  geom_smooth(method=lm)

######################################################################################
# LMER regression model
<<<<<<< HEAD
tstmod.8 <- glm(nevent.u.m2.min ~ (age.centerd+sex+group+thick.centerd)^2, data=alldata, family=poisson)

tstmod.7 <- update(tstmod.8, ~. -group:sex:age.centerd)
tstmod.6 <- update(tstmod.7, ~. -sex:age.centerd)
tstmod.5 <- update(tstmod.6, ~. -group:sex)
tstmod.4 <- update(tstmod.5, ~. -group:age.centerd)
tstmod.3 <- update(tstmod.4, ~. -sex)

tstmod.2 <- update(tstmod.3, ~. -age.centerd)
tstmod.1 <- update(tstmod.2, ~. -group)

anova(tstmod.1,tstmod.2,tstmod.3,tstmod.4,tstmod.5,tstmod.6,tstmod.7,tstmod.8, test="Chisq")

anova(tstmod.8, test="Chisq")
summary(tstmod.8, test="Chisq")
vif(tstmod.8)
=======
tstmod.Full3 <- glm(nevent.u.m2.min ~ (age.centerd+sex+group+thick.centerd)^3, data=alldata, family=poisson)

tstmod.ASG <- update(tstmod.Full3, ~. -group:sex:age.centerd)
tstmod.AST <- update(tstmod.Full3, ~. -group:sex:thick.centerd)
tstmod.AGT <- update(tstmod.Full3, ~. -group:age.centerd:thick.centerd)
tstmod.GST <- update(tstmod.Full3, ~. -group:sex:thick.centerd)
anova(tstmod.Full3,tstmod.ASG, test="Chisq")
anova(tstmod.Full3,tstmod.AST, test="Chisq")
anova(tstmod.Full3,tstmod.AGT, test="Chisq")
anova(tstmod.Full3,tstmod.GST, test="Chisq")

tstmod.Full2 <- glm(nevent.u.m2.min ~ (age.centerd+sex+group+thick.centerd)^2, data=alldata, family=poisson)

tstmod.SA <- update(tstmod.Full2, ~. -sex:age.centerd)
tstmod.GS <- update(tstmod.Full2, ~. -group:sex)
tstmod.GA <- update(tstmod.Full2, ~. -group:age.centerd)
tstmod.GT <- update(tstmod.Full2, ~. -sex:thick.centerd)
tstmod.ST <- update(tstmod.Full2, ~. -sex:thick.centerd)
tstmod.AT <- update(tstmod.Full2, ~. -age.centerd:thick.centerd)
anova(tstmod.Full2,tstmod.SA, test="Chisq")
anova(tstmod.Full2,tstmod.GS, test="Chisq")
anova(tstmod.Full2,tstmod.GA, test="Chisq")
anova(tstmod.Full2,tstmod.GT, test="Chisq")
anova(tstmod.Full2,tstmod.ST, test="Chisq")
anova(tstmod.Full2,tstmod.AT, test="Chisq")

tstmod.Full1 <- glm(nevent.u.m2.min ~ age.centerd+sex+group+thick.centerd, data=alldata, family=poisson)
tstmod.A <- update(tstmod.Full1, ~. -age.centerd)
tstmod.S <- update(tstmod.Full1, ~. -sex)
tstmod.G <- update(tstmod.Full1, ~. -group)
tstmod.T <- update(tstmod.Full1, ~. -thick.centerd)
anova(tstmod.Full1,tstmod.A, test="Chisq")
anova(tstmod.Full1,tstmod.S, test="Chisq")
anova(tstmod.Full1,tstmod.G, test="Chisq")
anova(tstmod.Full1,tstmod.T, test="Chisq")


anova(tstmod.Full3, test="Chisq")
summary(tstmod.Full3, test="Chisq")
vif(tstmod.Full)
>>>>>>> master





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