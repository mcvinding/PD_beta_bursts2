# PD beta burst statistics: analsyis of PSD output from FOOOF analysis
library(ggplot2)

# Load data
load(file='X://PD_longrest//groupanalysis//alldata_subj.Rdata')
load('C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj.Rdata')
alldata$age.centerd <- alldata$age-mean(alldata$age)
alldata$thick.centerd <- alldata$thick-mean(alldata$thick)

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
tstmod.a <- glm(a_intercept ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

# Plot winning model
new.dat <- alldata
new.dat$pred <- predict(tstmod2.8, re.form=NA)
ggplot(aes(x=age, y=a_intercept, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)

######################################################################################
# 1/f slope
tstmod.a <- glm(a_slope ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

# Plot winning model
new.dat <- alldata
new.dat$pred <- predict(tstmod2.8, re.form=NA)
ggplot(aes(x=age, y=a_slope, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)

######################################################################################
# Beta power
tstmod.a <- glm(log(beta_pw) ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

# Plot winning model
new.dat <- alldata
new.dat$pred <- exp(predict(tstmod2.8, re.form=NA))
ggplot(aes(x=age, y=beta_pw, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)

######################################################################################
# Beta peak freq
tstmod.a <- glm(beta_cf ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

######################################################################################
# Alpha power
tstmod.a <- glm(log(alpha_pw) ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

# Plot winning model
new.dat <- alldata[!is.na(alldata$alpha_bw),]
new.dat$pred <- exp(predict(tstmod2.8, re.form=NA))
ggplot(aes(x=age, y=alpha_pw, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)


######################################################################################
# Alpha peak freq
tstmod.a <- glm(alpha_cf ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

new.dat <- alldata[!is.na(alldata$alpha_bw),]
new.dat$pred <- predict(tstmod2.5, re.form=NA)
ggplot(aes(x=age, y=alpha_cf, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)


######################################################################################
# ~ Thickness
######################################################################################

######################################################################################
# 1/f intercept
tstmod.a <- glm(a_intercept ~ I(thick.centerd^2)*sex*group+thick.centerd*sex*group, data=alldata, family=gaussian)
tstmod2.1 <- update(tstmod.a, ~. -group:sex:I(thick.centerd^2))
tstmod2.2 <- update(tstmod2.1, ~. -sex:I(thick.centerd^2))
tstmod2.3 <- update(tstmod2.2, ~. -group:I(thick.centerd^2))
tstmod2.4 <- update(tstmod2.3, ~. -group:sex:thick.centerd)
tstmod2.5 <- update(tstmod2.4, ~. -sex:thick.centerd)
tstmod2.6 <- update(tstmod2.5, ~. -group:thick.centerd)
tstmod2.7 <- update(tstmod2.6, ~. -group:sex)
tstmod2.8 <- update(tstmod2.7, ~. -sex)
tstmod2.9 <- update(tstmod2.8, ~. -group)
tstmod2.10 <- update(tstmod2.9, ~. -I(thick.centerd^2))
tstmod2.11 <- update(tstmod2.10, ~. -thick.centerd)

anova(tstmod2.11,tstmod2.10,tstmod2.9,tstmod2.8,tstmod2.7,tstmod2.6,tstmod2.5,tstmod2.4,tstmod2.3,tstmod2.2,tstmod2.1,tstmod.a,
      test="Chisq")

######################################################################################
# 1/f slope
tstmod.a <- glm(a_slope ~ I(thick.centerd^2)*sex*group+thick.centerd*sex*group, data=alldata, family=gaussian)
tstmod2.1 <- update(tstmod.a, ~. -group:sex:I(thick.centerd^2))
tstmod2.2 <- update(tstmod2.1, ~. -sex:I(thick.centerd^2))
tstmod2.3 <- update(tstmod2.2, ~. -group:I(thick.centerd^2))
tstmod2.4 <- update(tstmod2.3, ~. -group:sex:thick.centerd)
tstmod2.5 <- update(tstmod2.4, ~. -sex:thick.centerd)
tstmod2.6 <- update(tstmod2.5, ~. -group:thick.centerd)
tstmod2.7 <- update(tstmod2.6, ~. -group:sex)
tstmod2.8 <- update(tstmod2.7, ~. -sex)
tstmod2.9 <- update(tstmod2.8, ~. -group)
tstmod2.10 <- update(tstmod2.9, ~. -I(thick.centerd^2))
tstmod2.11 <- update(tstmod2.10, ~. -thick.centerd)

anova(tstmod2.11,tstmod2.10,tstmod2.9,tstmod2.8,tstmod2.7,tstmod2.6,tstmod2.5,tstmod2.4,tstmod2.3,tstmod2.2,tstmod2.1,tstmod.a,
      test="Chisq")

######################################################################################
# Beta power
tstmod.a <- glm(log(beta_pw) ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

######################################################################################
# Beta peak freq
tstmod.a <- glm(beta_cf ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

######################################################################################
# Alpha power
tstmod.a <- glm(log(alpha_pw) ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

######################################################################################
# Alpha peak freq
tstmod.a <- glm(alpha_cf ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

t <- anova(tstmod2.11,tstmod2.10,tstmod2.9,tstmod2.8,tstmod2.7,tstmod2.6,tstmod2.5,tstmod2.4,tstmod2.3,tstmod2.2,tstmod2.1,tstmod.a,
           test="Chisq")
#END