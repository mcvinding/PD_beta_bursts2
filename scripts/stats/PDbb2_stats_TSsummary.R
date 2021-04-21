# PD beta burst statistics: analysis of skewness and kurtosis of beta time series
library(ggplot2)

# Load data
load(file='X://PD_longrest//groupanalysis//alldata_subj.Rdata')
alldata$age.centerd <- alldata$age-mean(alldata$age)

# NB!!!! FIND out what to do whith these subjects!
alldata$skw.u[alldata$skw.u>3] = NA
alldata$krt.u[alldata$krt.u>25] = NA

# Inspect
ggplot(aes(x=age, y=skw.u, color=group, shape=sex), data=alldata)+
  geom_point()+geom_smooth(method=lm)+
  theme_bw()+ggtitle('Mu+beta TS skewness ~ age')

ggplot(aes(x=age, y=krt.u, color=group, shape=sex), data=alldata)+
  geom_point()+geom_smooth(method=lm)+
  theme_bw()+ggtitle('Mu+beta TS kurtosis ~ age')

######################################################################################
# Skewness
tstmod.a <- glm(skw.u ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

new.dat <- alldata[!is.na(alldata$skw.u),]
new.dat$pred <- predict(tstmod2.4, re.form=NA)
ggplot(aes(x=age, y=skw.u, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)


######################################################################################
# Kurtosis
tstmod.a <- glm(krt.u ~ I(age.centerd^2)*sex*group+age.centerd*sex*group, data=alldata, family=gaussian)
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

new.dat <- alldata[!is.na(alldata$krt.u),]
new.dat$pred <- predict(tstmod2.10, re.form=NA)
ggplot(aes(x=age, y=krt.u, color=group, shape=sex), data=new.dat)+
  geom_point()+
  geom_line(aes(y = pred, linetype=sex), size = 1)

#END