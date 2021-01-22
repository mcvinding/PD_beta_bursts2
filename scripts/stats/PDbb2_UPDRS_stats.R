# UPDRS stats
library(ggplot2)

# Load data
load(file='X://PD_longrest//groupanalysis//clindata.Rdata')
#load(file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj.Rdata')

pd.data <- subset(alldata, group == 'patient')
pd.data <- subset(pd.data, !is.na(U.F7))
pd.data <- subset(pd.data, !is.na(U.F6))
pd.data <- subset(pd.data, !is.na(U.F3))
pd.data <- subset(pd.data, !is.na(U.F1))
pd.data$U.F45 <- pd.data$U.F4+pd.data$U.F5

# UPDRS ~ Neve
umod.u.m2 <- glm(UPDRS ~ nevent.u.m2.min * age * sex, data=pd.data, family=gaussian)
summary(umod.u.m2)

pd.data$pred <- predict(umod.u.m2, re.form=NA)
ggplot( aes(y=UPDRS, x=nevent.u.m2, color=sex), data=pd.data)+
  geom_point()+
  # geom_smooth(method='lm')+
  geom_line(aes(y = pred, linetype=sex), size = 1)
anova(umod.u.m2, test="Chisq")

######################################################################################
# UPDRS  Neve
nmod.u.0x <- glm(nevent.u.m2~age+sex, data=pd.data, family=poisson)

nmod.u.T <- glm(nevent.u.m2~UPDRS+age+sex, data=pd.data, family=poisson)
anova(nmod.u.0, nmod.u.T, test="Chisq")

umod.F1.Full <- glm(U.F1~nevent.u.m2.min+age+sex+thick.centerd+alpha_pw+beta_pw+a_slope+a_intercept+beta_cf+alpha_cf+
                      leneve+tueeve+maxeve,
                    data=pd.data, family=gaussian)
umod.F1.neve <- update(umod.F1.Full, .~. -nevent.u.m2.min)
umod.F1.age <- update(umod.F1.Full, .~. -age)
umod.F1.sex <- update(umod.F1.Full, .~. -sex)
umod.F1.thic <- update(umod.F1.Full, .~. -thick.centerd)
umod.F1.alpw <- update(umod.F1.Full, .~. -alpha_pw)
umod.F1.bepw <- update(umod.F1.Full, .~. -beta_pw)
umod.F1.slop <- update(umod.F1.Full, .~. -a_slope)
umod.F1.itcp <- update(umod.F1.Full, .~. -a_intercept)
umod.F1.becf <- update(umod.F1.Full, .~. -beta_cf)
umod.F1.alcf <- update(umod.F1.Full, .~. -alpha_cf)
umod.F1.lene <- update(umod.F1.Full, .~. -leneve)
umod.F1.tuee <- update(umod.F1.Full, .~. -tueeve)
umod.F1.maxe <- update(umod.F1.Full, .~. -maxeve)
anova(umod.F1.neve, umod.F1.Full, test="F")
anova(umod.F1.Full,umod.F1.age, test="F")
anova(umod.F1.Full,umod.F1.sex, test="F")
anova(umod.F1.Full,umod.F1.thic, test="F")
anova(umod.F1.Full,umod.F1.alpw, test="F")
anova(umod.F1.Full,umod.F1.bepw, test="F")
anova(umod.F1.Full,umod.F1.slop, test="F")
anova(umod.F1.Full,umod.F1.itcp, test="F")
anova(umod.F1.Full,umod.F1.becf, test="F")
anova(umod.F1.Full,umod.F1.alcf, test="F")
anova(umod.F1.Full,umod.F1.lene, test="F")
anova(umod.F1.Full,umod.F1.tuee, test="F")
anova(umod.F1.Full,umod.F1.maxe, test="F")

anova(umod.F1.Full, test="Chisq")



umod.F2.Full <- glm(U.F2~nevent.u.m2.min+age+sex+thick.centerd+alpha_pw+beta_pw+a_slope+a_intercept+beta_cf+alpha_cf+
                      leneve+tueeve+maxeve, data=pd.data, family=gaussian)
umod.F2.neve <- update(umod.F2.Full, .~. -nevent.u.m2.min)
umod.F2.age <- update(umod.F2.Full, .~. -age)
umod.F2.sex <- update(umod.F2.Full, .~. -sex)
umod.F2.thic <- update(umod.F2.Full, .~. -thick.centerd)
umod.F2.alpw <- update(umod.F2.Full, .~. -alpha_pw)
umod.F2.bepw <- update(umod.F2.Full, .~. -beta_pw)
umod.F2.slop <- update(umod.F2.Full, .~. -a_slope)
umod.F2.itcp <- update(umod.F2.Full, .~. -a_intercept)
umod.F2.becf <- update(umod.F2.Full, .~. -beta_cf)
umod.F2.alcf <- update(umod.F2.Full, .~. -alpha_cf)
umod.F2.lene <- update(umod.F2.Full, .~. -leneve)
umod.F2.tuee <- update(umod.F2.Full, .~. -tueeve)
umod.F2.maxe <- update(umod.F2.Full, .~. -maxeve)
anova(umod.F2.neve, umod.F2.Full, test="Chisq")
anova(umod.F2.age,umod.F2.Full, test="Chisq")
anova(umod.F2.sex,umod.F2.Full, test="Chisq")
anova(umod.F2.Full,umod.F2.thic, test="Chisq")
anova(umod.F2.Full,umod.F2.alpw, test="Chisq")
anova(umod.F2.Full,umod.F2.bepw, test="Chisq")
anova(umod.F2.Full,umod.F2.slop, test="Chisq")
anova(umod.F2.Full,umod.F2.itcp, test="Chisq")
anova(umod.F2.Full,umod.F2.becf, test="Chisq")
anova(umod.F2.Full,umod.F2.alcf, test="Chisq")
anova(umod.F2.Full,umod.F2.lene, test="Chisq")
anova(umod.F2.Full,umod.F2.tuee, test="Chisq")
anova(umod.F2.Full,umod.F2.maxe, test="Chisq")

anova(umod.F2.Full, test="F")



umod.F3.Full <- glm(U.F3~nevent.u.m2.min+age+sex+thick.centerd+alpha_pw+beta_pw+a_slope+a_intercept+beta_cf+alpha_cf+
                      leneve+tueeve+maxeve, data=pd.data, family=gaussian)
umod.F3.neve <- update(umod.F3.Full, .~. -nevent.u.m2.min)
umod.F3.age <- update(umod.F3.Full, .~. -age)
umod.F3.sex <- update(umod.F3.Full, .~. -sex)
umod.F3.thic <- update(umod.F3.Full, .~. -thick.centerd)
umod.F3.alpw <- update(umod.F3.Full, .~. -alpha_pw)
umod.F3.bepw <- update(umod.F3.Full, .~. -beta_pw)
umod.F3.slop <- update(umod.F3.Full, .~. -a_slope)
umod.F3.itcp <- update(umod.F3.Full, .~. -a_intercept)
umod.F3.becf <- update(umod.F3.Full, .~. -beta_cf)
umod.F3.alcf <- update(umod.F3.Full, .~. -alpha_cf)
umod.F3.lene <- update(umod.F3.Full, .~. -leneve)
umod.F3.tuee <- update(umod.F3.Full, .~. -tueeve)
umod.F3.maxe <- update(umod.F3.Full, .~. -maxeve)
anova(umod.F3.neve, umod.F3.Full, test="Chisq")
anova(umod.F3.age,umod.F3.Full, test="Chisq")
anova(umod.F3.sex,umod.F3.Full, test="Chisq")
anova(umod.F3.Full,umod.F3.thic, test="Chisq")
anova(umod.F3.Full,umod.F3.alpw, test="Chisq")
anova(umod.F3.Full,umod.F3.bepw, test="Chisq")
anova(umod.F3.Full,umod.F3.slop, test="Chisq")
anova(umod.F3.Full,umod.F3.itcp, test="Chisq")
anova(umod.F3.Full,umod.F3.becf, test="Chisq")
anova(umod.F3.Full,umod.F3.alcf, test="Chisq")
anova(umod.F3.Full,umod.F3.lene, test="Chisq")
anova(umod.F3.Full,umod.F3.tuee, test="Chisq")
anova(umod.F3.Full,umod.F3.maxe, test="Chisq")

anova(umod.F3.Full, test="F")

umod.F4.Full <- glm(U.F4~nevent.u.m2.min+age+sex+thick.centerd+alpha_pw+beta_pw+a_slope+a_intercept+beta_cf+alpha_cf+
                       leneve+tueeve+maxeve, data=pd.data, family=gaussian)
umod.F4.neve <- update(umod.F4.Full, .~. -nevent.u.m2.min)
umod.F4.age <- update(umod.F4.Full, .~. -age)
umod.F4.sex <- update(umod.F4.Full, .~. -sex)
umod.F4.thic <- update(umod.F4.Full, .~. -thick.centerd)
umod.F4.alpw <- update(umod.F4.Full, .~. -alpha_pw)
umod.F4.bepw <- update(umod.F4.Full, .~. -beta_pw)
umod.F4.slop <- update(umod.F4.Full, .~. -a_slope)
umod.F4.itcp <- update(umod.F4.Full, .~. -a_intercept)
umod.F4.becf <- update(umod.F4.Full, .~. -beta_cf)
umod.F4.alcf <- update(umod.F4.Full, .~. -alpha_cf)
umod.F4.lene <- update(umod.F4.Full, .~. -leneve)
umod.F4.tuee <- update(umod.F4.Full, .~. -tueeve)
umod.F4.maxe <- update(umod.F4.Full, .~. -maxeve)
anova(umod.F4.neve, umod.F4.Full, test="Chisq")
anova(umod.F4.age,umod.F4.Full, test="Chisq")
anova(umod.F4.sex,umod.F4.Full, test="Chisq")
anova(umod.F4.Full,umod.F4.thic, test="Chisq")
anova(umod.F4.Full,umod.F4.alpw, test="Chisq")
anova(umod.F4.Full,umod.F4.bepw, test="Chisq")
anova(umod.F4.Full,umod.F4.slop, test="Chisq")
anova(umod.F4.Full,umod.F4.itcp, test="Chisq")
anova(umod.F4.Full,umod.F4.becf, test="Chisq")
anova(umod.F4.Full,umod.F4.alcf, test="Chisq")
anova(umod.F4.Full,umod.F4.lene, test="Chisq")
anova(umod.F4.Full,umod.F4.tuee, test="Chisq")
anova(umod.F4.Full,umod.F4.maxe, test="Chisq")




umod.F45.Full <- glm(U.F45~nevent.u.m2.min+age+sex+thick.centerd+alpha_pw+beta_pw+a_slope+a_intercept+beta_cf+alpha_cf+
                      leneve+tueeve+maxeve, data=pd.data, family=gaussian)
umod.F45.neve <- update(umod.F45.Full, .~. -nevent.u.m2.min)
umod.F45.age <- update(umod.F45.Full, .~. -age)
umod.F45.sex <- update(umod.F45.Full, .~. -sex)
umod.F45.thic <- update(umod.F45.Full, .~. -thick.centerd)
umod.F45.alpw <- update(umod.F45.Full, .~. -alpha_pw)
umod.F45.bepw <- update(umod.F45.Full, .~. -beta_pw)
umod.F45.slop <- update(umod.F45.Full, .~. -a_slope)
umod.F45.itcp <- update(umod.F45.Full, .~. -a_intercept)
umod.F45.becf <- update(umod.F45.Full, .~. -beta_cf)
umod.F45.alcf <- update(umod.F45.Full, .~. -alpha_cf)
umod.F45.lene <- update(umod.F45.Full, .~. -leneve)
umod.F45.tuee <- update(umod.F45.Full, .~. -tueeve)
umod.F45.maxe <- update(umod.F45.Full, .~. -maxeve)
anova(umod.F45.neve, umod.F45.Full, test="Chisq")
anova(umod.F45.age,umod.F45.Full, test="Chisq")
anova(umod.F45.sex,umod.F45.Full, test="Chisq")
anova(umod.F45.Full,umod.F45.thic, test="Chisq")
anova(umod.F45.Full,umod.F45.alpw, test="Chisq")
anova(umod.F45.Full,umod.F45.bepw, test="Chisq")
anova(umod.F45.Full,umod.F45.slop, test="Chisq")
anova(umod.F45.Full,umod.F45.itcp, test="Chisq")
anova(umod.F45.Full,umod.F45.becf, test="Chisq")
anova(umod.F45.Full,umod.F45.alcf, test="Chisq")
anova(umod.F45.Full,umod.F45.lene, test="Chisq")
anova(umod.F45.Full,umod.F45.tuee, test="Chisq")
anova(umod.F45.Full,umod.F45.maxe, test="Chisq")

anova(umod.F45.Full, test="F")

umod.F6.Full <- glm(U.F6~nevent.u.m2.min+age+sex+thick.centerd+alpha_pw+beta_pw+a_slope+a_intercept+beta_cf+alpha_cf+
                       leneve+tueeve+maxeve, data=pd.data, family=gaussian)
umod.F6.neve <- update(umod.F6.Full, .~. -nevent.u.m2.min)
umod.F6.age <- update(umod.F6.Full, .~. -age)
umod.F6.sex <- update(umod.F6.Full, .~. -sex)
umod.F6.thic <- update(umod.F6.Full, .~. -thick.centerd)
umod.F6.alpw <- update(umod.F6.Full, .~. -alpha_pw)
umod.F6.bepw <- update(umod.F6.Full, .~. -beta_pw)
umod.F6.slop <- update(umod.F6.Full, .~. -a_slope)
umod.F6.itcp <- update(umod.F6.Full, .~. -a_intercept)
umod.F6.becf <- update(umod.F6.Full, .~. -beta_cf)
umod.F6.alcf <- update(umod.F6.Full, .~. -alpha_cf)
umod.F6.lene <- update(umod.F6.Full, .~. -leneve)
umod.F6.tuee <- update(umod.F6.Full, .~. -tueeve)
umod.F6.maxe <- update(umod.F6.Full, .~. -maxeve)
anova(umod.F6.neve, umod.F6.Full, test="Chisq")
anova(umod.F6.age,umod.F6.Full, test="Chisq")
anova(umod.F6.sex,umod.F6.Full, test="Chisq")
anova(umod.F6.Full,umod.F6.thic, test="Chisq")
anova(umod.F6.Full,umod.F6.alpw, test="Chisq")
anova(umod.F6.Full,umod.F6.bepw, test="Chisq")
anova(umod.F6.Full,umod.F6.slop, test="Chisq")
anova(umod.F6.Full,umod.F6.itcp, test="Chisq")
anova(umod.F6.Full,umod.F6.becf, test="Chisq")
anova(umod.F6.Full,umod.F6.alcf, test="Chisq")
anova(umod.F6.Full,umod.F6.lene, test="Chisq")
anova(umod.F6.Full,umod.F6.tuee, test="Chisq")
anova(umod.F6.Full,umod.F6.maxe, test="Chisq")

anova(umod.F6.Full, test="F")



umod.F7.Full <- glm(U.F7~nevent.u.m2.min+age+sex+thick.centerd+alpha_pw+beta_pw+a_slope+a_intercept+beta_cf+alpha_cf+
                      leneve+tueeve+maxeve, data=pd.data, family=gaussian)
umod.F7.neve <- update(umod.F7.Full, .~. -nevent.u.m2.min)
umod.F7.age <- update(umod.F7.Full, .~. -age)
umod.F7.sex <- update(umod.F7.Full, .~. -sex)
umod.F7.thic <- update(umod.F7.Full, .~. -thick.centerd)
umod.F7.alpw <- update(umod.F7.Full, .~. -alpha_pw)
umod.F7.bepw <- update(umod.F7.Full, .~. -beta_pw)
umod.F7.slop <- update(umod.F7.Full, .~. -a_slope)
umod.F7.itcp <- update(umod.F7.Full, .~. -a_intercept)
umod.F7.becf <- update(umod.F7.Full, .~. -beta_cf)
umod.F7.alcf <- update(umod.F7.Full, .~. -alpha_cf)
umod.F7.lene <- update(umod.F7.Full, .~. -leneve)
umod.F7.tuee <- update(umod.F7.Full, .~. -tueeve)
umod.F7.maxe <- update(umod.F7.Full, .~. -maxeve)
anova(umod.F7.neve, umod.F7.Full, test="Chisq")
anova(umod.F7.age,umod.F7.Full, test="Chisq")
anova(umod.F7.sex,umod.F7.Full, test="Chisq")
anova(umod.F7.Full,umod.F7.thic, test="Chisq")
anova(umod.F7.Full,umod.F7.alpw, test="Chisq")
anova(umod.F7.Full,umod.F7.bepw, test="Chisq")
anova(umod.F7.Full,umod.F7.slop, test="Chisq")
anova(umod.F7.Full,umod.F7.itcp, test="Chisq")
anova(umod.F7.Full,umod.F7.becf, test="Chisq")
anova(umod.F7.Full,umod.F7.alcf, test="Chisq")
anova(umod.F7.Full,umod.F7.lene, test="Chisq")
anova(umod.F7.Full,umod.F7.tuee, test="Chisq")
anova(umod.F7.Full,umod.F7.maxe, test="Chisq")

anova(umod.F7.Full, test="F")


umod.Full <- glm(UPDRS~nevent.u.m2.min+age+sex+thick.centerd+alpha_pw+beta_pw+a_slope+a_intercept+beta_cf+alpha_cf+
                      leneve+tueeve+maxeve, data=pd.data, family=gaussian)
umod.F7.neve <- update(umod.F7.Full, .~. -nevent.u.m2.min)
umod.F7.age <- update(umod.F7.Full, .~. -age)
umod.F7.sex <- update(umod.F7.Full, .~. -sex)
umod.F7.thic <- update(umod.F7.Full, .~. -thick.centerd)
umod.F7.alpw <- update(umod.F7.Full, .~. -alpha_pw)
umod.F7.bepw <- update(umod.F7.Full, .~. -beta_pw)
umod.F7.slop <- update(umod.F7.Full, .~. -a_slope)
umod.F7.itcp <- update(umod.F7.Full, .~. -a_intercept)
umod.F7.becf <- update(umod.F7.Full, .~. -beta_cf)
umod.F7.alcf <- update(umod.F7.Full, .~. -alpha_cf)
umod.F7.lene <- update(umod.F7.Full, .~. -leneve)
umod.F7.tuee <- update(umod.F7.Full, .~. -tueeve)
umod.F7.maxe <- update(umod.F7.Full, .~. -maxeve)
anova(umod.F7.neve, umod.F7.Full, test="Chisq")
anova(umod.F7.age,umod.F7.Full, test="Chisq")
anova(umod.F7.sex,umod.F7.Full, test="Chisq")
anova(umod.F7.Full,umod.F7.thic, test="Chisq")
anova(umod.F7.Full,umod.F7.alpw, test="Chisq")
anova(umod.F7.Full,umod.F7.bepw, test="Chisq")
anova(umod.F7.Full,umod.F7.slop, test="Chisq")
anova(umod.F7.Full,umod.F7.itcp, test="Chisq")
anova(umod.F7.Full,umod.F7.becf, test="Chisq")
anova(umod.F7.Full,umod.F7.alcf, test="Chisq")
anova(umod.F7.Full,umod.F7.lene, test="Chisq")
anova(umod.F7.Full,umod.F7.tuee, test="Chisq")
anova(umod.F7.Full,umod.F7.maxe, test="Chisq")

anova(umod.Full, test="F")
# nmod.u.F1 <- glm(nevent.u.m2~age+sex+U.F1, data=pd.data, family=poisson)
# nmod.u.0 <- update(nmod.u.F1, .~. -U.F1, subset=!(is.na(pd.data$U.F1)))
# anova(nmod.u.0, nmod.u.F1, test="Chisq")

###
nmod.u.F2 <- glm(U.F2~nevent.u.m2+age+sex, data=pd.data, family=gaussian)
nmod.u.F2.0 <- glm(U.F2~1+age+sex, data=pd.data, family=gaussian)
anova(nmod.u.F2.0, nmod.u.F2, test="Chisq")

nmod.u.F2 <- glm(nevent.u.m2~age+sex+U.F2, data=pd.data, family=poisson)
nmod.u.0 <- update(nmod.u.F2, .~. -U.F2, subset=!(is.na(pd.data$U.F2)))
anova(nmod.u.0, nmod.u.F2, test="Chisq")

###
nmod.u.F3 <- glm(nevent.u.m2~age+sex+U.F3, data=pd.data, family=poisson)
nmod.u.0 <- update(nmod.u.F3, .~. -U.F3, subset=!(is.na(pd.data$U.F3)))
anova(nmod.u.0, nmod.u.F3, test="Chisq")

nmod.u.F3 <- glm(U.F3~nevent.u.m2+age+sex, data=pd.data, family=gaussian)
nmod.u.F3.0 <- glm(U.F3~1+age+sex, data=pd.data, family=gaussian)
anova(nmod.u.F3.0, nmod.u.F3, test="Chisq")

###
nmod.u.F45 <- glm(U.F45~nevent.u.m2+age+sex, data=pd.data, family=gaussian)
nmod.u.F45.0 <- glm(U.F45~1+age+sex, data=pd.data, family=gaussian)
anova(nmod.u.F45.0, nmod.u.F45, test="Chisq")

nmod.u.F45 <- glm(nevent.u.m2~age+sex+U.F45, data=pd.data, family=poisson)
nmod.u.0 <- update(nmod.u.F45, .~. -U.F45, subset=!(is.na(pd.data$U.F45)))
anova(nmod.u.0, nmod.u.F45, test="Chisq")

###
nmod.u.F6 <- glm(U.F6~nevent.u.m2+age+sex, data=pd.data, family=gaussian)
nmod.u.F6.0 <- glm(U.F6~1+age+sex, data=pd.data, family=gaussian)
anova(nmod.u.F6.0, nmod.u.F6, test="Chisq")

nmod.u.F6 <- glm(nevent.u.m2~age+sex+U.F6, data=pd.data, family=poisson)
nmod.u.0 <- update(nmod.u.F6, .~. -U.F6, subset=!(is.na(pd.data$U.F6)))
anova(nmod.u.0, nmod.u.F6, test="Chisq")

###
nmod.u.F7 <- glm(U.F7~nevent.u.m2+age+sex, data=pd.data, family=gaussian)
nmod.u.F7.0 <- glm(U.F7~1+age+sex, data=pd.data, family=gaussian)
anova(nmod.u.F7.0, nmod.u.F7, test="Chisq")

nmod.u.F7 <- glm(nevent.u.m2~age+sex+U.F7, data=pd.data, family=poisson)
nmod.u.0 <- update(nmod.u.F7, .~. -U.F7, subset=!(is.na(pd.data$U.F7)))
anova(nmod.u.0, nmod.u.F7, test="Chisq")

############################################################################
library(brms)
bfnmod.u.F1 <- brm(U.F1~nevent.u.m2+age+sex, data=pd.data, family=poisson(), save_all_pars=TRUE)
bfnmod.u.F1.0 <- brm(U.F1~1+age+sex, data=pd.data, family=poisson, save_all_pars=TRUE)
bayes_factor(bfnmod.u.F1.0,nmod.u.F1)
hypothesis(nmod.u.F1, "nevent.u.m2<0")

nmod.u.F2 <- brm(U.F2~nevent.u.m2+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
nmod.u.F2.0 <- brm(U.F2~1+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
bayes_factor(nmod.u.F2.0,nmod.u.F2)

nmod.u.F3 <- brm(U.F3~nevent.u.m2+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
nmod.u.F3.0 <- brm(U.F3~1+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
bayes_factor(nmod.u.F3.0,nmod.u.F3)

nmod.u.F45 <- brm(U.F45~nevent.u.m2+age+sex, data=pd.data, family=poisson, save_all_pars=TRUE)
nmod.u.F45.0 <- brm(U.F45~1+age+sex, data=pd.data, family=poisson, save_all_pars=TRUE)
bayes_factor(nmod.u.F45.0,nmod.u.F45)
hypothesis(nmod.u.F45, "nevent.u.m2<0")


nmod.u.F6 <- brm(U.F6~nevent.u.m2+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
nmod.u.F6.0 <- brm(U.F6~1+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
bayes_factor(nmod.u.F45.0,nmod.u.F45)

nmod.u.F7 <- brm(U.F7~nevent.u.m2+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
nmod.u.F7.0 <- brm(U.F7~1+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
bayes_factor(nmod.u.F7.0,nmod.u.F7)

nmod.u.F45 <- brm(U.F45~nevent.u.m2+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
nmod.u.F45.0 <- brm(U.F45~1+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
bayes_factor(nmod.u.F45.0,nmod.u.F45)



hypothesis(nmod.u.F45, "nevent.u.m2<0")


nmod.u.F45 <- brm(U.F45~nevent.u.m2+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
nmod.u.F45.0 <- brm(U.F45~1+age+sex, data=pd.data, family=gaussian, save_all_pars=TRUE)
bayes_factor(nmod.u.F45.0,nmod.u.F45)



############################################################################
library(BayesFactor)
pd.data <- subset(pd.data, !is.na(U.F7))
pd.data <- subset(pd.data, !is.na(U.F6))
pd.data <- subset(pd.data, !is.na(U.F3))
pd.data <- subset(pd.data, !is.na(U.F1))

bf.0 <- lmBF(nevent.u.m2~age+sex, data=pd.data)
bf.u7 <- lmBF(nevent.u.m2~age+sex+U.F7, data=pd.data)
bf.u7/bf.0
bf.u6 <- lmBF(nevent.u.m2~age+sex+U.F6, data=pd.data)
bf.u6/bf.0
bf.u45 <- lmBF(U.F45~age+sex+nevent.u.m2, data=pd.data)
bf.u45.0 <- lmBF(U.F45~age+sex, data=pd.data)
bf.u45/bf.u45.0

bf.u45 <- lmBF(nevent.u.m2~U.F45+age+sex, data=pd.data)
bf.u45/bf.0

bf.u3 <- lmBF(nevent.u.m2~age+sex+U.F3, data=pd.data)
bf.u3/bf.0
bf.u2 <- lmBF(nevent.u.m2~age+sex+U.F2, data=pd.data)
bf.u2/bf.0
bf.u1 <- lmBF(nevent.u.m2~age+sex+U.F1, data=pd.data)
bf.u1/bf.0
bf.uT <- lmBF(nevent.u.m2~age+sex+UPDRS, data=pd.data)
bf.uT/bf.0

############################################################################
# summary(nmod.b.m2)
summary(nmod.u.m2)

pd.data$pred <- exp(predict(nmod.u.4, re.form=NA))
ggplot( aes(x=UPDRS, y=nevent.u.m2, color=sex), data=pd.data)+
  geom_point()+
  geom_smooth(mehtod='lm')+
  geom_line(aes(y = pred, linetype=sex), size = 1)

## SUBSCALES
nmod.F1 <- glm(nevent.u.m2 ~ U.F1 * age * sex, data=pd.data, family=poisson)
nmod.F2 <- glm(nevent.u.m2 ~ U.F2 * age * sex, data=pd.data, family=poisson)
nmod.F3 <- glm(nevent.u.m2 ~ U.F3 * age * sex, data=pd.data, family=poisson)
nmod.F4 <- glm(nevent.u.m2 ~ U.F4 * age * sex, data=pd.data, family=poisson)
nmod.F5 <- glm(nevent.u.m2 ~ U.F5 * age * sex, data=pd.data, family=poisson)
nmod.F6 <- glm(nevent.u.m2 ~ U.F6 * age * sex, data=pd.data, family=poisson)
nmod.F7 <- glm(nevent.u.m2 ~ U.F7 * age * sex, data=pd.data, family=poisson)
anova(nmod.F1, test="Chisq")
anova(nmod.F2, test="Chisq")
anova(nmod.F3, test="Chisq")
anova(nmod.F4, test="Chisq")
anova(nmod.F5, test="Chisq")
anova(nmod.F6, test="Chisq")
anova(nmod.F7, test="Chisq")


