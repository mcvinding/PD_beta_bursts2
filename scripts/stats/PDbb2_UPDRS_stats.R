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
# umod.b.m2 <- glm(UPDRS ~ nevent.b.m2 * age * sex, data=pd.data, family=gaussian)
umod.u.m2 <- glm(UPDRS ~ nevent.u.m2 * age * sex, data=pd.data, family=gaussian)
# summary(umod.b.m2)
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

# nmod.u.F1 <- glm(U.F1~nevent.u.m2+age+sex, data=pd.data, family=gaussian)
# nmod.u.F1.0 <- glm(U.F1~1+age+sex, data=pd.data, family=gaussian)

nmod.u.F1 <- glm(nevent.u.m2~age+sex+U.F1, data=pd.data, family=poisson)
nmod.u.0 <- update(nmod.u.F1, .~. -U.F1, subset=!(is.na(pd.data$U.F1)))
anova(nmod.u.0, nmod.u.F1, test="Chisq")

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


