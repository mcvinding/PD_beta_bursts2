# UPDRS stats
library(ggplot2)

# Load data
load(file='X://PD_longrest//groupanalysis//clindata.Rdata')
load('X://PD_longrest//groupanalysis//ndata_all.Rdata')


udata <- merge(ndata, clindata, by = c('subj','group','age','sex'))
pd.data <- subset(udata, group == 'patient')

# UPDRS ~ Neve
umod.b.m1 <- glm(UPDRS ~ nevent.b.m1 + age + sex, data=pd.data, family=poisson)
umod.b.m2 <- glm(UPDRS ~ nevent.b.m2 * age * sex, data=pd.data, family=poisson)
umod.b.pc <- glm(UPDRS ~ nevent.b.pc + age + sex, data=pd.data, family=poisson)
umod.u.m1 <- glm(UPDRS ~ nevent.u.m1 + age + sex, data=pd.data, family=poisson)
umod.u.m2 <- glm(UPDRS ~ nevent.u.m2 * age * sex, data=pd.data, family=poisson)
umod.u.pc <- glm(UPDRS ~ nevent.u.pc + age + sex, data=pd.data, family=poisson)
summary(umod.b.m1)
summary(umod.b.m2)
summary(umod.b.pc)
summary(umod.u.m1)
summary(umod.u.m2)
summary(umod.u.pc)


# Neve ~ UPDRS
nmod.b.m1 <- glm(nevent.b.m1 ~ UPDRS * age * sex, data=pd.data, family=poisson)
nmod.b.m2 <- glm(nevent.b.m2 ~ UPDRS * age * sex, data=pd.data, family=poisson)
nmod.b.pc <- glm(nevent.b.pc ~ UPDRS * age * sex, data=pd.data, family=poisson)
nmod.u.m1 <- glm(nevent.u.m1 ~ UPDRS * age * sex, data=pd.data, family=poisson)
nmod.u.m2 <- glm(nevent.u.m2 ~ UPDRS * age * sex, data=pd.data, family=poisson)
nmod.u.pc <- glm(nevent.u.pc ~ UPDRS * age * sex, data=pd.data, family=poisson)
summary(nmod.b.m1)
summary(nmod.b.m2)
summary(nmod.b.pc)
summary(nmod.u.m1)
summary(nmod.u.m2)
summary(nmod.u.pc)


pd.data$pred <- exp(predict(umod.b.m2, re.form=NA))

ggplot( aes(y=UPDRS, x=nevent.b.m2, color=sex), data=pd.data)+
  geom_point()+
  geom_smooth(method='lm')+
  geom_line(aes(y = pred, linetype=sex), size = 1)
anova(umod.b.m2, test="Chisq")



pd.data$pred <- exp(predict(umod.u.m2, re.form=NA))

ggplot( aes(y=UPDRS, x=nevent.u.m2, color=sex), data=pd.data)+
  geom_point()+
  geom_smooth(method='lm')+
  geom_line(aes(y = pred, linetype=sex), size = 1)
anova(umod.u.m2, test="Chisq")


## SUBSCALES

