# UPDRS stats
library(ggplot2)

# Load data
# load(file='X://PD_longrest//groupanalysis//clindata.Rdata')
load(file='C://Users//Mikkel//Documents//PDbb2//groupanalysis//alldata_subj.Rdata')

pd.data <- subset(alldata, group == 'patient')

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
# Neve ~ UPDRS
# nmod.b.m2 <- glm(nevent.b.m2 ~ UPDRS * age * sex, data=pd.data, family=poisson)
nmod.u.1 <- glm(nevent.u.m2 ~ UPDRS * age * sex, data=pd.data, family=poisson)
nmod.u.2 <- update(nmod.u.1, ~. -UPDRS:age:sex)
nmod.u.3 <- update(nmod.u.2, ~. -age:sex)
nmod.u.4 <- update(nmod.u.3, ~. -UPDRS:sex)
nmod.u.5 <- update(nmod.u.4, ~. -UPDRS:age)
nmod.u.6 <- update(nmod.u.5, ~. -sex)
nmod.u.7 <- update(nmod.u.6, ~. -age)
nmod.u.8 <- update(nmod.u.7, ~. -UPDRS)

anova(nmod.u.8,nmod.u.7,nmod.u.6,nmod.u.5,nmod.u.4,nmod.u.3,nmod.u.2,nmod.u.1, 
      test="Chisq")
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


