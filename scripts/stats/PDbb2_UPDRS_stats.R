# UPDRS stats
library(arm)

# Load data
wrkdir <- "X://PD_longrest//groupanalysis"
setwd(wrkdir)
load(file='X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
load(file='X://PD_longrest//groupanalysis//bbdata2.Rdata')

# Prepare data
bbsum <- aggregate(cbind(leneve, tueeve, maxeve)~subj, data=bbdata, FUN=median)
alldata <- merge(alldata, bbsum, by="subj")

pd.data <- subset(alldata, group == 'patient')
pd.data <- subset(pd.data, !is.na(U.F7))
pd.data <- subset(pd.data, !is.na(U.F6))
pd.data <- subset(pd.data, !is.na(U.F3))
pd.data <- subset(pd.data, !is.na(U.F1))
pd.data$U.F45 <- pd.data$U.F4+pd.data$U.F5

pd.data$log.leneve <- log(pd.data$leneve)
pd.data$log.tueeve <- log(pd.data$tueeve)
pd.data$log.maxeve <- log(pd.data$maxeve)

pd.data$beta_pwz <- zscore(pd.data$beta_pw)
pd.data$alpha_pwz <- zscore(pd.data$alpha_pw)

pd.data$beta_pwc <- (pd.data$beta_pw-min(pd.data$beta_pw))/(max(pd.data$beta_pw)-min(pd.data$beta_pw))*10
  

# Transform data to give intepretable
pd.data$nevent.u.m2 


######################################################################################
# MDS-UPDRS-III Factor 1
F1mod.x <- glm(U.F1 ~ nevent.u.m2.min+log.leneve+log.tueeve+log.maxeve+
                 age.centerd+sex+thick.centerd+
                 a_intercept+a_slope+alpha_pw+alpha_cf+beta_pwc+beta_cf, 
               data=pd.data, family=poisson)

summary(F1mod.x)
qqnorm(resid(F1mod.x))

qqline(resid(F1mod.x))

F1mod.neve <- update(F1mod.x, ~. -nevent.u.m2.min)
anova(F1mod.x, F1mod.neve, test="Chisq")

F1mod.leneve <- update(F1mod.x, ~. -log.leneve)
anova(F1mod.x, F1mod.leneve, test="Chisq")

F1mod.tueeve <- update(F1mod.x, ~. -log.tueeve)
anova(F1mod.x, F1mod.tueeve, test="Chisq")

F1mod.maxeve <- update(F1mod.x, ~. -log.maxeve)
anova(F1mod.x, F1mod.maxeve, test="Chisq")

F1mod.age <- update(F1mod.x, ~. -age.centerd)
anova(F1mod.x, F1mod.age, test="Chisq")

F1mod.sex <- update(F1mod.x, ~. -sex)
anova(F1mod.x, F1mod.sex, test="Chisq")

F1mod.thick <- update(F1mod.x, ~. -thick.centerd)
anova(F1mod.x, F1mod.thick, test="Chisq")

F1mod.intcpt <- update(F1mod.x, ~. -a_intercept)
anova(F1mod.x, F1mod.intcpt, test="Chisq")

F1mod.slope <- update(F1mod.x, ~. -a_slope)
anova(F1mod.x, F1mod.slope, test="Chisq")

F1mod.alpapw <- update(F1mod.x, ~. -alpha_pw)
anova(F1mod.x, F1mod.alpapw, test="Chisq")

F1mod.alpacf <- update(F1mod.x, ~. -alpha_cf)
anova(F1mod.x, F1mod.alpacf, test="Chisq")

F1mod.betapw <- update(F1mod.x, ~. -beta_pwc)
anova(F1mod.x, F1mod.betapw, test="Chisq")

F1mod.betacf <- update(F1mod.x, ~. -beta_cf)
anova(F1mod.x, F1mod.betacf, test="Chisq")

# SUMMARY
mod.sim <- sim(F1mod.x, n.sims=1000)
x1 <- summary(F1mod.x)$coefficients
x2 <- t(apply(coef(mod.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1[,1], x2)

(exp(sums)-1)*100

######################################################################################
# MDS-UPDRS-III Factor 2
F2mod.x <- glm(U.F2 ~ nevent.u.m2.min+log.leneve+log.tueeve+log.maxeve+
                 age.centerd+sex+thick.centerd+
                 a_intercept+a_slope+alpha_pw+alpha_cf+beta_pwz+beta_cf, 
               data=pd.data, family=poisson)

F2mod.neve <- update(F2mod.x, ~. -nevent.u.m2.min)
anova(F2mod.x, F2mod.neve, test="Chisq")

F2mod.leneve <- update(F2mod.x, ~. -log.leneve)
anova(F2mod.x, F2mod.leneve, test="Chisq")

F2mod.tueeve <- update(F2mod.x, ~. -log.tueeve)
anova(F2mod.x, F2mod.tueeve, test="Chisq")

F2mod.maxeve <- update(F2mod.x, ~. -log.maxeve)
anova(F2mod.x, F2mod.maxeve, test="Chisq")

F2mod.age <- update(F2mod.x, ~. -age.centerd)
anova(F2mod.x, F2mod.age, test="Chisq")

F2mod.sex <- update(F2mod.x, ~. -sex)
anova(F2mod.x, F2mod.sex, test="Chisq")

F2mod.thick <- update(F2mod.x, ~. -thick.centerd)
anova(F2mod.x, F2mod.thick, test="Chisq")

F2mod.intcpt <- update(F2mod.x, ~. -a_intercept)
anova(F2mod.x, F2mod.intcpt, test="Chisq")

F2mod.slope <- update(F2mod.x, ~. -a_slope)
anova(F2mod.x, F2mod.slope, test="Chisq")

F2mod.alpapw <- update(F2mod.x, ~. -alpha_pw)
anova(F2mod.x, F2mod.alpapw, test="Chisq")

F2mod.alpacf <- update(F2mod.x, ~. -alpha_cf)
anova(F2mod.x, F2mod.alpacf, test="Chisq")

F2mod.betapw <- update(F2mod.x, ~. -beta_pwz)
anova(F2mod.x, F2mod.betapw, test="Chisq")

F2mod.betacf <- update(F2mod.x, ~. -beta_cf)
anova(F2mod.x, F2mod.betacf, test="Chisq")

# SUMMARY
mod2.sim <- sim(F2mod.x, n.sims=1000)
x1 <- summary(F2mod.x)$coefficients
x2 <- t(apply(coef(mod2.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1[,1], x2)

(exp(sums)-1)*100


######################################################################################
# MDS-UPDRS-III Factor 3
F3mod.x <- glm(U.F3 ~ nevent.u.m2.min+log.leneve+log.tueeve+log.maxeve+
                 age.centerd+sex+thick.centerd+
                 a_intercept+a_slope+alpha_pw+alpha_cf+beta_pw+beta_cf, 
               data=pd.data, family=poisson)

F3mod.neve <- update(F3mod.x, .~. -nevent.u.m2.min)
anova(F3mod.x, F3mod.neve, test="Chisq")

F3mod.leneve <- update(F3mod.x, ~. -log.leneve)
anova(F3mod.x, F3mod.leneve, test="Chisq")

F3mod.tueeve <- update(F3mod.x, ~. -log.tueeve)
anova(F3mod.x, F3mod.tueeve, test="Chisq")

F3mod.maxeve <- update(F3mod.x, ~. -log.maxeve)
anova(F3mod.x, F3mod.maxeve, test="Chisq")

F3mod.age <- update(F3mod.x, ~. -age.centerd)
anova(F3mod.x, F3mod.age, test="Chisq")

F3mod.sex <- update(F3mod.x, ~. -sex)
anova(F3mod.x, F3mod.sex, test="Chisq")

F3mod.thick <- update(F3mod.x, ~. -thick.centerd)
anova(F3mod.x, F3mod.thick, test="Chisq")

F3mod.intcpt <- update(F3mod.x, ~. -a_intercept)
anova(F3mod.x, F3mod.intcpt, test="Chisq")

F3mod.slope <- update(F3mod.x, ~. -a_slope)
anova(F3mod.x, F3mod.slope, test="Chisq")

F3mod.alpapw <- update(F3mod.x, ~. -alpha_pw)
anova(F3mod.x, F3mod.alpapw, test="Chisq")

F3mod.alpacf <- update(F3mod.x, ~. -alpha_cf)
anova(F3mod.x, F3mod.alpacf, test="Chisq")

F3mod.betapw <- update(F3mod.x, ~. -beta_pw)
anova(F3mod.x, F3mod.betapw, test="Chisq")

F3mod.betacf <- update(F3mod.x, ~. -beta_cf)
anova(F3mod.x, F3mod.betacf, test="Chisq")

# SUMMARY
mod3.sim <- sim(F3mod.x, n.sims=1000)
x1 <- summary(F3mod.x)$coefficients
x2 <- t(apply(coef(mod3.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1[,1], x2)

(exp(sums)-1)*100

######################################################################################
# MDS-UPDRS-III Factor 4
# F4mod.x <- glm(U.F4 ~ nevent.u.m2.min+log.leneve+log.tueeve+log.maxeve+
#                  age.centerd+sex+thick.centerd+
#                  a_intercept+a_slope+alpha_pw+alpha_cf+beta_pw+beta_cf, 
#                data=pd.data, family=poisson)
# 
# F4mod.neve <- update(F4mod.x, .~. -nevent.u.m2.min)
# anova(F4mod.x, F4mod.neve, test="Chisq")
# 
# F4mod.leneve <- update(F4mod.x, ~. -log.leneve)
# anova(F4mod.x, F4mod.leneve, test="Chisq")
# 
# F4mod.tueeve <- update(F4mod.x, ~. -log.tueeve)
# anova(F4mod.x, F4mod.tueeve, test="Chisq")
# 
# F4mod.maxeve <- update(F4mod.x, ~. -log.maxeve)
# anova(F4mod.x, F4mod.maxeve, test="Chisq")
# 
# F4mod.age <- update(F4mod.x, ~. -age.centerd)
# anova(F4mod.x, F4mod.age, test="Chisq")
# 
# F4mod.sex <- update(F4mod.x, ~. -sex)
# anova(F4mod.x, F4mod.sex, test="Chisq")
# 
# F4mod.thick <- update(F4mod.x, ~. -thick.centerd)
# anova(F4mod.x, F4mod.thick, test="Chisq")
# 
# F4mod.intcpt <- update(F4mod.x, ~. -a_intercept)
# anova(F4mod.x, F4mod.intcpt, test="Chisq")
# 
# F4mod.slope <- update(F4mod.x, ~. -a_slope)
# anova(F4mod.x, F4mod.slope, test="Chisq")
# 
# F4mod.alpapw <- update(F4mod.x, ~. -alpha_pw)
# anova(F4mod.x, F4mod.alpapw, test="Chisq")
# 
# F4mod.alpacf <- update(F4mod.x, ~. -alpha_cf)
# anova(F4mod.x, F4mod.alpacf, test="Chisq")
# 
# F4mod.betapw <- update(F4mod.x, ~. -beta_pw)
# anova(F4mod.x, F4mod.betapw, test="Chisq")
# 
# F4mod.betacf <- update(F4mod.x, ~. -beta_cf)
# anova(F4mod.x, F4mod.betacf, test="Chisq")
# 
# ######################################################################################
# # MDS-UPDRS-III Factor 5
# F5mod.x <- glm(U.F5 ~ nevent.u.m2.min+log.leneve+log.tueeve+log.maxeve+
#                  age.centerd+sex+thick.centerd+
#                  a_intercept+a_slope+alpha_pw+alpha_cf+beta_pw+beta_cf, 
#                data=pd.data, family=poisson)
# 
# F5mod.neve <- update(F5mod.x, .~. -nevent.u.m2.min)
# anova(F5mod.x, F5mod.neve, test="Chisq")
# 
# F5mod.leneve <- update(F5mod.x, ~. -log.leneve)
# anova(F5mod.x, F5mod.leneve, test="Chisq")
# 
# F5mod.tueeve <- update(F5mod.x, ~. -log.tueeve)
# anova(F5mod.x, F5mod.tueeve, test="Chisq")
# 
# F5mod.maxeve <- update(F5mod.x, ~. -log.maxeve)
# anova(F5mod.x, F5mod.maxeve, test="Chisq")
# 
# F5mod.age <- update(F5mod.x, ~. -age.centerd)
# anova(F5mod.x, F5mod.age, test="Chisq")
# 
# F5mod.sex <- update(F5mod.x, ~. -sex)
# anova(F5mod.x, F5mod.sex, test="Chisq")
# 
# F5mod.thick <- update(F5mod.x, ~. -thick.centerd)
# anova(F5mod.x, F5mod.thick, test="Chisq")
# 
# F5mod.intcpt <- update(F5mod.x, ~. -a_intercept)
# anova(F5mod.x, F5mod.intcpt, test="Chisq")
# 
# F5mod.slope <- update(F5mod.x, ~. -a_slope)
# anova(F5mod.x, F5mod.slope, test="Chisq")
# 
# F5mod.alpapw <- update(F5mod.x, ~. -alpha_pw)
# anova(F5mod.x, F5mod.alpapw, test="Chisq")
# 
# F5mod.alpacf <- update(F5mod.x, ~. -alpha_cf)
# anova(F5mod.x, F5mod.alpacf, test="Chisq")
# 
# F5mod.betapw <- update(F5mod.x, ~. -beta_pw)
# anova(F5mod.x, F5mod.betapw, test="Chisq")
# 
# F5mod.betacf <- update(F5mod.x, ~. -beta_cf)
# anova(F5mod.x, F5mod.betacf, test="Chisq")

######################################################################################
# MDS-UPDRS-III Factor 4+5 combined
F45mod.x <- glm(U.F45 ~ nevent.u.m2.min+log.leneve+log.tueeve+log.maxeve+
                 age.centerd+sex+thick.centerd+
                 a_intercept+a_slope+alpha_pw+alpha_cf+beta_pw+beta_cf, 
                data=pd.data, family=poisson)

F45mod.neve <- update(F45mod.x, .~. -nevent.u.m2.min)
anova(F45mod.x, F45mod.neve, test="Chisq")

F45mod.leneve <- update(F45mod.x, ~. -log.leneve)
anova(F45mod.x, F45mod.leneve, test="Chisq")

F45mod.tueeve <- update(F45mod.x, ~. -log.tueeve)
anova(F45mod.x, F45mod.tueeve, test="Chisq")

F45mod.maxeve <- update(F45mod.x, ~. -log.maxeve)
anova(F45mod.x, F45mod.maxeve, test="Chisq")

F45mod.age <- update(F45mod.x, ~. -age.centerd)
anova(F45mod.x, F45mod.age, test="Chisq")

F45mod.sex <- update(F45mod.x, ~. -sex)
anova(F45mod.x, F45mod.sex, test="Chisq")

F45mod.thick <- update(F45mod.x, ~. -thick.centerd)
anova(F45mod.x, F45mod.thick, test="Chisq")

F45mod.intcpt <- update(F45mod.x, ~. -a_intercept)
anova(F45mod.x, F45mod.intcpt, test="Chisq")

F45mod.slope <- update(F45mod.x, ~. -a_slope)
anova(F45mod.x, F45mod.slope, test="Chisq")

F45mod.alpapw <- update(F45mod.x, ~. -alpha_pw)
anova(F45mod.x, F45mod.alpapw, test="Chisq")

F45mod.alpacf <- update(F45mod.x, ~. -alpha_cf)
anova(F45mod.x, F45mod.alpacf, test="Chisq")

F45mod.betapw <- update(F45mod.x, ~. -beta_pw)
anova(F45mod.x, F45mod.betapw, test="Chisq")

F45mod.betacf <- update(F45mod.x, ~. -beta_cf)
anova(F45mod.x, F45mod.betacf, test="Chisq")

# SUMMARY
mod45.sim <- sim(F45mod.x, n.sims=1000)
x1 <- summary(F45mod.x)$coefficients
x2 <- t(apply(coef(mod45.sim), 2, quantile, c(0.025, 0.975)))
sums <- cbind(x1[,1], x2)

(exp(sums)-1)*100

######################################################################################
# MDS-UPDRS-III Factor 6
F6mod.x <- glm(U.F6 ~ nevent.u.m2.min+log.leneve+log.tueeve+log.maxeve+
                 age.centerd+sex+thick.centerd+
                 a_intercept+a_slope+alpha_pw+alpha_cf+beta_pw+beta_cf, 
               data=pd.data, family=poisson)

F6mod.neve <- update(F6mod.x, .~. -nevent.u.m2.min)
anova(F6mod.x, F6mod.neve, test="Chisq")

F6mod.leneve <- update(F6mod.x, ~. -log.leneve)
anova(F6mod.x, F6mod.leneve, test="Chisq")

F6mod.tueeve <- update(F6mod.x, ~. -log.tueeve)
anova(F6mod.x, F6mod.tueeve, test="Chisq")

F6mod.maxeve <- update(F6mod.x, ~. -log.maxeve)
anova(F6mod.x, F6mod.maxeve, test="Chisq")

F6mod.age <- update(F6mod.x, ~. -age.centerd)
anova(F6mod.x, F6mod.age, test="Chisq")

F6mod.sex <- update(F6mod.x, ~. -sex)
anova(F6mod.x, F6mod.sex, test="Chisq")

F6mod.thick <- update(F6mod.x, ~. -thick.centerd)
anova(F6mod.x, F6mod.thick, test="Chisq")

F6mod.intcpt <- update(F6mod.x, ~. -a_intercept)
anova(F6mod.x, F6mod.intcpt, test="Chisq")

F6mod.slope <- update(F6mod.x, ~. -a_slope)
anova(F6mod.x, F6mod.slope, test="Chisq")

F6mod.alpapw <- update(F6mod.x, ~. -alpha_pw)
anova(F6mod.x, F6mod.alpapw, test="Chisq")

F6mod.alpacf <- update(F6mod.x, ~. -alpha_cf)
anova(F6mod.x, F6mod.alpacf, test="Chisq")

F6mod.betapw <- update(F6mod.x, ~. -beta_pw)
anova(F6mod.x, F6mod.betapw, test="Chisq")

F6mod.betacf <- update(F6mod.x, ~. -beta_cf)
anova(F6mod.x, F6mod.betacf, test="Chisq")

######################################################################################
# MDS-UPDRS-III Factor 7
F7mod.x <- glm(U.F7 ~ nevent.u.m2.min+log.leneve+log.tueeve+log.maxeve+
                 age.centerd+sex+thick.centerd+
                 a_intercept+a_slope+alpha_pw+alpha_cf+beta_pw+beta_cf, 
               data=pd.data, family=poisson)

F7mod.neve <- update(F7mod.x, .~. -nevent.u.m2.min)
anova(F7mod.x, F7mod.neve, test="Chisq")

F7mod.leneve <- update(F7mod.x, ~. -log.leneve)
anova(F7mod.x, F7mod.leneve, test="Chisq")

F7mod.tueeve <- update(F7mod.x, ~. -log.tueeve)
anova(F7mod.x, F7mod.tueeve, test="Chisq")

F7mod.maxeve <- update(F7mod.x, ~. -log.maxeve)
anova(F7mod.x, F7mod.maxeve, test="Chisq")

F7mod.age <- update(F7mod.x, ~. -age.centerd)
anova(F7mod.x, F7mod.age, test="Chisq")

F7mod.sex <- update(F7mod.x, ~. -sex)
anova(F7mod.x, F7mod.sex, test="Chisq")

F7mod.thick <- update(F7mod.x, ~. -thick.centerd)
anova(F7mod.x, F7mod.thick, test="Chisq")

F7mod.intcpt <- update(F7mod.x, ~. -a_intercept)
anova(F7mod.x, F7mod.intcpt, test="Chisq")

F7mod.slope <- update(F7mod.x, ~. -a_slope)
anova(F7mod.x, F7mod.slope, test="Chisq")

F7mod.alpapw <- update(F7mod.x, ~. -alpha_pw)
anova(F7mod.x, F7mod.alpapw, test="Chisq")

F7mod.alpacf <- update(F7mod.x, ~. -alpha_cf)
anova(F7mod.x, F7mod.alpacf, test="Chisq")

F7mod.betapw <- update(F7mod.x, ~. -beta_pw)
anova(F7mod.x, F7mod.betapw, test="Chisq")

F7mod.betacf <- update(F7mod.x, ~. -beta_cf)
anova(F7mod.x, F7mod.betacf, test="Chisq")



J <- 15
n <- J*(J+1)/2
group <- rep (1:J, 1:J)
mu.a <- 5
sigma.a <- 2
a <- rnorm (J, mu.a, sigma.a)
b <- -3
x <- rnorm (n, 2, 1)
sigma.y <- 6
y <- rnorm (n, a[group] + b*x, sigma.y)
u <- runif (J, 0, 3)
y123.dat <- cbind (y, x, group)

# Linear regression
x1 <- y123.dat[,2]
y1 <- y123.dat[,1]
M1 <- lm (y1 ~ x1)
display(M1)

M1.sim <- sim(F45mod.x)
coef.M1.sim <- coef(M1.sim)
sigma.M1.sim <- sigma.hat(M1.sim)

## to get the uncertainty for the simulated estimates
exp(coef(F45mod.x))-1
exp(apply(coef(M1.sim), 2, quantile, c(0.025, 0.975)))-1
quantile(sigma.hat(M1.sim))




#END