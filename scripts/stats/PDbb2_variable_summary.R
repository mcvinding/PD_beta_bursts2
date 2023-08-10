# Main outcome variable summaries

# Load data
load('alldata_subj2.Rdata')
load('bbdata2.Rdata')

bbsum <- aggregate(cbind(leneve.ms, tueeve.ms, maxeve)~subj, data=bbdata, FUN=median)
alldata <- merge(alldata, bbsum, by="subj")

######################################################################################
# Summaries by sex*group

sum.mean <- aggregate(cbind(nevent.u.m2.min, leneve.ms, tueeve.ms, maxeve, a_intercept, a_slope, beta_pw, beta_cf, alpha_pw, alpha_cf)~group*sex,
                       data=alldata, FUN=mean)
sum.sd <- aggregate(cbind(nevent.u.m2.min, leneve.ms, tueeve.ms, maxeve, a_intercept, a_slope, beta_pw, beta_cf, alpha_pw, alpha_cf)~group*sex,
                     data=alldata, FUN=sd)
sum.medi <- aggregate(cbind(nevent.u.m2.min, leneve.ms, tueeve.ms, maxeve, a_intercept, a_slope, beta_pw, beta_cf, alpha_pw, alpha_cf)~group*sex,
                       data=alldata, FUN=median)
sum.hdi <- aggregate(cbind(nevent.u.m2.min, leneve.ms, tueeve.ms, maxeve, a_intercept, a_slope, beta_pw, beta_cf, alpha_pw, alpha_cf)~group*sex,
                      data=alldata, FUN=quantile, probs=c(0.025, 0.975))

######################################################################################
# Compare FOOOF GOF (R^2)
t.test(rsqrd~group, data=alldata, alternative="two.sided")
summary(alldata$rsqrd)