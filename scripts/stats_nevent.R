# Analysis of N events
library(lme4)
library(arm)
library(ggplot2)
# Bayes?

## Load data
load('X://PD_longrest//groupanalysis//ndata.Rdata')

## Quick analysis

nmod <- lmer(nevent.min ~  I(age-mean(age)) * sex * group + I(thick-mean(thick))*I(age-mean(age))*group + (1|hemi) + (1|subj), data=ndata)
nmod.x3 <- lmer(nevent.min ~ I(age-mean(age)) * sex * group + (1|hemi) + (1|subj), data=ndata)

nplot <- ggplot(aes(x=thick, y=nevent, color=group, shape=sex), data=ndata)+
  geom_point()+
  geom_smooth(method=lm)
nplot

t.test(ndata$nevent.min[ndata$group=="Patient"], ndata$nevent.min[ndata$group=="Control"])
