# Plot results of regression models. Plot each model in a single panel figure (instad of all in one figure)
#
# Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., 
#  Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor 
#  rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease [Preprint]. 
#  medRxiv.org https://doi.org/10.1101/2021.06.27.21259592
#

library(sjPlot)
library(sjmisc)
library(lme4)
library(ggplot2)
library(plyr)
setwd('~/PD_longrest/output')
load(file='/home/mikkel/PD_longrest/groupanalysis/bbdata2.Rdata')

################################################################################
# Variable names
a.labels <- c("Group", "Age", "Sex", "Cortical thickness", "Age^2",
              "Group:Age", "Group:Sex", "Group:Cortical thickness", "Age:Sex", "Age:Cortical thickness", "Sex:Cortical thickness")
              # "Group:Age:Sex", "Group:Age:Cortical thickness", "Group:Sex:Cortical thickness", "Age:Sex:Cortical thickness")

a.labels <- rev(a.labels)

# de-center age
age.convert <- function(x) {
  lab = round(x + mean(bbdata$age))
  return(lab)
}

# thick.convert <- function(x) {
#   lab = round(x + mean(bbdata$thick))
#   return(lab)
# }

dpi <- 600
p1.h <- 50
p1.w <- 40
p2.h <- 40
p2.w <- 50

################################################################################
# 1/f intercept
################################################################################
load(file='mod_finter_2.RData')
finter.mod$data$sex <- revalue(finter.mod$data$sex, c("F"="Female", "M"="Male"))
finter.mod$data$group <- revalue(finter.mod$data$group, c("patient"="PD", "control"="HC"))

# Coefs plot
plt.finter <- plot_model(finter.mod, 
                         type = "est",
                         bpe = "mean",
                         bpe.color = "red",
                         prob.inner = .50,
                         prob.outer = .95,
                         show.p = TRUE,
                         p.style = "asterisk",
                         colors="bw",
                         vline.color = "gray",
                         axis.labels = a.labels,
                         show.values = TRUE,
                         title = "1/f offset"
) 
plt.finter <- plt.finter + scale_y_continuous(limits=c(-0.5, 0.5)) + theme_bw()
plt.finter
ggsave("modplt_finter.png", plt.finter, dpi=dpi,
       width=p1.w, height=p1.h, units="mm", scale=3)

# Age plot
p1 <- plot_model(finter.mod, 
           type = "pred",
           bpe = "median",
           terms = c("age.centerd","sex","group"),
           show.data = TRUE,
           ci.lvl = .95,
           title = " ",
           legend.title = ""
) 
ageplt.finter <- p1 + theme_bw() +
  ylab('1/f intercept (a.u.)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)+
  ylim(-2.25, -0)+
  scale_colour_discrete(labels=c("Female","Male"))
ageplt.finter

ggsave("ageplt_finter.png", ageplt.finter, dpi=dpi,
       width=p2.w, height=p2.h, units="mm", scale=3)

# # Thickz plot
# p2 <- plot_model(finter.mod, 
#                  type = "pred",
#                  bpe = "mean",
#                  terms = c("thickz","group","sex"),
#                  show.data = TRUE,
#                  ci.lvl = .95,
#                  color = c("red", "blue"),
#                  label= c("HC","PD"),
#                  title = "Cortical thickness",
#                  legend.title = "Group"
# ) 
# tckplt.finter <- p2 + theme_bw() +
#   ylab('1/f intercept (a.u.)') +
#   xlab("Cortical thickness (z-score)")
#   # scale_color_manual(labels=c("HC","PD"), values = c("red","blue", "red","blue"))
#   # scale_fill_manual(labels=c("HC","PD"), values = c("red", "blue")) +
# tckplt.finter
# 
# ggsave("tckplt_finter.png", tckplt.finter, dpi=900,
#        width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "*finter"))

################################################################################
# 1/f exponent
################################################################################
load(file='mod_fslope_2.RData')
fslope.mod$data$sex <- revalue(fslope.mod$data$sex, c("F"="Female", "M"="Male"))
fslope.mod$data$group <- revalue(fslope.mod$data$group, c("patient"="PD", "control"="HC"))

# Coef plot
plt.fslope <- plot_model(fslope.mod, 
                         type = "est",
                         bpe = "mean",
                         bpe.color = "red",
                         prob.inner = .50,
                         prob.outer = .95,
                         show.p = TRUE,
                         p.style = "asterisk",
                         colors="bw",
                         vline.color = "gray",
                         axis.labels = a.labels,
                         show.values = TRUE,
                         title = "1/f Exponent"
) 
plt.fslope <- plt.fslope + scale_y_continuous(limits=c(-0.2, 0.2)) + theme_bw()
plt.fslope
ggsave("modplt_fslope.png", plt.fslope, dpi=dpi,
       width=p1.w, height=p1.h, units="mm", scale=3)

# Age plot
p1 <- plot_model(fslope.mod, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","sex","group"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 title = "",
                 legend.title = ""
) 
ageplt.fslope <- p1 + theme_bw() +
  ylab('1/f exponent (a.u.)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)+
  ylim(0.3, 1.2)+
  scale_colour_discrete(labels=c("Female","Male"))
ageplt.fslope
ggsave("ageplt_fslope.png", ageplt.fslope, dpi=dpi,
       width=p2.w, height=p2.h, units="mm", scale=3)

# # Thickz plot
# p2 <- plot_model(fslope.mod, 
#                  type = "pred",
#                  bpe = "mean",
#                  terms = c("thickz","group","sex"),
#                  show.data = TRUE,
#                  ci.lvl = .95,
#                  color = c("red", "blue"),
#                  label= c("HC","PD"),
#                  title = "Cortical thickness",
#                  legend.title = "Group"
# ) 
# tckplt.fslope <- p2 + theme_bw() +
#   ylab('1/f intercept (a.u.)') +
#   xlab("Cortical thickness (z-score)")
# tckplt.fslope
# ggsave("tckplt_fslope.png", tckplt.fslope, dpi=900,
#        width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "*fslope"))

################################################################################
# Alpha power
################################################################################
load(file='mod_alphapw_2.RData')
alpha_pw.mod$data$sex <- revalue(alpha_pw.mod$data$sex, c("F"="Female", "M"="Male"))
alpha_pw.mod$data$group <- revalue(alpha_pw.mod$data$group, c("patient"="PD", "control"="HC"))

# Coef plot
plt.alp_pw <- plot_model(alpha_pw.mod, 
                         type = "est",
                         bpe = "mean",
                         bpe.color = "red",
                         prob.inner = .50,
                         prob.outer = .95,
                         show.p = TRUE,
                         p.style = "asterisk",
                         colors="bw",
                         vline.color = "gray",
                         axis.labels = a.labels,
                         show.values = TRUE,
                         title = "Alpha power"
) 
plt.alp_pw <- plt.alp_pw + scale_y_continuous(limits=c(-0.2, 0.2)) + theme_bw()
plt.alp_pw
ggsave("modplt_alpPow.png", plt.alp_pw, dpi=dpi,
       width=p1.w, height=p1.h, units="mm", scale=3)

# Age plot
p1 <- plot_model(alpha_pw.mod, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","sex","group"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 title = "",
                 legend.title = ""
) 
ageplt.alpPw <- p1 + theme_bw() +
  ylab('Power (a.u.)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)+
  ylim(0.1, 0.8)+
  scale_colour_discrete(labels=c("Female","Male"))
ageplt.alpPw
ggsave("ageplt_alpPow.png", ageplt.alpPw, dpi=dpi,
       width=p2.w, height=p2.h, units="mm", scale=3)

# # Thickz plot
# p2 <- plot_model(alpha_pw.mod,
#                  type = "pred",
#                  bpe = "mean",
#                  terms = c("thickz","group","sex"),
#                  show.data = TRUE,
#                  ci.lvl = .95,
#                  color = c("red", "blue"),
#                  label= c("HC","PD"),
#                  title = "Cortical thickness",
#                  legend.title = "Group"
# ) 
# tckplt.alpPw <- p2 + theme_bw() +
#   ylab('Power (a.u.)') +
#   xlab("Cortical thickness (z-score)")
# tckplt.alpPw
# ggsave("tckplt_alpPw.png", tckplt.alpPw, dpi=900,
#        width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "^alpha_pw"))

################################################################################
# Alpha peak frequency
################################################################################
load(file='mod_alphacf_2.RData')
alpha_cf.mod$data$sex <- revalue(alpha_cf.mod$data$sex, c("F"="Female", "M"="Male"))
alpha_cf.mod$data$group <- revalue(alpha_cf.mod$data$group, c("patient"="PD", "control"="HC"))

# Coef plot
plt.alp_cf <- plot_model(alpha_cf.mod, 
                         type = "est",
                         bpe = "mean",
                         bpe.color = "red",
                         prob.inner = .50,
                         prob.outer = .95,
                         show.p = TRUE,
                         p.style = "asterisk",
                         colors="bw",
                         vline.color = "gray",
                         axis.labels = a.labels,
                         show.values = TRUE,
                         title = "Alpha centre frequency"
) 
plt.alp_cf <- plt.alp_cf + scale_y_continuous(limits=c(-2, 2)) + theme_bw()
plt.alp_cf
ggsave("modplt_alpCF.png", plt.alp_cf, dpi=dpi,
       width=p1.w, height=p1.h, units="mm", scale=3)

# Age plot
p1 <- plot_model(alpha_cf.mod, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","sex","group"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 title = "",
                 legend.title = ""
) 
ageplt.alpCF <- p1 + theme_bw() +
  ylab('Frequency (Hz)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)+
  scale_colour_discrete(labels=c("Female","Male"))
ageplt.alpCF
ggsave("ageplt_alpCF.png", ageplt.alpCF, dpi=dpi,
       width=p2.w, height=p2.h, units="mm", scale=3)

# # thickz plot
# p2 <- plot_model(alpha_cf.mod,
#                  type = "pred",
#                  bpe = "mean",
#                  terms = c("thickz","group","sex"),
#                  show.data = TRUE,
#                  ci.lvl = .95,
#                  color = c("red", "blue"),
#                  label= c("HC","PD"),
#                  title = "Cortical thickness",
#                  legend.title = "Group"
# ) 
# tckplt.alpCF <- p2 + theme_bw() +
#   ylab('Frequency (Hz)') +
#   xlab("Cortical thickness (z-score)")
# tckplt.alpCF
# ggsave("tckplt_alpCf.png", tckplt.alpCF, dpi=900,
#        width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "^alpha_cf"))

################################################################################
# Beta power
################################################################################
load(file='mod_betapw_2.RData')
beta_pw.mod$data$sex <- revalue(beta_pw.mod$data$sex, c("F"="Female", "M"="Male"))
beta_pw.mod$data$group <- revalue(beta_pw.mod$data$group, c("patient"="PD", "control"="HC"))

# Coef plot
plt.bet_pw <- plot_model(beta_pw.mod,
                         type = "est",
                         bpe = "mean",
                         bpe.color = "red",
                         prob.inner = .50,
                         prob.outer = .95,
                         show.p = TRUE,
                         p.style = "asterisk",
                         colors="bw",
                         vline.color = "gray",
                         axis.labels = a.labels,
                         show.values = TRUE,
                         title = "Beta power"
) 
plt.bet_pw <- plt.bet_pw + scale_y_continuous(limits=c(-0.2, 0.2)) + theme_bw()
plt.bet_pw
ggsave("modplt_betPow.png", plt.bet_pw, dpi=dpi,
       width=p1.w, height=p1.h, units="mm", scale=3)

# Age plot
p1 <- plot_model(beta_pw.mod, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","sex","group"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 title = "",
                 legend.title = ""
) 
ageplt.betPw <- p1 + theme_bw() +
  ylab('Power (a.u.)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)+
  ylim(0.1, 0.7)+
  scale_colour_discrete(labels=c("Female","Male"))
ageplt.betPw
ggsave("ageplt_betPow.png", ageplt.betPw, dpi=dpi,
       width=p2.w, height=p2.h, units="mm", scale=3)

# # Thickz plot
# p2 <- plot_model(beta_pw.mod,
#                  type = "pred",
#                  bpe = "mean",
#                  terms = c("thickz","group","sex"),
#                  show.data = TRUE,
#                  ci.lvl = .95,
#                  color = c("red", "blue"),
#                  label= c("HC","PD"),
#                  title = "Cortical thickness",
#                  legend.title = "Group"
# ) 
# tckplt.betPw <- p2 + theme_bw() +
#   ylab('Power (a.u.)') +
#   xlab("Cortical thickness (z-score)")
# tckplt.betPw
# ggsave("tckplt_betPw.png", tckplt.betPw, dpi=900,
#        width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "*beta_pw"))

################################################################################
# Beta peak frequency
################################################################################
load(file='mod_betacf_2.RData')
beta_cf.mod$data$sex <- revalue(beta_cf.mod$data$sex, c("F"="Female", "M"="Male"))
beta_cf.mod$data$group <- revalue(beta_cf.mod$data$group, c("patient"="PD", "control"="HC"))

# Coef plot
plt.bet_cf <- plot_model(beta_cf.mod, 
                         type = "est",
                         bpe = "mean",
                         bpe.color = "red",
                         prob.inner = .50,
                         prob.outer = .95,
                         show.p = TRUE,
                         p.style = "asterisk",
                         colors="bw",
                         vline.color = "gray",
                         axis.labels = a.labels,
                         show.values = TRUE,
                         title = "Beta centre frequency"
) 
plt.bet_cf <- plt.bet_cf + scale_y_continuous(limits=c(-4.5, 4.5)) + theme_bw()
plt.bet_cf
ggsave("modplt_betCF.png", plt.bet_cf, dpi=dpi,
       width=p1.w, height=p1.h, units="mm", scale=3)

# Age plot
p1 <- plot_model(beta_cf.mod, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","sex","group"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 title = "",
                 legend.title = ""
) 
ageplt.betCF <- p1 + theme_bw() +
  ylab('Frequency (Hz)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)+
  ylim(10, 30)+
  scale_colour_discrete(labels=c("Female","Male"))
ageplt.betCF
ggsave("ageplt_betCF.png", ageplt.betCF, dpi=dpi,
       width=p2.w, height=p2.h, units="mm", scale=3)

# # thickz plot
# p2 <- plot_model(beta_cf.mod,
#                  type = "pred",
#                  bpe = "mean",
#                  terms = c("thickz","group","sex"),
#                  show.data = TRUE,
#                  ci.lvl = .95,
#                  color = c("red", "blue"),
#                  label= c("HC","PD"),
#                  title = "Cortical thickness",
#                  legend.title = "Group"
# ) 
# tckplt.betCF <- p2 + theme_bw() +
#   ylab('Frequency (Hz)') +
#   xlab("Cortical thickness (z-score)")
# tckplt.betCF
# ggsave("tckplt_betCF.png", tckplt.betCF, dpi=900,
#        width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "*betCf"))

################################################################################
# N events
################################################################################
load(file='mod_neveBF_2.RData')
mod_neve$data$sex <- revalue(mod_neve$data$sex, c("F"="Female", "M"="Male"))
mod_neve$data$group <- revalue(mod_neve$data$group, c("patient"="PD", "control"="HC"))

# Coef plot
plt.neve <- plot_model(mod_neve, 
                         type = "est",
                         bpe = "mean",
                         bpe.color = "red",
                         prob.inner = .50,
                         prob.outer = .95,
                         show.p = TRUE,
                         p.style = "asterisk",
                         colors="bw",
                         vline.color = "gray",
                         axis.labels = a.labels,
                         show.values = TRUE,
                         title = "Burst rate"
) 
plt.neve <- plt.neve + scale_y_continuous(limits=c(0.75, 1.25)) + theme_bw()
plt.neve
ggsave("modplt_neve.png", plt.neve, dpi=dpi,
       width=p1.w, height=p1.h, units="mm", scale=3)

# Age plot
p1 <- plot_model(mod_neve, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","sex","group"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 title = "",
                 legend.title = ""
) 
ageplt.neve <- p1 + theme_bw() +
  ylab('Burst/min') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)+
  ylim(10, 70)+
  scale_colour_discrete(labels=c("Female","Male"))
ageplt.neve
ggsave("ageplt_neve.png", ageplt.neve, dpi=dpi,
       width=p2.w, height=p2.h, units="mm", scale=3)

# # thickz plot
# p2 <- plot_model(mod_neve,
#                  type = "pred",
#                  bpe = "mean",
#                  terms = c("thickz","group","sex"),
#                  show.data = TRUE,
#                  ci.lvl = .95,
#                  color = c("red", "blue"),
#                  label= c("HC","PD"),
#                  title = "Cortical thickness",
#                  legend.title = "Group"
# ) 
# tckplt.neve <- p2 + theme_bw() +
#   ylab('Frequency (Hz)') +
#   xlab("Cortical thickness (z-score)")
# tckplt.neve
# ggsave("tckplt_betCF.png", tckplt.betCF, dpi=900,
#        width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "*neve"))

################################################################################
# Burst duration
################################################################################
load(file='lenmod_2.RData')
# lenmod$data$sex <- revalue(lenmod$data$sex, c("F"="Female", "M"="Male"))
# lenmod$data$group <- revalue(lenmod$data$group, c("patient"="PD", "control"="HC"))
lenmod$data <- data.frame(lenmod$data)

# Coef plot
plt.lene<- plot_model(lenmod, 
                       type = "est",
                       bpe = "mean",
                       bpe.color = "red",
                       prob.inner = .50,
                       prob.outer = .95,
                       show.p = TRUE,
                       p.style = "asterisk",
                       colors="bw",
                       vline.color = "gray",
                       axis.labels = a.labels,
                       show.values = TRUE,
                       title = "Burst duration"
) 
plt.lene <- plt.lene + scale_y_continuous(limits=c(-0.25, 0.25)) + theme_bw()
plt.lene
ggsave("modplt_lene.png", plt.lene, dpi=dpi,
       width=p1.w, height=p1.h, units="mm", scale=3)

# Age plot
p1 <- plot_model(lenmod,
                 type = "pred",
                 bpe = "mean",
                 terms = c('age.centerd',"sex","group"),
                 pred.type = "fe",
                 show.data = TRUE,
                 ci.lvl = .95,
                 title = "",
                 legend.title = ""
)   
ageplt.lene <- p1 + theme_bw() +
  ylab('log(s)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)+
  scale_colour_discrete(labels=c("Female","Male"))
  # scale_y_continuous(trans="exp", limits=c(-10, -1.5), breaks=log(c(0.5, 0.4, 0.3, 0.2, 0.1)), labels=c(0.5, 0.4, 0.3, 0.2, 0.1))
ageplt.lene
ggsave("ageplt_lene.png", ageplt.lene, dpi=dpi,
       width=p2.w, height=p2.h, units="mm", scale=3)

rm(list = ls(pattern = "*lene"))

################################################################################
# Burst interval
################################################################################
load(file='tuemod_2.RData')
tuemod$data <- data.frame(tuemod$data)

# Coef plot
plt.tue <- plot_model(tuemod, 
                      type = "est",
                      bpe = "mean",
                      bpe.color = "red",
                      prob.inner = .50,
                      prob.outer = .95,
                      show.p = TRUE,
                      p.style = "asterisk",
                      colors="bw",
                      vline.color = "gray",
                      axis.labels = a.labels,
                      show.values = TRUE,
                      title = "Burst inteval"
) 
plt.tue <- plt.tue + scale_y_continuous(limits=c(-0.4, 0.4)) + theme_bw()
plt.tue
ggsave("modplt_tue.png", plt.tue, dpi=dpi,
       width=p1.w, height=p1.h, units="mm", scale=3)

# Age plot
p1 <- plot_model(tuemod,
                 type = "pred",
                 bpe = "mean",
                 terms = c('age.centerd',"sex","group"),
                 pred.type = "fe",
                 show.data = TRUE,
                 ci.lvl = .95,
                 title = "",
                 legend.title = ""
)   
ageplt.tue<- p1 + theme_bw() +
  ylab('log(s)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)+
  scale_colour_discrete(labels=c("Female","Male"))
# scale_y_continuous(trans="exp", limits=c(-10, -1.5), breaks=log(c(0.5, 0.4, 0.3, 0.2, 0.1)), labels=c(0.5, 0.4, 0.3, 0.2, 0.1))
ageplt.tue
ggsave("ageplt_tue.png", ageplt.tue, dpi=dpi,
       width=p2.w, height=p2.h, units="mm", scale=3)

################################################################################
# Burst amplitude
################################################################################
load(file='maxmod_2.RData')
maxmod$data <- data.frame(maxmod$data)

# Coef plot
plt.max <- plot_model(maxmod, 
                      type = "est",
                      bpe = "mean",
                      bpe.color = "red",
                      prob.inner = .50,
                      prob.outer = .95,
                      show.p = TRUE,
                      p.style = "asterisk",
                      colors="bw",
                      vline.color = "gray",
                      axis.labels = a.labels,
                      show.values = TRUE,
                      title = "Burst amplitude"
) 
plt.max <- plt.max + scale_y_continuous(limits=c(-0.6, 0.6)) + theme_bw()
plt.max
ggsave("modplt_max.png", plt.max, dpi=dpi,
       width=40, height=50, units="mm", scale=3)

# Age plot
p1 <- plot_model(maxmod,
                 type = "pred",
                 bpe = "mean",
                 terms = c('age.centerd',"sex","group"),
                 pred.type = "fe",
                 show.data = TRUE,
                 ci.lvl = .95,
                 title = "",
                 legend.title = ""
)   
ageplt.max<- p1 + theme_bw() +
  ylab('log-power') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)+
  scale_colour_discrete(labels=c("Female","Male"))
# scale_y_continuous(trans="exp", limits=c(-10, -2), breaks=log(c(0.5, 0.4, 0.3, 0.2, 0.1)), labels=c(0.5, 0.4, 0.3, 0.2, 0.1))
ageplt.max
ggsave("ageplt_max.png", ageplt.max, dpi=dpi,
       width=p2.w, height=p2.h, units="mm", scale=3)

#END