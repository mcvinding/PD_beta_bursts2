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

# Variable names
a.labels <- c("Group", "Age", "Sex", "Cortical thickness", "Age^2",
              "Group:Age", "Group:Sex", "Group:Cortical thickness", "Age:Sex", "Age:Cortical thickness", "Sex:Cortical thickness")
              # "Group:Age:Sex", "Group:Age:Cortical thickness", "Group:Sex:Cortical thickness", "Age:Sex:Cortical thickness")

a.labels <- rev(a.labels)

age.convert <- function(x) {
  lab = round(x + mean(bbdata$age))
  return(lab)
}

# thick.convert <- function(x) {
#   lab = round(x + mean(bbdata$thick))
#   return(lab)
# }

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
ggsave("modplt_finter.png", plt.finter, dpi=900,
       width=60, height=80, units="mm", scale=2)

# Age plot
p1 <- plot_model(finter.mod, 
           type = "pred",
           bpe = "median",
           terms = c("age.centerd","group","sex"),
           show.data = TRUE,
           ci.lvl = .95,
           color = c("red", "blue"),
           label= c("HC","PD"),
           title = "Age",
           legend.title = "Group"
) 
ageplt.finter <- p1 + theme_bw() +
  ylab('1/f intercept (a.u.)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)
  # scale_color_manual(labels=c("HC","PD"), values = c("red", "blue")) +
ageplt.finter

ggsave("ageplt_finter.png", ageplt.finter, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

# Thickz plot
p2 <- plot_model(finter.mod, 
                 type = "pred",
                 bpe = "mean",
                 terms = c("thickz","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Cortical thickness",
                 legend.title = "Group"
) 
tckplt.finter <- p2 + theme_bw() +
  ylab('1/f intercept (a.u.)') +
  xlab("Cortical thickness (z-score)")
  # scale_color_manual(labels=c("HC","PD"), values = c("red","blue", "red","blue"))
  # scale_fill_manual(labels=c("HC","PD"), values = c("red", "blue")) +
tckplt.finter

ggsave("tckplt_finter.png", tckplt.finter, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "^finter"))

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
ggsave("modplt_fslope.png", plt.fslope, dpi=900,
       width=60, height=80, units="mm", scale=2)

# Age plot
p1 <- plot_model(fslope.mod, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Age",
                 legend.title = "Group"
) 
ageplt.fslope <- p1 + theme_bw() +
  ylab('1/f exponent (a.u.)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)
ageplt.fslope
ggsave("ageplt_fslope.png", ageplt.fslope, dpi=900,
       width=100, height=80, units="mm", scale=2)

# Thickz plot
p2 <- plot_model(fslope.mod, 
                 type = "pred",
                 bpe = "mean",
                 terms = c("thickz","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Cortical thickness",
                 legend.title = "Group"
) 
tckplt.fslope <- p2 + theme_bw() +
  ylab('1/f intercept (a.u.)') +
  xlab("Cortical thickness (z-score)")
tckplt.fslope
ggsave("tckplt_fslope.png", tckplt.fslope, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "^fslope"))

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
ggsave("modplt_alpPow.png", plt.alp_pw, dpi=900,
       width=60, height=80, units="mm", scale=2)

# Age plot
p1 <- plot_model(alpha_pw.mod, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Age",
                 legend.title = "Group"
) 
ageplt.alpPw <- p1 + theme_bw() +
  ylab('Power (a.u.)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)
ageplt.alpPw
ggsave("ageplt_alpPow.png", ageplt.alpPw, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

# Thickz plot
p2 <- plot_model(alpha_pw.mod,
                 type = "pred",
                 bpe = "mean",
                 terms = c("thickz","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Cortical thickness",
                 legend.title = "Group"
) 
tckplt.alpPw <- p2 + theme_bw() +
  ylab('Power (a.u.)') +
  xlab("Cortical thickness (z-score)")
tckplt.alpPw
ggsave("tckplt_alpPw.png", tckplt.alpPw, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

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
ggsave("modplt_alpCF.png", plt.alp_cf, dpi=900,
       width=60, height=80, units="mm", scale=2)

# Age plot
p1 <- plot_model(alpha_cf.mod, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Age",
                 legend.title = "Group"
) 
ageplt.alpCF <- p1 + theme_bw() +
  ylab('Frequency (Hz)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)
ageplt.alpCF
ggsave("ageplt_alpCF.png", ageplt.alpCF, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

# thickz plot
p2 <- plot_model(alpha_cf.mod,
                 type = "pred",
                 bpe = "mean",
                 terms = c("thickz","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Cortical thickness",
                 legend.title = "Group"
) 
tckplt.alpCF <- p2 + theme_bw() +
  ylab('Frequency (Hz)') +
  xlab("Cortical thickness (z-score)")
tckplt.alpCF
ggsave("tckplt_alpCf.png", tckplt.alpCF, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

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
ggsave("modplt_betPow.png", plt.bet_pw, dpi=900,
       width=60, height=80, units="mm", scale=2)

# Age plot
p1 <- plot_model(beta_pw.mod, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Age",
                 legend.title = "Group"
) 
ageplt.betPw <- p1 + theme_bw() +
  ylab('Power (a.u.)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)
ageplt.betPw
ggsave("ageplt_betPow.png", ageplt.betPw, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

# Thickz plot
p2 <- plot_model(beta_pw.mod,
                 type = "pred",
                 bpe = "mean",
                 terms = c("thickz","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Cortical thickness",
                 legend.title = "Group"
) 
tckplt.betPw <- p2 + theme_bw() +
  ylab('Power (a.u.)') +
  xlab("Cortical thickness (z-score)")
tckplt.betPw
ggsave("tckplt_betPw.png", tckplt.betPw, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "^beta_pw"))

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
                         title = "Alpha centre frequency"
) 
plt.bet_cf <- plt.bet_cf + scale_y_continuous(limits=c(-4.5, 4.5)) + theme_bw()
plt.bet_cf
ggsave("modplt_betCF.png", plt.bet_cf, dpi=900,
       width=60, height=80, units="mm", scale=2)

# Age plot
p1 <- plot_model(beta_cf.mod, 
                 type = "pred",
                 bpe = "median",
                 terms = c("age.centerd","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Age",
                 legend.title = "Group"
) 
ageplt.betCF <- p1 + theme_bw() +
  ylab('Frequency (Hz)') +
  scale_x_continuous(name="Age (years)", breaks=seq(-20, 20, by=10), labels=age.convert)
ageplt.betCF
ggsave("ageplt_betCF.png", ageplt.betCF, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

# thickz plot
p2 <- plot_model(beta_cf.mod,
                 type = "pred",
                 bpe = "mean",
                 terms = c("thickz","group","sex"),
                 show.data = TRUE,
                 ci.lvl = .95,
                 color = c("red", "blue"),
                 label= c("HC","PD"),
                 title = "Cortical thickness",
                 legend.title = "Group"
) 
tckplt.betCF <- p2 + theme_bw() +
  ylab('Frequency (Hz)') +
  xlab("Cortical thickness (z-score)")
tckplt.betCF
ggsave("tckplt_betCF.png", tckplt.betCF, dpi=900,
       width=100, height=80, units="mm", scale=1.5)

rm(list = ls(pattern = "^beta_cf"))

################################################################################
## Plot burst stats
# Load models
load(file='mod_neveBF.RData')
load(file='lenmod.RData')
load(file='tuemod.RData')
load(file='maxmod.RData')

mods <- list(mod.neve.Full2, lenmod, tuemod, maxmod)
m.labels <- c("Rate", "Length", "Interval", "Amplitude")
mods <- rev(mods)
m.labels <- rev(m.labels)

# Plot
eve.plt <- plot_models(mod.neve.ST, grid = TRUE, std.est = "std", transform = NULL, ci.lvl = 0.95,
            axis.labels = a.labels, m.labels = m.labels, vline.color = "gray", show.intercept = FALSE,
            show.values = TRUE, show.p = FALSE, p.shape = FALSE, auto.label = TRUE, show.legend = FALSE,
            digits = 3, value.size=3, colors="bw")
eve.plt <- eve.plt + scale_y_continuous(limits=c(-1.25, 1.4)) + theme_bw()
eve.plt <- eve.plt + scale_y_continuous(limits=c(-0.25, 0.25)) + theme_bw()
eve.plt

# USE THIS PLOT!!!1
# and different inner/outer probabilities of the HDI
plot_model(mod.neve.Full2_bf, 
           type = "est",
  bpe = "mean",
  # bpe.style = "dot",
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
) + scale_y_continuous(limits=c(0.7, 1.3)) + theme_bw()
                       
ggsave("modplt_eve.png", eve.plt, dpi=900,
       width=100, height=50, units="mm", scale=3)

plot_model(mod.neve.Full2_bf, 
           type = "pred",
           bpe = "mean",
           # bpe.style = "dot",
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

## Plot PSD stats
# Load models
load(file='mod_finter.RData')
load(file='mod_fslope.RData')
load(file='mod_betapw.RData')
load(file='mod_betacf.RData')
load(file='mod_alphapw.RData')
load(file='mod_alphacf.RData')

mods <- list(beta_pw.Full3,
             beta_cf.Full3,
             alpha_pw.Full3,
             alpha_cf.Full3,
             finter.Full3,
             fslope.Full3)
m.labels <- c("Beta power", "Beta centre freq.", "Alpha power", "Alpha centre freq.", "1/f offset", "1/f exponent")
mods <- rev(mods)
m.labels <- rev(m.labels)

psd.plt <- plot_models(mods, grid = TRUE, std.est = "std2", transform = NULL,
                       axis.labels = a.labels, m.labels = m.labels, vline.color = "gray",
                       show.values = TRUE, show.p = FALSE, p.shape = FALSE, auto.label = FALSE, show.legend = FALSE,
                       digits = 3, value.size=3, colors="bw")
psd.plt <- psd.plt + scale_y_continuous(limits=c(-1.25, 1.5)) + theme_bw()
psd.plt

ggsave("modplt_psd.png", psd.plt, dpi=900,
       width=140, height=50, units="mm", scale=3)

## PLOT ANALYSIS 2

# Load models
load(file="updrsmods.RData")

mods <- list(F1mod.x, F2mod.x, F3mod.x, F45mod.x, F6mod.x, F7mod.x)
mods <- rev(mods)

a.labels <- c("Rate", "Length", "Interval", "Amplitude",
              "1/f offset", "1/f exponent", "Alpha power", "Alpha centre freq.", "Beta power", "Beta centre freq.",
              "Age", "Sex", "Cortical thickness")
a.labels <- rev(a.labels)
m.labels <- c("Midline function", 
              "Rest tremor", 
              "Rigidity", 
              "Bradykinesia upper limb", 
              "Postural and kinetic tremor", 
              "Bradykinesia lower limb")
m.labels <- rev(m.labels)

updrs.plt <- plot_models(mods, grid = TRUE, std.est = NULL, transform = NULL, axis.lim=c(-2,2),
                         axis.labels = a.labels, m.labels = m.labels, vline.color = "gray",
                         show.values = TRUE, show.p = FALSE, p.shape = FALSE, auto.label = FALSE, show.legend = FALSE,
                         digits = 3, value.size=3, colors="bw")
updrs.plt <- updrs.plt + scale_y_continuous(limits=c(-2.5, 2)) + theme_bw()
updrs.plt

updrs.plt <- plot_models(F45mod.x, show.intercept = FALSE, rm.terms = 'b_Intercept')
updrs.plt


, grid = TRUE, std.est = NULL, transform = NULL, axis.lim=c(-2,2),
                         axis.labels = a.labels, m.labels = m.labels, vline.color = "gray",
                         show.values = TRUE, show.p = FALSE, p.shape = FALSE, auto.label = FALSE, show.legend = FALSE,
                         digits = 3, value.size=3, colors="bw")
updrs.plt


ggsave("modplt_updrs.png", updrs.plt, dpi=900,
       width=140, height=50, units="mm", scale=3)

