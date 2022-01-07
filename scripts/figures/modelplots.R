# Plot results of regression models
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
setwd('X://PD_longrest//output')

tab_model(mod.neve.Full3, transform = NULL)
tab_model(lenmod, transform = NULL)

## PLOT ANALYSIS 1
# Variable names
a.labels <- c("Group", "Age", "Sex", "Cortical thickness", "Age^2",
              "Group:Age", "Group:Sex", "Group:Cortical thickness", "Age:Sex", "Age:Cortical thickness", "Sex:Cortical thickness",
              "Group:Age:Sex", "Group:Age:Cortical thickness", "Group:Sex:Cortical thickness", "Age:Sex:Cortical thickness"
)
a.labels <- rev(a.labels)

## Plot burst stats
# Load models
load(file='mod_neveBF.RData')
load(file='lenmod.RData')
load(file='tuemod.RData')
load(file='maxmod.RData')

mods <- list(mod.neve.Full3, lenmod, tuemod, maxmod)
m.labels <- c("Rate", "Length", "Interval", "Amplitude")
mods <- rev(mods)
m.labels <- rev(m.labels)

# Plot
eve.plt <- plot_models(mods, grid = TRUE, std.est = "std", transform = NULL,
            axis.labels = a.labels, m.labels = m.labels, vline.color = "gray",
            show.values = TRUE, show.p = FALSE, p.shape = FALSE, auto.label = FALSE, show.legend = FALSE,
            digits = 3, value.size=3, colors="bw")
eve.plt <- eve.plt + scale_y_continuous(limits=c(-1.25, 1.4)) + theme_bw()
eve.plt

ggsave("modplt_eve.png", eve.plt, dpi=900,
       width=100, height=50, units="mm", scale=3)

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

ggsave("modplt_updrs.png", updrs.plt, dpi=900,
       width=140, height=50, units="mm", scale=3)

