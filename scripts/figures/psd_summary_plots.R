# Plot summary of PSD features from FOOOF analysis
#
# Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M., 
#  Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor 
#  rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease [Preprint]. 
#  medRxiv.org https://doi.org/10.1101/2021.06.27.21259592
#

library(ggplot2)

# load source code
devtools::source_gist("2a1bb0133ff568cbe28d", filename = "geom_flat_violin.R") # sourced from github "dgrtwo/geom_flat_violin.R"

# Load data
load(file='/home/mikkel/PD_longrest/groupanalysis/alldata_subj2.Rdata')
setwd('~/PD_longrest/figures')

# 1/f intercept/offset plot
icp.plt <- ggplot(alldata, aes(x = group, y = a_intercept)) + 
  geom_flat_violin(scale = "count", trim = FALSE, aes(fill=group), alpha=0.6) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", position = position_nudge(0.05)) + 
  geom_dotplot(binaxis = "y", dotsize = 0.75, stackdir = "down", position = position_nudge(-0.025)) + 
  labs(x = "", y = "", title="1/f offset") +
  scale_x_discrete(labels=c("HC","PD")) +
  scale_fill_manual(values=c('red','blue'), guide=FALSE)+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
icp.plt
ggsave("hist_oof_intercept.jpg", plot=icp.plt,
       device="jpg", units="mm", width=30, height=30, dpi=600, scale=4)

# 1/f exponent
exp.plt <- ggplot(alldata, aes(x = group, y = a_slope)) + 
  geom_flat_violin(scale = "count", trim = FALSE, aes(fill=group), alpha=0.6) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", position = position_nudge(0.05)) + 
  geom_dotplot(binaxis = "y", dotsize = 0.75, stackdir = "down",  position = position_nudge(-0.025)) + 
  labs(x = "", y = "", title="1/f exponent") +
  scale_x_discrete(labels=c("HC","PD")) +
  scale_fill_manual(values=c('red','blue'), guide=FALSE)+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
exp.plt
ggsave("hist_oof_exp.jpeg", plot=exp.plt,
       device="jpg", units="mm", width=30, height=30, dpi=600, scale=4)

# 1/f beta power
bpow.plt <- ggplot(alldata, aes(x = group, y = beta_pw)) + 
  geom_flat_violin(scale = "count", trim = FALSE, aes(fill=group), alpha=0.6) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", position = position_nudge(0.05)) + 
  geom_dotplot(binaxis = "y", dotsize = 0.75, stackdir = "down",  position = position_nudge(-0.025)) + 
  labs(x = "", y = "", title="Beta power") +
  scale_x_discrete(labels=c("HC","PD")) +
  scale_fill_manual(values=c('red','blue'), guide=FALSE)+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
bpow.plt
ggsave("hist_oof_bpow.jpeg", plot=bpow.plt,
       device="jpg", units="mm", width=30, height=30, dpi=600, scale=4)

# 1/f alpha power
apow.plt <- ggplot(alldata, aes(x = group, y = alpha_pw)) + 
  geom_flat_violin(scale = "count", trim = FALSE, aes(fill=group), alpha=0.6) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", position = position_nudge(0.05)) + 
  geom_dotplot(binaxis = "y", dotsize = 0.75, stackdir = "down",  position = position_nudge(-0.025)) + 
  labs(x = "", y = "", title="Alpha power") +
  scale_x_discrete(labels=c("HC","PD")) +
  scale_fill_manual(values=c('red','blue'), guide=FALSE)+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
apow.plt
ggsave("hist_oof_apow.jpeg", plot=apow.plt,
       device="jpg", units="mm", width=30, height=30, dpi=600, scale=4)

# 1/f beta center freq.
bcf.plt <- ggplot(alldata, aes(x = group, y = beta_cf)) + 
  geom_flat_violin(scale = "count", trim = FALSE, aes(fill=group), alpha=0.6) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", position = position_nudge(0.05)) + 
  geom_dotplot(binaxis = "y", dotsize = 0.75, stackdir = "down",  position = position_nudge(-0.025)) + 
  labs(x = "", y = "", title="Beta center frequency") +
  scale_x_discrete(labels=c("HC","PD")) +
  scale_fill_manual(values=c('red','blue'), guide=FALSE)+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
bcf.plt
ggsave("hist_oof_bcf.jpeg", plot=bcf.plt,
       device="jpg", units="mm", width=30, height=30, dpi=600, scale=4)

# 1/f alpha center freq.
acf.plt <- ggplot(alldata, aes(x = group, y = alpha_cf)) + 
  geom_flat_violin(scale = "count", trim = FALSE, aes(fill=group), alpha=0.6) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", position = position_nudge(0.05)) + 
  geom_dotplot(binaxis = "y", dotsize = 0.75, stackdir = "down",  position = position_nudge(-0.025)) + 
  labs(x = "", y = "", title="Alpha center frequency") +
  scale_x_discrete(labels=c("HC","PD")) +
  scale_fill_manual(values=c('red','blue'), guide=FALSE)+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
acf.plt
ggsave("hist_oof_acf.jpeg", plot=acf.plt,
       device="jpg", units="mm", width=30, height=30, dpi=600, scale=4)
