#Plots
library(merTools)
library(ggplot2)
library(tibble)
library(lme4)

## Load data
outdir <- 'X://PD_longrest//output'
figdir <- 'X://PD_longrest//figures'

## make general plot function
main_plot <- function(mod, new.dat, old.dat, yname=NA, link="identity") {
  
  if (is.na(yname)){
    yname <- names(mod$model)[1]
  }
  plt.dat <- new.dat
  
  sims <- sim(mod, n.sims = 1000)
  if (class(sims) == "sim.merMod"){
    fs <- fixef(sims)
  } else {
    fs <- coef(sims)
  }
  
  Xmat <- model.matrix(~ ~ (group+age.centerd+sex+thickz)^3, data = new.dat)
  fitmat <- matrix(ncol = nrow(fs), nrow = nrow(newavg))
  for (i in 1:nrow(fs)) {
    fitmat[, i] <- Xmat %*% as.matrix(fs)[i, ]
  }
  
  plt.dat$lwr <- apply(fitmat, 1, quantile, prob = 0.05)
  plt.dat$median <- apply(fitmat, 1, quantile, prob = 0.5)
  plt.dat$upr <- apply(fitmat, 1, quantile, prob = 0.95) 
  
  if (link=="log"){
    plt.dat$lwr <- exp(plt.dat$lwr)
    plt.dat$median <- exp(plt.dat$median)
    plt.dat$upr <- exp(plt.dat$upr)
  }
  
  plt <- ggplot(aes_string(x="age", y=yname, color="group", shape="sex"), data=old.dat)+
    geom_point(size=2)+
    geom_line(aes(y=median, linetype=sex), size=1, data=plt.dat)+
    geom_ribbon(data=plt.dat, aes(ymin=lwr, ymax=upr, x=age, fill=group, shape=sex),
                alpha = 0.07, inherit.aes = FALSE)+
    scale_color_manual(values=c("red", "blue"), labels=c("HC","PD"))+
    scale_fill_manual(values=c("red", "blue"), labels=c("HC","PD"))+
    theme_bw()+ 
    theme(legend.position = c(.02, .98),
          legend.justification = c("left", "top"),
          legend.box.just = "left",
          legend.margin = margin(1, 1, 1, 1),
          legend.title = element_text(face="bold"),
          legend.spacing = unit(0, "lines"),
          # legend.direction = "horziontal",
          text = element_text(size = 11),
          plot.title = element_text(size=18, vjust=1.5, face="bold", lineheight = NULL, hjust=0.5),
          axis.title = element_text(size = 13, vjust = .5, face="bold"),
          axis.text = element_text(face="bold", size=12),
          panel.grid.major = element_line(size=.5,colour="grey"),
          panel.grid.minor = element_blank())+
    xlab("Age (years)")+
    labs(color  = "Group", fill="Group", linetype = "Sex", shape = "Sex")+
    guides(colour = guide_legend(order = 1), 
           shape = guide_legend(order = 2),
           fill = guide_legend(order = 1),
           linetype = guide_legend(order = 2))
    
    # facet_wrap(~group)
  # plt
  return(plt)
}

###################
## Read models
setwd(outdir)
load(file='mod_finter.RData')
load(file='mod_fslope.RData')
load(file='mod_betapw.RData')
load(file='mod_betacf.RData')
load(file='mod_alphapw.RData')
load(file='mod_alphacf.RData')
load(file='lenmod.RData')
load(file='tuemod.RData')
load(file='maxmod.RData')
load(file='mod_neve.RData')

## Read data
load(file='X://PD_longrest//groupanalysis//alldata_subj2.Rdata')
load(file='X://PD_longrest//groupanalysis//bbdata2.Rdata')

## Make dummy data
agespan <- seq(min(alldata$age.centerd), max(alldata$age.centerd), 0.01)
agespan2 <- seq(min(alldata$age),max(alldata$age), 0.01)

new.dat <- data.frame(age.centerd=rep(agespan,4), 
                      group=as.factor(rep(c("patient","control"),each=length(agespan)*2)), 
                      sex=as.factor(rep(rep(c("M","F"), each=length(agespan)), 2)),
                      thickz=0,
                      age=rep(agespan2, 4))

bbsum <- aggregate(cbind(leneve, tueeve, maxeve)~subj, data=bbdata, FUN=median)
alldata <- merge(alldata, bbsum, by="subj")

### Plot
plt_fint <- main_plot(finter.Full3, new.dat, alldata) + 
  ylab('1/f intercept (a.u.)') + ggtitle('1/f Intercept')
plt_fslo <- main_plot(fslope.Full3, new.dat, alldata) +
  ylab('1/f exponent (a.u.)') + ggtitle('1/f Exponent')
plt_apow <- main_plot(alpha_pw.Full3, new.dat, alldata) +
  ylab('Power (a.u.)') + ggtitle('Alpha Power')
plt_acfq <- main_plot(alpha_cf.Full3, new.dat, alldata) +
  ylab('Frequency (Hz)') + ggtitle('Alpha Peak Frequency')
plt_bpow <- main_plot(beta_pw.Full3, new.dat, alldata) +
  ylab('Power (a.u.)') + ggtitle('Beta Power')
plt_bcfq <- main_plot(beta_cf.Full3, new.dat, alldata) +
  ylab('Frequency (Hz)') + ggtitle('Beta Peak Frequency')
plt_lene <- main_plot(lenmod, new.dat, alldata, yname="leneve", link = "log") +
  ylab('Time (ms)') + ggtitle('Burst Duration')
plt_maxe <- main_plot(maxmod, new.dat, alldata, yname="maxeve", link = "log") +
  ylab('Amplitude') + ggtitle('Burst Amplitude')
plt_tuee <- main_plot(tuemod, new.dat, alldata, yname="tueeve", link = "log") +
  ylab('Time (ms)') + ggtitle('Time Between Bursts')
plt_neve <- main_plot(mod.neve.Full3, new.dat, alldata, link = "log") +
  ylab('Burst per minute') + ggtitle('Burst Rate')

# Inspect
plt_fint
plt_fslo
plt_apow
plt_acfq
plt_bpow
plt_bcfq
plt_lene
plt_maxe
plt_tuee
plt_neve

## Save
setwd(outdir)
ggsave("plt_fint.jpeg", plot=plt_fint, device="jpeg", units="mm", width=60, height=50, dpi=600, scale=3)
ggsave("plt_fslo.jpeg", plot=plt_fslo, device="jpeg", units="mm", width=60, height=50, dpi=600, scale=3)
ggsave("plt_apow.jpeg", plot=plt_apow, device="jpeg", units="mm", width=60, height=50, dpi=600, scale=3)
ggsave("plt_acfq.jpeg", plot=plt_acfq, device="jpeg", units="mm", width=60, height=50, dpi=600, scale=3)
ggsave("plt_bpow.jpeg", plot=plt_bpow, device="jpeg", units="mm", width=60, height=50, dpi=600, scale=3)
ggsave("plt_bcfq.jpeg", plot=plt_bcfq, device="jpeg", units="mm", width=60, height=50, dpi=600, scale=3)
ggsave("plt_lene.jpeg", plot=plt_lene, device="jpeg", units="mm", width=60, height=50, dpi=600, scale=3)
ggsave("plt_maxe.jpeg", plot=plt_maxe, device="jpeg", units="mm", width=60, height=50, dpi=600, scale=3)
ggsave("plt_tuee.jpeg", plot=plt_tuee, device="jpeg", units="mm", width=60, height=50, dpi=600, scale=3)
ggsave("plt_neve.jpeg", plot=plt_neve, device="jpeg", units="mm", width=60, height=50, dpi=600, scale=3)

#END