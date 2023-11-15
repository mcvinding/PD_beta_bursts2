# Plot summary of burst features
#
# Vinding, M. C., Eriksson, A., Low, C. M. T., Waldthaler, J., Ferreira, D., Ingvar, M.,  Svenningsson, P., & Lundqvist, D. (2021). Different features of the cortical sensorimotor rhythms are uniquely linked to the severity of specific symptoms in Parkinson's disease [Preprint]. medRxiv.org https://doi.org/10.1101/2021.06.27.21259592
#
library(ggplot2)

# Load data
load(file='/home/mikkel/PD_longrest/groupanalysis/alldata_subj2.Rdata')
load(file='/home/mikkel/PD_longrest/groupanalysis/bbdata2.Rdata')
setwd('~/PD_longrest/figures')

# Plot burst rate histogram
mean.br <- aggregate(alldata$nevent.u.m2, list(alldata$group), mean)

br.hist <- ggplot(alldata, aes(x=nevent.u.m2, color=group, fill=group)) +
  geom_histogram(binwidth=8, alpha=0.5) +
  # geom_vline(data=mean.br, aes(xintercept=x, color=Group.1), linetype="dashed")+
  theme_bw()+
  scale_fill_manual(values=c('red','blue'), labels=c("HC","PD")) +
  scale_color_manual(values=c('red','blue'), guide="none") +
  labs(x = "Bursts/min", y = "", title="Burst rate") +
  theme(legend.position = c(.9, .9),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1)),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
br.hist
ggsave("hist_bbrate.jpg", plot=br.hist, 
       device="jpg", units="mm", width=40, height=30, dpi=600, scale=4)

# Burst amplitudes
mean1.ba <- aggregate(bbdata$maxeve, list(bbdata$subj, bbdata$group), mean)
mean2.ba <- aggregate(mean1.ba$x, list(mean1.ba$Group.2), mean)

ba.hist <- ggplot(bbdata, aes(x=maxeve, color=group)) +
  geom_density(bw=0.15, size=1, show.legend=TRUE) +
  geom_vline(data=mean2.ba, aes(xintercept=x, color=Group.1), linetype="dashed", alpha=0.8, show.legend=FALSE) +
  theme_bw()+
  scale_color_manual(values=c('red','blue'), labels=c("HC","PD")) +
  labs(x = "Burst amp. (F)", y = "Density", title="Burst amplitude") +
  theme(legend.position = c(.9, .9),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1)),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  xlim(c(0.7, 7))
ba.hist
ggsave("hist_bamp.jpg", plot=ba.hist, 
       device="jpg", units="mm", width=40, height=30, dpi=600, scale=4)

# Burst duration
mean1.bd <- aggregate(bbdata$leneve.ms, list(bbdata$subj, bbdata$group), median)
mean2.bd <- aggregate(mean1.bd$x, list(mean1.bd$Group.2), mean)

bd.hist <- ggplot(bbdata, aes(x=leneve.ms, color=group)) +
  geom_density(bw=3, size=1, show.legend=TRUE) +
  geom_vline(data=mean2.bd, aes(xintercept=x, color=Group.1), linetype="dashed", alpha=0.8, show.legend=FALSE) +
  theme_bw()+
  scale_color_manual(values=c('red','blue'), labels=c("HC","PD")) +
  labs(x = "Duration (ms)", y = "Density", title="Burst duration") +
  theme(legend.position = c(.9, .9),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1)),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(c(0.7, 250))
bd.hist
ggsave("hist_bdur.jpg", plot=bd.hist, 
       device="jpg", units="mm", width=40, height=30, dpi=600, scale=4)

# Burst interval
mean1.bi <- aggregate(bbdata$tueeve.ms, list(bbdata$subj, bbdata$group), median)
mean2.bi <- aggregate(mean1.bi$x, list(mean1.bi$Group.2), mean)

bi.hist <- ggplot(bbdata, aes(x=tueeve.ms, color=group)) +
  geom_density(bw=50, size=1, show.legend=TRUE) +
  geom_vline(data=mean2.bi, aes(xintercept=x, color=Group.1), linetype="dashed", alpha=0.8, show.legend=FALSE) +
  theme_bw()+
  scale_color_manual(values=c('red','blue'), labels=c("HC","PD")) +
  labs(x = "Duration (ms)", y = "Density", title="Burst interval") +
  theme(legend.position = c(.9, .9),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1)),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(c(0, 4500))
bi.hist
ggsave("hist_bint.jpg", plot=bi.hist, 
       device="jpg", units="mm", width=40, height=30, dpi=600, scale=4)


ggplot(bbdata, aes(x=group, y=maxeve)) +
  geom_flat_violin(scale = "count", trim = FALSE, aes(fill=group), alpha=0.6)



l.den <- ggplot(leneve.data, aes(x=eve.len.ms, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(bw=2.5,alpha=0.6) +
  # geom_line(stat="density",data=pp.len.long, aes(x=pp, fill=NA), bw=2.5, linetype = "dashed", size=0.5, color='black') +
  facet_wrap(~label) + 
  xlim(0,200)+
  labs(x='Duration (ms)', y = "Density")+
  theme_bw()+
  scale_fill_manual(values=c('red','blue'), guide=FALSE)+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", lineheight=12, size=12),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, lineheight=12, face="bold", size=16),
        axis.title = element_text(face="bold", lineheight=11, size=12),
        axis.text = element_text(color="black", lineheight=9, size=10),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle('Beta event duration')
l.den