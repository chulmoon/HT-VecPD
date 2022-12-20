library(tidyverse)
library(reshape2)

#################################################################
# Method comparison
#################################################################

# FPR
## Persistence landscape
load("test_persistence_landscape_fpr.Rdata")
fprpl = rowSums(perz<0.05)/500
## Kernel
fprpwgk = c(0.042,0.048,0.072,0.032)
## Two stage
load("test_two_stage_constant_fpr.Rdata")
fprts = rowMeans(pvals[,,5])
## Persistence diagram
load("test_persistence_diagram_fpr.Rdata")
fprpd = rowSums(perz<0.05)/500

fprres=data.frame(cbind(noise=seq(0.05,0.2,by=0.05),fprts,fprpd,fprpl,fprpwgk))

colnames(fprres)=c("noise","Two-stage","PD", "PL","Kernel")
fprdf=reshape2::melt(fprres, id.var = "noise")
ggplot(fprdf,aes(x=noise,y=value)) +
  geom_point(size=2.5,aes(shape=variable,group=variable,color=variable))+
  geom_line(size=1.05,aes(linetype=variable,group=variable,color=variable)) +
  labs(x="Noise Level",y="False Positive Rate")+
  guides(shape=guide_legend("Method"),linetype=guide_legend("Method"),color=guide_legend("Method"))+
  theme_bw()+
  ylim(0,0.1)+
  theme(legend.position=c(0.85,0.80))

# Power
# PL
load("test_persistence_landscape_power.Rdata")
powerpl=rowSums(perz<0.05)/500
# PWGK
powerpwgk = c(0.974,0.604,0.278,0.184)
# two stage
load("test_two_stage_constant_power.Rdata")
powerts = rowMeans(pvals[,,5])
# PD
load("/test_persistence_diagram_power.Rdata")
powerpd=rowSums(perz<0.05)/500

meanres=data.frame(cbind(noise=seq(0.05,0.2,by=0.05),powerts,powerpd,powerpl,powerpwgk))

colnames(meanres)=c("noise","Two-stage","PD", "PL","Kernel")
meandf=reshape2::melt(meanres, id.var = "noise")

ggplot(meandf,aes(x=noise,y=value)) +
  geom_point(size=2.5,aes(shape=variable,group=variable,color=variable))+
  geom_line(size=1.05,aes(linetype=variable,group=variable,color=variable)) +
  labs(x="Noise Level",y="Power")+
  guides(shape=guide_legend("Method"),linetype=guide_legend("Method"),color=guide_legend("Method"))+
  theme_bw()+
  theme(legend.position=c(0.85,0.80))

#################################################################
# Threshold comparison
#################################################################
# power

load("test_two_stage_constant_power.Rdata")
threspower=cbind(colMeans(pvals[1,,]),colMeans(pvals[2,,]),colMeans(pvals[3,,]),colMeans(pvals[4,,]))
load("test_two_stage_constant_nofilter_power.Rdata")
threspower=data.frame(rbind(rowMeans(pvals[,,1]),threspower))

threspower=cbind(c("No pre-filtering","C=0%","C=20%","C=40%","C=60%","C=80%"),threspower)
colnames(threspower)=c("filter","0.05","0.1","0.15", "0.2")
threspowerdf=reshape2::melt(threspower, id.var = "filter")
threspowerdf$variable=as.numeric(as.character(threspowerdf$variable))
threspowerdf$filter=factor(threspowerdf$filter,levels=c("C=80%","C=60%","C=40%","C=20%","C=0%","No pre-filtering"))

ggplot(threspowerdf,aes(x=variable,y=value)) +
  geom_point(size=2.5,aes(shape=filter,group=filter,color=filter))+
  geom_line(size=1.05,aes(linetype=filter,group=filter,color=filter)) +
  labs(x="Noise Level",y="Power")+
  guides(shape=guide_legend("Filtering"),linetype=guide_legend("Filtering"),color=guide_legend("Filtering"))+
  theme_bw()+
  theme(legend.position=c(0.805,0.76))

# fpr
load("test_two_stage_constant_fpr.Rdata")
thresfpr=cbind(colMeans(pvals[1,,]),colMeans(pvals[2,,]),colMeans(pvals[3,,]),colMeans(pvals[4,,]))
load("test_two_stage_constant_nofilter_fpr.Rdata")
thresfpr=data.frame(rbind(rowMeans(pvals[,,1]),thresfpr))

thresfpr=cbind(c("No pre-filtering","C=0%","C=20%","C=40%","C=60%","C=80%"),thresfpr)
colnames(thresfpr)=c("filter","0.05","0.1","0.15", "0.2")
thresfprdf=reshape2::melt(thresfpr, id.var = "filter")
thresfprdf$variable=as.numeric(as.character(thresfprdf$variable))
thresfprdf$filter=factor(thresfprdf$filter,levels=c("C=80%","C=60%","C=40%","C=20%","C=0%","No pre-filtering"))

ggplot(thresfprdf,aes(x=variable,y=value)) +
  geom_point(size=2.5,aes(shape=filter,group=filter,color=filter))+
  geom_line(size=1.05,aes(linetype=filter,group=filter,color=filter)) +
  labs(x="Noise Level",y="False Positive Rate")+
  guides(shape=guide_legend("Filtering"),linetype=guide_legend("Filtering"),color=guide_legend("Filtering"))+
  theme_bw()+
  ylim(0,0.1)+
  theme(legend.position=c(0.805,0.76))



#################################################################
# weight comparison
#################################################################

# FPR
# constant
load("test_two_stage_constant_fpr.Rdata")
fprconstant60=rowMeans(pvals[,,4])
fprconstant80=rowMeans(pvals[,,5])
# arctan
load("test_two_stage_arctan_fpr.Rdata")
fprarctan60=rowMeans(pvals[,,4])
fprarctan80=rowMeans(pvals[,,5])
# linear
load("test_two_stage_linear_fpr.Rdata")
fprlinear60=rowMeans(pvals[,,4])
fprlinear80=rowMeans(pvals[,,5])

fprres=data.frame(cbind(noise=seq(0.05,0.2,by=0.05),fprconstant60, fprconstant80,fprarctan60, fprarctan80,fprlinear60,fprlinear80))

colnames(fprres)=c("noise","Constant (C=60%)","Constant (C=80%)","Arctangent (C=60%)", "Arctangent (C=80%)", "Linear (C=60%)","Linear (C=80%)")
fprdf=reshape2::melt(fprres, id.var = "noise")

fprdf$Weight=rep(c("Constant","Arctangent","Linear"),each=8)
fprdf$Threshold=rep(rep(c("C=60%","C=80%"),each=4),3)

ggplot(fprdf,aes(x=noise,y=value,group=interaction(Weight,Threshold))) +
  geom_point(size=2.5,aes(shape=Threshold,color=Weight))+
  geom_line(size=1.05,aes(linetype=Threshold,color=Weight)) +
  scale_colour_discrete(guide = guide_legend(override.aes = list(linetype = c(1,1,1), shape = c(NA,NA,NA)))) +
  labs(x="Noise Level",y="False Positive Rate")+
  theme_bw()+
  ylim(0,0.1)+
  theme(legend.position=c(0.735,0.86),legend.box = "horizontal")


# power
# constant
load("test_two_stage_constant_power.Rdata")
powerconstant60=rowMeans(pvals[,,4])
powerconstant80=rowMeans(pvals[,,5])
# arctan
load("test_two_stage_arctan_power.Rdata")
powerarctan60=rowMeans(pvals[,,4])
powerarctan80=rowMeans(pvals[,,5])
# linear
load("test_two_stage_linear_power.Rdata")
powerlinear60=rowMeans(pvals[,,4])
powerlinear80=rowMeans(pvals[,,5])

powerres=data.frame(cbind(noise=seq(0.05,0.2,by=0.05),powerconstant60, powerconstant80,powerarctan60, powerarctan80,powerlinear60,powerlinear80))

colnames(powerres)=c("noise","Constant (C=60%)","Constant (C=80%)","Arctangent (C=60%)", "Arctangent (C=80%)", "Linear (C=60%)","Linear (C=80%)")
powerdf=reshape2::melt(powerres, id.var = "noise")

powerdf$Weight=rep(c("Constant","Arctangent","Linear"),each=8)
powerdf$Threshold=rep(rep(c("C=60%","C=80%"),each=4),3)

ggplot(powerdf,aes(x=noise,y=value,group=interaction(Weight,Threshold))) +
  geom_point(size=2.5,aes(shape=Threshold,color=Weight))+
  geom_line(size=1.05,aes(linetype=Threshold,color=Weight)) +
  scale_colour_discrete(guide = guide_legend(override.aes = list(linetype = c(1,1,1),
                                                                 shape = c(NA,NA,NA)) )) +
  labs(x="Noise Level",y="Power")+
  theme_bw()+
  theme(legend.position=c(0.735,0.86),legend.box = "horizontal")


#################################################################
# adjustment comparison
#################################################################

# power
# BH
load("test_two_stage_constant_power.Rdata")
power_constant_bh_60=rowMeans(pvals[,,4])
power_constant_bh_80=rowMeans(pvals[,,5])
# BY
load("test_two_stage_constant_by_power.Rdata")
power_constant_by_60=rowMeans(pvals[,,4])
power_constant_by_80=rowMeans(pvals[,,5])
# q-value
load("test_two_stage_constant_q_power.Rdata")
power_constant_q_60=rowMeans(pvals[,,4])
power_constant_q_80=rowMeans(pvals[,,5])

powerres=data.frame(cbind(noise=seq(0.05,0.2,by=0.05),
                        power_constant_bh_60, power_constant_bh_80,
                        power_constant_by_60, power_constant_by_80,
                        power_constant_q_60, power_constant_q_80))

colnames(powerres)=c("noise","BH (C=60%)","BH (C=80%)","BY (C=60%)", "BY (C=80%)", "q-value (C=60%)","q-value (C=80%)")
powerdf=reshape2::melt(powerres, id.var = "noise")

powerdf$Method=rep(c("BH","BY","q-value"),each=8)
powerdf$Threshold=rep(rep(c("C=60%","C=80%"),each=4),3)

ggplot(powerdf,aes(x=noise,y=value,group=interaction(Method,Threshold))) +
  geom_point(size=2.5,aes(shape=Threshold,color=Method))+
  geom_line(size=1.05,aes(linetype=Threshold,color=Method)) +
  scale_colour_discrete(guide = guide_legend(override.aes = list(linetype = c(1,1,1),
                                                                 shape = c(NA,NA,NA)
  ) )) +
  labs(x="Noise Level",y="Power")+
  theme_bw()+
  theme(legend.position=c(0.765,0.86),legend.box = "horizontal")


# FPR
# BH
load("test_two_stage_constant_fpr.Rdata")
fpr_constant_bh_60=rowMeans(pvals[,,4])
fpr_constant_bh_80=rowMeans(pvals[,,5])
# BY
load("test_two_stage_constant_by_fpr.Rdata")
fpr_constant_by_60=rowMeans(pvals[,,4])
fpr_constant_by_80=rowMeans(pvals[,,5])
# q-value
load("test_two_stage_constant_q_fpr.Rdata")
fpr_constant_q_60=rowMeans(pvals[,,4])
fpr_constant_q_80=rowMeans(pvals[,,5])

fprres=data.frame(cbind(noise=seq(0.05,0.2,by=0.05),
                          fpr_constant_bh_60, fpr_constant_bh_80,
                          fpr_constant_by_60, fpr_constant_by_80,
                          fpr_constant_q_60, fpr_constant_q_80))

colnames(fprres)=c("noise","BH (C=60%)","BH (C=80%)","BY (C=60%)", "BY (C=80%)", "q-value (C=60%)","q-value (C=80%)")
fprdf=reshape2::melt(fprres, id.var = "noise")

fprdf$Method=rep(c("BH","BY","q-value"),each=8)
fprdf$Threshold=rep(rep(c("C=60%","C=80%"),each=4),3)

fprplot=ggplot(fprdf,aes(x=noise,y=value,group=interaction(Method,Threshold))) +
  geom_point(size=2.5,aes(shape=Threshold,color=Method))+
  geom_line(size=1.05,aes(linetype=Threshold,color=Method)) +
  scale_colour_discrete(guide = guide_legend(override.aes = list(linetype = c(1,1,1),
                                                                 shape = c(NA,NA,NA)) )) +
  labs(x="Noise Level",y="False Positive Rate")+
  theme_bw()+
  ylim(0,0.1)+
  theme(legend.position=c(0.765,0.86),legend.box = "horizontal")
