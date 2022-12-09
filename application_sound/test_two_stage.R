library(nonlinearTseries)
library(TDAstats)
library(tidyverse)
library(TDA)
library(seewave)

set.seed(1)

source("functions.R")

##########################################
# music data
##########################################

# flute
flute_pd <- list()
Flute_full_note <- read.csv("Flute A full time series.csv")
for (i in 1:60) {
  # sample time series
  ai <- sample(10000:37999,1,replace = FALSE, prob = NULL)
  bi <- ai + 1000
  sample_sequence_ai <- Flute_full_note[ai:bi,]
  
  # Takens embedding
  flute_A4_matrix <- data.matrix(sample_sequence_ai)
  tak <- buildTakens(flute_A4_matrix,2,3)
  
  # compute PH
  flute_pd[[i]] <- calculate_homology(tak,return_df = TRUE) 
}

# clarinet
clarinet_pd <- list()
Clarinet_full_note <- read.csv("Clarinet A full time series.csv")
for (i in 1:60) {
  ai2 <- sample(10000:95000,1,replace = FALSE, prob = NULL)
  bi2 <- ai2 + 1000
  sample_sequence_ai2 <- Clarinet_full_note[ai2:bi2,]
  
  # Takens embedding
  Clarinet_A4_matrix <- data.matrix(sample_sequence_ai2)
  tak2 <- buildTakens(Clarinet_A4_matrix,2,3)
  
  clarinet_pd[[i]] <- calculate_homology(tak2,return_df = TRUE) 
}

##########################################
# plots
##########################################

set.seed(9)

# sample time series
ai <- sample(10000:37999,1,replace = FALSE, prob = NULL)
bi <- ai + 1000
sample_sequence_ai <- Flute_full_note[ai:bi,]
sample_sequence_ai2 <- Clarinet_full_note[ai:bi,]


npoint=1001
plotdata=data.frame(times = 50*rep(ai:(ai+npoint-1)/1000000,2),
                    y=c(sample_sequence_ai2[1:npoint],sample_sequence_ai[1:npoint]),
                    Instrument=c(rep("Clarinet",npoint),rep("Flute",npoint)))
ggplot(plotdata, aes(x=times,y=y,color=Instrument,linetype=Instrument)) +
  geom_line()+
  labs(x="Second",y="")+
  scale_color_manual(values = c("blue","red"), 
                     labels = c( "Clarinet", "Flute"))+
  scale_linetype_manual(values = c(1,5), 
                        labels = c("Clarinet", "Flute"))+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         legend.position = "top", 
                         legend.background=element_rect(colour = "white"),
                         legend.title = element_text(size=15),
                         legend.text = element_text(size=15))

# spectrum
flute_spec=meanspec(sample_sequence_ai,44100,dB = 'max0')
clarinet_spec=meanspec(sample_sequence_ai2,44100,dB = 'max0')
specdata=data.frame(times = rep(flute_spec[,1],2),
                    y=c(clarinet_spec[,2],flute_spec[,2]),
                    Instrument=c(rep("Clarinet",nrow(flute_spec)),rep("Flute",nrow(flute_spec))))

ggplot(specdata, aes(x=times,y=y,color=Instrument,linetype=Instrument)) +
  geom_line()+
  labs(x="Frequency (kHz)",y="Amplitude (dB)")+
  scale_color_manual(values = c("blue","red"), 
                     labels = c( "Clarinet", "Flute"))+
  scale_linetype_manual(values = c(1,5), 
                        labels = c("Clarinet", "Flute"))+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "top", 
        legend.background=element_rect(colour = "white"),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

# flute
x<- data.matrix(sample_sequence_ai)
a<- buildTakens(x,2,3)
ggplot(data.frame(a,col="a"))+
  geom_point(aes(x=X1,y=X2,color=col))+
  labs(x=expression(Z[t]),y=expression(Z[t+3]))+
  scale_color_manual(values = c("red"))+
  theme_minimal()+
  theme(legend.position = "none")+
  ylim(-0.4,0.4)+
  xlim(-0.4,0.4)

# clarinet
x<- data.matrix(sample_sequence_ai2)
a<- buildTakens(x,2,3)
ggplot(data.frame(a,col="a"))+
  geom_point(aes(x=X1,y=X2,color=col))+
  labs(x=expression(Z[t]),y=expression(Z[t+3]))+
  scale_color_manual(values = c("blue"))+
  theme_minimal()+
  theme(legend.position = "none")+
  ylim(-0.4,0.4)+
  xlim(-0.4,0.4)

########################################################
# persistence image
########################################################
res=40
rangex=rangey=mm=c(0,0.25)
h=0.025

piflute_arctan = lapply(1:60, pimain, flute_pd,
                 rangex=rangex, rangey=rangey, 
                 nbins=res, h=h, weight="arctan")
piclarinet_arctan = lapply(1:60, pimain, clarinet_pd,
                    rangex=rangex, rangey=rangey,
                    nbins=res, h=h, weight="arctan")

piflute=array(NA,dim=c(res,res,60))
piclarinet=array(NA,dim=c(res,res,60))

for (pp in 1:60){
  piflute[,,pp]=piflute_arctan[[pp]] # flute
  piclarinet[,,pp]=piclarinet_arctan[[pp]] # clarinet
}

########################################################
# two-stage test
########################################################
cc=0.8

# between flute
group1 = piflute[,,21:40]
group2 = piflute[,,41:60]
testlist = testpi(group1,group2,res)
testfun(testlist,mm,res=res,cc=cc) # p-value

# between clarinet
group1 = piclarinet[,,21:40]
group2 = piclarinet[,,41:60]
testlist = testpi(group1,group2,res)
testfun(testlist,mm,res=res,cc=cc) # p-value

# between clarinet and flute
group1 = piflute[,,1:20]
group2 = piclarinet[,,21:40]
testlist = testpi(group1,group2,res)
testfun(testlist,mm,res=res,cc=cc) # p-value


########################################################
# plot
########################################################

sdoverall=testlist[[1]]
tpval=testlist[[2]]

nr=nrow(tpval)
nc=ncol(tpval)

cunit=(mm[2]-mm[1])/(2*res)
gridx=seq(mm[1]+cunit,mm[2]-cunit,length.out = res)
gridy=seq(mm[2]-cunit,mm[1]+cunit,length.out = res)

# x and y location
xgrid = 1:res
ygrid = 1:res
ind=expand.grid(xgrid,ygrid)
# p value
pval=as.vector(tpval)
# overall variance
oversd=as.vector(sdoverall)

# make data frame
df.pval=data.frame(idx=ind[,1],idy=ind[,2],pval,oversd)

# remove the upper triangle pixels
df.pval=df.pval %>%
  mutate(ind=((idx-idy)>=0))

df.pval.rm=df.pval[(df.pval$idx-df.pval$idy)>=0,]
npix=nrow(df.pval.rm)

# filter 
## sd
df.pval.rm$sdrank= (rank(df.pval.rm$oversd)/npix)
df.pval.sd=df.pval.rm[df.pval.rm$sdrank>cc,]
## sd
df.pval.sd$BH=p.adjust(df.pval.sd$pval,method=c("BH"))

df.new.sd=df.pval %>%
  select(idx,idy,ind) %>%
  left_join(df.pval.sd,by=c("idx"="idx","idy"="idy"))


final=df.new.sd %>%
  mutate(BH_new=case_when(
    BH<0.001 ~ " <0.001",
    BH>=0.001 & BH<0.005 ~ "0.001 to 0.005",
    BH>=0.005 & BH<0.01 ~ "0.005 to 0.01",
    BH<0.05 & BH>=0.01 ~ "0.01 to 0.05",
    (BH<0.1 & BH>=0.05) ~ "0.05 to 0.1",
    (BH<2 & BH>0.1) ~ "0.1 to 1",
    ind.x==TRUE ~ "filtered",
    TRUE ~ " ")
  )

df=data.frame(x=rep(gridx,each=res),y=rep(gridy,res),z=final$BH_new)
ggplot(df,aes(x,y,fill=z))+
  geom_raster()+
  xlim(range[1]-0.01,range[2]+0.01)+
  ylim(range[1]-0.01,range[2]+0.01)+
  scale_fill_manual(values=c("white","royalblue4","wheat"))+
  theme_minimal()+
  labs(x="",y="",fill="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = c(0.85,0.75))
