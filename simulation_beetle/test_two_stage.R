source("functions.R")

# import persistence diagrams
load("./data_persistence_diagram.Rdata")

##########################################################
# Persistence diagram and persistence images
##########################################################
# parameters for persistence images
res=40
rangex=rangey=mm=c(0,100)
h=3
nset=2000*2
ndat=2*nset

pilist_arctan = lapply(1:ndat,pimain,pdlist,
                       rangex=rangex,rangey=rangey,
                       nbins=res,h=h,weight="arctan")

pilist_linear = lapply(1:ndat,pimain,pdlist,
                       rangex=rangex,rangey=rangey,
                       nbins=res,h=h,weight="linear")

pilist_constant = lapply(1:ndat,pimain,pdlist,
                         rangex=rangex,rangey=rangey,
                         nbins=res,h=h,weight="constant")

pione_constant=pione_arctan=pione_linear=array(NA,dim=c(res,res,ndat))


for (pp in 1:ndat){
  pione_constant[,,pp]=pilist_constant[[pp]]
  pione_arctan[,,pp]=pilist_arctan[[pp]]
  pione_linear[,,pp]=pilist_linear[[pp]]
}

##########################################################
# Two-stage tests
##########################################################

################# arctan weight ##########################

# between stable and aperiodic regimes
cc=0.8
nrep=100
nset=20
ntotalset=4000
pvals=rep(NA,nrep)
for (rr in 1:nrep){
  group1 = pione_arctan[,,rr*(1:nset)]
  group2 = pione_arctan[,,rr*(1:nset)+ntotalset]
  testlist = ttestpi(group1,group2,res)
  pvals[rr]= testfunpi(testlist,rangex,rangey,cc,res)
}
sum(pvals<0.05)/nrep

# between aperiodic regimes
cc=0.8
nrep=100
nset=20
ntotalset=2000
pvals=rep(NA,nrep)
for (rr in 1:nrep){
  group1 = pione_arctan[,,rr*(1:nset)]
  group2 = pione_arctan[,,rr*(1:nset)+ntotalset]
  testlist = ttestpi(group1,group2,res)
  pvals[rr]= testfunpi(testlist,rangex,rangey,cc,res)
}
sum(pvals<0.05)/nrep

# between stable regimes
cc=0.8
nrep=100
nset=20
ntotalset=2000
pvals=rep(NA,nrep)
for (rr in 1:nrep){
  group1 = pione_arctan[,,4000+rr*(1:nset)]
  group2 = pione_arctan[,,4000+rr*(1:nset)+ntotalset]
  testlist = ttestpi(group1,group2,res)
  pvals[rr]= testfunpi(testlist,rangex,rangey,cc,res)
}
sum(pvals<0.05)/nrep

##############################################################
# plot
##############################################################
cc=0.8
nrep=100
nset=20
ntotalset=4000
pvals=rep(NA,nrep)
for (rr in 1:1){
  group1 = pione_arctan[,,rr*(1:nset)]
  group2 = pione_arctan[,,rr*(1:nset)+ntotalset]
}

testlist = ttestpi(group1,group2,res)
sdoverall=testlist[[1]]
tpval=testlist[[2]]

# x and y location
rangex=mm
rangey=mm
range=mm

cunit=(range[2]-range[1])/(2*res)
gridx=seq(range[1]+cunit,range[2]-cunit,length.out = res)
gridy=seq(range[2]-cunit,range[1]+cunit,length.out = res)

# x and y location
xgrid = 1:res
ygrid = 1:res
ind=expand.grid(xgrid,ygrid)
# p value
pval=as.vector(tpval)
# overall mean
overmean=as.vector(meanoverall)
# overall variance
oversd=as.vector(sdoverall)

# make data frame
df.pval=data.frame(idx=ind[,1],idy=ind[,2],pval,overmean,oversd)

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
  scale_fill_manual(values=c("white","royalblue4","royalblue3","royalblue","royalblue2","royalblue1","skyblue","wheat"))+
  theme_minimal()+
  labs(x="",y="",fill="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = c(0.86,0.80), legend.background=element_rect(colour = "white"))
