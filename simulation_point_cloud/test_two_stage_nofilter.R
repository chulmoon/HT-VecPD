library(tidyverse)

# constant
load("data_pi_onecircle_constant.Rdata")
load("data_pi_twocircle_constant.Rdata")

onepi=onepi_constant
twopi=twopi_constant

# settings
sig=seq(0.05,0.20,length.out = 4)
npc=10
nset=500
nsig=length(sig)

n1=npc*nset
n2=npc*nset
range=c(0,2)
res=40
alpha=0.05

#######################################################################
# power
#######################################################################

# take vectors
onepix=list()
twopix=list()
for (ii in 1:nsig) {
  onepix[[ii]]=list()
  twopix[[ii]]=list()
  for (jj in 1:nset) {
    onepix[[ii]][[jj]]=array(NA, dim = c(res,res,npc))
    twopix[[ii]][[jj]]=array(NA, dim = c(res,res,npc))
    for (kk in 1:npc){
      onepix[[ii]][[jj]][,,kk]=onepi[[ii]][[(jj-1)*npc+kk]]
      twopix[[ii]][[jj]][,,kk]=twopi[[ii]][[(jj-1)*npc+kk]]
    }
  }
}


filterfun = function(ind,ii,jj,onepix,twopix,alpha){
  xx=ind[1]
  yy=ind[2]
  sdres=sd( c(onepix[[ii]][[jj]][xx,yy,],twopix[[ii]][[jj]][xx,yy,]) )
  tpvalres=t.test(onepix[[ii]][[jj]][xx,yy,],twopix[[ii]][[jj]][xx,yy,], 
                  alternative=c("two.sided"), conf.level = (1-alpha))$p.value
  return(c(sdres, tpvalres))
}

xgrid = 1:res
ygrid = 1:res
ind=expand.grid(xgrid,ygrid)

sdoverall=array(NA,dim=c(nsig,nset,res*res))
tpval=array(NA,dim=c(nsig,nset,res*res))

# filter using overall sd
for (ii in 1:nsig) {
  for (jj in 1:nset) {
    if(jj%%500==0) print(jj)
    res=t(apply(ind,1,filterfun,ii,jj,onepix,twopix,alpha))
    sdoverall[ii,jj,]=res[,1]
    tpval[ii,jj,]=res[,2]
  }
}


# filtering thresholds
C=c(0)
# create BH correction and BY correction
pvals=array(NA,dim=c(nsig,nset,length(C)))

for (ii in 1:nsig) {
  for (jj in 1:nset) {
    # make data frame
    df.pval=data.frame(idx=ind[,1],idy=ind[,2],pval=tpval[ii,jj,],oversd=sdoverall[ii,jj,])
    
    # filter 
    df.pval.rm=df.pval
    npix=nrow(df.pval.rm)
    
    for (cc in 1:length(C)) {
      df.pval.rm$sdrank= (rank(df.pval.rm$oversd)/npix)
      df.pval.sd=df.pval.rm[df.pval.rm$sdrank>C[cc],]
      pvals[ii,jj,cc]=(min(p.adjust(df.pval.sd$pval,method=c("BH")))<alpha)*1
    }
  }
}

save(pvals,file="test_two_stage_constant_nofilter_power.Rdata")


#######################################################################
# fpr
#######################################################################

# filter using overall sd
sdoverall=array(NA,dim=c(nsig,nset,res,res))
tpval=array(NA,dim=c(nsig,nset,res,res))

# combinations
cbn=combn(1:nset,2)

for (ii in 1:nsig) {
  ind=t(cbn[,sample(1:ncol(cbn),nset)])
  for (jj in 1:nset) {
    for (xx in 1:res) {
      for (yy in 1:res) {
        sdoverall[ii,jj,xx,yy]=sd( c(twopix[[ii]][[ind[jj,1]]][xx,yy,],twopix[[ii]][[ind[jj,2]]][xx,yy,]) )
        
        tpval[ii,jj,xx,yy]=t.test(twopix[[ii]][[ind[jj,1]]][xx,yy,],twopix[[ii]][[ind[jj,2]]][xx,yy,], 
                                  alternative=c("two.sided"), conf.level = (1-alpha))$p.value
      }
    }
  }
}

# filtering thresholds
C=c(0)
# create BH correction and BY correction
pvals=array(NA,dim=c(nsig,nset,length(C)))

for (ii in 1:nsig) {
  for (jj in 1:nset) {
    nr=nrow(tpval[ii,jj,,])
    nc=ncol(tpval[ii,jj,,])
    # x and y location
    idx=rep(1:nr,each=nc)
    idy=rep(1:nc,nr)
    
    # p value
    pval=as.vector(tpval[ii,jj,,])
    # overall sd
    oversd=as.vector(sdoverall[ii,jj,,])
    
    # make data frame
    df.pval=data.frame(idx,idy,pval,oversd)
    
    # remove the upper triangle pixels
    # filter 
    df.pval.rm=df.pval
    npix=nrow(df.pval.rm)
    
    for (cc in 1:length(C)) {
      df.pval.rm$sdrank= (rank(df.pval.rm$oversd)/npix)
      df.pval.sd=df.pval.rm[df.pval.rm$sdrank>C[cc],]
      pvals[ii,jj,cc]=(min(p.adjust(df.pval.sd$pval,method=c("BH")))<alpha)*1
    }
  }
}

save(pvals,file="test_two_stage_constant_nofilter_fpr.Rdata")
