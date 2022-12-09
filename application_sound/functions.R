library(nonlinearTseries)
library(TDAstats)
library(tidyverse)
library(TDA)

# arctan weight
arctanw=function(pd,Cdist,ddist) {
  arctan=atan(Cdist*(pd$death)^ddist)
  return(arctan)
}

# linear weight 
linearw=function(pd) {
  linear=pd$death
  return(linear)
}

pimain = function(ii,pdlist, rangex, rangey, nbins, h, weight){
  pd=pdlist[[ii]]
  pd1 = pd %>%  
    dplyr::filter(dimension==1) %>%
    dplyr::select(-dimension)
  pd1$death = pd1$death - pd1$birth
  
  if (weight=="arctan") wgt = arctanw(pd1,0.5,0.5)
  else if (weight=="linear") wgt = linearw(pd1)
  else if (weight=="constant") wgt = rep(1,nrow(pd1))
  else stop("Select either arctan, linear, or constant weights.")
  
  pi1=pers_image(pd=pd1, rangex=rangex, rangey=rangey, wgt=wgt, nbins=res, h=h)
  
  return(pi1)
}

# create persistence image
pers_image = function(pd, rangex, rangey, wgt, nbins, h){
  
  x_grid = seq(rangex[1], rangex[2], length.out = nbins+1)
  x_lower = x_grid[-(nbins+1)]
  x_upper = x_grid[-1]
  
  y_grid = seq(rangey[2], rangey[1], length.out = nbins+1)
  y_lower = y_grid[-1]
  y_upper = y_grid[-(nbins+1)]
  
  PSurface = function(point_weight){
    
    x = point_weight[1]
    y = point_weight[2]
    weight = point_weight[3]
    
    out1 = pnorm(x_upper, mean = x, sd = h) - pnorm(x_lower, mean = x, sd = h)
    out2 = pnorm(y_upper, mean = y, sd = h) - pnorm(y_lower, mean = y, sd = h)
    
    return(out2 %o% out1 * weight)
    
  }
  
  point_weight = cbind(pd,wgt)
  Psurf_mat = apply(point_weight, 1, PSurface)
  out = apply(Psurf_mat, 1, sum)
  
  return(matrix(out, nrow = nbins))
  
}


# compute filter statistics and p-values
testpi = function(group1,group2,res) {
  sdoverall=array(NA,dim=c(res,res))
  tpval=array(NA,dim=c(res,res))
  
  for (xx in 1:res) {
    for (yy in 1:res) {
      sdoverall[xx,yy]=sd(c(group1[xx,yy,],group2[xx,yy,]))
      tpval[xx,yy]=t.test(group1[xx,yy,],group2[xx,yy,], 
                          alternative=c("two.sided"))$p.value
    }
  }
  return(list(sdoverall,tpval))
}

# filtering and multiple testing adjustment
testfun = function(testlist,mm,res,cc=0.8) {
  sdoverall=testlist[[1]]
  tpval=testlist[[2]]
  
  nr=nrow(tpval)
  nc=ncol(tpval)
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
  
  return(min(df.new.sd$BH,na.rm=T))
}

