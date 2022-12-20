library(tidyverse)
library(TDA)
library(parallel)
library(gtools)

###############################################################
# Functions for persistence image
###############################################################

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

pers.image = function(pd, rangex, rangey, wgt, nbins, h){
  
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

###############################################################
# Functions for two-stage
###############################################################

ts_main_fpr = function(twopi,sig,npc,nset,range,res,alpha){
  nsig=length(sig)
  n1=npc*nset
  n2=npc*nset
  
  # take vectors
  twopix=list()
  for (ii in 1:nsig) {
    twopix[[ii]]=list()
    for (jj in 1:nset) {
      twopix[[ii]][[jj]]=array(NA, dim = c(res,res,npc))
      for (kk in 1:npc){
        twopix[[ii]][[jj]][,,kk]=twopi[[ii]][[(jj-1)*npc+kk]]
      }
    }
  }
  
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
  C=c(0,0.2,0.4,0.6,0.8)
  
  # Adjust p-values
  pvals=array(NA,dim=c(nsig,nset,length(C)))
  
  for (ii in 1:nsig) {
    for (jj in 1:nset) {
      xgrid = 1:res
      ygrid = 1:res
      ind=expand.grid(xgrid,ygrid)
      
      # make data frame
      df.pval=data.frame(idx=ind[,1],idy=ind[,2],pval=as.vector(tpval[ii,jj,,]),oversd=as.vector(sdoverall[ii,jj,,]))
      
      # remove the upper triangle pixels
      df.pval=df.pval %>%
        mutate(ind=((idx-idy)>=0))
      # filter 
      df.pval.rm=df.pval[(df.pval$idx-df.pval$idy)>=0,]
      npix=nrow(df.pval.rm)
      
      for (cc in 1:length(C)) {
        df.pval.rm$sdrank= (rank(df.pval.rm$oversd)/npix)
        df.pval.sd=df.pval.rm[df.pval.rm$sdrank>C[cc],]
        pvals[ii,jj,cc]=(min(p.adjust(df.pval.sd$pval,method=c("BH")))<alpha)*1
      }
    }
  }
  return(pvals)
}



ts_main_power = function(onepi,twopi,sig,npc,nset,range,res,alpha){
  nsig=length(sig)
  n1=npc*nset
  n2=npc*nset
  
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
  C=c(0,0.2,0.4,0.6,0.8)
  
  # Adjust p-values
  pvals=array(NA,dim=c(nsig,nset,length(C)))
  
  for (ii in 1:nsig) {
    for (jj in 1:nset) {
      # make data frame
      df.pval=data.frame(idx=ind[,1],idy=ind[,2],pval=tpval[ii,jj,],oversd=sdoverall[ii,jj,])
      
      # remove the upper triangle pixels
      df.pval=df.pval %>%
        mutate(ind=((idx-idy)>=0))
      # filter 
      df.pval.rm=df.pval[(df.pval$idx-df.pval$idy)>=0,]
      npix=nrow(df.pval.rm)
      
      for (cc in 1:length(C)) {
        df.pval.rm$sdrank= (rank(df.pval.rm$oversd)/npix)
        df.pval.sd=df.pval.rm[df.pval.rm$sdrank>C[cc],]
        pvals[ii,jj,cc]=(min(p.adjust(df.pval.sd$pval,method=c("BH")))<alpha)*1
      }
    }
  }
  return(pvals)
}

###############################################################
# Functions for Robins and Turner
###############################################################

# test Robins and Turner
parrttest=function(kk,mergeddat,totallabel,shufflelabel,npc,nset) {
  # 45 combinations for 10 npc
  pd1=mergeddat[shufflelabel[kk,]]
  pd2=mergeddat[totallabel[-shufflelabel[kk,]]]
  
  jointlabel=gtools::combinations(npc,2)
  perf=rep(NA,nrow(jointlabel))
  for( ii in 1:nrow(jointlabel)) {
    perf[ii]=TDA::wasserstein(pd1[[jointlabel[ii,1]]],pd1[[jointlabel[ii,2]]])+TDA::wasserstein(pd2[[jointlabel[ii,1]]],pd2[[jointlabel[ii,2]]])
  }
  return(mean(perf))
}

# test Robins and Turner
rttest=function(pd1,pd2,npc,nset) {
  jointlabel=combinations(npc,2)
  perf=rep(NA,nrow(jointlabel))
  for (ii in 1:nrow(jointlabel)) {
    perf[ii]=wasserstein(pd1[[jointlabel[ii,1]]],pd1[[jointlabel[ii,2]]])+
      wasserstein(pd2[[jointlabel[ii,1]]],pd2[[jointlabel[ii,2]]])
  }
  return(mean(perf))
}



###############################################################
# Functions for persistence landscape
# Modified from https://people.clas.ufl.edu/peterbubenik/files/intro_tda.txt
###############################################################

landscape.matrix.from.list <- function(PL.list){
  n <- length(PL.list)
  m <- ncol(PL.list[[1]])
  max.depth <- integer(n)
  for (i in 1:n)
    max.depth[i] <- nrow(PL.list[[i]])
  K <- max(max.depth)
  PL.matrix <- matrix(0, nrow = n, ncol = m*K)
  for (i in 1:n)
    for (j in 1:max.depth[i])
      PL.matrix[i,(1+(j-1)*m):(j*m)] <- PL.list[[i]][j,]
  return(PL.matrix)
}

euclidean.distance <- function(a, b) sqrt(sum((a - b)^2))

permutation.test <- function(M1 ,M2, num.repeats = 5000){
  # append zeros if necessary so that the matrices have the same number of columns
  num.columns <- max(ncol(M1),ncol(M2))
  M1 <- cbind(M1, matrix(0,nrow=nrow(M1),ncol=num.columns-ncol(M1)))
  M2 <- cbind(M2, matrix(0,nrow=nrow(M2),ncol=num.columns-ncol(M2)))
  t.obs <- euclidean.distance(colMeans(M1),colMeans(M2))
  k <- dim(M1)[1]
  M <- rbind(M1,M2)
  n <- dim(M)[1]
  count <- 0
  for (i in 1:num.repeats){
    permutation <- sample(1:n)
    t <- euclidean.distance(colMeans(M[permutation[1:k],]),colMeans(M[permutation[(k+1):n],]))
    if (t >= t.obs)
      count <- count + 1
  }
  return(count/num.repeats)  
}
