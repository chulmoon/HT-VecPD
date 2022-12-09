library(TDA)
library(pracma)
library(mvtnorm)
library(tidyverse)
library(RColorBrewer)
library(gtools)
library(parallel)

###############################################################
# Functions for generating persistence image
###############################################################

bpdist=function(pd){
  pd=cbind(pd[,1],pd[,2]-pd[,1])
  return(pd)
}

pwgk=function(pd,C,d) {
  pwgk=atan(C*(pd[,2])^d)
  return(pwgk)
}

gaussian.bump=function(pd,gridbirth,griddeath,sig,weight,cluster){
  grid=meshgrid(gridbirth,griddeath)
  gridx=as.vector(grid$X)
  gridy=as.vector(grid$Y)
  ZZ=array(NA,dim=c(res,res,dim(pd)[1]))
  for (ii in 1:dim(pd)[1]) {
    AA=weight[ii]*parApply(cluster,cbind(gridx,gridy),1,pmvnorm,lower=-Inf,
                           mean=c(pd[ii,1],pd[ii,2]),sigma=sig*diag(2))
    AA=matrix(AA,length(gridbirth),length(griddeath))
    ZZ[,,ii]=(-AA[-1,-1]-AA[-length(gridbirth),-length(griddeath)]
              +AA[-1, -length(griddeath)]+AA[-length(gridbirth),-1])
  }
  return(apply(ZZ,c(1,2),sum))
}

makePI=function(pd,res,range,weight,sig=NULL,cluster) {
  pd=bpdist(pd)
  if (is.null(sig)==T) sig=((range[2]-range[1])/(res+1))*1.5
  gridbirth=seq(range[1],range[2],length.out=res+1)
  griddeath=seq(range[2]-range[1],0,length.out=res+1)
  weightres=do.call(weight,list(pd,0.5,0.5))
  return(gaussian.bump(pd,gridbirth,griddeath,sig,weightres,cluster))
}


###############################################################
# Functions for two-stage
###############################################################

# compute p-values for t-test
ttestpi = function(group1,group2,res,cc=0.5) {
  # select dataset for vectorized PD
  sdoverall=array(NA,dim=c(res,res))
  tpval=array(NA,dim=c(res,res))
  
  for (xx in 1:res) {
    for (yy in 1:res) {
      sdoverall[xx,yy]=sd(c(group1[xx,yy,],group2[xx,yy,]))
      tpval[xx,yy]=t.test(group1[xx,yy,],group2[xx,yy,], 
                          alternative=c("two.sided"))$p.value
    }
  }

  # x and y location
  idx=rep(1:res,res)
  idy=rep(1:res,each=res)
  # p value
  pval=as.vector(tpval)
  # overall variance
  oversd=as.vector(sdoverall)
  
  # make data frame
  df.pval=data.frame(idx,idy,pval,oversd)
  
  # remove the upper triangle pixels
  df.pval.rm=df.pval[(idx-idy)>=0,]
  npix=nrow(df.pval.rm)
  # filter 
  ## sd
  df.pval.rm$sdrank= (rank(df.pval.rm$oversd)/npix)
  df.pval.sd=df.pval.rm[df.pval.rm$sdrank>cc,]
  
  # BH correction
  ## sd
  df.pval.sd$BH=p.adjust(df.pval.sd$pval,method=c("BH"))
  
  df.new.sd=df.pval %>%
    select(idx,idy) %>%
    left_join(df.pval.sd,by=c("idx"="idx","idy"="idy"))
  
  return(min(df.new.sd$BH,na.rm=T))
}



###############################################################
# Functions for Robins and Turner
###############################################################

# parallel version of Robins and Turner
parrttest=function(kk,mergeddat,shufflelabel,totalset,dim) {
  print(kk)
  totallabel=1:(2*totalset)
  # shuffle the labels
  pd1=mergeddat[as.numeric(shufflelabel[kk,])]
  pd2=mergeddat[totallabel[-as.numeric(shufflelabel[kk,])]]
  
  # do the permutation tests
  jointlabel=gtools::combinations(totalset,2)
  perf=rep(0,nrow(jointlabel))
  jointloss=function(ii,pd1,pd2,dim,jointlabel){
    TDA::wasserstein(pd1[[jointlabel[ii,1]]],pd1[[jointlabel[ii,2]]],dimension=dim)+TDA::wasserstein(pd2[[jointlabel[ii,1]]],pd2[[jointlabel[ii,2]]],dimension=dim)
  }
  perf=purrr::map_dbl(1:nrow(jointlabel),jointloss,pd1,pd2,dim,jointlabel)
  return(mean(perf))
}

# for a single reference test: test Robins and Turner
rttest=function(pd1,pd2,totalset,dim) {
  jointlabel=gtools::combinations(totalset,2)
  perf=rep(0,nrow(jointlabel))
  jointloss=function(ii,pd1,pd2,dim,jointlabel){
    TDA::wasserstein(pd1[[jointlabel[ii,1]]],pd1[[jointlabel[ii,2]]],dimension=dim)+
      TDA::wasserstein(pd2[[jointlabel[ii,1]]],pd2[[jointlabel[ii,2]]],dimension=dim)
  }
  perf=purrr::map_dbl(1:nrow(jointlabel),jointloss,pd1,pd2,dim,jointlabel)
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

testfunpl = function(ii,M,n,k){
  euclidean.distance <- function(a, b) sqrt(sum((a - b)^2))
  permutation <- sample(1:n)
  t <- euclidean.distance(colMeans(M[permutation[1:k],]),colMeans(M[permutation[(k+1):n],]))
  return(t)
}

permutation.test <- function(M1 ,M2, num.repeats = 5000, cl){
  euclidean.distance <- function(a, b) sqrt(sum((a - b)^2))
  # append zeros if necessary so that the matrices have the same number of columns
  num.columns <- max(ncol(M1),ncol(M2))
  M1 <- cbind(M1, matrix(0,nrow=nrow(M1),ncol=num.columns-ncol(M1)))
  M2 <- cbind(M2, matrix(0,nrow=nrow(M2),ncol=num.columns-ncol(M2)))
  t.obs <- euclidean.distance(colMeans(M1),colMeans(M2))
  k <- dim(M1)[1]
  M <- rbind(M1,M2)
  n <- dim(M)[1]
  
  tvec=parApply(cl,as.matrix(1:num.repeats),1,testfun,M,n,k)
  return(sum(tvec>=t.obs)/num.repeats)  
}