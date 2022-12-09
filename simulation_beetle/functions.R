library(tidyverse)
library(TDA)
library(parallel)
library(dtw)

###############################################################
# Functions for two-stage
###############################################################

# arctangent weight
arctanw=function(pd,Cdist,ddist) {
  arctan=atan(Cdist*(pd$death)^ddist)
  return(arctan)
}

# linear weight 
linearw=function(pd) {
  linear=pd$death
  return(linear)
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

# function for running persistence image
pimain = function(ii,pdlist, rangex, rangey, nbins, h, weight){
  pd=pdlist[[ii]]
  pddf=data.frame(matrix(as.vector(pd),nrow(pd),ncol(pd)))
  colnames(pddf)=c("dimension","birth","death")
  pd1 = pddf %>%  
    dplyr::filter(dimension==1) %>%
    dplyr::select(-dimension)
  pd1$death = pd1$death - pd1$birth
  
  if (weight=="arctan") wgt = arctanw(pd1,0.5,0.5)
  else if (weight=="linear") wgt = linearw(pd1)
  else if (weight=="constant") wgt = rep(1,nrow(pd1))
  else stop("Select either arctan, linear, or constant weights.")
  
  pi1=pers_image(pd=pd1, rangex=rangex, rangey=rangey, wgt=wgt, nbins=nbins, h=h)
  return(pi1)
}

# compute p-values for t-test
ttestpi = function(group1,group2,res) {
  
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

# compute p-values for ANOVA
anovapi = function(group1,group2,group3,res) {
  
  sdoverall=array(NA,dim=c(res,res))
  anovapval=array(NA,dim=c(res,res))
  
  for (xx in 1:res) {
    for (yy in 1:res) {
      sdoverall[xx,yy]=sd(c(group1[xx,yy,],group2[xx,yy,],group3[xx,yy,]))
      data_df=data.frame(val=c(group1[xx,yy,],group2[xx,yy,],group3[xx,yy,]),
                         group=rep(c(1,2,3),each=length(group1[xx,yy,])))
      anova_one_way = aov(val~group, data = data_df)
      anova_summary=summary(anova_one_way)
      anovapval[xx,yy]=anova_summary[[1]][["Pr(>F)"]][1]
    }
  }
  
  return(list(sdoverall,anovapval))
}

testfunpi = function(testlist,rangex,rangey,cc,res) {
  sdoverall=testlist[[1]]
  tpval=testlist[[2]]
  
  # x and y location
  cunit=(rangex[2]-rangex[1])/(2*res)
  gridx=seq(rangex[1]+cunit,rangex[2]-cunit,length.out = res)
  gridy=seq(rangey[2]-cunit,rangey[1]+cunit,length.out = res)
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
  
  # BH correction
  ## sd
  df.pval.sd$BH=p.adjust(df.pval.sd$pval,method=c("BH"))
  
  df.new.sd=df.pval %>%
    select(idx,idy,ind) %>%
    left_join(df.pval.sd,by=c("idx"="idx","idy"="idy"))
  
  return(min(df.new.sd$BH,na.rm=T))
}


test_results = function(pilist0,pilist1,mmzero,mmone,cc=0.6){
  
  # parameters
  res0=mmzero[2]-mmzero[1]+1
  res1=mmone[2]-mmone[1]+1
  
  # dimension zero
  pi0group1 = pi0group2 = pi0group3 = pi0group4 = array(
    dim=c(mmzero[2]-mmzero[1]+1,mmzero[2]-mmzero[1]+1,50))
  # dimension one
  pi1group1 = pi1group2 = pi1group3 = pi1group4 = array(
    dim=c(mmone[2]-mmone[1]+1,mmone[2]-mmone[1]+1,50) )
  
  for (ii in 1:50){
    pi0group1[,,ii] = pilist0[[ii]] # group 1
    pi0group2[,,ii] = pilist0[[ii+50]] # group 2
    pi0group3[,,ii] = pilist0[[ii+100]] # group 3
    pi0group4[,,ii] = pilist0[[ii+150]] # group 4
    
    pi1group1[,,ii] = pilist1[[ii]]
    pi1group2[,,ii] = pilist1[[ii+50]]
    pi1group3[,,ii] = pilist1[[ii+100]]
    pi1group4[,,ii] = pilist1[[ii+150]]
  }
  
  # scenario 1
  ## p values, unfiltered
  sc1_arctan_dim0=ttestpi(pi0group1,pi0group2,res0) # dimension zero
  sc1_arctan_dim1=ttestpi(pi1group1,pi1group2,res1) # dimension one
  ## p values, filtered and adjusted
  pval10=testfun(sc1_arctan_dim0,mmzero,mmzero+abs(min(mmzero)),cc)
  pval11=testfun(sc1_arctan_dim1,mmone,mmone+abs(min(mmone)),cc)
  
  # scenario 2
  ## p values, unfiltered
  sc2_arctan_dim0=ttestpi(pi0group1,pi0group3,res0) # dimension zero
  sc2_arctan_dim1=ttestpi(pi1group1,pi1group3,res1) # dimension one
  ## p values, filtered and adjusted
  pval20=testfun(sc2_arctan_dim0,mmzero,mmzero+abs(min(mmzero)),cc)
  pval21=testfun(sc2_arctan_dim1,mmone,mmone+abs(min(mmone)),cc)
  
  # scenario 3
  ## p values, unfiltered
  sc3_arctan_dim0=ttestpi(pi0group1,pi0group4,res0) # dimension zero
  sc3_arctan_dim1=ttestpi(pi1group1,pi1group4,res1) # dimension one
  ## p values, filtered and adjusted
  pval30=testfun(sc3_arctan_dim0,mmzero,mmzero+abs(min(mmzero)),cc)
  pval31=testfun(sc3_arctan_dim1,mmone,mmone+abs(min(mmone)),cc)
  
  # scenario 4
  ## p values, unfiltered
  sc4_arctan_dim0=anovapi(pi0group1,pi0group3,pi0group4,res0) # dimension zero
  sc4_arctan_dim1=anovapi(pi1group1,pi1group3,pi1group4,res1) # dimension one
  ## p values, filtered and adjusted
  pval40=testfun(sc4_arctan_dim0,mmzero,mmzero+abs(min(mmzero)),cc)
  pval41=testfun(sc4_arctan_dim1,mmone,mmone+abs(min(mmone)),cc)
  
  return(c(pval10,pval11,pval20,pval21,pval30,pval31,pval40,pval41))
}

###############################################################
# Functions for Robins and Turner
###############################################################

# test Robins and Turner
rttest=function(pd1,pd2,nset) {
  jointlabel=t(combn(1:nset,2))
  perf=rep(NA,nrow(jointlabel))
  for (ii in 1:nrow(jointlabel)) {
    perf[ii]=TDA::wasserstein(pdlist1[[jointlabel[ii,1]]],pdlist1[[jointlabel[ii,2]]])+
      TDA::wasserstein(pdlist2[[jointlabel[ii,1]]],pdlist2[[jointlabel[ii,2]]])
  }
  return(mean(perf))
}

# test Robins and Turner
parrttest=function(kk,pdlist,rr,nset,ntotalset) {
  pdlist1 = list()
  pdlist2 = list()
  shufflelabel=sample(1:(2*nset),(2*nset))
  for (ss in 1:nset){
    pdlist1[[ss]]=pdlist[[ (rr-1)*nset + shufflelabel[ss] ]]
    pdlist2[[ss]]=pdlist[[ (rr-1)*nset + ntotalset + shufflelabel[ss] ]]
  }
  jointlabel=t(combn(1:nset,2))
  perf=rep(NA,nrow(jointlabel))
  for( ii in 1:nrow(jointlabel)) {
    perf[ii]=TDA::wasserstein(pdlist1[[jointlabel[ii,1]]],pdlist1[[jointlabel[ii,2]]])+
      TDA::wasserstein(pdlist2[[jointlabel[ii,1]]],pdlist2[[jointlabel[ii,2]]])
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


###############################################################
# Functions for DTW
###############################################################


# baseline of the permutation test
baseline_dtw=function(ts1,ts2,nts) {
  jointlabel=t(combn(1:nts,2))
  perf=rep(NA,nrow(jointlabel))
  for (ii in 1:nrow(jointlabel)) {
    perf[ii]=dtw::dtw(ts1[jointlabel[ii,1],],ts1[jointlabel[ii,2],])$distance+
      dtw::dtw(ts2[jointlabel[ii,1],],ts2[jointlabel[ii,2],])$distance
  }
  return(mean(perf))
}

# permutation test
permtest_dtw=function(kk,tsdata,rr,nts,nset,ntotalset) {
  shufflelabel=sample(1:(2*nset),(2*nset))
  ts1 = tsdata[shufflelabel[1:nset],]
  ts2 = tsdata[shufflelabel[(nset+1):(2*nset)],]
  jointlabel=t(combn(1:nts,2))
  perf=rep(NA,nrow(jointlabel))
  for( ii in 1:nrow(jointlabel)) {
    perf[ii]=dtw::dtw(ts1[jointlabel[ii,1],],ts1[jointlabel[ii,2],])$distance+
      dtw::dtw(ts2[jointlabel[ii,1],],ts2[jointlabel[ii,2],])$distance
  }
  return(mean(perf))
}
