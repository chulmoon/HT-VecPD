source("functions.R")

totalset=27

# persistence diagrams
fpath="./data/"
flist=list.files(fpath)
totalpd=list()
for (ii in 1:length(flist)){
  fname=paste(fpath,flist[ii],sep="")
  pd=read.table(fname)
  pd=pd %>%
    rename(dimension=V1,Birth=V2,Death=V3)
  totalpd[[ii]]=as.matrix(pd)
}

load("permulabel.Rdata")

cl=makeCluster(detectCores()-1)

######################################################
# F42B and F42C
######################################################
group1=totalpd[1:27]
group2=totalpd[(1*27+1):(2*27)]

# dimension zero
permun=500
dim=0
# initial joint loss
basef=rttest(group1,group2,totalset,dim)
# permutation test
mergeddat=append(group1,group2)
kk=as.matrix(1:permun)
perf_pd=parApply(cl,kk,1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(perf_pd<basef)/permun # p-value

# dimension one
permun=100
dim=1
# initial joint loss
basef=rttest(group1,group2,totalset,dim)
# permutation test
mergeddat=append(group1,group2)
kk=as.matrix(1:permun)
perf_pd=parApply(cl,kk,1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(perf_pd<basef)/permun # p-value

# dimension two
permun=500
dim=2
# initial joint loss
basef=rttest(group1,group2,totalset,dim)
# permutation test
mergeddat=append(group1,group2)
kk=as.matrix(1:permun)
perf_pd=parApply(cl,kk,1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(perf_pd<basef)/permun # p-value

######################################################
# LV60A and LV60C
######################################################
group1=totalpd[(2*27+1):(3*27)]
group2=totalpd[(3*27+1):(4*27)]

# dimension zero
permun=500
dim=0
# initial joint loss
basef=rttest(group1,group2,totalset,dim)
# permutation test
mergeddat=append(group1,group2)
kk=as.matrix(1:permun)
perf_pd=parApply(cl,kk,1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(perf_pd<basef)/permun # p-value

# dimension one
permun=100
dim=1
# initial joint loss
basef=rttest(group1,group2,totalset,dim)
# permutation test
mergeddat=append(group1,group2)
kk=as.matrix(1:permun)
perf_pd=parApply(cl,kk,1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(perf_pd<basef)/permun # p-value

# dimension two
permun=500
dim=2
# initial joint loss
basef=rttest(group1,group2,totalset,dim)
# permutation test
mergeddat=append(group1,group2)
kk=as.matrix(1:permun)
perf_pd=parApply(cl,kk,1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(perf_pd<basef)/permun # p-value

######################################################
# F42B and LV60A
######################################################
group1=totalpd[1:27]
group2=totalpd[(2*27+1):(3*27)]

# dimension zero
permun=500
dim=0
# initial joint loss
basef=rttest(group1,group2,totalset,dim)
# permutation test
mergeddat=append(group1,group2)
kk=as.matrix(1:permun)
perf_pd=parApply(cl,kk,1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(perf_pd<basef)/permun # p-value

# dimension one
permun=100
dim=1
# initial joint loss
basef=rttest(group1,group2,totalset,dim)
# permutation test
mergeddat=append(group1,group2)
kk=as.matrix(1:permun)
perf_pd=parApply(cl,kk,1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(perf_pd<basef)/permun # p-value

# dimension two
permun=500
dim=2
# initial joint loss
basef=rttest(group1,group2,totalset,dim)
# permutation test
mergeddat=append(group1,group2)
kk=as.matrix(1:permun)
perf_pd=parApply(cl,kk,1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(perf_pd<basef)/permun # p-value

stopCluster(cl)