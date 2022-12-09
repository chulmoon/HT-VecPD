source("functions.R")

# import persistence diagrams
load("./data_persistence_diagram.Rdata")

set.seed(75)

cl = makeCluster(detectCores()-1)

# create persistence landscape
PL.list=list()
for (ii in 1:length(pdlist)){
  pd=pdlist[[ii]]
  PL.list[[ii]]=t(landscape(pd,dimension=1,KK=1:50,tseq=seq(0,100,length=500)))
}

np=500 # number of permutation

#########################################################
# between aperiodic and stable regimes
#########################################################
nrep=100 # number of sets
nset=20 # number of data in one set
ntotalset=4000
pvals=rep(NA,nrep)
for (rr in 1:nrep){
  PL.list1 = list()
  PL.list2 = list()
  for (ss in 1:nset){
    PL.list1[[ss]]=PL.list[[(rr-1)*nset+ss]]
    PL.list2[[ss]]=PL.list[[(rr-1)*nset+ss+ntotalset]]
  }
  PL.mat1 = landscape.matrix.from.list(PL.list1)
  PL.mat2 = landscape.matrix.from.list(PL.list2)
  pvals[rr]=permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)
}
sum(pvals<0.05)/nrep # power

#########################################################
# between aperiodic regimes
#########################################################
nrep=100
nset=20
ntotalset=2000
pvals=rep(NA,nrep)
for (rr in 1:nrep){
  PL.list1 = list()
  PL.list2 = list()
  for (ss in 1:nset){
    PL.list1[[ss]]=PL.list[[(rr-1)*nset+ss]]
    PL.list2[[ss]]=PL.list[[(rr-1)*nset+ss+ntotalset]]
  }
  PL.mat1 = landscape.matrix.from.list(PL.list1)
  PL.mat2 = landscape.matrix.from.list(PL.list2)
  pvals[rr]=permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)
}
sum(pvals<0.05)/nrep # false positive rate

#########################################################
# between stable regimes
#########################################################
nrep=100
nset=20
ntotalset=2000
pvals=rep(NA,nrep)
for (rr in 1:nrep){
  PL.list1 = list()
  PL.list2 = list()
  for (ss in 1:nset){
    PL.list1[[ss]]=PL.list[[(rr-1)*nset+ss+4000]]
    PL.list2[[ss]]=PL.list[[(rr-1)*nset+ss+4000+ntotalset]]
  }
  PL.mat1 = landscape.matrix.from.list(PL.list1)
  PL.mat2 = landscape.matrix.from.list(PL.list2)
  pvals[rr]=permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)
}
sum(pvals<0.05)/nrep # false positive rate

stopCluster(cl)