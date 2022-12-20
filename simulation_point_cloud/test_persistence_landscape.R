source("functions.R")

load("data_onecircle.Rdata")
load("data_twocircle.Rdata")

############################################################
# Power
############################################################
set.seed(28)

sig=seq(0.05,0.20,length.out = 4) # noise
nsig=length(sig)

npc=10 # number of point cloud in one set
nset=500 # number of sets
testn=1000 # N_p

perz=matrix(NA,nsig,nset)

for (ii in 1:nsig) {
  print(ii)
  for (jj in 1:nset) {
    PL.list1=list()
    PL.list2=list()
    for (kk in 1:npc){
      PL.list1[[kk]]=t(landscape(onepd[[ii]][[npc*(jj-1)+kk]],dimension=1,KK=1:10,tseq=seq(0,2,length=100)))
      PL.list2[[kk]]=t(landscape(twopd[[ii]][[npc*(jj-1)+kk]],dimension=1,KK=1:10,tseq=seq(0,2,length=100)))
    }
    PL.mat1 = landscape.matrix.from.list(PL.list1)
    PL.mat2 = landscape.matrix.from.list(PL.list2)
    perz[ii,jj]=permutation.test(PL.mat1,PL.mat2,num.repeats=testn)
  }
}

save(perz, file="test_persistence_landscape_power.Rdata")

############################################################
# False positive rate
############################################################
# combinations
cbn=combn(1:nset,2)
perz=matrix(NA,nsig,nset)

for (ii in 1:nsig) {
  ind=t(cbn[,sample(1:ncol(cbn),nset)])
  for (jj in 1:nset) {
    PL.list1=list()
    PL.list2=list()
    for (kk in 1:npc){
      PL.list1[[kk]]=t(landscape(twopd[[ii]][[ npc*(ind[jj,1]-1)+kk ]],dimension=1,KK=1:10))
      PL.list2[[kk]]=t(landscape(twopd[[ii]][[ npc*(ind[jj,2]-1)+kk ]],dimension=1,KK=1:10))
    }
    PL.mat1 = landscape.matrix.from.list(PL.list1)
    PL.mat2 = landscape.matrix.from.list(PL.list2)
    perz[ii,jj]=permutation.test(PL.mat1,PL.mat2,num.repeats=testn)
  }
}

save(perz, file="test_persistence_landscape_fpr.Rdata")
