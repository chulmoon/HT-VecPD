source("functions.R")

# import persistence diagrams
load("./data_persistence_diagram.Rdata")

#########################################################
# between aperiodic and stable regimes
#########################################################
set.seed(6)

np=100 # number of permutation

nrep=100 # number of sets
nset=20 # number of data in one set
ntotalset=4000

basef=rep(NA,nrep)
perf=matrix(NA,nrep,np)

cl = makeCluster(detectCores()-1)

for (rr in 1:nrep) {
  # base
  pdlist1 = list()
  pdlist2 = list()
  for (ss in 1:nset){
    pdlist1[[ss]]=pdlist[[(rr-1)*nset + ss]]
    pdlist2[[ss]]=pdlist[[(rr-1)*nset  +ntotalset+ ss]]
  }
  basef[rr]=rttest(pdlist1,pdlist2,nset)
  
  # permutation
  kk=as.matrix(1:np)
  perf[rr,]=parApply(cl,kk,1,parrttest,pdlist,rr,nset,ntotalset)
}

stopCluster(cl)


perz=rep(NA,nrep)
for (ii in 1:nrep) {
  perz[ii]=sum(perf[ii,]<basef[ii])/np
}


#########################################################
# between aperiodic regimes
#########################################################

set.seed(8)

nset=20
np=100 # number of permutation

nrep=100 # number of sets
nset=20 # number of data in one set
ntotalset=2000

basef=rep(NA,nrep)
perf=matrix(NA,nrep,np)

cl = makeCluster(detectCores()-1)

for (rr in 1:nrep) {
  # base
  pdlist1 = list()
  pdlist2 = list()
  for (ss in 1:nset){
    pdlist1[[ss]]=pdlist[[(rr-1)*nset + ss]]
    pdlist2[[ss]]=pdlist[[(rr-1)*nset  +ntotalset+ ss]]
  }
  basef[rr]=rttest(pdlist1,pdlist2,nset)
  
  # permutation
  kk=as.matrix(1:np)
  perf[rr,]=parApply(cl,kk,1,parrttest,pdlist,rr,nset,ntotalset)
}

stopCluster(cl)

perz=rep(NA,nrep)
for (ii in 1:nrep) {
  perz[ii]=sum(perf[ii,]<basef[ii])/np
}


#########################################################
# between stable regimes
#########################################################
set.seed(7)

np=100 # number of permutation

nrep=100 # number of sets
nset=20 # number of data in one set
ntotalset=2000


basef=rep(NA,nrep)
perf=matrix(NA,nrep,np)

cl = makeCluster(detectCores()-1)

for (rr in 1:nrep) {
  # base
  pdlist1 = list()
  pdlist2 = list()
  for (ss in 1:nset){
    pdlist1[[ss]]=pdlist[[(rr-1)*nset + ss+4000]]
    pdlist2[[ss]]=pdlist[[(rr-1)*nset  +ntotalset+ ss+4000]]
  }
  basef[rr]=rttest(pdlist1,pdlist2,nset)
  
  # permutation
  kk=as.matrix(1:np)
  perf[rr,]=parApply(cl,kk,1,parrttest,pdlist,rr,nset,ntotalset)
}

stopCluster(cl)

perz=rep(NA,nrep)
for (ii in 1:nrep) {
  perz[ii]=sum(perf[ii,]<basef[ii])/np
}

