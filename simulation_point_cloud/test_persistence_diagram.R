source("functions.R")

load("data_onecircle.Rdata")
load("data_twocircle.Rdata")

#####################################################################
# False positive rate
#####################################################################
sig=seq(0.05,0.20,length.out = 4) # noise
nsig=length(sig)

npc=10 # number of point cloud in one set
nset=500 # number of sets
testn=1000 # N_p

shufflelabel=gtools::combinations(2*npc,npc) %>%
  as.data.frame() %>%
  sample_n(testn) %>%
  as.matrix()
totallabel=1:(2*npc)

basef=matrix(NA,nsig,nset)
perf=array(NA,c(nsig,nset,testn))


cl = makeCluster(detectCores()-1)

set.seed(17)
# combinations
cbn=combn(1:nset,2)

for (ii in 1:nsig) {
  ind=t(cbn[,sample(1:ncol(cbn),nset)])
  for (jj in 1:nset) {
    orgpd1=twopd[[ii]][(npc*(ind[jj,1]-1)+1):(npc*ind[jj,1])]
    orgpd2=twopd[[ii]][(npc*(ind[jj,2]-1)+1):(npc*ind[jj,2])]
    basef[ii,jj]=rttest(orgpd1,orgpd2,npc,nset)
    mergeddat=append(orgpd1,orgpd2)
    kk=as.matrix(1:testn)
    perf[ii,jj,]=parApply(cl,kk,1,parrttest,mergeddat,totallabel,shufflelabel,npc,nset)
  }
}

stopCluster(cl)

perz=matrix(NA,nsig,nset)
for (ii in 1:nsig) {
  for (jj in 1:nset) {
    perz[ii,jj]=sum(perf[ii,jj,]<basef[ii,jj])/testn
  }
}

save(perz, file="test_persistence_diagram_fpr.Rdata")

#####################################################################
# Power
#####################################################################
set.seed(16)

shufflelabel=gtools::combinations(2*npc,npc) %>%
  as.data.frame() %>%
  sample_n(testn) %>%
  as.matrix()
totallabel=1:(2*npc)

basef=matrix(NA,nsig,nset)
perf=array(NA,c(nsig,nset,testn))

cl = makeCluster(detectCores()-1)

for (ii in 1:nsig) {
  for (jj in 1:nset) {
    orgpd1=onepd[[ii]][(npc*(jj-1)+1):(npc*jj)]
    orgpd2=twopd[[ii]][(npc*(jj-1)+1):(npc*jj)]
    basef[ii,jj]=rttest(orgpd1,orgpd2,npc,nset)
    mergeddat=append(orgpd1,orgpd2)
    kk=as.matrix(1:testn)
    perf[ii,jj,]=parApply(cl,kk,1,parrttest,mergeddat,totallabel,shufflelabel1,npc,nset)
  }
}

stopCluster(cl)

perz=matrix(NA,nsig,nset)
for (ii in 1:nsig) {
  for (jj in 1:nset) {
    perz[ii,jj]=sum(perf[ii,jj,]<basef[ii,jj])/testn
  }
}

save(perz, file="test_persistence_diagram_power.Rdata")
