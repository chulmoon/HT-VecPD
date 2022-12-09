source("functions.R")

# import persistence diagrams
load("./data_time_series.Rdata")

#########################################################
# between aperiodic and stable regimes
#########################################################
set.seed(6)

nts=20
np=500 # number of permutation
nrep=100 # number of cases
nset=20 # number of data in one set
ntotalset=4000

basef=rep(NA,nrep)
perf=matrix(NA,nrep,np)

for (rr in 1:nrep) {
  print(rr)
  # base
  ts1 = Time_Series_Data[((rr-1)*nts+1):(rr*nts),]
  ts2 = Time_Series_Data[(ntotalset+(rr-1)*nts+1):(ntotalset+rr*nts),]

  basef[rr]=baseline(ts1,ts2,nts)
  
  # permutation
  tsdata = rbind_dtw(ts1,ts2)
  perf[rr,]=apply(as.matrix(1:np),1,permtest_dtw,tsdata,rr,nts,nset,ntotalset)
}

perz=rep(NA,nrep)
for (ii in 1:nrep) {
  perz[ii]=sum(perf[ii,]<basef[ii])/np
}

sum(perz<0.05)/100 # power

#########################################################
# between aperiodic regimes
#########################################################

basef=rep(NA,nrep)
perf=matrix(NA,nrep,np)

for (rr in 1:nrep) {
  print(rr)
  # base
  ts1 = Time_Series_Data[((rr-1)*nts+1):(rr*nts),]
  ts2 = Time_Series_Data[(ntotalset/2+(rr-1)*nts+1):(ntotalset/2+rr*nts),]
  
  basef[rr]=baseline_dtw(ts1,ts2,nts)
  
  # permutation
  tsdata = rbind(ts1,ts2)
  perf[rr,]=apply(as.matrix(1:np),1,permtest_dtw,tsdata,rr,nts,nset,ntotalset)
}


perz=rep(NA,nrep)
for (ii in 1:nrep) {
  perz[ii]=sum(perf[ii,]<basef[ii])/np
}

sum(perz<0.05)/100 # false discovery rate


#########################################################
# between stable regimes
#########################################################

basef=rep(NA,nrep)
perf=matrix(NA,nrep,np)

for (rr in 1:nrep) {
  print(rr)
  # base
  ts1 = Time_Series_Data[(ntotalset+(rr-1)*nts+1):(ntotalset+rr*nts),]
  ts2 = Time_Series_Data[(ntotalset+ntotalset/2+(rr-1)*nts+1):(ntotalset+ntotalset/2+rr*nts),]
  
  basef[rr]=baseline_dtw(ts1,ts2,nts)
  
  # permutation
  tsdata = rbind(ts1,ts2)
  perf[rr,]=apply(as.matrix(1:np),1,permtest_dtw,tsdata,rr,nts,nset,ntotalset)
}


perz=rep(NA,nrep)
for (ii in 1:nrep) {
  perz[ii]=sum(perf[ii,]<basef[ii])/np
}

sum(perz<0.05)/100 # false discovery rate
