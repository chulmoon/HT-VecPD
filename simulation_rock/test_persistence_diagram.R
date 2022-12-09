load("./rock_simulation_pd.Rdata")
totalset=50 # number of images in a group
np=100 # number of permutations
cl=makeCluster(detectCores()-1)
clusterEvalQ(cl, library(tidyverse))

###################### Scenario 1 ######################

# dimension 0
set.seed(46)
dim=0
pd1=pdzero[1:50]
pd2=pdzero[51:100]
mergeddat=append(pd1,pd2)

## initial joint loss
basef0=rttest(pd1,pd2,totalset,dim)

## shuffled labels
shufflelabel=matrix(NA,permun,totalset)
for (ii in 1:100){
  shufflelabel[ii,]=sample(1:(2*totalset),totalset)
}

## permutation test
perf_pd0=parApply(cl,as.matrix(1:np),1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(basef0>perf_pd0)/np # p-value


# dimension 1
set.seed(46)
dim=1
pd1=pdone[1:50]
pd2=pdone[51:100]
mergeddat=append(pd1,pd2)

## initial joint loss
basef1=rttest(pd1,pd2,totalset,dim)

## shuffled labels
shufflelabel=matrix(NA,permun,totalset)
for (ii in 1:100){
  shufflelabel[ii,]=sample(1:(2*totalset),totalset)
}

## permutation test
perf_pd1=parApply(cl,as.matrix(1:np),1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(basef1>perf_pd1)/np # p-value

###################### Scenario 2 ######################

# dimension 0
set.seed(47)
dim=0
pd1=pdzero[1:50]
pd2=pdzero[101:150]
mergeddat=append(pd1,pd2)

## initial joint loss
basef0=rttest(pd1,pd2,totalset,dim)

## shuffled labels
shufflelabel=matrix(NA,permun,totalset)
for (ii in 1:100){
  shufflelabel[ii,]=sample(1:(2*totalset),totalset)
}

## permutation test
perf_pd0=parApply(cl,as.matrix(1:np),1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(basef0>perf_pd0)/np # p-value


# dimension 1
set.seed(47)
dim=1
pd1=pdone[1:50]
pd2=pdone[101:150]
mergeddat=append(pd1,pd2)

## initial joint loss
basef1=rttest(pd1,pd2,totalset,dim)

## shuffled labels
shufflelabel=matrix(NA,permun,totalset)
for (ii in 1:100){
  shufflelabel[ii,]=sample(1:(2*totalset),totalset)
}

## permutation test
perf_pd1=parApply(cl,as.matrix(1:np),1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(basef1>perf_pd1)/np # p-value

###################### Scenario 3 ######################

# dimension 0
set.seed(48)
dim=0
pd1=pdzero[1:50]
pd2=pdzero[151:200]
mergeddat=append(pd1,pd2)

## initial joint loss
basef0=rttest(pd1,pd2,totalset,dim)

## shuffled labels
shufflelabel=matrix(NA,permun,totalset)
for (ii in 1:100){
  shufflelabel[ii,]=sample(1:(2*totalset),totalset)
}

## permutation test
perf_pd0=parApply(cl,as.matrix(1:np),1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(basef0>perf_pd0)/np # p-value


# dimension 1
set.seed(48)
dim=1
pd1=pdone[1:50]
pd2=pdone[151:200]
mergeddat=append(pd1,pd2)

## initial joint loss
basef1=rttest(pd1,pd2,totalset,dim)

## shuffled labels
shufflelabel=matrix(NA,permun,totalset)
for (ii in 1:100){
  shufflelabel[ii,]=sample(1:(2*totalset),totalset)
}

## permutation test
perf_pd1=parApply(cl,as.matrix(1:np),1,parrttest,mergeddat,shufflelabel,totalset,dim)
sum(basef1>perf_pd1)/np # p-value

stopCluster(cl)