source("functions.R")

# import persistence diagrams
load("./rock_simulation_pd.Rdata")

set.seed(83)

cl = makeCluster(detectCores()-1)

# number of permutation N_p
np=1000

PL.list1=list()
PL.list2=list()

###################### Scenario 1 ######################

# dimension 0
for (kk in 1:50){
  pdtemp1 = pdzero[[kk]]
  pdtemp1$dimension=0
  pdtemp1= pdtemp1[,c(3,1,2)]
  colnames(pdtemp1)=c("dimension","Birth","Death")
  PL.list1[[kk]]=t(landscape(pdtemp1,dimension=0,KK=1:150,tseq=seq(mmzero[1],mmzero[2],length=500) ))
  pdtemp2 = pdzero[[kk+50]]
  pdtemp2$dimension=0
  pdtemp2= pdtemp2[,c(3,1,2)]
  colnames(pdtemp2)=c("dimension","Birth","Death")
  PL.list2[[kk]]=t(landscape(pdtemp2,dimension=0,KK=1:150,tseq=seq(mmzero[1],mmzero[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

# dimension 1
for (kk in 1:50){
  pdtemp1 = pdone[[kk]]
  pdtemp1$dimension=1
  pdtemp1= pdtemp1[,c(3,1,2)]
  colnames(pdtemp1)=c("dimension","Birth","Death")
  PL.list1[[kk]]=t(landscape(pdtemp1,dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
  pdtemp2 = pdone[[kk+50]]
  pdtemp2$dimension=1
  pdtemp2= pdtemp2[,c(3,1,2)]
  colnames(pdtemp2)=c("dimension","Birth","Death")
  PL.list2[[kk]]=t(landscape(pdtemp2,dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)


###################### Scenario 2 ######################

# dimension 0
for (kk in 1:50){
  pdtemp1 = pdzero[[kk]]
  pdtemp1$dimension=0
  pdtemp1= pdtemp1[,c(3,1,2)]
  colnames(pdtemp1)=c("dimension","Birth","Death")
  PL.list1[[kk]]=t(landscape(pdtemp1,dimension=0,KK=1:150,tseq=seq(mmzero[1],mmzero[2],length=500) ))
  pdtemp2 = pdzero[[kk+100]]
  pdtemp2$dimension=0
  pdtemp2= pdtemp2[,c(3,1,2)]
  colnames(pdtemp2)=c("dimension","Birth","Death")
  PL.list2[[kk]]=t(landscape(pdtemp2,dimension=0,KK=1:150,tseq=seq(mmzero[1],mmzero[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

# dimension 1
for (kk in 1:50){
  pdtemp1 = pdone[[kk]]
  pdtemp1$dimension=1
  pdtemp1= pdtemp1[,c(3,1,2)]
  colnames(pdtemp1)=c("dimension","Birth","Death")
  PL.list1[[kk]]=t(landscape(pdtemp1,dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
  pdtemp2 = pdone[[kk+100]]
  pdtemp2$dimension=1
  pdtemp2= pdtemp2[,c(3,1,2)]
  colnames(pdtemp2)=c("dimension","Birth","Death")
  PL.list2[[kk]]=t(landscape(pdtemp2,dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)


###################### Scenario 3 ######################

# dimension 0
for (kk in 1:50){
  pdtemp1 = pdzero[[kk]]
  pdtemp1$dimension=0
  pdtemp1= pdtemp1[,c(3,1,2)]
  colnames(pdtemp1)=c("dimension","Birth","Death")
  PL.list1[[kk]]=t(landscape(pdtemp1,dimension=0,KK=1:150,tseq=seq(mmzero[1],mmzero[2],length=500) ))
  pdtemp2 = pdzero[[kk+150]]
  pdtemp2$dimension=0
  pdtemp2= pdtemp2[,c(3,1,2)]
  colnames(pdtemp2)=c("dimension","Birth","Death")
  PL.list2[[kk]]=t(landscape(pdtemp2,dimension=0,KK=1:150,tseq=seq(mmzero[1],mmzero[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

# dimension 1
for (kk in 1:50){
  pdtemp1 = pdone[[kk]]
  pdtemp1$dimension=1
  pdtemp1= pdtemp1[,c(3,1,2)]
  colnames(pdtemp1)=c("dimension","Birth","Death")
  PL.list1[[kk]]=t(landscape(pdtemp1,dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
  pdtemp2 = pdone[[kk+150]]
  pdtemp2$dimension=1
  pdtemp2= pdtemp2[,c(3,1,2)]
  colnames(pdtemp2)=c("dimension","Birth","Death")
  PL.list2[[kk]]=t(landscape(pdtemp2,dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

stopCluster(cl)