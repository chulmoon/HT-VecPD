source("functions.R")

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

# range
mmzero=c(-15,14)
mmone=c(-12,27) 
mmtwo=c(-2,28)

# number of permutation
np=1000

cl = makeCluster(detectCores()-1)

PL.list1=list()
PL.list2=list()

######################################################
# F42B and F42C
######################################################
group1=totalpd[(0*27+1):(1*27)]
group2=totalpd[(1*27+1):(2*27)]

## dim 0
for (kk in 1:27){
  PL.list1[[kk]]=t(landscape(group1[[kk]][,1:3],dimension=0,KK=1:300,tseq=seq(mmzero[1],mmzero[2],length=500) ))
  PL.list2[[kk]]=t(landscape(group2[[kk]][,1:3],dimension=0,KK=1:300,tseq=seq(mmzero[1],mmzero[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

## dim 1
for (kk in 1:27){
  PL.list1[[kk]]=t(landscape(group1[[kk]][,1:3],dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
  PL.list2[[kk]]=t(landscape(group2[[kk]][,1:3],dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

## dim 2
for (kk in 1:27){
  PL.list1[[kk]]=t(landscape(group1[[kk]][,1:3],dimension=2,KK=1:300,tseq=seq(mmtwo[1],mmtwo[2],length=500) ))
  PL.list2[[kk]]=t(landscape(group2[[kk]][,1:3],dimension=2,KK=1:300,tseq=seq(mmtwo[1],mmtwo[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

######################################################
# LV60A and LV60C
######################################################
group1=totalpd[(2*27+1):(3*27)]
group2=totalpd[(3*27+1):(4*27)]

## dim 0
for (kk in 1:27){
  PL.list1[[kk]]=t(landscape(group1[[kk]][,1:3],dimension=0,KK=1:300,tseq=seq(mmzero[1],mmzero[2],length=500) ))
  PL.list2[[kk]]=t(landscape(group2[[kk]][,1:3],dimension=0,KK=1:300,tseq=seq(mmzero[1],mmzero[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

## dim 1
for (kk in 1:27){
  PL.list1[[kk]]=t(landscape(group1[[kk]][,1:3],dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
  PL.list2[[kk]]=t(landscape(group2[[kk]][,1:3],dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

## dim 2
for (kk in 1:27){
  PL.list1[[kk]]=t(landscape(group1[[kk]][,1:3],dimension=2,KK=1:300,tseq=seq(mmtwo[1],mmtwo[2],length=500) ))
  PL.list2[[kk]]=t(landscape(group2[[kk]][,1:3],dimension=2,KK=1:300,tseq=seq(mmtwo[1],mmtwo[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

######################################################
# F42B and LV60A
######################################################
group1=totalpd[(0*27+1):(1*27)]
group2=totalpd[(2*27+1):(3*27)]

## dim 0
for (kk in 1:27){
  PL.list1[[kk]]=t(landscape(group1[[kk]][,1:3],dimension=0,KK=1:300,tseq=seq(mmzero[1],mmzero[2],length=500) ))
  PL.list2[[kk]]=t(landscape(group2[[kk]][,1:3],dimension=0,KK=1:300,tseq=seq(mmzero[1],mmzero[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

## dim 1
for (kk in 1:27){
  PL.list1[[kk]]=t(landscape(group1[[kk]][,1:3],dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
  PL.list2[[kk]]=t(landscape(group2[[kk]][,1:3],dimension=1,KK=1:300,tseq=seq(mmone[1],mmone[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

## dim 2
for (kk in 1:27){
  PL.list1[[kk]]=t(landscape(group1[[kk]][,1:3],dimension=2,KK=1:300,tseq=seq(mmtwo[1],mmtwo[2],length=500) ))
  PL.list2[[kk]]=t(landscape(group2[[kk]][,1:3],dimension=2,KK=1:300,tseq=seq(mmtwo[1],mmtwo[2],length=500) ))
}
PL.mat1 = landscape.matrix.from.list(PL.list1)
PL.mat2 = landscape.matrix.from.list(PL.list2)
permutation.test(PL.mat1,PL.mat2,num.repeats=np,cl)

stopCluster(cl)