source("functions.R")

# filepath
fpath="./data/"
flist=list.files(fpath)

cl = makeCluster(detectCores()-1)

mmzero=c(-15,14) # range
reszero=(mmzero[2]-mmzero[1])
pizero=array(NA,c(reszero,reszero,length(flist)))

mmone=c(-12,27) # range
resone=(mmone[2]-mmone[1])
pione=array(NA,c(resone,resone,length(flist)))

mmtwo=c(-2,28) # range
restwo=(mmtwo[2]-mmtwo[1])
pitwo=array(NA,c(restwo,restwo,length(flist)))

for (ii in 1:length(flist)){
  fname=paste(fpath,flist[ii],sep="")
  pd=read.table(fname)
  pdzero=pd%>%
    filter(pd$V1==0) %>%
    select(V2,V3)
  pizero[,,ii]=makePI(pdzero,reszero,mmzero,weight="pwgk",sig=NULL,cluster=cl)
  pdone=pd%>%
    filter(pd$V1==1) %>%
    select(V2,V3)
  pione[,,ii]=makePI(pdone,resone,mmone,weight="pwgk",sig=NULL,cluster=cl)
  pdtwo=pd%>%
    filter(pd$V1==2) %>%
    select(V2,V3)
  pitwo[,,ii]=makePI(pdtwo,restwo,mmtwo,weight="pwgk",sig=NULL,cluster=cl)
}
stopCluster(cl)

save(pizero, mmzero, file="rock_dim0_arctan.Rdata")
save(pione, mmone, file="rock_dim1_arctan.Rdata")
save(pitwo, mmtwo, file="rock_dim2_arctan.Rdata")
