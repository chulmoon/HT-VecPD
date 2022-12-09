source("functions.R")

#######################################################
# persistence diagram
#######################################################
# filepath
fpath="./data/"
flist=list.files(fpath,pattern="+pd")

# find minimum and maximum for each dimension
totalpdzero=c()
totalpdone=c()

# aggregated pd
for (ii in 1:length(flist)){
  fname=paste(fpath,flist[ii],sep="")
  pd=read.table(fname)
  colnames(pd)=c("dimension","birth","death")
  pd0=pd %>%
    filter(dimension==0,death<Inf) %>%
    mutate(death=death-birth) %>%
    select(birth,death)
  pd1=pd%>%
    filter(dimension==1,death<Inf)%>%
    mutate(death=death-birth) %>%
    select(birth,death)
  totalpdzero=rbind(totalpdzero,pd0)
  totalpdone=rbind(totalpdone,pd1)
}

# range of pd
mmzero=c(floor(min(totalpdzero[,c(1,2)])),
         ceiling(max(totalpdzero[,c(1,2)])))
mmone=c(floor(min(totalpdone[,c(1,2)])),ceiling(max(totalpdone[,c(1,2)])))

# save pd
pdzero=list()
pdone=list()

for (ii in 1:length(flist)){
  fname=paste(fpath,flist[ii],sep="")
  pd=read.table(fname)
  colnames(pd)=c("dimension","birth","death")
  pd0=pd %>%
    filter(dimension==0,death<Inf) %>%
    select(birth,death)
  pd1=pd%>%
    filter(dimension==1,death<Inf) %>%
    select(birth,death)
  
  pdzero[[ii]]=pd0
  pdone[[ii]]=pd1
}

# save pd
save(mmzero,mmone,pdone,pdzero,file = "./rock_simulation_pd.Rdata")