source("functions.R")

# settings
res=40 # resolution 40 by 40
rangex=c(0,2)
rangey=c(0,2)
h=0.5

# one circle
## load data
load("data_onecircle.Rdata")

onepi_constant = list()
onepi_linear = list()
onepi_arctan = list()

for (ii in 1:4){
  onepi_constant[[ii]] = list()
  onepi_linear[[ii]] = list()
  onepi_arctan[[ii]] = list()
  for (jj in 1:5000){
    pd=onepd[[ii]][[jj]]
    pdmat=data.frame(matrix(as.vector(pd),nrow(pd),3))
    colnames(pdmat)=c("dimension","birth","death") 
    pd1 = pdmat %>%  
      dplyr::filter(dimension==1) %>%
      dplyr::select(-dimension)
    pd1$death = pd1$death - pd1$birth
    
    # weights
    arctanweight = arctanw(pd1,0.5,0.5)
    linearweight = linearw(pd1)
    
    # compute persistence image
    pi1_constant=pers.image(pd=pd1, rangex=rangex, rangey=rangey, wgt=rep(1,nrow(pd1)), nbins=res, h=h)
    onepi_constant[[ii]][[jj]] = pi1_constant
    pi1_arctan=pers.image(pd=pd1, rangex=rangex, rangey=rangey, wgt=arctanweight, nbins=res, h=h)
    onepi_arctan[[ii]][[jj]] = pi1_arctan
    pi1_linear=pers.image(pd=pd1, rangex=rangex, rangey=rangey, wgt=linearweight, nbins=res, h=h)
    onepi_linear[[ii]][[jj]] = pi1_linear
  }
}

## save persistence images
save(onepi_constant, file=paste0("data_pi_onecircle_constant.Rdata"))
save(onepi_arctan, file=paste0("data_pi_onecircle_arctan.Rdata"))
save(onepi_linear, file=paste0("data_pi_onecircle_linear.Rdata"))


# two circle
## load data
load("data_twocircle.Rdata")

twopi_constant = list()
twopi_linear = list()
twopi_arctan = list()

for (ii in 1:4){
  twopi_constant[[ii]] = list()
  twopi_linear[[ii]] = list()
  twopi_arctan[[ii]] = list()
  for (jj in 1:5000){
    pd=twopd[[ii]][[jj]]
    pdmat=data.frame(matrix(as.vector(pd),nrow(pd),3))
    colnames(pdmat)=c("dimension","birth","death")
    pd1 = pdmat %>%  
      dplyr::filter(dimension==1) %>%
      dplyr::select(-dimension)
    pd1$death = pd1$death - pd1$birth
    
    # weights
    arctanweight = arctanw(pd1,0.5,0.5)
    linearweight = linearw(pd1)
    
    # compute persistence image
    pi1_constant=pers.image(pd=pd1, rangex=rangex, rangey=rangey, wgt=rep(1,nrow(pd1)), nbins=res, h=h)
    twopi_constant[[ii]][[jj]] = pi1_constant
    pi1_arctan=pers.image(pd=pd1, rangex=rangex, rangey=rangey, wgt=arctanweight, nbins=res, h=h)
    twopi_arctan[[ii]][[jj]] = pi1_arctan
    pi1_linear=pers.image(pd=pd1, rangex=rangex, rangey=rangey, wgt=linearweight, nbins=res, h=h)
    twopi_linear[[ii]][[jj]] = pi1_linear
  }
}

## save persistence images
save(twopi_constant, file=paste0("data_pi_twocircle_constant.Rdata"))
save(twopi_arctan, file=paste0("data_pi_twocircle_arctan.Rdata"))
save(twopi_linear, file=paste0("data_pi_twocircle_linear.Rdata"))
