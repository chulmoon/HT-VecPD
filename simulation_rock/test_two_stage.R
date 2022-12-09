source("functions.R")

# import persistence diagrams
load("./rock_simulation_pd.Rdata")

#######################################################
# persistence image
#######################################################
# parameters for persistence image
res0=mmzero[2]-mmzero[1]+1
res1=mmone[2]-mmone[1]+1
h=3 # smoothing parameter

# generate persistence images
## arctangent weight
pilist0_arctan = lapply(1:length(flist),pimain,pdzero,
                 rangex=mmzero, rangey=mmzero+abs(min(mmzero)), 
                 nbins=res0, h=h, weight="arctan")
pilist1_arctan = lapply(1:length(flist),pimain,pdone,rangex=mmone, 
                 rangey=mmone+abs(min(mmone)), 
                 nbins=res1, h=h, weight="arctan")

## linear weight
pilist0_linear = lapply(1:length(flist),pimain,pdzero,
                        rangex=mmzero, rangey=mmzero+abs(min(mmzero)), 
                        nbins=res0, h=h, weight="linear")
pilist1_linear = lapply(1:length(flist),pimain,pdone,rangex=mmone, 
                        rangey=mmone+abs(min(mmone)), 
                        nbins=res1, h=h, weight="linear")

## constant weight
pilist0_constant = lapply(1:length(flist),pimain,pdzero,
                        rangex=mmzero, rangey=mmzero+abs(min(mmzero)), 
                        nbins=res0, h=h, weight="constant")
pilist1_constant = lapply(1:length(flist),pimain,pdone,rangex=mmone, 
                        rangey=mmone+abs(min(mmone)), 
                        nbins=res1, h=h, weight="constant")

#######################################################
# test
#######################################################
## arctangent
test_results(pilist0_arctan,pilist1_arctan,mmzero,mmone,cc=0.6)
## linear
test_results(pilist0_linear,pilist1_linear,mmzero,mmone,cc=0.6)
## constant
test_results(pilist0_constant,pilist1_constant,mmzero,mmone,cc=0.6)
