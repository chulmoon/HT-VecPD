source("functions.R")

#######################################################
# dimension zero
#######################################################
# persistence images
load("rock_dim0_arctan.Rdata")

result=pizero
range=mmzero

res=dim(result)[1]

# test between F42B and F42C
group1=result[,,(1:27)]
group2=result[,,((1*27+1):(2*27))]
ttestpi(group1,group2,res,cc=0.5)

# test between LV60A and LV60C
group1=result[,,((2*27+1):(3*27))]
group2=result[,,((3*27+1):(4*27))]
ttestpi(group1,group2,res,cc=0.5)
 
# test between F42B and LV60A
group1=result[,,(1:27)]
group2=result[,,((2*27+1):(3*27))]
ttestpi(group1,group2,res,cc=0.5)

#######################################################
# dimension one
#######################################################
# persistence images
load("rock_dim1_arctan.Rdata")

result=pione
range=mmone

res=dim(result)[1]

# test between F42B and F42C
group1=result[,,(1:27)]
group2=result[,,((1*27+1):(2*27))]
ttestpi(group1,group2,res,cc=0.5)

# test between LV60A and LV60C
group1=result[,,((2*27+1):(3*27))]
group2=result[,,((3*27+1):(4*27))]
ttestpi(group1,group2,res,cc=0.5)

# test between F42B and LV60A
group1=result[,,(1:27)]
group2=result[,,((2*27+1):(3*27))]
ttestpi(group1,group2,res,cc=0.5)

#######################################################
# dimension two
#######################################################
# persistence images
load("rock_dim2_arctan.Rdata")

result=pitwo
range=mmtwo

res=dim(result)[1]

# test between F42B and F42C
group1=result[,,(1:27)]
group2=result[,,((1*27+1):(2*27))]
ttestpi(group1,group2,res,cc=0.5)

# test between LV60A and LV60C
group1=result[,,((2*27+1):(3*27))]
group2=result[,,((3*27+1):(4*27))]
ttestpi(group1,group2,res,cc=0.5)

# test between F42B and LV60A
group1=result[,,(1:27)]
group2=result[,,((2*27+1):(3*27))]
ttestpi(group1,group2,res,cc=0.5)
