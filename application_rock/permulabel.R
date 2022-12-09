totalset=27

# number of permutations for RT test
permun=500

shufflelabel=data.frame(matrix(NA,permun,totalset))
for (ii in 1:permun) {
	shufflelabel[ii,]=sample(1:(2*totalset),totalset,replace=F)
}

totallabel=1:(2*totalset)

save(shufflelabel,totallabel,file="permulabel.Rdata")
