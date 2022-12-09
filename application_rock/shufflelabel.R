library(TDA)
library(pracma)
library(mvtnorm)
library(tidyverse)
library(RColorBrewer)
library(gtools)
library(parallel)

# number of combinations of 20-subsets
testn=100

totalset=27
selectedset=20

# label for F42A
shufflelabel1=combinations(totalset,selectedset) %>%
	as.data.frame() %>%
	sample_n(testn) %>%
	as.matrix()

# label for F42B
shufflelabel2=combinations(totalset,selectedset) %>%
	as.data.frame() %>%
	sample_n(testn) %>%
	as.matrix()
shufflelabel2=shufflelabel2+totalset


# label for F42C
shufflelabel3=combinations(totalset,selectedset) %>%
	as.data.frame() %>%
	sample_n(testn) %>%
	as.matrix()
shufflelabel3=shufflelabel3+totalset*2

# label for LV60A
shufflelabel4=combinations(totalset,selectedset) %>%
	as.data.frame() %>%
	sample_n(testn) %>%
	as.matrix()
shufflelabel4=shufflelabel4+totalset*3

# label for LV60C
shufflelabel5=combinations(totalset,selectedset) %>%
	as.data.frame() %>%
	sample_n(testn) %>%
	as.matrix()
shufflelabel5=shufflelabel5+totalset*4

# save labels for rock combinations
save(shufflelabel1,shufflelabel2,shufflelabel3,shufflelabel4,shufflelabel5,file="C:/Users/Chul/Box Sync/r_Papers/HypothesisTesting/code/shufflelabel.Rdata")
