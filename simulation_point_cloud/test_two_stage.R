source("functions.R")

# settings
sig=seq(0.05,0.20,length.out = 4)
npc=10
nset=500
range=c(0,2)
res=40
alpha=0.05


# arctan
load("data_pi_onecircle_arctan.Rdata")
load("data_pi_twocircle_arctan.Rdata")
## FPR
pvals_arctan_fpr=ts_main_fpr(twopi_arctan,sig,npc,nset,range,res,alpha)
save(pvals_arctan_fpr,file="test_two_stage_arctan_fpr.Rdata")
## Power
pvals_arctan_power=ts_main_power(onepi_arctan,twopi_arctan,sig,npc,nset,range,res,alpha)
save(pvals_arctan_power,file="test_two_stage_arctan_power.Rdata")


# constant
load("data_pi_onecircle_constant.Rdata")
load("data_pi_twocircle_constant.Rdata")
## FPR
pvals_constant_fpr=ts_main_fpr(twopi_constant,sig,npc,nset,range,res,alpha)
save(pvals_constant_fpr,file="test_two_stage_constant_fpr.Rdata")
## Power
pvals_constant_power=ts_main_power(onepi_constant,twopi_constant,sig,npc,nset,range,res,alpha)
save(pvals_constant_power,file="test_two_stage_constant_power.Rdata")


# linear
load("data_pi_onecircle_linear.Rdata")
load("data_pi_twocircle_linear.Rdata")
## FPR
pvals_linear_fpr=ts_main_fpr(twopi_linear,sig,npc,nset,range,res,alpha)
save(pvals_linear_fpr,file="test_two_stage_linear_fpr.Rdata")
## Power
pvals_linear_power=ts_main_power(onepi_linear,twopi_linear,sig,npc,nset,range,res,alpha)
save(pvals_linear_power,file="test_two_stage_linear_power.Rdata")



