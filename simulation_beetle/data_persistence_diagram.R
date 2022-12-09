library(deSolve)
library(tidyverse)
library(TDA)

##########################################################
# Data
##########################################################
set.seed(5)
nset=2000*2
ntime=120
Time_Series_Data = matrix(NA,2*nset,ntime+1)

# aperiodic
for (i in 1:nset){
  parameters2 <- c(b = 7.48, c_ea = 0.009, c_pa = 0.004, c_el = 0.012, u_p = 0, u_l = 0.267, u_a = 0.96)
  state2 <- c(L = sample(2:100,1), P = sample(2:100,1), A = sample(2:100,1))
  beetles2<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      sig=0.1
      E1t = rnorm(1,0,sig)
      E2t = rnorm(1,0,sig)
      E3t = rnorm(1,0,sig)
      L1 <- b * A * exp((-c_el * L) - (c_ea * A)+E1t)
      P1 <- L * (1 - u_l) * exp(E2t)
      A1 <- (P * exp(-c_pa * A) + A * (1 - u_a)) * exp(E3t)
      list(c(L1, P1, A1))
    }) 
  }
  Aperiodic_i <- ode(y = state2, times = seq(0, ntime), func = beetles2, parms = parameters2, method = "iteration")
  Time_Series_Data[i,] <- as.vector(Aperiodic_i[,4])
}

# stable
for (i in 1:nset){
  parameters1 <- c(b = 7.48, c_ea = 0.009, c_pa = 0.004, c_el = 0.012, u_p = 0, u_l = 0.267, u_a = 0.73)
  state1 <- c(L = sample(2:100,1), P = sample(2:100,1), A = sample(2:100,1))
  beetles1<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      sig=0.1
      E1t = rnorm(1,0,sig)
      E2t = rnorm(1,0,sig)
      E3t = rnorm(1,0,sig)
      L1 <- b * A * exp((-c_el * L) - (c_ea * A)+E1t)
      P1 <- L * (1 - u_l) * exp(E2t)
      A1 <- (P * exp(-c_pa * A) + A * (1 - u_a)) * exp(E3t) 
      list(c(L1, P1, A1))
    }) 
  }
  Stable_i <- ode(y = state1, times = seq(0, ntime), func = beetles1, parms = parameters1, method = "iteration")
  Time_Series_Data[nset+i,] <- as.vector(Stable_i[,4])
}

save(Time_Series_Data, file="data_time_series.Rdata")

##########################################################
# Persistence diagram
##########################################################
pdlist=list()
for (i in 1:(2*nset)){
  if (i %% 100 == 0) print(i)
  x<- data.matrix(Time_Series_Data[i,])
  a<- buildTakens(x,2,3)
  pdlist[[i]] <- ripsDiag(a,maxdimension=1,maxscale=150)$diag
}

ndat=nset*2
maxdeath=rep(NA,ndat)
for (ii in 1:ndat){
  pd=data.frame(matrix(as.vector(pdlist[[ii]]),
                       nrow=nrow(pdlist[[ii]]),
                       ncol=ncol(pdlist[[ii]])))
  colnames(pd) = c("dimension","birth","death")
  pd1=pd %>%
    filter(dimension==1) 
  maxdeath[ii]=max(pd1$death)
}
summary(maxdeath)

save(pdlist, file="data_beetle.Rdata")

##########################################################
# plots
##########################################################
plotdata=data.frame(weeks = 2*rep(1:50,2),
                    y=c(Time_Series_Data[4,1:50],Time_Series_Data[nset+4,1:50]),
                    Regime=c(rep("Aperiodic",50),rep("Stable",50)))
ggplot(plotdata, aes(x=weeks,y=y,color=Regime,linetype=Regime)) +
  geom_line()+
  labs(x="Weeks",y="Adult Population")+
  scale_color_manual(values = c("red","blue"), 
                     labels = c("Aperiodic", "Stable"))+
  scale_linetype_manual(values = c(1,5), 
                        labels = c("Aperiodic", "Stable"))+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = "top", 
        legend.background=element_rect(colour = "white"))


plotdata=data.frame(weeks = 2*rep(1:50,2),
                    y=c(Time_Series_Data[1,1:50],Time_Series_Data[nset+1,1:50]),
                    Regime=c(rep("Aperiodic",50),rep("Stable",50)))
ggplot(plotdata, aes(x=weeks,y=y,color=Regime,linetype=Regime)) +
  geom_line()+
  labs(x="Weeks",y="Adult Population")+
  scale_color_manual(values = c("red","blue"), 
                     labels = c("Aperiodic", "Stable"))+
  scale_linetype_manual(values = c(1,5), 
                        labels = c("Aperiodic", "Stable"))+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "top", 
        legend.background=element_rect(colour = "white"),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))+
  ylim(0,150)


# aperiodic
x<- data.matrix(Time_Series_Data[1,])
a<- buildTakens(x,2,3)
ggplot(data.frame(a,col="a"))+
  geom_point(aes(x=X1,y=X2,color=col))+
  labs(x=expression(Z[t]),y=expression(Z[t+3]))+
  scale_color_manual(values = c("red"))+
  theme_minimal()+
  theme(legend.position = "none")+
  ylim(0,160)+
  xlim(0,160)

# stable
x<- data.matrix(Time_Series_Data[nset+1,])
a<- buildTakens(x,2,3)
ggplot(data.frame(a,col="a"))+
  geom_point(aes(x=X1,y=X2,color=col))+
  labs(x=expression(Z[t]),y=expression(Z[t+3]))+
  scale_color_manual(values = c("blue"))+
  theme_minimal()+
  theme(legend.position = "none")+
  ylim(0,160)+
  xlim(0,160)
