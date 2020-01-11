virtual_eco_integ <- function(P.TL1=11195.44, TE=0.1, Temp=15){

#1) Set Tl vectors and Production of TL 1
seq_tl<-c(1, seq(2,7,.1));length(seq_tl)
seq_tl.2up<-seq(2,7,.1);length(seq_tl.2up)    # generate a sequence of TL without TL1

#2) calculate production, flow and mean flow of TL =>2  

# for TL=2
P<-c()
P[1]<-P.TL1
flow_mean<-c()
flow_mean[1]<-P.TL1
flow<-c()
flow[1]<-P.TL1


# I am not back calculating the production of trophic level 2. I will use the 
# flow of trophic level 1 to calculate the flow at trophic level 2.
flow[2] <- flow[1] * TE^1

# for TL>2
for (tl in 2:(length(seq_tl)-1)) { # tl [-1] because we start the iterative process at tl=2.1
  # production flow
  flow[tl+1] <- flow[tl] * TE^0.1
  for (tl2 in 2:length(seq_tl)){
  # mean production flow
  flow_mean[tl2] <- (1/0.1) * flow[tl2] * (1-TE^0.1)/(-log(TE))
  # production
  P[tl2] <- flow_mean[tl2] * 0.1
}}

#3) calculate biomass
B.TL1<-P.TL1/kin(1,Temp);length(B.TL1)
B.TL2up<-P[-1]/kin(seq_tl.2up,Temp);length(B.TL2up)

#4) Calculate kinetics
Kinet<-kin(tl = seq_tl, sst = Temp)

#5) Fishing vector
Fish_mort_acc=rep(0,52)                                      # cero fishing
Fish_mort=rep(0,52)

# 6) Y_tot
Y_tot=rep(0,52)

#7) Calculate fishing_loss rate and accessible fishing loss rate
F_loss <- rep(0,52)
F_loss_acc <- rep(0,52)

#8) paste, save and plot results

B<-c(B.TL1, B.TL2up);length(B)
res<-as.data.frame(cbind(seq_tl, P,flow, flow_mean, B, Kinet, Y_tot=Y_tot, F_loss=F_loss,
                         F_loss_acc=F_loss_acc, Fish_mort=Fish_mort, 
                         Fish_mort_acc=Fish_mort_acc))
return(res)
}
