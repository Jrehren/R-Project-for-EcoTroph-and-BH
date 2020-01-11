dat_frame_integrate_access_adap <- function(VE,asymptote,TL50,slope, Temp=15, TE=0.1){
  
  
  # 1.) SELECTIVITY ----
  
  seq_tl<-c(1, seq(2,7,.1));length(seq_tl)                # Trophic level vector 
  
  seq_tl.2up<-seq(2,7,.1);length(seq_tl.2up)              # generate a sequence 
  # of TL without TL1
  
  accessibility<-s_selec(seq_tl.2up, asymptote=asymptote, TL50 = TL50, 
                         slope=slope)                        # parameters are set 
  # according to Gascuels 
  # excel sheet
  
  accessibility<-c(0,accessibility)                   # accesssible from TL 2
  names(accessibility)<-seq_tl
  
  
  # 2.) CALCULATE ACCESSIBLE KINETICS ----
  
  # We first calculate which of the trophic levels are fully (>=0.99) and which 
  # are not fully accessible
  seq_tl_full_expl<-names(which(accessibility>=0.99))
  seq_tl_notfull_expl<-names(which(accessibility<0.99))
  
  # For those that are fully accessible we will use equation 05 from Gascuel et al. 2008
  # to calculate kinetics and for those that are not fully accessible we will use
  # equation 04
  Kin_acc_notfull_expl<-kin_acc(tl = as.numeric(seq_tl_notfull_expl[-1]), sst = Temp)
  Kin_acc_full_expl<-kin(tl = as.numeric(seq_tl_full_expl), sst = Temp)
  
  # We bind them together with the kinetics of the first TL
  Kin_acc<-c(Kin_acc_notfull_expl, Kin_acc_full_expl)
  
  # Because now we have a plunge between those trophic levels that are fully
  # accessible and those that are not fully accessible, we will smooth the  
  # trophic distribution of accessible kinetics
  if (accessibility[2]<0.99){    # I only want to smooth, if I have differences
    # in the accessibility
    # LOESS SMOOTHING
    kin_acc_log<-log(Kin_acc)
    lo <- loess.smooth(seq_tl[-1], kin_acc_log,span=1/3, evaluation = 51)
    kin_acc_lo<-exp(lo$y);length(kin_acc_lo)
    # Because of the smoothing I might have some values for the accessible kinetics
    # that are above the kinetics of the total biomass and that should not be the
    # case. So I test which trophic levels have higher accessible kinetics than
    # total kinetics and than I set it to be as the total kinetics
    pos<-which(kin_acc_lo>VE$Kinet[-1])
    pos_ch<-seq(pos[1],51, 1)
    kin_acc_lo[c(pos_ch)]<-VE$Kinet[c(pos_ch+1)];length(kin_acc_lo)
    tl1_kin<-VE$Kinet[1]*as.numeric(accessibility[1])
    # bind all
    Kin_acc_fin<-c(tl1_kin, kin_acc_lo)
  } else {                           # here I don't have any differences in 
    # accessibility so I dont want to run a 
    # loess smoothing
    tl1_kin<-VE$Kinet[1]*as.numeric(accessibility[1])
    # bind all
    Kin_acc_fin<-c(tl1_kin, Kin_acc)
  }
  
  
  # 3.) CALCULATE THE ACCESSIBLE FLOW AND MEAN FLOW ----
  # Because we 1) calculate B* by P*/K*, and 2) we want to relate the accessibility 
  # to the biomass and not to the production or the biomass flow, we first must 
  # find an accessible biomass flow vector that leads to the desired B* shape. 
  # We, however, know that B_acc = B*S or P_acc=P*S*(K*/K). The only thing
  # we do not know is the relationship between P_acc and flow_acc and therefore, 
  # we loop over it to find the solution. 
  
  ## A) I set the condition and starting values of the loop ----
  # This is the condition
  sc <- 2;length(sc)
  
  # Here, I set the starting values. 
  
  B_acc<-accessibility*VE$B   # Here, I set out how I want my accessible biomass vector
  # to look like
  P_acc_desired <- VE$P*accessibility*(Kin_acc_fin/VE$Kinet) # this is the accessible P given
  # the desired accessible B
  diff <- rep(1,52)        # This is the vector which will be used to adjust the 
  # flow at each iteration. At the beginning it is 1 later
  # it will be P_acc_desired/P_acc
  
  # we calculate a starting vector for the flow
  flow_acc<-VE$flow*accessibility*(Kin_acc_fin/VE$Kinet) ;length(flow_acc)
  
  flow_acc_temp<-flow_acc          # and save it to a temporal flow
  
  ## B) We run the loop to find the actual accessible flow that gives us the ----
  #     desired B* or P* vector ----
  
  while (sc!=0){
    
    # Here, I am adjusting the accessible flow to deliver the desired accessible P. 
    
    flow_acc<-flow_acc_temp*diff
    
    # N_loss and N_loss_acc 
    seq_Tl_diff<-c(1,rep(0.1, 50));length(seq_Tl_diff)
    N_loss <- c(log(VE$flow[-length(VE$flow)]/VE$flow[-1])/seq_Tl_diff - VE$F_loss[-length(VE$F_loss)], log(VE$flow[51]/VE$flow[52])/0.1 - VE$F_loss[51]);length(N_loss)  
    
    N_loss_acc <- c(log(flow_acc[-length(flow_acc)]/flow_acc[-1])/seq_Tl_diff - VE$F_loss_acc[-length(VE$F_loss_acc)], log(flow_acc[51]/flow_acc[52])/0.1 - VE$F_loss_acc[51])
    
    ## Now I have to recalculate from the production my P and flow mean
    
    flow_mean_acc <- c()
    P_acc <- c()
    
    # TL=1
    flow_mean_acc[1]<-flow_acc[1]
    P_acc[1]<-flow_acc[1]
    
    # mean flow TL>1
    flow_mean_acc <- flow_acc*((1-(exp(-N_loss_acc*0.1)))/(N_loss_acc*0.1))
    
    # mean_flow_acc TL>1
    P_acc<-flow_mean_acc*0.1
    
    # I calculate the difference between the desired and the resulting accessible P
    diff<-P_acc_desired/P_acc;diff
    
    # storing the flow as the temporal flow
    flow_acc_temp<-flow_acc
    
    # We only calculate the condition on TL2:52, because there is no accessible B
    #  on the TL 1
    sc<-round((sum(P_acc_desired[-1]-P_acc[-1]))^2, digits = 24)
    
  }
  
  # 4.) CALCULATE THE ACCESSIBLE BIOMASS ----
  B_acc<-P_acc/Kin_acc_fin;length(B_acc)
  
  
  # 5.) Fishing loss rate ----
  
  F_loss<-VE$Fish_mort*VE$Kinet
  F_loss_acc<-VE$Fish_mort_acc*Kin_acc_fin
  
  
  # 6.) MERGE AND CREATE THE DATA INPUT FOR ET.DIAGNOSE  ----
  
  ## A) ET_Main
  ET_Main<-data.frame(B=VE$B, B_acc=B_acc, P=VE$P, P_acc=P_acc, flow=VE$flow, 
                      flow_acc=flow_acc, flow_mean=VE$flow_mean, 
                      flow_mean_acc=flow_mean_acc, Kin=VE$Kinet, 
                      Kin_acc=Kin_acc_fin, N_loss=N_loss, Fish_mort=VE$Fish_mort, 
                      Fish_mort_acc=VE$Fish_mort_acc, Selec=accessibility, 
                      N_loss_acc=N_loss_acc, F_loss=F_loss, 
                      F_loss_acc=F_loss_acc, row.names=seq_tl)
  
  
  ## B) Merge data ET.Main and catch to create the input list
  dat<-list(ET_Main=ET_Main)
  
  return(dat)
}
