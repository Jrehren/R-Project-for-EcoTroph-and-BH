mf.diagnosis.adapted=function(x,ET_Main,TL_out,fleet,n.fleet,B.Input,Beta,TopD,FormD,TLpred,n.TL,range.TLpred,lim.high.TL,range.highTL,Fish_mort, Fish_mort_acc){
  
  #---------------------------------- LOOP -------------------------------------
  
  #  Setting the conditions of the while loop. SC is the ratio btw. output ---
  #  kinetics of each iteration of the while loop and the temporal kinetics.
  #  SC2= same for production.----
  
  x[['SC']] <- 3
  x[['SC2']] <- 3
  
  It=0
  while (!((x[['SC2']] < 1e-6) && (x[['SC']] < 1e-6))) {
    
  # 1.) CALCULATING BIOMASS (Gascuel et al., 2011; Equation A8) ----
    x[['BIOM_MF']] <- x[['Prod_MF']]/x[['Kin_MF']] 
    x[['BIOM_MF']] <- x[['BIOM_MF']][-1] ##sum of BIOM_MF for TL>=2 in the biomass input control equation
    x[['BIOM_MF_acc']] <- x[['Prod_MF_acc']]/x[['Kin_MF_acc']] 
    
  # 2.) CALCULATING PRODUCTION ----
   ## A.) CALCULATING BIOMASS INPUT CONTROL (Gascuel, 2005; Equation 6.2) ----
    x[['Prod_MF_TMP']] <- x[['Prod_MF']]
    if(B.Input){
      x[['Prod_MF']][1] <- (1 - Beta) * ET_Main[1, "P"]+
      Beta * ET_Main[1, "P"] * (sum(x[['BIOM_MF']])/sum(ET_Main[-1,'B']))
    # As production of the first trophic level has changed, so does the flow
    # and mean_flow. Because the production at trophic level 1 is also the flow
    # and therewith the mean_flow at trophic level 1. I just have to change
    # the flow and mean_flow according to the changes in production
    x[["flow"]][1]<-x[['Prod_MF']][1]
    x[["mean_flow"]][1]<-x[['Prod_MF']][1]
    }
    
   ## B.) CALCULATING PRODUCTION FOR TL>1 (Gascuel et al., 2011; Equation A4 & A6 & A9) ----
    
    # I need to recalculate TL2, because if the biomass input 
    # changes it will effect the production/flow at TL=1 
    # and therewith the flow at TL=2
    
    # Gascuel et al. 2011; Equation A4 (The general flow equation)
    x[["flow"]][2]<- x[["flow"]][1]*exp(-(ET_Main[(1), "N_loss"]+
                      (-log(1-x[["Fish_mort"]][1]/x[['Kin_MF']][1])))*
                                (TL_out[2] - TL_out[1]))
    
    # Huberts implementation
    # x[["flow"]][2] <- x[["flow"]][1]*
    #   exp(-ET_Main[(1), "N_loss"]*1)*(1-(x[["Fish_mort"]][1]/x[['Kin_MF']][1]))
    
    # here we start to calculate the flow of trophic levels from 2.1 on and we
    # also recalculate the mean_flow and the production from trophic level 2 on
    
    for (compteur in 2:(n.TL-1)) {
      # flow : Gascuel et al. 2011; Equation A4 (The general flow eq.)
      x[["flow"]][compteur+1] <- x[["flow"]][compteur] * exp(-(ET_Main[(compteur), "N_loss"]+
                                (-log(1-x[["Fish_mort"]][compteur]/x[['Kin_MF']][compteur])))*
                                (TL_out[compteur+1] - TL_out[compteur]))
      
      # Huberts implementation
      # x[["flow"]][compteur+1] <- x[["flow"]][compteur]*
      #     exp(-ET_Main[(compteur), "N_loss"]*0.1)*(1-(x[["Fish_mort"]][compteur]/x[['Kin_MF']][compteur]))
         
      for (compteur2 in 2:n.TL){
        # mean flow : Gascuel et al. 2011; Equation A6 (Integrated version of the flow eq.)
       
        # Huberts code
        # x[["flow_mean"]][compteur2] <- (1/0.1) * x[["flow"]][compteur2] *
        #   (exp(-ET_Main[(compteur2), "N_loss"]*0.1)*(1-(x[["Fish_mort"]][compteur2]/x[['Kin_MF']][compteur2]))-1)/
        #   ((-ET_Main[(compteur2), "N_loss"]) + (1/0.1)*log(1-(x[["Fish_mort"]][compteur2]/x[['Kin_MF']][compteur2])))
        # 
        x[["flow_mean"]][compteur2] <- x[["flow"]][compteur2]*
                                       (1-(exp(-(ET_Main[(compteur2), "N_loss"]+
                                      (log(1-x[["Fish_mort"]][compteur2]/x[['Kin_MF']][compteur2])))*
                                      (0.1))))/
                                      ((ET_Main[(compteur2), "N_loss"]+
                                          (log(1-x[["Fish_mort"]][compteur2]/x[['Kin_MF']][compteur2])))*
                                         (0.1))

        # production : Gascuel et al., 2011; Equation A9 (Production at a given TL)
        x[['Prod_MF']][compteur2] <- x[["flow_mean"]][compteur2] * 0.1
      }
    }
    
   ## C.) CALCULATING ACCESSIBLE PRODUCTION ----
      # We must adjust the accessible flow/production, given the changes in the
      # biomass input that is directly effecting the first trophic level. HOwever,
      # if we have no fishing at trophic level one, we also have no accessible
      # production/flow. So in this case we adjust the accessible flow of TL=2

    x[['flow_acc']][2] <- ET_Main[2, "flow_acc"] * x[['flow']][2]/ET_Main[2, "flow"]

      
      ### CALCULATE THE PRODUCTION FROM TL>1
      # here we start to calculate the flow of trophic levels from 2.1 on and we
      # also recalculate the mean_flow and the production from trophic level 2 on 
    
         for (compteur in 2:(n.TL-1)){
       # Accessible flow : Gascuel et al. 2011; Equation A4 (The general flow eq.)
        x[["flow_acc"]][compteur+1] <- x[["flow_acc"]][compteur]*
          exp(-(ET_Main[(compteur), "N_loss_acc"]+
                  (-log(1-x[["Fish_mort_acc"]][compteur]/x[['Kin_MF_acc']][compteur])))*
                (TL_out[compteur+1] - TL_out[compteur]))

        # Huberts code       
        # x[["flow_acc"]][compteur+1] <- x[["flow_acc"]][compteur]*
        #   exp(-ET_Main[(compteur), "N_loss_acc"]*0.1)*(1-(x[["Fish_mort_acc"]][compteur]/x[['Kin_MF_acc']][compteur]))
        # 
        
        for (compteur2 in 2:n.TL){
        # mean accessible flow : Gascuel et al. 2011; Equation A6 (Integrated version of the flow eq.)
          # Huberts code
          # x[["flow_mean_acc"]][compteur2] <- (1/0.1) * x[["flow_acc"]][compteur2] *
          #   ((exp((-ET_Main[(compteur2), "N_loss_acc"])*0.1))*
          # (1-(x[["Fish_mort_acc"]][compteur2]/x[['Kin_MF_acc']][compteur2]))-1)/
          #   (-ET_Main[(compteur2), "N_loss_acc"]+(1/0.1)*
          #      log(1-(x[["Fish_mort_acc"]][compteur2]/x[['Kin_MF_acc']][compteur2])))

          
          x[["flow_mean_acc"]][compteur2] <- x[["flow_acc"]][compteur2] *
            (1-(exp(-(ET_Main[(compteur2), "N_loss_acc"]+
              (log(1-x[["Fish_mort_acc"]][compteur2]/x[['Kin_MF_acc']][compteur2])))*
                      (0.1))))/
            ((ET_Main[(compteur2), "N_loss_acc"]+
                (log(1-x[["Fish_mort_acc"]][compteur2]/x[['Kin_MF_acc']][compteur2])))*
               (0.1))


        # production : Gascuel et al., 2011; Equation A9 (Production at a given TL)
          x[['Prod_MF_acc']][compteur2] <- x[["flow_mean_acc"]][compteur2] * 0.1  
        }
      }
    
  # 3.) RECALCULATING BIOMASS (Gascuel et al., 2011; Equation A8)-----  
    x[['BIOM_MF']] <- x[['Prod_MF']]/x[['Kin_MF']] 
    x[['BIOM_MF_acc']] <- x[['Prod_MF_acc']]/x[['Kin_MF_acc']]
    
  # 4.) CALCULATING KINETICS FROM TOP-DOWN EQUATION (Gascuel et al., 2011; Equation A13)----
   ## A.) FOR TL=1----
    x[['Kin_MF']][1] <- (ET_Main[1, "Kin"] - ET_Main[1, "Fish_mort"]) * (1 + TopD[1]* 
                       (sum(x[['BIOM_MF']][TL_out[TL_out>=2&TL_out<=2.3]])^FormD[1] - 
                        sum(ET_Main[TL_out[TL_out>=2&TL_out<=2.3], "B"])^FormD[1])/
                       (sum(ET_Main[TL_out[TL_out>=2&TL_out<=2.3], "B"])^FormD[1]))+ 
                        x[["Fish_mort"]][1]
    
    
   ## B.) FOR TL 2:38 ----
    x[['Kin_MF']][2:lim.high.TL] <-sapply(2:lim.high.TL,a13.eq,ET_Main,
                                          x[['BIOM_MF']], x[["Fish_mort"]], TopD,
                                          FormD,range.TLpred)
    
    x[['Kin_MF_acc']][2:lim.high.TL] <-sapply(2:lim.high.TL,a13.eq.ac,ET_Main,
                                              x[['BIOM_MF']], x[["Fish_mort_acc"]],
                                              TopD,FormD,range.TLpred) 
    
    
   ## C.) FOR TL >38 (Linear regression analysis to approximate the kinetics for high TL) ----
    if (sum(TopD)==0 || sum(x[["mf"]])==0) # I test if there is top-down control &
    {                                      # if it is not the virgin state
      for (compteur in 39:52){             # If there is no top-down control &
                                           # it is the reference state, I do not
                                           # use the linear regression for calculation
        x[['Kin_MF']][compteur] <- (ET_Main$Kin[compteur]- 
                                   ET_Main$Fish_mort[compteur])+
                                   x[["Fish_mort"]][compteur]
        
        x[['Kin_MF_acc']][compteur] <- (ET_Main$Kin_acc[compteur]- 
                                        ET_Main$Fish_mort_acc[compteur])+
                                        x[["Fish_mort_acc"]][compteur]
      }}
    
    else {                                 # else means there is top-down control &
                                           # so we have to apply the regression analysis
      x[['Kin_MF']][(lim.high.TL+1):n.TL]=sapply((lim.high.TL+1):n.TL,regPB,
                                                x[['Kin_MF']],TL_out,range.highTL)
      x[['Kin_MF_acc']][(lim.high.TL+1):n.TL]=sapply((lim.high.TL+1):n.TL,regPB.ac,
                                                 x[['Kin_MF_acc']],TL_out,range.highTL)
    }
    
  # 5.) CALCULATING THE DIFFERENCE BETWEEN TEMPORAL PRODUCTION AND KINETICS AND----
    #    OUTPUT OF THE ITERATION+1 ----
    # x[['SC2']] <- round((sum(x[['Prod_MF']]) - sum(x[['Prod_MF_TMP']])) * 1E3)
    x[['SC']] <- abs(1-(sum(x[['Kin_MF']])/sum(x[['TEMP_Kin']])))
    x[['SC2']] <- abs(1-(sum(x[['Prod_MF']])/sum(x[['Prod_MF_TMP']])))
    # x[['SC']] <- round(sum(x[['Kin_MF']])/sum(x[['TEMP_Kin']]),3)
    
    It=It+1
    
    x[['TEMP_Kin']] <- x[['Kin_MF']]
    
  }  
  
  #---------------------------- END OF THE ITERATIONS --------------------------
  
  #-------------------------- CALCULATING OTHER OUTPUTS ------------------------
  
  # 1.) EXPLOITATION RATE AND FISHING LOSS RATE----
  x[['E']]=x[['Fish_mort']]/x[['Kin_MF']]
  x[['E_acc']]=x[['Fish_mort_acc']]/x[['Kin_MF_acc']]
  
  x[['F_loss']]=-log(1-x[['E']])
  x[['F_loss_acc']]=-log(1-x[['E_acc']])
  
  # 2.) TOTAL BIOMASS FOR EACH MULTIPLIER ----
  TOT_B <-sum(x[['BIOM_MF']][-1]) #without TL=1 
  TOT_B_acc <- sum(x[['BIOM_MF_acc']][-1])#without TL=1
  
  # 3.) TOTAL PRODUCTION FOR EACH MULTIPLIER ----
  Pred_B <- sum(x[['BIOM_MF']][as.numeric(names(TL_out[TL_out>=TLpred]))])
  TOT_P <- sum(x[['Prod_MF']])  
  TOT_P_acc <- sum(x[['Prod_MF_acc']]) 
  Pred_P <- sum(x[['Prod_MF']][as.numeric(names(TL_out[TL_out>=TLpred]))])
  
  # 4.) CATCHES PER TROPHIC LEVEL AND FLEET ----
  Catches <-x[['BIOM_MF_acc']][-1]*x[['Fish_mort_acc']][-1]

  # 5.) TOTAL CATCHES PER TROPHIC LEVEL ----
  Y <- sum(Catches)
  Pred_Y <- sum(Catches[as.numeric(names(Catches))%in%TL_out[TL_out>=TLpred]],na.rm=T)
  
  # 6.) TOTAL BIOMASS AND TOTAL PRODUCTION RELATIVE TO REFERENCE STATE ----
  R_TOT_B <- TOT_B/sum(ET_Main[-1,'B'])
  R_TOT_B_acc <- TOT_B_acc/sum(ET_Main[-1,'B_acc'])
  R_Pred_B <- Pred_B/sum(ET_Main[as.numeric(names(TL_out[TL_out>=TLpred])),'B'],na.rm=T)
  R_TOT_P <- TOT_P/sum(ET_Main[,'P'])
  R_TOT_P_acc <- TOT_P_acc/sum(ET_Main[,'P_acc'])
  R_Pred_P <- Pred_P/sum(ET_Main[as.numeric(names(TL_out[TL_out>=TLpred])),'P'],na.rm=T)
  
  # 7.) TOTAL CATCHES PER TROPHIC LEVEL ----
    x[['Catches.tot']]=c(0,Catches)
  names(x[['Catches.tot']])=TL_out
  
  
  # 8.) TROPHIC LEVEL OF BIOMASS AND CATCHES PER MULTIPLIER ----
  TL_TOT_B <- sum(x[['BIOM_MF']][-1]*TL_out[-1])/TOT_B
  TL_TOT_B_acc <-sum(x[['BIOM_MF_acc']][-1]*TL_out[-1])/TOT_B_acc 
  
    TL_Y<-sum(Catches*TL_out[-1])/Y

  
  # 9.) TOTAL CATCHES PER TROPHIC LEVEL AND FLEET ----
  
  x[['Catches']]=list()
  for(i in 1:n.fleet){
    x[['Catches']][[fleet[i]]]=c(0,Catches*x[['mf']][[i]]* Fish_mort[[i]][-1]/x[['Fish_mort']][-1])
    names(x[['Catches']][[fleet[i]]])=TL_out
  }
  


  #----------------------------- END OF CALCULATIONS ---------------------------
  
  #------------------------------- BINDING TO LIST -----------------------------
  ET_Main_diagnose<- list(TOT_B=TOT_B,TOT_B_acc=TOT_B_acc,Pred_B=Pred_B, 
                          TOT_P=TOT_P,TOT_P_acc=TOT_P_acc,Pred_P=Pred_P,Y=Y,
                          Pred_Y=Pred_Y,R_TOT_B=R_TOT_B,R_TOT_B_acc=R_TOT_B_acc,
                          R_Pred_B=R_Pred_B,R_TOT_P=R_TOT_P,
                          R_TOT_P_acc=R_TOT_P_acc,R_Pred_P=R_Pred_P,
                          TL_TOT_B=TL_TOT_B, TL_TOT_B_acc=TL_TOT_B_acc,TL_Y=TL_Y)
  
  
  x[['ET_Main_diagnose']]=ET_Main_diagnose
  return(x)
}

