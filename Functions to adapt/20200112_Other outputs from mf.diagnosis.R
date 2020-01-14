################################################################################
############################ CALCULATING OTHER OUTPUTS #########################
################################################################################

# 1.) Explanation ----
# Here, I calculate other outputs used in the EcoTroph analysis. This is the 
# ET.Main data frame. I externalize this, because the original code runs with a 
# lot of if clauses to check for the amount of fleets, etc. There are also a lot
# of parameters that might not be used by Hubert and will only slow down the analysis

ET.summary=function(dat, dat.ETMain, TLpred=3.5, fleet_names="catch.1", 
                    Fish_mort, Fish_mort_acc){

# Some variables used in the calculations

  seq_tl<-c(1,seq(2,7, 0.1))
  TLpredpos<-which(seq_tl==TLpred)
  
TL_out <- as.numeric(rownames(dat.ETMain))
fleet=fleet_names ; n.fleet=length(fleet)
  
# 2.) TOTAL BIOMASS FOR EACH MULTIPLIER ----
TOT_B <-sum(dat[['BIOM_MF']][-1]) #without TL=1 
TOT_B_acc <- sum(dat[['BIOM_MF_acc']][-1])#without TL=1

# 3.) TOTAL PRODUCTION FOR EACH MULTIPLIER ----
Pred_B <- sum(dat[['BIOM_MF']][TLpredpos:52])
TOT_P <- sum(dat[['Prod_MF']])  
TOT_P_acc <- sum(dat[['Prod_MF_acc']][-1]) 
Pred_P <- sum(dat[['Prod_MF']][TLpredpos:52])

# 4.) CATCHES PER TROPHIC LEVEL AND FLEET ----
Catches <-dat[['BIOM_MF_acc']][-1]*dat[['Fish_mort_acc']][-1]

# 5.) TOTAL CATCHES PER TROPHIC LEVEL ----
Y <- sum(Catches)
Pred_Y <- sum(Catches[TLpredpos:52],na.rm=T)

# 6.) TOTAL BIOMASS AND TOTAL PRODUCTION RELATIVE TO REFERENCE STATE ----
R_TOT_B <- TOT_B/sum(dat.ETMain[-1,'B'])
R_TOT_B_acc <- TOT_B_acc/sum(dat.ETMain[-1,'B_acc'])
R_Pred_B <- Pred_B/sum(dat.ETMain[TLpredpos:52,'B'],na.rm=T)
R_TOT_P <- TOT_P/sum(dat.ETMain[,'P'])
R_TOT_P_acc <- TOT_P_acc/sum(dat.ETMain[,'P_acc'][-1])
R_Pred_P <- Pred_P/sum(dat.ETMain[TLpredpos:52,'P'],na.rm=T)

# 7.) TOTAL CATCHES PER TROPHIC LEVEL ----
dat[['Catches.tot']]=c(0,Catches)
names(dat[['Catches.tot']])=TL_out


# 8.) TROPHIC LEVEL OF BIOMASS AND CATCHES PER MULTIPLIER ----
TL_TOT_B <- sum(dat[['BIOM_MF']][-1]*TL_out[-1])/TOT_B
TL_TOT_B_acc <-sum(dat[['BIOM_MF_acc']][-1]*TL_out[-1])/TOT_B_acc 

TL_Y<-sum(Catches*TL_out[-1])/Y


# 9.) TOTAL CATCHES PER TROPHIC LEVEL AND FLEET ----

dat[['Catches']]=list()
for(i in 1:n.fleet){
  dat[['Catches']][[fleet[i]]]=c(0,Catches*dat[['mf']][[i]]* Fish_mort[[i]][-1]/dat[['Fish_mort']][-1])
  names(dat[['Catches']][[fleet[i]]])=TL_out
}



#----------------------------- END OF CALCULATIONS ---------------------------

#------------------------------- BINDING TO LIST -----------------------------
ET_Main_diagnose<- data.frame(TOT_B=TOT_B,TOT_B_acc=TOT_B_acc,Pred_B=Pred_B, 
                        TOT_P=TOT_P,TOT_P_acc=TOT_P_acc,Pred_P=Pred_P,Y=Y,
                        Pred_Y=Pred_Y,R_TOT_B=R_TOT_B,R_TOT_B_acc=R_TOT_B_acc,
                        R_Pred_B=R_Pred_B,R_TOT_P=R_TOT_P,
                        R_TOT_P_acc=R_TOT_P_acc,R_Pred_P=R_Pred_P,
                        TL_TOT_B=TL_TOT_B, TL_TOT_B_acc=TL_TOT_B_acc,TL_Y=TL_Y)


return(ET_Main_diagnose)

}
