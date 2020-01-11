# ---------------------------------------------------------------------
# ET biomass modelling with catch
# ---------------------------------------------------------------------
rm(list=ls())
#-------------------------------------------------------------------------------
library(RPostgreSQL)
library(parallel)
library(reshape2)
library(tidyverse) 
# library(car)
library(EcoTroph)
# ---------------------------------------------------------------------
######################## FUNCTIONS AND PARAMETER DEFINITION ####################
# ------------------------------------------------------------------------------
# Speed of biomass flows: TL + SST Gascuel 2008
# ------------------------------------------------------------------------------
kin_f <- function(tl,sst){return(20.19*(tl^(-3.26))*exp(.041*sst))}

# XXXXXXX could we put this in a seperate function, in order to separate
#         the functions from the analysis?

# -------------------------------------------------------------------
# top down function
# -------------------------------------------------------------------
# dt_base<-dt_fixed$`2000`
# a_td<-0.4
# b_td<-0.5
# bio<-dt_base$bio 
# bio_ref<-dt_base$bio_v 
# tl<-20
# dim(dt_base)
# tl <= 5.6 (=37 dans seq_tl)
topd_classic <- function(a_td,b_td,bio,bio_ref,tl) {
  topd_para <- a_td * (sum(bio[(tl+8):(tl+13)])^b_td - sum(bio_ref[(tl+8):(tl+13)])^b_td) / (sum(bio_ref[(tl+8):(tl+13)])^b_td)
  return(topd_para)
}
# tl for high trophic level  tl > 5.6 (=37 dans seq_tl)
topd_high_tl <- function(a_td,b_td,bio,bio_ref,tl) {
  para_56<-topd_classic(a_td,b_td,bio,bio_ref,37)
  para_69<-0
  return((-para_56/13)*tl + (50*para_56)/13)
}
# topd_high_tl(a_td,b_td,bio,bio_ref,50) 
# topd_classic(a_td,b_td,bio,bio_ref,30)
# cell_no<-10612
# dt_input<-select_data(cell_no)
# dt_base<-dt_input[[56]]
# res<-lapply(dt_input[56],core_function)

# XXXXXXXXX he does create a top-down control function seperately, which is 
#           probably the sub.mf function. We need to save this seperately 

################################################################################
# ------------------------------------------------------------------------------
#################################### CORE FUNCTION #############################
# ------------------------------------------------------------------------------
# base function to calculate production and biomass
core_function <-function(dt_base) {
  # top down parameters
  a_td<-0.4
  b_td<-0.5
  
  # XXXXXXX here he defines the parameters, that is something that we have to 
  #         homogenize
  
  iteration=0
  # seq_tl
  seq_tl <- seq(2,7,by=0.1)
  dt_add <-as.data.frame(matrix(nrow = nrow(dt_base), ncol = 15))
  colnames(dt_add)<-c("flow","flow_mean","flow_mean_v" ,"prod","prod_v", "bio","bio_v", "kin", 
                      "f_thau_ref", "f_thau","flow_optm","prod_optm","bio_optm","kin_optm","supported")
  dt_base<-cbind(dt_base,dt_add)
  
  # XXXXXX Is that the final data base he is using in the diagnosis run?
  
  # indicate if production support or not catch 
  # 0 = prodcution support catches
  # 1 = catch > prod => prod = 0
  # 2 = catch>catc ==> bio = 0; prod=0
  dt_base$supported<-0
  
  # initialization total at TL=2 
  dt_base$flow[1] <- dt_base$flow2[1]
  # Firt iteration
  # for (tl_ite in seq(tl_r,1,-1)) {
  
  # Step 1: Virgin production flow
  for (tl_ite in 1:(length(seq_tl)-1)) { # tl [-1] because we start the iterative process at tl=2.1
    # production flow
    dt_base$flow[tl_ite+1] <- dt_base$flow[tl_ite] * dt_base$te[tl_ite]^0.1
    # mean production flow
    dt_base$flow_mean[tl_ite] <- (1/0.1) * dt_base$flow[tl_ite] * (1-dt_base$te[tl_ite]^0.1)/(-log(dt_base$te[tl_ite]))
    # we store virgin production virgin mean production flow
    dt_base$flow_mean_v[tl_ite] <- dt_base$flow_mean[tl_ite]
    # production
    dt_base$prod[tl_ite] <- dt_base$flow_mean[tl_ite] * 0.1
    # we store virgin production
    dt_base$prod_v[tl_ite] <- dt_base$prod[tl_ite]
    # biomass = flow/kin
    dt_base$bio[tl_ite] <- dt_base$flow_mean[tl_ite] * 0.1 / dt_base$kin_v[tl_ite]
    # we store virgin biomass
    dt_base$bio_v[tl_ite] <- dt_base$bio[tl_ite]
    # we evualate if catch < biomass otherwise biomass in the next trophic levels are null 
    if ((dt_base$catch[tl_ite])>(dt_base$bio[tl_ite])) {
      dt_base$flow[tl_ite:(length(seq_tl)-1)]<-0
      dt_base$flow_mean[tl_ite:(length(seq_tl)-1)]<-0
      dt_base$prod[tl_ite:(length(seq_tl)-1)]<-0
      dt_base$bio[tl_ite:(length(seq_tl)-1)]<-0
      dt_base$kin[tl_ite:(length(seq_tl)-1)]<-0
      dt_base$f_thau[tl_ite:(length(seq_tl)-1)]<-NA
      # to track this situation an cancel the other steps
      dt_base$supported<-2
    } else {
      # Fihsing mortality reference: catch/biomass
      dt_base$f_thau_ref[tl_ite] <- dt_base$catch[tl_ite]/dt_base$bio[tl_ite]
      # Fishing mortality
      dt_base$f_thau[tl_ite] <- dt_base$f_thau_ref[tl_ite]
      # Kinetics = kin_ref - Fishing morta_ref + Fishing mortality
      dt_base$kin[tl_ite] <- dt_base$kin_v[tl_ite] - dt_base$f_thau_ref[tl_ite] + dt_base$f_thau[tl_ite]
    }
  }
  
  # XXXXXXXXXX Step 1 is the virgin ecosystem flow and should be outside of the 
  #            iterative run, as in my case I have a different way of calculating
  #            virgin ecosystem. Or at least we need to check what is overlapping
  #            here. 
  # --------------------------------------------------------------
  # if biomass does not support catch, in this case everything is fucked
  if (dt_base$supported[1]!=2) {
    
    
    # Step 2: we recalculate the flow using (1-Y/P) and B depends on this new flow and kin virgin
    for (tl_ite in 1:(length(seq_tl)-1)) { # tl [-1] because we start the iterative process at tl=2.1
      if ((dt_base$catch[tl_ite])>(dt_base$prod[tl_ite])) {
        dt_base$flow[tl_ite:(length(seq_tl)-1)]<-0
        dt_base$flow_mean[tl_ite:(length(seq_tl)-1)]<-0
        dt_base$prod[tl_ite:(length(seq_tl)-1)]<-0
        dt_base$bio[tl_ite:(length(seq_tl)-1)]<-0
        dt_base$kin[tl_ite:(length(seq_tl)-1)]<-0
        dt_base$f_thau[tl_ite:(length(seq_tl)-1)]<-NA
        dt_base$supported[tl_ite:(length(seq_tl)-1)]<-1
      } else {
        # biomass flow
        dt_base$flow[tl_ite+1] <- dt_base$flow[tl_ite] * dt_base$te[tl_ite]^0.1 * 
          (1-(dt_base$catch[tl_ite]/dt_base$prod[tl_ite]))
        # mean biomass flow
        dt_base$flow_mean[tl_ite] <- (1/0.1) * dt_base$flow[tl_ite] *
          ((dt_base$te[tl_ite]^0.1)*(1-(dt_base$catch[tl_ite]/dt_base$prod[tl_ite]))-1)/
          (log(dt_base$te[tl_ite]) + (1/0.1)*log(1-(dt_base$catch[tl_ite]/dt_base$prod[tl_ite])))
        # production
        dt_base$prod[tl_ite] <- dt_base$flow_mean[tl_ite] * 0.1
        # biomass
        dt_base$bio[tl_ite] <- (dt_base$flow_mean[tl_ite] * 0.1) / dt_base$kin[tl_ite]
        # Fishing mortality
        dt_base$f_thau[tl_ite] <- dt_base$catch[tl_ite]/dt_base$bio[tl_ite]
        # Kinetics = (kin_ref - Fishing morta_ref) * (1 + topdown) + Fishing mortality
        # top down function before tl = 5.6 (line number 37)
        # top down effect as described in Gascuel et al., 2011
        if (tl_ite<=37) {
          dt_base$kin[tl_ite] <- (dt_base$kin_v[tl_ite] - dt_base$f_thau_ref[tl_ite]) *
            (1+topd_classic(a_td,b_td,dt_base$bio,dt_base$bio_v,tl_ite)) +
            dt_base$f_thau[tl_ite]
        } else {
          # top down function after tl = 5.6 (line number 37) - 
          #linear decrease from topdown (tl=5.6) at tl=5.6 to 0 at tl=6.9
          dt_base$kin[tl_ite] <- (dt_base$kin_v[tl_ite] - dt_base$f_thau_ref[tl_ite]) *
            (1+topd_high_tl(a_td,b_td,dt_base$bio,dt_base$bio_v,tl_ite)) +
            dt_base$f_thau[tl_ite]
        }
      }
    }
    # --------------------------------------------------------------
    # Step 3: we recalculate the flow using (1-Y/P) and B depends on this new flow and on kin modfied
    for (tl_ite in 1:(length(seq_tl)-1)) { # tl [-1] because we start the iterative process at tl=2.1
      if ((dt_base$catch[tl_ite])>(dt_base$prod[tl_ite])) {
        dt_base$flow[tl_ite:(length(seq_tl)-1)]<-0
        dt_base$flow_mean[tl_ite:(length(seq_tl)-1)]<-0
        dt_base$prod[tl_ite:(length(seq_tl)-1)]<-0
        dt_base$bio[tl_ite:(length(seq_tl)-1)]<-0
        dt_base$kin[tl_ite:(length(seq_tl)-1)]<-0
        dt_base$f_thau[tl_ite:(length(seq_tl)-1)]<-NA
        dt_base$supported[tl_ite:(length(seq_tl)-1)]<-1
      } else {
        # biomass flow
        dt_base$flow[tl_ite+1] <- dt_base$flow[tl_ite] * dt_base$te[tl_ite]^0.1 * 
          (1-(dt_base$catch[tl_ite]/dt_base$prod[tl_ite]))
        # mean biomass flow
        dt_base$flow_mean[tl_ite] <- (1/0.1) * dt_base$flow[tl_ite] *
          ((dt_base$te[tl_ite]^0.1)*(1-(dt_base$catch[tl_ite]/dt_base$prod[tl_ite]))-1)/
          (log(dt_base$te[tl_ite]) + (1/0.1)*log(1-(dt_base$catch[tl_ite]/dt_base$prod[tl_ite])))
        # production
        dt_base$prod[tl_ite] <- dt_base$flow_mean[tl_ite] * 0.1
        # biomass
        dt_base$bio[tl_ite] <- (dt_base$flow_mean[tl_ite] * 0.1) / dt_base$kin[tl_ite]
        # Fishing mortality
        dt_base$f_thau[tl_ite] <- dt_base$catch[tl_ite]/dt_base$bio[tl_ite]
        # Kinetics = (kin_ref - Fishing morta_ref) * (1 + topdown) + Fishing mortality
        # top down function before tl = 5.6 (line number 37)
        # top down effect as described in Gascuel et al., 2011
        if (tl_ite<=37) {
          dt_base$kin[tl_ite] <- (dt_base$kin_v[tl_ite] - dt_base$f_thau_ref[tl_ite]) *
            (1+topd_classic(a_td,b_td,dt_base$bio,dt_base$bio_v,tl_ite)) +
            dt_base$f_thau[tl_ite]
        } else {
          # top down function after tl = 5.6 (line number 37) - 
          #linear decrease from topdown (tl=5.6) at tl=5.6 to 0 at tl=6.9
          dt_base$kin[tl_ite] <- (dt_base$kin_v[tl_ite] - dt_base$f_thau_ref[tl_ite]) *
            (1+topd_high_tl(a_td,b_td,dt_base$bio,dt_base$bio_v,tl_ite)) +
            dt_base$f_thau[tl_ite]
        }
      }
      # --------------------------------------------------------------
      dt_base$bio_optm[tl_ite] <- dt_base$bio[tl_ite]
      dt_base$prod_optm[tl_ite] <- dt_base$prod[tl_ite]
      dt_base$flow_optm[tl_ite] <- dt_base$flow_mean[tl_ite]
      dt_base$kin_optm[tl_ite] <- dt_base$kin[tl_ite]
    }
    
    mini_v1=1
    mini_v2=1
    mini_v3=1
    mini_v4=1
    iteration=0
    supported=0
    while (!((mini_v1<1e-6) && (mini_v2<1e-6) && (mini_v3<1e-6) && (mini_v4<1e-6)) & (supported<50)) {
      for (tl_ite in 1:(length(seq_tl)-1)) {
        dt_base$bio_optm[tl_ite] <- dt_base$bio[tl_ite]
        dt_base$prod_optm[tl_ite] <- dt_base$prod[tl_ite]
        dt_base$flow_optm[tl_ite] <- dt_base$flow_mean[tl_ite]
        dt_base$kin_optm[tl_ite] <- dt_base$kin[tl_ite]
        if ((dt_base$catch[tl_ite])>(dt_base$prod[tl_ite])) {
          dt_base$flow[tl_ite:(length(seq_tl)-1)]<-0
          dt_base$flow_mean[tl_ite:(length(seq_tl)-1)]<-0
          dt_base$prod[tl_ite:(length(seq_tl)-1)]<-0
          dt_base$bio[tl_ite:(length(seq_tl)-1)]<-0
          dt_base$kin[tl_ite:(length(seq_tl)-1)]<-0
          dt_base$f_thau[tl_ite:(length(seq_tl)-1)]<-NA
          dt_base$supported[tl_ite:(length(seq_tl)-1)]<-1
        } else {
          # biomass flow
          dt_base$flow[tl_ite+1] <- dt_base$flow[tl_ite] * dt_base$te[tl_ite]^0.1 * 
            (1-(dt_base$catch[tl_ite]/dt_base$prod[tl_ite]))
          # mean biomass flow
          dt_base$flow_mean[tl_ite] <- (1/0.1) * dt_base$flow[tl_ite] *
            ((dt_base$te[tl_ite]^0.1)*(1-(dt_base$catch[tl_ite]/dt_base$prod[tl_ite]))-1)/
            (log(dt_base$te[tl_ite]) + (1/0.1)*log(1-(dt_base$catch[tl_ite]/dt_base$prod[tl_ite])))
          # production
          dt_base$prod[tl_ite] <- dt_base$flow_mean[tl_ite] * 0.1
          # biomass
          dt_base$bio[tl_ite] <- (dt_base$flow_mean[tl_ite] * 0.1) / dt_base$kin[tl_ite]
          # Fishing mortality
          dt_base$f_thau[tl_ite] <- dt_base$catch[tl_ite]/dt_base$bio[tl_ite]
          # Kinetics = (kin_ref - Fishing morta_ref) * (1 + topdown) + Fishing mortality
          if (tl_ite<=37) {
            dt_base$kin[tl_ite] <- (dt_base$kin_v[tl_ite] - dt_base$f_thau_ref[tl_ite]) *
              (1+topd_classic(a_td,b_td,dt_base$bio,dt_base$bio_v,tl_ite)) +
              dt_base$f_thau[tl_ite]
          } else {
            dt_base$kin[tl_ite] <- (dt_base$kin_v[tl_ite] - dt_base$f_thau_ref[tl_ite]) *
              (1+topd_high_tl(a_td,b_td,dt_base$bio,dt_base$bio_v,tl_ite)) +
              dt_base$f_thau[tl_ite]
          }
        }
        # --------------------------------------------------------------
      }
      mini_v1<-abs(1-(sum(dt_base$bio_optm[-51])/sum(dt_base$bio[-51])))
      mini_v2<-abs(1-(sum(dt_base$kin_optm[-51])/sum(dt_base$kin[-51])))
      mini_v3<-abs(1-(sum(dt_base$flow_optm[-51])/sum(dt_base$flow_mean[-51])))
      mini_v4<-abs(1-(sum(dt_base$prod_optm[-51])/sum(dt_base$prod[-51])))
      supported<-sum(dt_base$supported[-51]) # if production does not support catch
      iteration=iteration+1
    }
# XXXXXXXXXXXX Somehow he manages to not put the while loop at the beginning
#              I am not sure how this works in fact. 
    
  } 
  # XXXXXXX He needs to have two steps one which is, when production supports the
  #        catch and one if it does not do this. It seems that the iteration process
  #        is not including the construction of the virgin ecosystem. So we could
  #        seperate this. 
  
  # --------------------------------------------------------------------
  # result<-list(dt_base,iteration)
  # result<-dt_base
  # return(result)
  dt_base$iteration<-iteration
  result<-as.data.frame(dt_base[-51,c("cell_no","year","tl","catch","prod","prod_v","bio","bio_v","supported","f_thau","iteration")])
  colnames(result)<-c("cell_no","year","tl","catch","production","production_v","biomass","biomass_v","supported","f_morta","iteration")
  return(result)
}