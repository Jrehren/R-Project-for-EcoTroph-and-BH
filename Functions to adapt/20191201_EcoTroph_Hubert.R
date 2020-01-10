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
# Speed of biomass flows: TL + SST Gascuel 2008
kin_f <- function(tl,sst){return(20.19*(tl^(-3.26))*exp(.041*sst))}
# ---------------------------------------------------------------------
# model du pontavice et al., in review
source("model_te.R") 
# ---------------------------------------------------------------------
# Data
data_npp<-"model_ubc.pp_gfdl_esm2m_85"
data_sst<-"model_ubc.sst_gfdl_esm2m_85"
# --------------------------------------------------------------------
# linking mt database
drv <- dbDriver("PostgreSQL")
liaisons2 <- dbConnect(drv, dbname="environnement",host="halieut.agrocampus-ouest.fr",port=5432,user="hubert",password="dup098!")
liaisons1 <- dbConnect(drv, dbname="world_mapping",host="sirs.agrocampus-ouest.fr",port=5432,user="hubert",password="dup098!")
# ---------------------------------------------------------------------
cell_world <- as.data.frame(dbGetQuery(liaisons2,"select distinct cell_no from geo.world_grid1x1_v2 inner join geo.biogeo3 using(cell_no)")) 
# ---------------------------------------------------------------------
# output model_te
TE <- function(biome2,sst) {
  if (biome2 == "polar"){b=0; d=0; intercept=(coef(model_te)["(Intercept)"]);a=(coef(model_te)["sst"])}
  if (biome2 == "temperate"){b=(coef(model_te)["eco_typetemperate"]); d=(coef(model_te)["sst:eco_typetemperate"]); intercept=(coef(model_te)["(Intercept)"]);a=(coef(model_te)["sst"])}
  if (biome2 == "tropical"){b=(coef(model_te)["eco_typetropical"]); d=(coef(model_te)["sst:eco_typetropical"]); intercept=(coef(model_te)["(Intercept)"]);a=(coef(model_te)["sst"])}
  if (biome2 == "upwelling"){b=(coef(model_te)["eco_typeupwelling"]); d=(coef(model_te)["sst:eco_typeupwelling"]); intercept=(coef(model_te)["(Intercept)"]);a=(coef(model_te)["sst"])}
  return((exp((intercept + a*sst + b + d*sst)+0.5*var(residuals(model_te)))))
}
# -------------------------------------------------------------------
# Build catch trophic spectrum for a given year
# Required data: data_catch, data_clim, info cell
data_spectrum <- function(x,data_catch,data_clim,info_cell) {
  tl<-seq(2,7,.1)
  # catch
  data_catch_y <- filter(data_catch, year == as.numeric(x)) %>%
    select(taxonkey, tl, catch) %>%
    rename(group_name=taxonkey, TL=tl, catches=catch)%>%
    mutate(group_name=as.character(group_name))
  data_catch_smooth <- create.smooth(data_catch_y)
  spectrum_catch <- Transpose(data_catch_smooth, data_catch_y, column="catches")
  spectrum_data<-data.frame(tl=tl,catch=apply(spectrum_catch,1,sum)[2:52])
  # temperature data
  data_clim_y <- filter(data_clim, year == as.numeric(x))
  # Kin
  spectrum_data<-data.frame(spectrum_data,kin_v=kin_f(tl,data_clim_y$sst))
  # TE
  te_htl<-TE(info_cell$eco_type, data_clim_y$sst)
  te_ltl<-te_tl23(tl[1:11],data_clim_y$te_meso,te_htl)
  # te_ltl<-te_tl23(tl[1:11],0.1,te_htl)
  te_tot<-c(te_ltl,rep(te_htl,length(tl[12:51])))
  spectrum_data<-data.frame(spectrum_data,te=te_tot)
  # flow at tl 2
  spectrum_data<-data.frame(spectrum_data,flow2=rep(data_clim_y$flow2,length(tl)))
  # add year and cell_no
  spectrum_data<-data.frame(spectrum_data,flow2=rep(data_clim_y$flow2,length(tl)))
  spectrum_data<-data.frame(cell_no=rep(data_clim_y$cell_no[1],length(tl)),
                            year=rep(x,length(tl)),
                            spectrum_data)
  return(spectrum_data)
}
# -------------------------------------------------------------------
# function to calculate te between 2 and 3 as a linear relation from TE of the lower tl and TE of the upper tl
te_tl23<- function(x,TE2,TE3){
  a<- TE3 - TE2
  b<- 3*TE2 - 2*TE3
  return(a*x + b)
}
# te_tl23(2.9,0.8,0.1)
# --------------------------------------------------------------------
# --------------------------------------------------------------------
# cell_no<-cell_world[25600,1]
# cell_no <- 14933 # Celtic sea
# cell_no <- 17803
# --------------------------------------------------------------------
select_data<-function(cell_no) {
  # PP in g.an-1 (data in mole.m-2.s-1 --> x 12 g.mol-1 (masse molaire du carbone) 
  #  1e+6 : gram => ton
  #  x 365*24*3600 (day*hour*sec)
  # surface is in m^2
  # x * 9 (carbone ==> biomass)
  # TOTAL = 365*24*3600*9*12
  # SST in K -> °C 
  # linking mt database
  data_clim<-as.data.frame(dbGetQuery(liaisons2,paste(
    "select distinct cell_no, avg(sst) as sst, avg(pp) as pp, year, surface_st
    from ",data_npp," A inner join ",data_sst," B using (cell_no,year)
    inner join geo.world_grid1x1_v2 C using (cell_no)
    where cell_no = ",cell_no," and pp>=0 and sst>-100 and (year between 1950 and 2010)
    group by year, C.surface_st, cell_no",sep="")))
  data_catch<-as.data.frame(dbGetQuery(liaisons1,paste(
    "select distinct y_year as year, spe_code as taxonkey, tl,  catch
    from catch.catch_by_cell_1x1 A inner join ref.tl B on(A.spe_code=B.taxonkey)
    where cell_no = ",cell_no," and (y_year between 1950 and 2010)",sep="")))
  # ---------------------
  if (nrow(data_clim)>0) {
    # if no catch SAU ==> No fisheries ==> table with catch=0
    if (nrow(data_catch)==0) {
      data_catch<-data.frame(year=data_clim$year,taxonkey=100,tl=3,catch=0)
    } else {
      # if no catch SAU for some years ==> No fisheries for some uear ==> table with catch=0 for some years
      if (length(unique(data_clim$year)) != length(unique(data_catch$year))) {
        data_catch_year<-data.frame(year=data_clim$year)
        data_catch<-left_join(data_catch_year,data_catch,by="year")
        data_catch[is.na(data_catch$taxonkey),"taxonkey"]<-100
        data_catch[is.na(data_catch$tl),"tl"]<-3
        data_catch[is.na(data_catch$catch),"catch"]<-0
        data_catch<-as.data.frame(data_catch)
      }
    }
    data_clim$pp <- data_clim$pp * 365 * 24 * 3600 * 9 * 12 * data_clim$surface_st / 1e+6 
    data_clim$sst <- data_clim$sst- 273.15
    # info on cell
    info_cell<-as.data.frame(dbGetQuery(liaisons2,paste(
      "select distinct cell_no, eco_type from geo.biogeo3
      where cell_no = ",cell_no,sep="")))
    data_clim<-data_clim[order(data_clim$year),]
    # ------------------------------
    # TE planktonic food web
    # P2=P1*0.1 # because by convention p2/p1 = 0.1
    TE_data <-as.data.frame(dbGetQuery(liaisons1,paste(
      "select distinct (te_meso) as te_meso,year from ltl.ltl_data4
      where cell_no = ",cell_no,sep="")))
    if (nrow(TE_data)>0) {
      data_clim<-inner_join(data_clim,TE_data,by="year")
      # P2
      data_clim$flow2=data_clim$pp*data_clim$te_meso
      # ---------------------------------------------------------------------
      # Gather all the fixed data at each trophic level for all year over the considered time period
      dt_fixed <- lapply(data_clim$year, data_spectrum,data_catch=data_catch,data_clim=data_clim,info_cell=info_cell)
      names(dt_fixed)<-data_clim$year
      rm(data_catch)
      rm(data_clim)
      rm(info_cell)
      return(dt_fixed)
    }
  }
  }
# select_data(14933)
# test function data spectrum
# data_spectrum(1950)
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
# --------------------------------------------------------------------
# base function to calculate production and biomass
core_function <-function(dt_base) {
  # top down parameters
  a_td<-0.4
  b_td<-0.5
  iteration=0
  # seq_tl
  seq_tl <- seq(2,7,by=0.1)
  dt_add <-as.data.frame(matrix(,nrow = nrow(dt_base), ncol = 15))
  colnames(dt_add)<-c("flow","flow_mean","flow_mean_v" ,"prod","prod_v", "bio","bio_v", "kin", 
                      "f_thau_ref", "f_thau","flow_optm","prod_optm","bio_optm","kin_optm","supported")
  dt_base<-cbind(dt_base,dt_add)
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
  } 
  # --------------------------------------------------------------------
  # result<-list(dt_base,iteration)
  # result<-dt_base
  # return(result)
  dt_base$iteration<-iteration
  result<-as.data.frame(dt_base[-51,c("cell_no","year","tl","catch","prod","prod_v","bio","bio_v","supported","f_thau","iteration")])
  colnames(result)<-c("cell_no","year","tl","catch","production","production_v","biomass","biomass_v","supported","f_morta","iteration")
  return(result)
}
# --------------------------------------------------------------------
cell_fao28 <- as.data.frame(dbGetQuery(liaisons1,"select distinct cell_no from geo.world_grid1x1_v2  inner join geo.biogeo3 using(cell_no) inner join geo.cross_fao_cells using (cell_no) inner join geo.fao_major using (f_code) where gid=28"))
# dim(cell_fao28)
# cell_no<-10612
run_et<-function(cell_no) {
  dt_input<-select_data(cell_no)
  res<-lapply(dt_input,core_function)
  if (length(res)>0) {
    res <- do.call("rbind", res)
    # return(res)
    dbWriteTable(liaisons1,c("et_paper2",name_table),res, row.names=FALSE, append=TRUE)
  }
}
# run_et(17803)
# create table to store data
# output 
name_table <- "rcp85_gfdl_catch_fao28"
dbGetQuery(liaisons1,paste("create table et_paper2.",name_table," (
      cell_no integer, 
      year integer, 
      tl numeric,
      catch numeric, 
      production numeric, 
      production_v numeric, 
      biomass numeric, 
      biomass_v numeric, 
      supported integer, 
      f_morta numeric, 
      iteration integer)",sep=""))
# system.time({
# run_et(11719) 
#   })
# system.time({
#   for (i in cell_fao28[1:160,]) {
#     print(i)
#     print(system.time({run_et(i)}))
#   }
# })
# -------------------------------------------------------------------------------
# # num_cores <- detectCores()
system.time({
  cl <- makeCluster(20) 
  clusterEvalQ(cl, {
    library(RPostgreSQL)
    library(parallel)
    library(reshape2)
    library(tidyverse) 
    library(EcoTroph)
    # linking mt database
    drv <- dbDriver("PostgreSQL")
    liaisons2 <- dbConnect(drv, dbname="environnement",host="halieut.agrocampus-ouest.fr",port=5432,user="hubert",password="dup098!")
    liaisons1 <- dbConnect(drv, dbname="world_mapping",host="sirs.agrocampus-ouest.fr",port=5432,user="hubert",password="dup098!")
  })
  clusterExport(cl, c('cell_fao28','select_data','TE','kin_f','core_function','run_et','data_spectrum','model_te',
                      'topd_high_tl','topd_classic','te_tl23','name_table','data_sst','data_npp'))
  values <-cell_fao28[,1] #
  result <- parLapply(cl,values,run_et)   # paralel execution
  clusterEvalQ(cl, {
    dbDisconnect(liaisons2)
    dbDisconnect(liaisons1)
  })
  stopCluster(cl)
})
# -------------------------------------------------------------------------------