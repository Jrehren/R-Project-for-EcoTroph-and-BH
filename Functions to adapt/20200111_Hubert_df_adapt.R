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

# XXXXXXX could we put this in a seperate function, in order to separate
#         the functions from the analysis?

# ---------------------------------------------------------------------
# model du pontavice et al., in review
source("model_te.R") 
# ---------------------------------------------------------------------
# Data
data_npp<-"model_ubc.pp_gfdl_esm2m_85"
data_sst<-"model_ubc.sst_gfdl_esm2m_85"

# XXXXXXXX So does this mean that he is calculating the transfer efficiency and 
#          biomass in this script from his net primary production? 

# --------------------------------------------------------------------
# linking mt database
drv <- dbDriver("PostgreSQL")
liaisons2 <- dbConnect(drv, dbname="environnement",host="halieut.agrocampus-ouest.fr",port=5432,user="hubert",password="dup098!")
liaisons1 <- dbConnect(drv, dbname="world_mapping",host="sirs.agrocampus-ouest.fr",port=5432,user="hubert",password="dup098!")

# XXXXXXX For me this is all in a seperate script in which he generates his input
#         Table

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

# XXXXXXXXXX All of this seems to be for generating the basis for the data frame
#            for the analysis
