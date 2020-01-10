################################################################################
############################# SETTING THE SCENE ################################
################################################################################

# -------------------- EXPLANATION AND INDEXING OF SCENARIOS -------------------
# 1.) Explanation of analysis ----
# This is the collection of analysis for my first France post-doc paper. The
# equations here used are mainly based on the equations from Gascuel et al., 2011. 
# The code used comes from the EcoTroph package, but has been modified to 1.) 
# be able to run the code for my purpose, 2) optimize and correct some few 
# things in the code, 3) correct an inconsistency coming from the production equation. 
# To tackle the latter point, we changed the mathematical implementation of the 
# fishing loss rate in the production equation. To further improve the model we
# used the proper continous implementation as presented in Gascuel et al., 2011, 
# instead of using an approximation. All corrections made to the source code
# can be found in the doc file: "00_protocol_adapting the source code.docx".
# We adapted the code to be run on a virgin ecosystem with a starting fishing 
# mortality of 0. 
# While originally I run many different scenarios with fixed fishing mortalities, 
# compensated harvest scenarios and others, here only the scenarios are run,
# which are presented in the paper. The scenarios I used in the analysis are 
# listed below. 
# 2.) Indexing the different scenarios (Balanced harvest = BH) ----
# 1. F=c*K, selec=1, TLst=2        -->   bhm12     --> BH F~P/B
# 2. F=c*P, selec=1, TLst=2        -->   bhm12_P   --> BH F~P
# 3. F=c*K, selec=s-shape          -->   bumss     --> BH F~P/B with selectivity
# 4. Searching the best F selec=1  -->   alt       --> F vector with no impact on
#                                                      the biomass structure (AH)
# 5. Searching the best F          -->   alt_cc    --> AH with selectivity
#    selec=s-shape

# Scenarios were we always varied the intensity of fishing and the trophic level
# at first catch (at which trophic level are we starting to fish fully).

# 7. Fishing, c fixed, High-low,   -->  fehtla_acc   --> Multiple runs of TL50 and asymptote
#    TL50 + asymptote + accessibility                with the "normal" logistic function + on top of it we run an accessibility function
# 8. Fishing, c fixed, High-low,   -->  feltla_acc   --> Multiple runs of TL50 and asymptote + on top of it we run an accessibility function
#    TL50 + asymptote + accessibility  
# 9. Fishing, c fixed, High-low,   -->  ffhtla_acc   --> Multiple runs of TL50 and asymptote + on top of it we run an accessibility function
#    TL50 + asymptote + accessibility  

#--------------------------- LOAD PACKAGES & FUNCTIONS -------------------------
# 1.) CLEAR SCREEN, SET DIRECTORY ----
ls()
rm(list=ls())

# saving the old par
# op<-par(no.readonly=TRUE)
# par(op)

dir<-setwd("C:/Users/Ursel/ownCloud/France_Post-doc/Data analysis/Paper I_Balanced harvest with EcoTroph/Analysis of the logarithmic solution")

# 2.) LOAD PACKAGES ----
library(viridis)
library(dplyr)
library(DescTools)
library(tidyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(mosaic)

# 3.) LOAD FUNCTIONS ----
source("Functions/20181120_Functions for ET.Main.R")
source("Functions/sub.mf_cran.R")
source("Functions/Logarithmic solution/20191130_Data frame function_integrate_access adapted.R")
source("Functions/Logarithmic solution/20190411_Virtual ecosystem function_continous implementation.R")
source("Functions/20190214_Dynamic fishing mort_function.R")
source("Functions/20190214_Fixed fishing mort_function.R")
source("Functions/20190805_Performance measures.R")

# loading the new generalized logistic function for the run ofmultiple fishing
# distributions
source("Functions/20190409_Generalized logistic function.R")

# these are the new files for the alternative harvest calcualtions
source("Functions/20190307_Alternative harvest function.R")
source("Functions/20190306_et.diagnosis_adapted_with f_mort_alt.R")
# Setting F ~ P 
source("Functions/20190701_Dynamic fishing mort_F ~ production.R")
# loading the new solution
source("Functions/Logarithmic solution/20191014_mf.diagnosis_adapted_log_solution_die zweite.R")
source("Functions/Logarithmic solution/20191014_et.diagnosis_adapted_log_solution_integrate.R")
# for the alternative harvest function
source("Functions/Logarithmic solution/20191022_Alternative harvest function_log_solution_adapting the top-down control_integrated.R")
source("Functions/Logarithmic solution/20191016_et.diagnosis_adapted_with f_mort_alt_log_solution.integrated.R")

################################################################################
#################################### ANALYSIS ##################################
################################################################################

# ------------------------------ RUN SCENARIOS ---------------------------------
# 1.) CREATE A VIRGIN ECOSYSTEM ---- 

VE<-virtual_eco_integ(TE = 0.1)

# 2.) BHM12 : BALANCED HARVEST F~K ---- 

 ## A.) CREATE A DATA FRAME (ET_Main) used in the create.diagnosis function ----

dat_bhm12<-dat_frame_integrate_access_adap(VE=VE, asymptote = 2, TL50 = 2, slope = 0)

 ## B.) CREATE A FISHING MORTALITY STARTING VECTOR ----

F_morts_bhm12<-F_mort_dyn(dat = dat_bhm12, const=0.1)

 ## C.) RUN CREATE.DIAGNOSIS ----

df_bhm12<-create.ETdiagnosis.logsol.integrated(data=dat_bhm12, 
                                   Mul_eff=seq(0,10, 0.5), 
                                   Fish_mort = F_morts_bhm12$Fish_mort,
                                   Fish_mort_acc = F_morts_bhm12$Fish_mort_acc, 
                                   TopD = 0.4)


# 3.) BHM12_P : BALANCED HARVEST F~P ---- 
 ## A.) CREATE A DATA FRAME (ET_Main) used in the create.diagnosis function ----  

dat_bhm12_P<-dat_frame_integrate_access_adap(VE=VE, asymptote = 2, TL50 = 2, slope = 0)

 ## B.) CREATE A FISHING MORTALITY STARTING VECTOR ----

# the constant is set so that trophic level 2 has the same fishing mortality 
# as in the variant of F~K
# with a 0.0039 constant the fishing mortalities have a similar value at TL 2
# with a value of 0.008 the total fishing mortalities are similar
# with a value of 0.005 the fishing loss rates of TL 2 are similar

F_morts_bhm12_P<-F_mort_prod(dat = dat_bhm12_P, const=0.0039)

# Compare the fishing mortality vectors:

F_morts_bhm12_P$Fish_mort$catch.1[2]/F_morts_bhm12$Fish_mort$catch.1[2]

 ## C.) RUN CREATE.DIAGNOSIS ----

df_bhm12_P<-create.ETdiagnosis.logsol.integrated(data=dat_bhm12_P, 
                                   Mul_eff=seq(0,10,0.5), 
                                   Fish_mort = F_morts_bhm12_P$Fish_mort, 
                                   Fish_mort_acc = F_morts_bhm12_P$Fish_mort_acc, 
                                   TopD = 0.4)

# 4.) BUMSS : BALANCED HARVEST F~K WITH ACCESSIBILITY ---- 
 ## A.) CREATE A DATA FRAME (ET_Main) used in the create.diagnosis function ----

dat_bumss<-dat_frame_integrate_access_adap(VE=VE, asymptote = 1, TL50 = 2.5, slope = 4.84)

 ## B.) CREATE A FISHING MORTALITY STARTING VECTOR ----

F_morts_bumss<-F_mort_dyn(dat=dat_bumss, const=0.1)

 ## C.) RUN CREATE.DIAGNOSIS ----

df_bumss<-create.ETdiagnosis.logsol.integrated(data=dat_bumss, 
                                  Mul_eff=seq(0,10,0.5), 
                                  Fish_mort = F_morts_bumss$Fish_mort, 
                                  Fish_mort_acc = F_morts_bumss$Fish_mort_acc, 
                                  TopD = 0.4)


# 5.) BUMSS_P : BALANCED HARVEST F~P WITH ACCESSIBILITY ----
 ## A.) CREATE A DATA FRAME (ET_Main) used in the create.diagnosis function ----

dat_bumss_P<-dat_frame_integrate_access_adap(VE=VE, asymptote = 1, TL50 = 2.5, slope = 4.84)

 ## B.) CREATE A FISHING MORTALITY STARTING VECTOR ----

F_morts_bumss_P<-F_mort_prod(dat = dat_bumss_P, const=0.0383)
F_morts_bumss_P$Fish_mort_acc$catch.1[2]
F_morts_bumss$Fish_mort_acc$catch.1[2]
 ## C.) RUN CREATE.DIAGNOSIS ----

df_bumss_P<-create.ETdiagnosis.logsol.integrated(data=dat_bumss_P, 
                                  Mul_eff=seq(0,10,0.5), 
                                  Fish_mort = F_morts_bumss_P$Fish_mort, 
                                  Fish_mort_acc = F_morts_bumss_P$Fish_mort_acc)


# 6.) AH : ALTERNATIVE HARVEST SCENARIO ----
 ## A.) CREATE DATA FRAME ----

dat_alt<-dat_frame_integrate_access_adap(VE=VE, asymptote = 2, TL50 = 2, slope = 0)

 ## B.) CREATE A FISHING MORTALITY STARTING VECTOR ----

# the constant is set so that trophic level 2 has the same fishing mortality 
# as in the variant of F~K & F~P
# with a 0.0039 constant the fishing mortalities have a similar value at TL 2
# with a value of 0.008 the total fishing mortalities are similar
# with a value of 0.005 the fishing loss rates of TL 2 are similar

# F_morts_alt<-F_mort_dyn(dat = dat_alt, const=0.1)
F_morts_alt<-F_mort_prod(dat = dat_alt, const=0.0033)
F_morts_alt$Fish_mort_acc$catch.1[2]
df_bhm12$`1`$Fish_mort_acc[2]
dat_alt$ET_Main$Kin/dat_alt$ET_Main$Kin_acc
 ## C.) RUN CREATE.DIAGNOSIS ----

df_alt <- create.ETdiagnosis.alt_logsol.integrated(data=dat_alt, 
          Mul_eff=seq(0,10,0.5), Fish_mort = F_morts_alt$Fish_mort, 
          Fish_mort_acc = F_morts_alt$Fish_mort_acc, TopD = 0.4)


# comparing the fishing mortalities with the BH impelementation, and adjusting 
# if necessary

round(df_alt$`1`$Fish_mort_acc[2], digits = 3)
round(df_bhm12$`1`$Fish_mort_acc[2], digits = 3)
round(df_bhm12_P$`1`$Fish_mort_acc[2], digits = 3)

round(df_alt$`2`$Fish_mort_acc[2], digits = 2)
round(df_bhm12$`2`$Fish_mort_acc[2], digits = 2)
round(df_bhm12_P$`2`$Fish_mort_acc[2], digits = 2)

round(df_alt$`5`$Fish_mort_acc[2], digits = 2)
round(df_bhm12$`5`$Fish_mort_acc[2], digits = 2)
round(df_bhm12_P$`5`$Fish_mort_acc[2], digits = 2)

round(df_alt$`10`$Fish_mort_acc[2], digits = 4)
round(df_bhm12$`10`$Fish_mort_acc[2], digits = 4)
round(df_bhm12_P$`10`$Fish_mort_acc[2], digits = 4)


# 7.) AH_acc : ALTERNATIVE HARVEST SCENARIO WITH ACCESSIBILITY ----
 ## A) CREATE DATA FRAME ----

dat_alt_acc<-dat_frame_integrate_access_adap(VE=VE, asymptote = 1, TL50 = 2.5, 
                                             slope = 4.84)

 ## B.) CREATE A FISHING MORTALITY STARTING VECTOR ----

# I do not find an optimum solution for the F-vector with accessibility, 
# but I first take out the optimized F-vector (for each mE) from a run 
# without accessibility. These optimized F-vectors are then used to calculate
# accessibile F-vectors when following the s-shaped accessibility assumption 
# the F and accessible F vectors are then used in the normal et.diagnosis analysis 

# create the nested list of fishing mortalities from your f_vectors taken from 
# df_alt

Fish_mort_acc_altacc=c()
for (i in seq_along(df_alt)){
  Fish_mort_acc_altacc[[i]]=list(catch.1=df_alt[[i]]$Fish_mort_acc)
}
names(Fish_mort_acc_altacc)<-names(df_alt)

# use this list of F to calculate the accessible F
F_new_acc<-c()
Fish_mort_alt_acc<-c()
for (i in names(df_alt)){
  F_new_acc[[i]]<-c(0,(Fish_mort_acc_altacc[[i]]$catch.1[-1]*
                         dat_alt_acc$ET_Main$Selec[-1]))
  Fish_mort_alt_acc[[i]]=list(catch.1=F_new_acc[[i]])
}

 ## C.) RUN CREATE.DIAGNOSIS ----

df_alt_acc<-mapply(create.ETdiagnosis.logsol.integrated, 
                   Fish_mort_acc=Fish_mort_acc_altacc,
                   Fish_mort=Fish_mort_alt_acc, 
                   MoreArgs = list(data= dat_alt_acc, Mul_eff=c(0,1)), 
                   SIMPLIFY = F, USE.NAMES = T)



# ----------------------- CALCULATE PERFORMANCE INDICATORS ---------------------
# ----
bhm12_perform<-lapply(df_bhm12,function (x) perform(data=x, VE=df_bhm12$`0`))
bhm12_perform_P<-lapply(df_bhm12_P,function (x) perform(data=x, VE=df_bhm12_P$`0`))
bumss_perform<-lapply(df_bumss,function (x) perform(data=x, VE=df_bumss$`0`))
bumss_perform_P<-lapply(df_bumss_P,function (x) perform(data=x, VE=df_bumss_P$`0`))
alt_perform<-lapply(df_alt,function (x) perform(data=x, VE=df_alt$`0`))

df_alt_acc_sub<-c()
for (i in names(df_alt_acc)){
  df_alt_acc_sub[[i]]<-df_alt_acc[[i]]$`1`
}
alt_acc_perform<-lapply(df_alt_acc_sub,
                        function (x) perform(data=x, VE=df_alt_acc$`1`$`0`))


# -------------------- PART III : MULITPLE FISHING SCENARIOS -------------------
# 1.) USING A FIXED EXPLOITATION WITH ACCESSIBILITY ----
 ## A) HIGH TL FISHED HARD WITH ACCESSIBILITY F~P/B (fehtla_acc)----
# setting the vector for TL50 and the asymptote and defining the seq_tl

seq_tl<-c(1,seq(2,7,0.1));length(seq_tl)
seq_tl_2up<-c(seq(2,7,0.1));length(seq_tl_2up)
TL50_vec_fehtla_acc=seq(1,5, 0.1);length(TL50_vec_fehtla_acc)
names(TL50_vec_fehtla_acc)<-TL50_vec_fehtla_acc

K_vec_fehtla_acc=seq(0,1,0.1);length(K_vec_fehtla_acc)    # asymptote
names(K_vec_fehtla_acc)<-K_vec_fehtla_acc


# creating a new ET_Main for each combination of TL50 and asymptote value
expand_fehtla_acc<-expand.grid(TL50_vec_fehtla_acc,K_vec_fehtla_acc);dim(expand_fehtla_acc)

names(expand_fehtla_acc)<-c("TL50","Asymptote")
head(expand_fehtla_acc)
expand_fehtla_acc$name <- paste(expand_fehtla_acc$TL50,expand_fehtla_acc$Asymptote)

length(df_alt)
# with this new combination of TL50 and asymptotes, I create multiple vectors of
# exploitation 
exploit_vec_fehtla_acc<-mapply(FUN = gl_selec, TL50=expand_fehtla_acc[,1], 
                               K=expand_fehtla_acc[,2],
                               MoreArgs = list(x=seq_tl_2up, slope=4.84, p=0, A=0), SIMPLIFY = F)
length(exploit_vec_fehtla_acc)
exploit_vec_fehtla_acc$`0`[2:37]
names(exploit_vec_fehtla_acc)<-expand_fehtla_acc$name
names(exploit_vec_fehtla_acc)

# here I visualize the exploitation pattern 
cl_fehtla_acc<-viridis(length(exploit_vec_fehtla_acc), direction = -1)
min=2
max=37
plot(seq_tl[min:max],exploit_vec_fehtla_acc$`2 0`[min:max], type="n", ylim=c(0,1))
for (i in seq_along(exploit_vec_fehtla_acc)){
  lines(seq_tl[min:max], exploit_vec_fehtla_acc[[i]][min:max], col=cl_fehtla_acc[i])
}

# Here, I calculate my data frame. Because we do not assume differences in the
# accessibility per trophic level, as this is a theoretical exploration, I use
# an accessibility of 1 for all, except trophic level 1

dat_fehtla_acc<-dat_frame_integrate_access_adap(VE=VE, asymptote = 1, TL50 = 2.5, slope = 4.84)

# now I use the multiple vectors of exploitation to calculate the fishing
# mortality over trophic level, using the F_mort_dyn function 

F_morts_fehtla_acc<-lapply(exploit_vec_fehtla_acc,
                           function (a) F_mort_dyn(const = a, dat = dat_fehtla_acc))
exploit_vec_fehtla_acc$`1.2 0`
length(F_morts_fehtla_acc)
str(F_morts_fehtla_acc)


# Because, a) I have now a list of 17 (depending on the TL50 vector) each with
# two vectors: fishing mortality and accessible fishing mortality and b) I want to loop
# over each list but taking either of the vectors and c) I don't know how to do this
# in the simple mapply command, I use a for loop to seperate these two vectors. 

F_mort_fehtla_acc<-c()
F_mort_acc_fehtla_acc<-c()
for (i in seq_along(F_morts_fehtla_acc)){
  F_mort_fehtla_acc[i]<-F_morts_fehtla_acc[[i]]["Fish_mort"]
  F_mort_acc_fehtla_acc[i]<-F_morts_fehtla_acc[[i]]["Fish_mort_acc"]
}
names(F_mort_fehtla_acc)<-expand_fehtla_acc$name;names(F_mort_fehtla_acc)
names(F_mort_acc_fehtla_acc)<-expand_fehtla_acc$name
length(F_mort_fehtla_acc)
# Now I can run the et diagnosis:
df_fehtla_acc<-mapply(create.ETdiagnosis.logsol.integrated, 
                      Fish_mort_acc=F_mort_acc_fehtla_acc, 
                      Fish_mort=F_mort_fehtla_acc, 
                      MoreArgs = list(data= dat_fehtla_acc, Mul_eff=c(0,1)), 
                      SIMPLIFY = F, USE.NAMES = T)
df_fehtla_acc_test<-df_fehtla_acc
# Now I calculate the performance measure for all runs and remove the last 
# performance measure from the lists
length(df_fehtla_acc)
perform_fehtla_acc<-c()
perform_fehtla_acc_sub<-c()
for (i in seq_along(df_fehtla_acc)){
  perform_fehtla_acc[[i]]<-perform(data=df_fehtla_acc[[i]][[2]], 
                                   VE=df_fehtla_acc[[i]][[1]])
  perform_fehtla_acc_sub[[i]]<-perform_fehtla_acc[[i]][1:10]
}
names(perform_fehtla_acc_sub)<-expand_fehtla_acc$name
names(perform_fehtla_acc_sub)
length(perform_fehtla_acc_sub)
length(perform_fehtla_acc_sub$`1 0`)
# I collapse the performance measures to a data frame
# as it is easier to work with a data frame later on

perform_fehtla_acc_dat <- data.frame(matrix(unlist(perform_fehtla_acc_sub), nrow=451 
                                            ,byrow = T
),
stringsAsFactors=F)
dim(perform_fehtla_acc_dat)
colnames(perform_fehtla_acc_dat)<-names(perform_fehtla_acc[[1]][1:10])
rownames(perform_fehtla_acc_dat)<-expand_fehtla_acc$name
names(perform_fehtla_acc_dat)
# add the expand.grid values to reshape each performance measure as a matrix for
# plotting
perform_fehtla_acc_dat_new<-cbind(expand_fehtla_acc[,1:2], perform_fehtla_acc_dat)
length(perform_fehtla_acc_dat_new)
names(perform_fehtla_acc_dat_new)
# I need 8 matrices with values arranged so that the first value of the rowname 
# becomes the rownames of the new matrix and the second name becomes the columnnames

fehtla_acc_mat<-c()
for (i in names(perform_fehtla_acc_dat_new[,3:10])){
  fehtla_acc_mat[[i]]<-acast(perform_fehtla_acc_dat_new,Asymptote ~ TL50, 
                             value.var=i)
}


 ## B) HIGH TL FISHED HARD WITH ACCESSIBILITY AH (fehtla_acc_AH) ----

# setting the vector for TL50 and the asymptote and defining the seq_tl

seq_tl<-c(1,seq(2,7,0.1));length(seq_tl)
seq_tl_2up<-c(seq(2,7,0.1));length(seq_tl_2up)
TL50_vec_fehtla_acc_Ah=seq(1,5, 0.1);length(TL50_vec_fehtla_acc_Ah)
names(TL50_vec_fehtla_acc_Ah)<-TL50_vec_fehtla_acc_Ah

# Now I create selectivity vectors, that are used to create the actual fishing
# mortalities.
accessibility_vec<-lapply(TL50_vec_fehtla_acc_Ah, 
                          function (i) s_selec(TL50=i, seq_tl_2up, asymptote=1,
                                               slope=4.84))

# The fishing intensity (The asymptote) is calculated differently:
# We take the accessible fishing mortality vectors and use them directly as 
# fishing intensities and use the TL50 values to calculate new fishing mortalites

dat_alt<-dat_frame_integrate_access_adap(VE=VE, asymptote = 2, TL50 = 2, slope = 0)
dat_alt$ET_Main$Selec

# F_morts_alt<-F_mort_dyn(dat = dat_alt, const=0.043)
F_morts_alt<-F_mort_prod(dat = dat_alt, const=0.0033)

# we only want 11 intensities so we use only 1 step
df_alt_formultirun<-create.ETdiagnosis.alt_logsol.integrated(data=dat_alt, 
                                                  Mul_eff=seq(0,10,1), Fish_mort = F_morts_alt$Fish_mort, 
                                                  Fish_mort_acc = F_morts_alt$Fish_mort_acc, TopD = 0.4)
names(df_alt_formultirun)
length(df_alt_formultirun)
Mul_eff=seq(0,10,1)
E_vec<-Mul_eff*0.046
# store the accessible fishing mortalitiy vector in an object

F_mort_acc_AH<-lapply(1:11, function(i) df_alt_formultirun[[i]][["Fish_mort_acc"]][2:52])
length(F_mort_acc_AH)
names(F_mort_acc_AH)<-names(df_alt_formultirun)
length(F_mort_acc_AH$`1`)

# I calculate the accessible fishing mortality given the different TL50 strategies

F_mort_acc_fehtla_acc_AH1<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[1]])
F_mort_acc_fehtla_acc_AH2<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[2]])
F_mort_acc_fehtla_acc_AH3<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[3]])
F_mort_acc_fehtla_acc_AH4<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[4]])
F_mort_acc_fehtla_acc_AH5<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[5]])
F_mort_acc_fehtla_acc_AH6<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[6]])
F_mort_acc_fehtla_acc_AH7<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[7]])
F_mort_acc_fehtla_acc_AH8<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[8]])
F_mort_acc_fehtla_acc_AH9<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[9]])
F_mort_acc_fehtla_acc_AH10<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[10]])
F_mort_acc_fehtla_acc_AH11<-lapply(accessibility_vec,function (a) a*F_mort_acc_AH[[11]])


F_mort_acc_fehtla_acc_AH_all<-list("1"=F_mort_acc_fehtla_acc_AH1, 
                                   "2"=F_mort_acc_fehtla_acc_AH2, 
                                   "3"=F_mort_acc_fehtla_acc_AH3, "4"=F_mort_acc_fehtla_acc_AH4, 
                                   "5"=F_mort_acc_fehtla_acc_AH5, "6"=F_mort_acc_fehtla_acc_AH6, 
                                   "7"=F_mort_acc_fehtla_acc_AH7, "8"=F_mort_acc_fehtla_acc_AH8, 
                                   "9"=F_mort_acc_fehtla_acc_AH9, "10"=F_mort_acc_fehtla_acc_AH10, 
                                   "11"=F_mort_acc_fehtla_acc_AH11)  
length(F_mort_acc_fehtla_acc_AH_all)
F_mort_acc_fehtla_acc_AH_all_flat <- do.call(c, F_mort_acc_fehtla_acc_AH_all)
names(F_mort_acc_fehtla_acc_AH_all_flat)<- names(df_fehtla_acc)

length(F_mort_acc_fehtla_acc_AH_all_flat)

# I have the final accessible fishing mortality vectors now and I want to
# calculate the total fishing mortality given the assumed biological selectivity 
# with a TL50 of 2.5.

# I first must add the F value (0) of the first trophic level to get an F-vector of 52

F_mort_acc_fehtla_acc_AH_all_flat_com<-c()
for (i in names(F_mort_acc_fehtla_acc_AH_all_flat)){
  F_mort_acc_fehtla_acc_AH_all_flat_com[[i]]<-list(catch.1=c(0,F_mort_acc_fehtla_acc_AH_all_flat[[i]]))
}

# Now I multiply the accessible fishing mortality with the biological selectivity
# to get the total fishing mortality. 

# Here, I calculate my data frame, given that we have a biological accessibility

dat_fehtla_acc_AH<-dat_frame_integrate_access_adap(VE=VE, asymptote = 1, TL50 = 2.5, slope = 4.84)

F_mort_acc_fehtla_acc_AH_tot_fishmort<-c()
for (i in names(F_mort_acc_fehtla_acc_AH_all_flat_com)){
  F_mort_acc_fehtla_acc_AH_tot_fishmort[[i]]<-list(catch.1=(F_mort_acc_fehtla_acc_AH_all_flat_com[[i]]$catch.1*
                                                              dat_fehtla_acc_AH$ET_Main$Selec))
}
length(dat_fehtla_acc_AH)

# Now I have the two fishing mortalities seperated in two different lists
# and can directly loop over them
df_fehtla_acc_AH<-mapply(create.ETdiagnosis.logsol.integrated, 
                         Fish_mort_acc=F_mort_acc_fehtla_acc_AH_all_flat_com, 
                         Fish_mort=F_mort_acc_fehtla_acc_AH_tot_fishmort, 
                         MoreArgs = list(data= dat_fehtla_acc_AH, Mul_eff=c(0,1)), 
                         SIMPLIFY = F, USE.NAMES = T)
length(df_fehtla_acc_AH)



# Now I calculate the performance measure for all runs and remove the last 
# performance measure from the lists
length(df_fehtla_acc_AH)
perform_fehtla_acc_AH<-c()
perform_fehtla_acc_sub_AH<-c()
for (i in names(df_fehtla_acc_AH)){
  perform_fehtla_acc_AH[[i]]<-perform(data=df_fehtla_acc_AH[[i]][[2]], 
                                      VE=df_fehtla_acc_AH$`1 0`$`0`)
  perform_fehtla_acc_sub_AH[[i]]<-perform_fehtla_acc_AH[[i]][1:10]
}


names(perform_fehtla_acc_sub_AH)
names(perform_fehtla_acc_AH)
length(perform_fehtla_acc_sub_AH$`1 0`)
names(perform_fehtla_acc_AH)
length(perform_fehtla_acc_sub_AH$`1.1`)
perform_fehtla_acc_sub_AH$`1 0`$d
# I collapse the performance measures to a data frame
# as it is easier to work with a data frame later on

perform_fehtla_acc_dat_AH <- data.frame(matrix(unlist(perform_fehtla_acc_sub_AH),
                                               nrow=451,byrow = T),
                                        stringsAsFactors=F)
dim(perform_fehtla_acc_dat_AH)

names(perform_fehtla_acc_dat_AH)
length(perform_fehtla_acc_dat_AH)
colnames(perform_fehtla_acc_dat_AH)<-names(perform_fehtla_acc_AH[[1]][1:10])
rownames(perform_fehtla_acc_dat_AH)<-names(df_fehtla_acc_AH)

# add the expand.grid values to reshape each performance measure as a matrix for
# plotting

# I need to create name vectors for my asyptote and my TL50 because I need to 
# add it to later seperate it into lists
Asymptote=c(0:10);length(Asymptote)
nam_vec<-expand.grid(TL50_vec_fehtla_acc_Ah, Asymptote)
names(nam_vec)<-c("TL50", "Asymptote")
head(names(df_fehtla_acc_AH))
head(nam_vec)


perform_fehtla_acc_dat_new_AH<-cbind(nam_vec, perform_fehtla_acc_dat_AH)
names(perform_fehtla_acc_dat_new_AH)
dim(perform_fehtla_acc_dat_new_AH)


# I need 8 matrices with values arranged so that the first value of the rowname 
# becomes the rownames of the new matrix and the second name becomes the columnnames

fehtla_acc_mat_AH<-c()
for (i in names(perform_fehtla_acc_dat_new_AH[,3:10])){
  fehtla_acc_mat_AH[[i]]<-acast(perform_fehtla_acc_dat_new_AH,
                                Asymptote ~ TL50,value.var=i)
}
fehtla_acc_mat_AH$d



# --------------------------- RESULTS OF THE SCENARIOS -------------------------

# 1.) SCENARIOS WITHOUT ACCESSIBILITY ----
 ## A.) BIND THE SCENARIOS IN ONE DATA FRAME TO LOOP OVER ----
# I first have to remove one level from the accessible alternative harvest 
# data frame, because it is a nested list.

df_alt_acc_adapt<-list()
for (i in names(df_alt_acc)){
df_alt_acc_adapt[[i]]<-df_alt_acc[[i]]$`1`
}

scenarios<-list(bhm12=df_bhm12, bhm12_P=df_bhm12_P, bumss=df_bumss, 
                bumss_P=df_bumss_P, AH=df_alt, AH_acc=df_alt_acc_adapt)

# Names
names<-c()          # storing the effort multipliers in a list to use in legends
for (i in names(df_bhm12)){
  names[i]<-c(df_bhm12[[i]][["mf"]])
}
names<-names[-2]            # Because AH does not contain the mE=0.5, I remove it
# from names to be able to loop over 20 mE instead of 21


 ## B.) EXTRACT SOME VARIABLES FROM THE DATA FRAMES ----
   ### i.) Extracting total catch & TL catch ----
total_catch<-list()
total_catch2.1up<-list()
C_tot2.7<-list()
C_tot3<-list()
C_tot3.5<-list()
C_tot4<-list()
C_tot4.5<-list()
C_tot5<-list()
C_tot5.5<-list()
C_tot7<-list()

for (j in names(scenarios)){
  for (i in names(names)){
    total_catch[[j]][i]<-sum(scenarios[[j]][[i]][["Catches.tot"]])
    total_catch2.1up[[j]][i]<-sum(scenarios[[j]][[i]]$Catches.tot[3:52]) 
    C_tot2.7[[j]][i]<-sum(scenarios[[j]][[i]]$Catches.tot[9])
    C_tot3[[j]][i]<-sum(scenarios[[j]][[i]]$Catches.tot[12])
    C_tot3.5[[j]][i]<-sum(scenarios[[j]][[i]]$Catches.tot[17])
    C_tot4[[j]][i]<-sum(scenarios[[j]][[i]]$Catches.tot[22])
    C_tot4.5[[j]][i]<-sum(scenarios[[j]][[i]]$Catches.tot[27])
    C_tot5[[j]][i]<-sum(scenarios[[j]][[i]]$Catches.tot[32])
    C_tot5.5[[j]][i]<-sum(scenarios[[j]][[i]]$Catches.tot[37])
    C_tot7[[j]][i]<-sum(scenarios[[j]][[i]]$Catches.tot[52])
  }
}
 ## C.) SET GLOBAL VARIABLES (colors, etc.) ----

E_acc<-names*0.1            # setting the x-axis for the catch plots
seq_tl<-c(1, seq(2,7,0.1))  # setting the x-acis for the biomass, production plots

cl1 <- viridis(length(names), direction = -1) # setting the colors
scenario_nam<-names(scenarios)
 ## D.) VISUALIZE OUTPUTS ----

pdf("Results/Figures/Scenarios_acc adapted.pdf")

for (j in seq_along(scenarios)){
   ### i.) BIOMASS ----
par(mar=c(4,4,3,4), mfrow=c(1,1))
plot(seq_tl[-1], scenarios[[j]]$`0`$BIOM_MF[-1], type="n", 
     ylab = "Biomass", xlab = "Trophic class", log="y", main= scenario_nam[j])
for (ii in seq_along(names)){
  i<-scenarios[[j]][[ii]]
  lines(seq_tl[-1], scenarios[[j]][[ii]][["BIOM_MF"]][-1], col=cl1[ii])
}

   ### ii.) ACCESSIBLE BIOMASS ----
par(mar=c(4,4,3,4), mfrow=c(1,1))
plot(seq_tl[-1], scenarios[[j]]$`0`$BIOM_MF_acc[-1], type="n", 
     ylab = "Accessible Biomass", xlab = "Trophic class", log="y", 
     main= scenario_nam[j])
for (ii in seq_along(names)){
  i<-scenarios[[j]][[ii]]
  lines(seq_tl[-1], scenarios[[j]][[ii]][["BIOM_MF_acc"]][-1], col=cl1[ii])
}
   ### iii.) CATCHES ----
par(mar=c(4,4,3,4), mfrow=c(1,1))
plot(seq_tl, scenarios[[j]]$`10`$Catches.tot, type="n", 
     ylab = "Catch", xlab = "Trophic class", main = scenario_nam[[j]])
for (ii in seq_along(names)){
  i<-scenarios[[j]][[ii]]
  lines(seq_tl, scenarios[[j]][[ii]][["Catches.tot"]], col=cl1[ii])
}

   ### iv.) KINETICS ----
par(mar=c(4,4,3,4), mfrow=c(1,1))
plot(seq_tl[-1], scenarios[[j]]$`10`$Kin_MF[-1], type="n", 
     ylab = "Kinetics", xlab = "Trophic class", main = scenario_nam[[j]])
for (ii in seq_along(names)){
  i<-scenarios[[j]][[ii]]
  lines(seq_tl[-1], scenarios[[j]][[ii]][["Kin_MF"]][-1], col=cl1[ii])
}

   ### v.) PRODUCTION ----
par(mar=c(4,4,3,4), mfrow=c(1,1))
plot(seq_tl[-1], scenarios[[j]]$`0`$Prod_MF[-1], type="n",
     ylab = "Production", xlab = "Trophic class", log = "y", 
     main = scenario_nam[[j]])
for (ii in seq_along(names)){
  i<-scenarios[[j]][[ii]]
  lines(seq_tl[-1], scenarios[[j]][[ii]][["Prod_MF"]][-1], col=cl1[ii])
}

   ### vi.) ACCESSIBLE PRODUCTION ----
par(mar=c(4,4,3,4), mfrow=c(1,1))
plot(seq_tl[-1], scenarios[[j]]$`0`$Prod_MF_acc[-1], type="n",
     ylab = "Accessible Production", xlab = "Trophic class", log = "y", 
     main = scenario_nam[[j]])
for (ii in seq_along(names)){
  i<-scenarios[[j]][[ii]]
  lines(seq_tl[-1], scenarios[[j]][[ii]][["Prod_MF_acc"]][-1], col=cl1[ii])
}

   ### vii.) TOTAL YIELD VS TOTAL FISHING MORTALITY ----
par(mar=c(4,4,3,4), mfrow=c(1,1))
plot(E_acc, total_catch[[j]], pch=19, xlab = "Starting exploitation", 
     ylab="Total catch", main=scenario_nam[j])

par(mar=c(4,4,3,4), mfrow=c(1,1))
plot(E_acc, total_catch2.1up[[j]], pch=19, xlab = "Starting exploitation", 
     ylab="Total catch from TL 2.1 up", main=scenario_nam[j])


   ### ix.) CATCH OF DIFFERENT TROPHIC LEVELS ----

par(mar=c(4,4,3,4), mfrow=c(1,1))
plot(E_acc, C_tot2.7[[j]], pch=19, xlab = "Starting exploitation", 
     ylab="Catch", type="n", main=scenario_nam[j])
lines(E_acc, C_tot2.7[[j]])
lines(E_acc, C_tot3[[j]], col="red")
lines(E_acc, C_tot3.5[[j]], col="blue")
legend("topright", legend = c("TL 2.7", "TL 3", "TL 3.5"), pch=19, 
       col = c("black", "red", "blue"))


par(mar=c(4,4,3,4), mfrow=c(1,1))
plot(E_acc, C_tot4.5[[j]], pch=19, xlab = "Starting exploitation", 
     ylab="Catch", type="n", main=scenario_nam[j])
lines(E_acc, C_tot4.5[[j]], col="red")
lines(E_acc, C_tot5[[j]], col="blue")
lines(E_acc, C_tot5.5[[j]], col="green")
legend("topright", legend = c("TL 4.5", "TL 5", "TL 5.5"), pch=19, 
       col = c("black","red", "blue"))


#_______________________________________________________________________________
#_______________________________________________________________________________


# ----
}
dev.off()




#--------------------------- SENSITIVITY ANALYSIS ------------------------------
# 1.) TOP-DOWN CONTROL ----

   ### i.) Preparing the vector to loop over and the data frame----
# create a vector of top-down control
topdown=seq(0,0.9,0.1);length(topdown)
names(topdown)<-topdown

formd=seq(0,0.9,0.1);length(formd)
names(formd)<-formd

# creating a new ET_Main for each combination of formD and topD value
exp_topd<-expand.grid(formd,topdown);dim(exp_topd)
names(exp_topd)<-c("FormD","TopD")
exp_topd$name <- paste(exp_topd$FormD,exp_topd$TopD)
dim(exp_topd)


   ### ii.) Run create Etdiagnosis with the vector----
# Now I can run the et diagnosis:

df_topd_alpha<-mapply(create.ETdiagnosis.logsol.integrated, FormD=exp_topd[,1], 
                      TopD=exp_topd[,2],
                      MoreArgs = list(data=dat_bhm12, 
                                      Fish_mort=F_morts_bhm12$Fish_mort,
                                      Fish_mort_acc=F_morts_bhm12$Fish_mort_acc,
                                      Mul_eff=seq(0,3,0.5)), 
                      SIMPLIFY = F, USE.NAMES = T)
length(df_topd_alpha)


df_topd_alpha_P<-mapply(create.ETdiagnosis.logsol.integrated, FormD=exp_topd[,1], 
                        TopD=exp_topd[,2],
                        MoreArgs = list(data=dat_bhm12_P, 
                                        Fish_mort=F_morts_bhm12_P$Fish_mort,
                                        Fish_mort_acc=F_morts_bhm12_P$Fish_mort_acc,
                                        Mul_eff=seq(0,10,0.5)), 
                        SIMPLIFY = F, USE.NAMES = T)
df_topd_alpha_P$`0.1 0`[[7]]$mf
# naming the data frame 
names(df_topd_alpha)<-exp_topd$name
df_topd_alpha$`0 0`$`1.5`$Fish_mort
names(df_topd_alpha_P)<-exp_topd$name
df_topd_alpha_P$`0 0`[[2]]$mf
length(df_topd_alpha$`0 0`)

sum(abs((df_topd_alpha$`0.1 0.1`$`3`$BIOM_MF[2:37]/sum(df_topd_alpha$`0.1 0.1`$`3`$BIOM_MF[2:37]))-
  (df_topd_alpha$`0.1 0.1`$`0`$BIOM_MF[2:37]/sum(df_topd_alpha$`0.1 0.1`$`0`$BIOM_MF[2:37]))))

sum(abs((df_topd_alpha$`0.9 0.9`$`3`$BIOM_MF[2:37]/sum(df_topd_alpha$`0.9 0.9`$`3`$BIOM_MF[2:37]))-
  (df_topd_alpha$`0.9 0.9`$`0`$BIOM_MF[2:37]/sum(df_topd_alpha$`0.9 0.9`$`0`$BIOM_MF[2:37]))))

sum(abs((df_topd_alpha_P$`0.1 0.1`$`3`$BIOM_MF[2:37]/sum(df_topd_alpha_P$`0.1 0.1`$`3`$BIOM_MF[2:37]))-
          (df_topd_alpha_P$`0.1 0.1`$`0`$BIOM_MF[2:37]/sum(df_topd_alpha_P$`0.1 0.1`$`0`$BIOM_MF[2:37]))))

sum(abs((df_topd_alpha_P$`0.9 0.9`$`3`$BIOM_MF[2:37]/sum(df_topd_alpha_P$`0.9 0.9`$`3`$BIOM_MF[2:37]))-
          (df_topd_alpha_P$`0.9 0.9`$`0`$BIOM_MF[2:37]/sum(df_topd_alpha_P$`0.9 0.9`$`0`$BIOM_MF[2:37]))))

   ### iii.) Calculate the performance----
# Now I calculate the performance measure for all runs and remove the last 
# performance measure from the lists

perform_topd_alpha<-c()
perform_topd_alpha_sub<-c()
for (i in seq_along(df_topd_alpha)){
  perform_topd_alpha[[i]]<-perform(data=df_topd_alpha[[i]][[7]], 
                                   VE=df_topd_alpha$`0.1 0`$`0`)
  perform_topd_alpha_sub[[i]]<-perform_topd_alpha[[i]][1:10]
}
names(perform_topd_alpha_sub)<-exp_topd$name;names(perform_topd_alpha_sub)
length(perform_topd_alpha_sub)
perform_topd_alpha_sub$`0.1 0.1`$d
perform_topd_alpha_sub$`0.9 0.9`$d

perform_topd_alpha_P<-c()
perform_topd_alpha_sub_P<-c()
for (i in seq_along(df_topd_alpha_P)){
  perform_topd_alpha_P[[i]]<-perform(data=df_topd_alpha_P[[i]][[7]], 
                                     VE=df_topd_alpha_P$`0.1 0`$`0`)
  perform_topd_alpha_sub_P[[i]]<-perform_topd_alpha_P[[i]][1:10]
}
names(perform_topd_alpha_sub_P)<-exp_topd$name;names(perform_topd_alpha_sub_P)
length(perform_topd_alpha_sub_P[[1]])
length(perform_topd_alpha_sub_P)
names(perform_topd_alpha_sub_P[[1]])
perform_topd_alpha_sub_P$`0.1 0.1`$d
perform_topd_alpha_sub_P$`0.9 0.9`$d


# I collapse the performance measures to a data frame
# as it is easier to work with a data frame later on

perform_topd_alpha_dat <- data.frame(matrix(unlist(perform_topd_alpha_sub), 
                                            nrow=100, byrow = T), stringsAsFactors=F)

colnames(perform_topd_alpha_dat)<-names(perform_topd_alpha[[1]][1:10])
rownames(perform_topd_alpha_dat)<-exp_topd$name
perform_topd_alpha_dat

perform_topd_alpha_dat_P <- data.frame(matrix(unlist(perform_topd_alpha_sub_P), 
                                              nrow=100, byrow = T), stringsAsFactors=F)
dim(perform_topd_alpha_dat_P)

colnames(perform_topd_alpha_dat_P)<-names(perform_topd_alpha_P[[1]][1:10])
rownames(perform_topd_alpha_dat_P)<-exp_topd$name
perform_topd_alpha_dat_P

# add the expand.grid values to reshape each performance measure as a matrix for
# plotting
perform_topd_alpha_dat_new<-cbind(exp_topd[,1:2], perform_topd_alpha_dat)

perform_topd_alpha_dat_new_P<-cbind(exp_topd[,1:2], perform_topd_alpha_dat_P)

# adding the relative predator and forage fish catch to the matrix
perform_topd_alpha_dat_new<- perform_topd_alpha_dat_new %>% 
  mutate(rel.ypred=y.pred/y) %>% 
  mutate(rel.y.forage=y.forage/y)
perform_topd_alpha_dat_new$rel.y.forage

perform_topd_alpha_dat_new_P<- perform_topd_alpha_dat_new_P %>% 
  mutate(rel.ypred=y.pred/y) %>% 
  mutate(rel.y.forage=y.forage/y)
perform_topd_alpha_dat_new_P$rel.y.forage

# 2.) TRANSFER EFFICIENCY ----
 ## A) Set the transfer efficiency vector ----
TE=seq(0.05,0.2, 0.05);length(TE)           # setting the Transfer Efficiency vector
names(TE)<-TE
 ## B) RUN SENSITIVITY ANALYSIS WITH DIFFERENT TE INTENSITIES F~K ----
   ### i.) Create virtual ecosystems ----

# First I create different Virtual Ecosystems using different TE's
VE_bhm12_TE<-lapply(TE,function (a) virtual_eco_integ(TE=a, Temp = 15))

   ### ii.) Create the data frames ---- 

# Then I use these different VE'S to create different data frames
dat_bhm12_TE<-lapply(VE_bhm12_TE,function (a) dat_frame_integrate_access_adap(VE=a, 
                                                                  asymptote = 2, 
                                                        TL50 = 2, slope = 0))

   ### iii.) Create the fishing vectors ---- 
# Create a fishing mortality vector. Since Kinetics are the same, I can only use
# one data frame to calculate the fishing mortality vector. In other words the 
# vector is the same for all TE values. 
F_morts_bhm12_TE<-F_mort_dyn(dat = dat_bhm12_TE$`0.05`, const=0.1)

   ### iv.) Run et.diagnosis ---- 
df_bhm12_TE<-lapply(dat_bhm12_TE, 
                    function (a) create.ETdiagnosis.logsol.integrated(data = a,
                                    Fish_mort=F_morts_bhm12_TE$Fish_mort,
                                    Fish_mort_acc=F_morts_bhm12_TE$Fish_mort_acc,
                                    Mul_eff=seq(0,10,0.5)))

   ### v.) Check the results ----
cl<-viridis(length(df_bhm12_TE$`0.1`))
plot(df_bhm12_TE$`0.1`$`0`$BIOM_MF[-1], type="n")
for (i in seq_along(df_bhm12_TE$`0.1`)){
  lines(df_bhm12_TE$`0.1`[[i]]$BIOM_MF[-1], col=cl[i])
  lines(df_bhm12[[i]]$BIOM_MF[-1], col=cl[i], lty=2)
}
diff<-c()
for (i in seq_along(df_bhm12_TE$`0.1`)){
  diff[[i]]<-df_bhm12_TE$`0.1`[[i]]$BIOM_MF/df_bhm12[[i]]$BIOM_MF
}

 ## C) RUN SENSITIVITY ANALYSIS WITH DIFFERENT TE INTENSITIES F~P ----
   ### i.) Create virtual ecosystems ----
# First I create different Virtual Ecosystems using different TE's
VE_bhm12_TE_P<-lapply(TE,function (a) virtual_eco_integ(TE=a, Temp = 15))
   ### ii.) Create the data frames ----
# Then I use these different VE'S to create different data frames
dat_bhm12_TE_P<-lapply(VE_bhm12_TE_P,function (a) dat_frame_integrate_access_adap(VE=a, 
                                                                  asymptote = 2, 
                                                        TL50 = 2, slope = 0))
length(dat_bhm12_TE_P)

   ### iii.) Create the fishing vectors ----
# Then I use the different data frames to create different fishing mortality vectors
# Because now the production has changed, I also have to change the 
con<-c(0.00806, 0.0039, 0.002548, 0.001885)
F_morts_bhm12_P<-mapply(F_mort_prod, dat=dat_bhm12_TE, const=con, 
                        SIMPLIFY = F, USE.NAMES = T)

round(F_morts_bhm12_P$`0.05`$Fish_mort_acc$catch.1[10]/dat_bhm12_TE$`0.05`$ET_Main$Kin_acc[10], digits = 4)
round(F_morts_bhm12_P$`0.1`$Fish_mort_acc$catch.1[10]/dat_bhm12_TE$`0.1`$ET_Main$Kin_acc[10], digits=4)
round(F_morts_bhm12_P$`0.15`$Fish_mort_acc$catch.1[10]/dat_bhm12_TE$`0.15`$ET_Main$Kin_acc[10], digits=4)
round(F_morts_bhm12_P$`0.2`$Fish_mort_acc$catch.1[10]/dat_bhm12_TE$`0.2`$ET_Main$Kin_acc[10], digits=4)

F_mort_TE<-c()
F_mort_acc_TE<-c()
for (i in seq_along(F_morts_bhm12_P)){
  F_mort_TE[i]<-F_morts_bhm12_P[[i]]["Fish_mort"]
  F_mort_acc_TE[i]<-F_morts_bhm12_P[[i]]["Fish_mort_acc"]
}
names(F_mort_TE)<-names(TE)
names(F_mort_acc_TE)<-names(TE)

   ### iv.) Run et.diagnosis ----
df_bhm12_TE_P_2<-mapply(create.ETdiagnosis.logsol.integrated, 
                      data=dat_bhm12_TE_P, Fish_mort=F_mort_TE,
                      Fish_mort_acc=F_mort_acc_TE,
                      MoreArgs = list(Mul_eff=seq(0,10,0.5)), 
                      SIMPLIFY = F, USE.NAMES = T)


   ### v.) Check the results ----
cl<-viridis(length(df_bhm12_TE_P$`0.1`))
plot(df_bhm12_TE_P$`0.1`$`0`$BIOM_MF[-1], type="n")
for (i in seq_along(df_bhm12_TE_P$`0.1`)){
  lines(df_bhm12_TE_P$`0.1`[[i]]$BIOM_MF[-1], col=cl[i])
  lines(df_bhm12_P[[i]]$BIOM_MF[-1], col=cl[i], lty=2)
}
diff<-c()
for (i in seq_along(df_bhm12_TE_P$`0.1`)){
  diff[[i]]<-df_bhm12_TE_P$`0.1`[[i]]$BIOM_MF/df_bhm12_P[[i]]$BIOM_MF
}


 ## D) PERFORMANCE MEASURE ----
# Now I calculate the performance measure for all runs and remove the last 
# performance measure from the lists
# Calculating the performance measures
TE_perform<-c()
TE_perform_P<-c()
for (i in names(df_bhm12_TE)){
  for (j in names(df_bhm12_TE$`0.05`)){
    TE_perform[[i]][[j]]<-perform(data=df_bhm12_TE[[i]][[j]], 
                                  VE=df_bhm12_TE[[i]][[1]])
    TE_perform_P[[i]][[j]]<-perform(data=df_bhm12_TE_P[[i]][[j]], 
                                    VE=df_bhm12_TE_P[[i]][[1]])
  }}

# 3.) MIDPOINTS ----
 ## A) CALCULATING DIFFERENT MIDPOINTS AND THE DATA FRAMES----

   ### i.) Calculate different midpoints----
midpoints_vec=seq(1,5, 0.5);length(midpoints_vec)
names(midpoints_vec)<-midpoints_vec

   ### ii.) Calculate the Data frame ----
dat_bumss_midpoints<-lapply(midpoints_vec,
                            function (a) dat_frame_integrate_access_adap(TL50=a, 
                                                                VE=VE, slope=4.84,
                                                                 asymptote=1))
length(dat_bumss_midpoints)
dat_bumss_midpoints$`1`$ET_Main$Selec

   ### iii.) Calculate fishing mortalities for F~P/B----
F_morts_midpoints<-lapply(dat_bumss_midpoints,
                          function (a) F_mort_dyn(dat = a, const = 0.1))

F_mort_midpoints<-c()
F_mort_midpoints_acc<-c()
for (i in seq_along(F_morts_midpoints)){
  F_mort_midpoints[i]<-F_morts_midpoints[[i]]["Fish_mort"]
  F_mort_midpoints_acc[i]<-F_morts_midpoints[[i]]["Fish_mort_acc"]
}
names(F_mort_midpoints)<-midpoints_vec

   ### iv.) Calculate fishing mortalities for AH----

# Create the nested list of fishing mortalities from your f_vectors taken from df_alt
# Fish_mort_altacc=c()
Fish_mort_acc_altacc=c()
for (i in seq_along(df_alt)){
  # Fish_mort_altacc[[i]]=list(catch.1=df_alt[[i]]$Fish_mort)
  Fish_mort_acc_altacc[[i]]=list(catch.1=df_alt[[i]]$Fish_mort_acc)
}
# names(Fish_mort_altacc)<-names(df_alt)
names(Fish_mort_acc_altacc)<-names(df_alt)

# Calculate the fishing mortality under the different midpoints 

Fnew.midpoints.alt=function(Fish_mort_acc,dat){
  F[1]<-0
  F[2:52]<-Fish_mort_acc$catch.1[-1]*dat$ET_Main$Selec[-1]
  F_new_acc=list(catch.1=F)
}

F_mort_midpoints_alt<-c()
for (i in names(dat_bumss_midpoints)){
  F_mort_midpoints_alt[[i]]<-lapply(Fish_mort_acc_altacc,
                                    function (a) Fnew.midpoints.alt(Fish_mort_acc = a, 
                                                          dat=dat_bumss_midpoints[[i]]))
}
F_mort_midpoints_alt$`3`$`10`


 ## B) RUN THE ANALYSIS----
   ### i.) F~P/B ----
df_bumss_midpoints<-mapply(create.ETdiagnosis.logsol.integrated, 
                           data = dat_bumss_midpoints, 
                           Fish_mort=F_mort_midpoints, 
                           Fish_mort_acc=F_mort_midpoints_acc,
                           MoreArgs = list(Mul_eff=seq(0,10,0.5)), 
                           SIMPLIFY = F, USE.NAMES = T)

# Check
plot(dat_bumss_midpoints$`2`$ET_Main$Selec)
lines(dat_bumss_midpoints$`2.5`$ET_Main$Selec)
lines(dat_bumss_midpoints$`3`$ET_Main$Selec)
lines(dat_bumss_midpoints$`3.5`$ET_Main$Selec)
lines(dat_bumss_midpoints$`4.5`$ET_Main$Selec)
lines(dat_bumss_midpoints$`5.5`$ET_Main$Selec)

plot(df_bumss_midpoints$`2`$`1`$E, ylim = c(0,0.11))
lines(df_bumss_midpoints$`2.5`$`1`$E)
lines(df_bumss_midpoints$`3`$`1`$E)
lines(df_bumss_midpoints$`3.5`$`1`$E)
lines(df_bumss_midpoints$`4.5`$`1`$E)
lines(df_bumss_midpoints$`5.5`$`1`$E)

   ### ii.) AH ----

df_alt_midpoints<-c()

for (i in names(dat_bumss_midpoints)){
  df_alt_midpoints[[i]]<-mapply(create.ETdiagnosis.logsol.integrated, 
                                Fish_mort_acc=Fish_mort_acc_altacc,
                                Fish_mort=F_mort_midpoints_alt[[i]], 
                                MoreArgs = list(data= dat_bumss_midpoints[[i]], 
                                                Mul_eff=c(1)), 
                                SIMPLIFY = F, USE.NAMES = T)
}

# Check
plot(dat_bumss_midpoints$`2`$ET_Main$Selec)
lines(dat_bumss_midpoints$`2.5`$ET_Main$Selec)
lines(dat_bumss_midpoints$`3`$ET_Main$Selec)
lines(dat_bumss_midpoints$`3.5`$ET_Main$Selec)
lines(dat_bumss_midpoints$`4.5`$ET_Main$Selec)
lines(dat_bumss_midpoints$`5.5`$ET_Main$Selec)

plot(df_alt_midpoints$`1.5`$`1`$`1`$Fish_mort)
lines(df_alt_midpoints$`1`$`1`$`1`$Fish_mort)
lines(df_alt_midpoints$`2.5`$`1`$`1`$Fish_mort)
lines(df_alt$`1`$Fish_mort, col="red")
lines(df_alt_midpoints$`3`$`1`$`1`$Fish_mort)
lines(df_alt_midpoints$`3.5`$`1`$`1`$Fish_mort)
lines(df_alt_midpoints$`4.5`$`1`$`1`$Fish_mort)
lines(df_alt_midpoints$`5.5`$`1`$`1`$Fish_mort)

 ## C) CALCULATE PERFORMANCE MEASURES ----
   ### i.) F~P/B ----
midpoints_performance<-lapply(df_bumss_midpoints, 
                              function(a) lapply(a,
                                    function(a) perform(data=a,
                                                VE=df_bumss_midpoints$`1`$`0`)))
df_bumss_midpoints$`1`$`1`

   ### ii.) AH ----
# removing the effort multiplier 0
df_alt_midpoints_sub<-lapply(df_alt_midpoints,
                             function(a) lapply(a, function(a) a[[1]]))

# Running the performance function
midpoints_performance_alt<-lapply(df_alt_midpoints_sub, 
                                  function(a) lapply(a, 
                                    function(a) perform(data=a,
                                        VE=df_alt_midpoints_sub$`1`$`0`)))


# 4.) MULTI-SCENARIO-RUNS WITHOUT ACCESSIBILITY ----
 ## A) HIGH TL FISHED HARD WITHOUT ACCESSIBILITY F~P/B (fehtla)
   ### i.) Setting the vectors to loop over ----
# setting the vector for TL50 and the asymptote and defining the seq_tl

seq_tl<-c(1,seq(2,7,0.1));length(seq_tl)
seq_tl_2up<-c(seq(2,7,0.1));length(seq_tl_2up)
TL50_vec=seq(1,5, 0.1);length(TL50_vec)
names(TL50_vec)<-TL50_vec

K_vec=seq(0,1,0.1);length(K_vec)    # asymptote
names(K_vec)<-K_vec


# creating a new ET_Main for each combination of TL50 and asymptote value
expand_fehtla<-expand.grid(TL50_vec,K_vec);dim(expand_fehtla)

names(expand_fehtla)<-c("TL50","Asymptote")
expand_fehtla$name <- paste(expand_fehtla$TL50,expand_fehtla$Asymptote)


# with this new combination of TL50 and asymptotes, I create multiple vectors of
# exploitation 
exploit_vec<-mapply(FUN = gl_selec, TL50=expand_fehtla[,1], K=expand_fehtla[,2],
                    MoreArgs = list(x=seq_tl_2up, slope=4.84, p=0, A=0), SIMPLIFY = F)

names(exploit_vec)<-expand_fehtla$name
names(exploit_vec)

   ### ii.) Calculate the data frames ----
# Here, I calculate my data frame. Because we do not assume differences in the
# accessibility per trophic level, as this is a theoretical exploration, I use
# an accessibility of 1 for all, except trophic level 1

dat_fehtla<-dat_frame_integrate_access_adap(VE=VE, asymptote = 2, TL50 = 2, slope = 0)


   ### iii.) Calculate the fishing mortality vector ----
# now I use the multiple vectors of exploitation to calculate the fishing
# mortality over trophic level, using the F_mort_dyn function 

F_morts_fehtla<-lapply(exploit_vec,function (a) F_mort_dyn(const = a, dat = dat_fehtla))
length(F_morts_fehtla)
str(F_morts_fehtla)


# Because, a) I have now a list of 17 (depending on the TL50 vector) each with
# two vectors: fishing mortality and accessible fishing mortality and b) I want to loop
# over each list but taking either of the vectors and c) I don't know how to do this
# in the simple mapply command, I use a for loop to seperate these two vectors. 

F_mort_fehtla<-c()
F_mort_acc_fehtla<-c()
for (i in seq_along(F_morts_fehtla)){
  F_mort_fehtla[i]<-F_morts_fehtla[[i]]["Fish_mort"]
  F_mort_acc_fehtla[i]<-F_morts_fehtla[[i]]["Fish_mort_acc"]
}
names(F_mort_fehtla)<-expand_fehtla$name
names(F_mort_acc_fehtla)<-expand_fehtla$name


   ### iv.) Run et.diagnosis ----
# Now I can run the et diagnosis:
df_fehtla<-mapply(create.ETdiagnosis.logsol.integrated, Fish_mort_acc=F_mort_acc_fehtla, 
                  Fish_mort=F_mort_fehtla, 
                  MoreArgs = list(data= dat_fehtla, Mul_eff=c(0,1)), 
                  SIMPLIFY = F, USE.NAMES = T)

   ### v.) Calculate performance measures ----
# Now I calculate the performance measure for all runs and remove the last 
# performance measure from the lists
   
perform_fehtla<-c()
perform_fehtla_sub<-c()
for (i in seq_along(df_fehtla)){
  perform_fehtla[[i]]<-perform(data=df_fehtla[[i]][[2]], VE=df_fehtla[[i]][[1]])
  perform_fehtla_sub[[i]]<-perform_fehtla[[i]][1:10]
}
names(perform_fehtla_sub)<-expand_fehtla$name
names(perform_fehtla_sub)
length(perform_fehtla_sub)


# I collapse the performance measures to a data frame
# as it is easier to work with a data frame later on

perform_fehtla_dat <- data.frame(matrix(unlist(perform_fehtla_sub), nrow=451 
                                        ,byrow = T
),
stringsAsFactors=F)
dim(perform_fehtla_dat)

colnames(perform_fehtla_dat)<-names(perform_fehtla[[1]][1:10])
rownames(perform_fehtla_dat)<-expand_fehtla$name
names(perform_fehtla_dat)
# add the expand.grid values to reshape each performance measure as a matrix for
# plotting
perform_fehtla_dat_new<-cbind(expand_fehtla[,1:2], perform_fehtla_dat)

# I need 8 matrices with values arranged so that the first value of the rowname 
# becomes the rownames of the new matrix and the second name becomes the columnnames

fehtla_mat<-c()
for (i in names(perform_fehtla_dat_new[,3:10])){
  fehtla_mat[[i]]<-acast(perform_fehtla_dat_new,Asymptote ~ TL50, 
                         value.var=i)
}


 ## B) HIGH TL FISHED HARD AH (fehtla_AH)----
   ### i.) Setting the vectors to loop over ----
# setting the vector for TL50 and the asymptote and defining the seq_tl

seq_tl<-c(1,seq(2,7,0.1));length(seq_tl)
seq_tl_2up<-c(seq(2,7,0.1));length(seq_tl_2up)
TL50_vec_fehtla_acc_Ah=seq(1,5, 0.1);length(TL50_vec_fehtla_acc_Ah)
names(TL50_vec_fehtla_acc_Ah)<-TL50_vec_fehtla_acc_Ah


   ### ii.) Run the baseline data frame used for the extraction of the ----
     #       fishing mortalities ----

# The fishing intensity (The asymptote) is calculated differently:
# We take the accessible fishing mortality vectors and use them directly as 
# fishing intensities and use the TL50 values to calculate new fishing mortalites

dat_alt<-dat_frame_integrate_access_adap(VE=VE, asymptote = 2, TL50 = 2, slope = 0)
dat_alt$ET_Main$Selec

# F_morts_alt<-F_mort_dyn(dat = dat_alt, const=0.043)
F_morts_alt<-F_mort_prod(dat = dat_alt, const=0.0033)

# we only want 11 intensities so we use only 1 step
df_alt_formultirun<-create.ETdiagnosis.alt_logsol.integrated(data=dat_alt, 
                                      Mul_eff=seq(0,10,1), 
                                      Fish_mort = F_morts_alt$Fish_mort, 
                                      Fish_mort_acc = F_morts_alt$Fish_mort_acc)
names(df_alt_formultirun$`1`)

   ### iii.) Extract the accessible fishing mortalities and calculate the ----
     #      total F vectors ----

# Now I create selectivity vectors, that are used to create the actual fishing
# # mortalities.
accessibility_vec<-lapply(TL50_vec_fehtla_acc_Ah,
                          function (i) s_selec(TL50=i, seq_tl_2up, asymptote=1,
                                               slope=4.84))

# store the accessible fishing mortalitiy vector in an object

F_mort_AH<-lapply(1:11, function(i) df_alt_formultirun[[i]][["Fish_mort_acc"]][2:52])
length(F_mort_AH)
names(F_mort_AH)<-names(df_alt_formultirun)
length(F_mort_AH$`1`)

# I calculate the accessible fishing mortality given the different TL50 strategies

F_mort_fehtla_acc_AH1<-lapply(accessibility_vec,function (a) a*F_mort_AH[[1]])
F_mort_fehtla_acc_AH2<-lapply(accessibility_vec,function (a) a*F_mort_AH[[2]])
F_mort_fehtla_acc_AH3<-lapply(accessibility_vec,function (a) a*F_mort_AH[[3]])
F_mort_fehtla_acc_AH4<-lapply(accessibility_vec,function (a) a*F_mort_AH[[4]])
F_mort_fehtla_acc_AH5<-lapply(accessibility_vec,function (a) a*F_mort_AH[[5]])
F_mort_fehtla_acc_AH6<-lapply(accessibility_vec,function (a) a*F_mort_AH[[6]])
F_mort_fehtla_acc_AH7<-lapply(accessibility_vec,function (a) a*F_mort_AH[[7]])
F_mort_fehtla_acc_AH8<-lapply(accessibility_vec,function (a) a*F_mort_AH[[8]])
F_mort_fehtla_acc_AH9<-lapply(accessibility_vec,function (a) a*F_mort_AH[[9]])
F_mort_fehtla_acc_AH10<-lapply(accessibility_vec,function (a) a*F_mort_AH[[10]])
F_mort_fehtla_acc_AH11<-lapply(accessibility_vec,function (a) a*F_mort_AH[[11]])


F_mort_acc_fehtla_AH_all<-list("1"=F_mort_fehtla_acc_AH1, 
                               "2"=F_mort_fehtla_acc_AH2, 
                               "3"=F_mort_fehtla_acc_AH3, "4"=F_mort_fehtla_acc_AH4, 
                               "5"=F_mort_fehtla_acc_AH5, "6"=F_mort_fehtla_acc_AH6, 
                               "7"=F_mort_fehtla_acc_AH7, "8"=F_mort_fehtla_acc_AH8, 
                               "9"=F_mort_fehtla_acc_AH9, "10"=F_mort_fehtla_acc_AH10, 
                               "11"=F_mort_fehtla_acc_AH11)  
length(F_mort_acc_fehtla_AH_all$`11`)
F_mort_acc_fehtla_AH_all <- do.call(c, F_mort_acc_fehtla_AH_all)
names(F_mort_acc_fehtla_AH_all)<- names(df_fehtla)

length(F_mort_acc_fehtla_AH_all)

# # I have the final accessible fishing mortality vectors now and I want to
# # calculate the total fishing mortality given the assumed biological selectivity
# # with a TL50 of 2.5.
#
# # I first must add the F value (0) of the first trophic level to get an F-vector of 52

F_mort_acc_fehtla_AH_all_flat_com<-c()
for (i in names(F_mort_acc_fehtla_AH_all)){
  F_mort_acc_fehtla_AH_all_flat_com[[i]]<-list(catch.1=c(0,F_mort_acc_fehtla_AH_all[[i]]))
}
length(F_mort_acc_fehtla_AH_all_flat_com)
# Now I multiply the accessible fishing mortality with the biological selectivity
# to get the total fishing mortality. 

   ### iv.) Calculate the data frame ----
# Here, I calculate my data frame with FULL ACCESSIBILITY

dat_fehtla_AH<-dat_frame_integrate_access_adap(VE=VE, asymptote = 2, TL50 = 2, slope = 0)

F_mort_acc_fehtla_AH_tot_fishmort<-c()
for (i in names(F_mort_acc_fehtla_AH_all_flat_com)){
  F_mort_acc_fehtla_AH_tot_fishmort[[i]]<-list(catch.1=(F_mort_acc_fehtla_AH_all_flat_com[[i]]$catch.1*
                                                          dat_fehtla_AH$ET_Main$Selec))
}
length(F_mort_acc_fehtla_AH_tot_fishmort)

   ### v.) Run et.diagnosis ----
# Now I have the two fishing mortalities seperated in two different lists
# and can directly loop over them
df_fehtla_AH<-mapply(create.ETdiagnosis.logsol.integrated, 
                     Fish_mort_acc=F_mort_acc_fehtla_AH_all_flat_com, 
                     Fish_mort=F_mort_acc_fehtla_AH_tot_fishmort, 
                     MoreArgs = list(data= dat_fehtla_AH, Mul_eff=c(0,1)), 
                     SIMPLIFY = F, USE.NAMES = T)


names(df_fehtla_AH)

   ### vi.) Calculate performance measures ----
# Now I calculate the performance measure for all runs and remove the last 
# performance measure from the lists

perform_fehtla_AH<-c()
perform_fehtla_sub_AH<-c()
for (i in names(df_fehtla_AH)){
  perform_fehtla_AH[[i]]<-perform(data=df_fehtla_AH[[i]][[2]], 
                                  VE=df_fehtla_AH$`1 0`$`0`)
  perform_fehtla_sub_AH[[i]]<-perform_fehtla_AH[[i]][1:10]
}

names(perform_fehtla_sub_AH)
names(perform_fehtla_sub_AH)
length(perform_fehtla_sub_AH)
length(perform_fehtla_sub_AH$`1.1`)
perform_fehtla_sub_AH$`1 0`$d


# I collapse the performance measures to a data frame
# as it is easier to work with a data frame later on

perform_fehtla_dat_AH <- data.frame(matrix(unlist(perform_fehtla_sub_AH),
                                           nrow=451,byrow = T),
                                    stringsAsFactors=F)

colnames(perform_fehtla_dat_AH)<-names(perform_fehtla_AH[[1]][1:10])
rownames(perform_fehtla_dat_AH)<-names(df_fehtla_AH)
names(perform_fehtla_dat_AH)

# add the expand.grid values to reshape each performance measure as a matrix for
# plotting

# I need to create name vectors for my asyptote and my TL50 because I need to 
# add it to later seperate it into lists
Asymptote=c(0:10);length(Asymptote)
nam_vec<-expand.grid(TL50_vec_fehtla_acc_Ah, Asymptote)
names(nam_vec)<-c("TL50", "Asymptote")
head(names(df_fehtla_AH))
head(nam_vec)


perform_fehtla_dat_new_AH<-cbind(nam_vec, perform_fehtla_dat_AH)
names(perform_fehtla_dat_new_AH)
dim(perform_fehtla_dat_new_AH)



# --------------------------- MANUSCRIPT MAIN FIGURES --------------------------
# 0.) BINDING THE SCENARIOS ----
# We remove one level of the df_al_acc data frame. To have the same levels of 
# nested objects as in the other data frames. 
df_alt_acc_adapt<-list()
for (i in names(df_alt_acc)){
  df_alt_acc_adapt[[i]]<-df_alt_acc[[i]]$`1`
}

scenarios<-list(bhm12=df_bhm12, bhm12_P=df_bhm12_P, bumss=df_bumss, 
                bumss_P=df_bumss_P, AH=df_alt, AH_acc=df_alt_acc_adapt)

# Making a sublist for Figure 4:
scenarios_sub<-list(bhm12=df_bhm12, AH=df_alt)

# Making a sublist for Part II:
scenarios_sub_acc<-list(bumss=df_bumss, AH=df_alt_acc_adapt)

# Names
names<-c()          # storing the effort multipliers in a list to use in legends
for (i in names(df_bhm12)){
  names[i]<-c(df_bhm12[[i]][["mf"]])
}
names_alt<-c()          # storing the effort multipliers in a list to use in legends
for (i in names(df_alt)){
  names_alt[i]<-c(df_alt[[i]][["mf"]])
}

# names<-names[-2]            # Because AH does not contain the mE=0.5, I remove it
# from names to be able to loop over 20 mE instead of 21

# 1.) METHODS ----
 ## A.) SETTING PLOTTING PARAMETER ----
   ### a.) General parameter ----
seq_tl=c(1,seq(2,7,0.1))
min=2                      # trophic level range
max=32
seq_tl[32]

   ### b.) Setting parameter for the multi-scenario F's ----
# BH F~P/B
# we first extract all combinations with an asymptote of 0.4 of an  
F_mort_acc_fehtla_acc_0.4<-F_mort_acc_fehtla_acc[185:230]
names(F_mort_acc_fehtla_acc_0.4)<-seq(1,5.5,0.1)
length(F_mort_acc_fehtla_acc_0.4)

# we now extract onyl five TL50's
F_mort_acc_fehtla_acc_0.4_sub<-F_mort_acc_fehtla_acc_0.4[c(1,11,21,31,41)]

# setting the colors
cl<-viridis(length(F_mort_acc_fehtla_acc_0.4_sub), direction = -1)

 ## B.) ACCESSIBILITY & FISHING MORTALITY GRAPH ----

tiff("Results/Figures/Integrated implementation/Methods/Accessibility+Fishing mortality for multiruns.tif", 
     res=300, compression = "lzw", pointsize = 9,
     height=60, width=169, units="mm",
     family = "sans")

   ### a.) Setting par ----

par(mar = c(4, 4, 2, 1), # Dist' from plot to side of page
    mfrow=c(1,2),
    mgp = c(2.2, 0.7, 0), # Dist' plot to label
    # mgp = c(2, 0.5, 0.3), # Dist' plot to label
    # las = 1, # Rotate y-axis text
    tck = -.02 # Reduce tick length
    # xaxs = "i", yaxs = "i", # Remove plot padding
    # family="serif"
) 
  
   ### b.) plotting ----
plot(seq_tl[min:max], dat_bumss$ET_Main$Selec[min:max], type="n", 
     ylim=c(0,1.05), xlab="Trophic level", yaxp  = c(0, 1, 5), 
     yaxs="i", xaxs="i", las=1, ylab="",
     # ylab=expression("Exploitatble biomass fraction [S"[tau]*"]"), 
     bg="grey"
     #, panel.first = abline(h = seq(0, 1, 0.2), col = "grey80")
     )
mtext(text = expression("Exploitatble"), side=2, 
      las=3, line=3) 
mtext(text = expression("biomass fraction [S"[tau]*"]"), side=2, 
      las=3, line=2) 
lines(seq_tl[min:max], dat_bumss$ET_Main$Selec[min:max], col="black", lwd=1.6)
# adding the 50 % and 90% lines
text(2.75,0.45,labels = "50%")
text(3.25,0.9, labels = "90%")
lines(x=c(2.5,2.5), y=c(0,0.5), lty=2, lwd=1)
lines(x=c(3,3), y=c(0,0.9), lty=2, lwd=1)
# mtext('(a)', side = 3, outer=TRUE, adj = 0.04, line=-1.65, cex=1.1)


# Fishing mortalities of the multi-scenarios
plot(seq_tl[2:max],F_mort_acc_fehtla_acc_0.4_sub$`1`$catch.1[2:max], type="n", 
     ylim=c(0,0.6),  yaxp  = c(0, 0.6, 3),yaxs="i", xaxs="i", las=1,
     xlab="Trophic level",
     ylab="Fishing mortality"
)

for (i in seq_along(F_mort_acc_fehtla_acc_0.4_sub)){
  lines(seq_tl[2:max],F_mort_acc_fehtla_acc_0.4_sub[[i]]$catch.1[2:max], col=cl[i])
}
# mtext('(b)', side = 3, outer=TRUE, adj = 0.56, line=-1.65, cex=1.1)


#----
dev.off()

# 2.) RESULTS ----
 ## PART I: BALANCED HARVEST DOES NOT MAINTAIN ECOSYSTEM STRUCTURE. THE IDEAL -----
#           SOLUTION IS TO EXPLOIT LOWER TROPHIC LEVELS HARDER THAN HIGHER TL'S ----
   ### SETTING PLOTTING PARAMETER ----

seq_tl=c(1,seq(2,7,0.1))
min=2                      # trophic level range
max=32

lmi<-0                     # y limitis of the axis
lma<-2

# colours for the plot
cl<-viridis(length(names), direction = -1)         # mE
cl1<-viridis(length(scenarios), direction = -1)    # scenarios

# Plot Numbering
nam<-c(expression("BH"[P/B]), expression("BSH"[]))
ynam<-c("Biomass [log]", "")
ynam2<-c("Catch", "")

   ### FIG.3 : Combined figure of the exploitation pattern and the and biomass ---- 
#              trophic spectrum along the food web ----

tiff("Results/Figures/Integrated implementation/Part I_BH does not maintain BTS, but AH/B and F of BH and AH.tif", 
     res=300, compression = "lzw", pointsize = 9,
     height=60, width=169, units="mm",
     family = "sans")

      #### a.) Setting par ----

par(mar = c(4, 4, 2, 1), # Dist' from plot to side of page
    mfrow=c(1,2),
    mgp = c(2.7, 0.7, 0), # Dist' plot to label
    # mgp = c(2, 0.5, 0.3), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.02 # Reduce tick length
    # xaxs = "i", yaxs = "i", # Remove plot padding
    # family="serif"
) 

      #### b.) Fishing pattern ----
plot(seq_tl[min:max], df_bhm12$`5`$Fish_mort_acc[min:max], type="n", 
     xlab="Trophic level", ylab="Fishing mortality", bg="grey", 
     ylim=c(0,2),
     yaxs="i", xaxs="i",
     main="")

lines(seq_tl[min:max], df_bhm12$`5`$Fish_mort_acc[min:max], col="black", lty=1)
lines(seq_tl[min:max], df_bhm12_P$`5`$Fish_mort_acc[min:max], col="black", lty=2)
lines(seq_tl[min:max], df_alt$`5`$Fish_mort_acc[min:max], col="grey", lty=1)


      #### c.) Biomass trophic spectrum ----

plot(seq_tl[min:max], scenarios$bhm12$`0`$BIOM_MF[min:max], type="n", log = "y", 
     xlab="Trophic level", ylab="Biomass [log]"
     ,ylim=c(0.1,100), yaxs="i", xaxs="i"
)

lines(seq_tl[min:max], df_bhm12$`5`$BIOM_MF[min:max], col="black", lty=1)
lines(seq_tl[min:max], df_bhm12$`0`$BIOM_MF[min:max], col="red", lty=1)
lines(seq_tl[min:max], df_bhm12_P$`5`$BIOM_MF[min:max], col="black", lty=2)
lines(seq_tl[min:max], df_alt$`5`$BIOM_MF[min:max], col="grey", lty=1)


# ----

dev.off()

   ### FIG: 4 : RUNNING THE NORMAL BST AND, CST PLOT WITH DIFFERENT mE's----

tiff("Results/Figures/Integrated implementation/Part I_BH does not maintain BTS, but AH/B & C of BH and AH.tif", 
     res=300, compression = "lzw", height=120, width=169, units="mm", 
     family = "sans", pointsize = 9)

      #### a.) Setting par ----
par(mar = c(0,0,0,0), # Dist' from plot to side of page
    mfcol=c(2,2),
    mgp = c(1, 0.7, 0), # Dist' plot to label
    mai = c(0.2, 0.2, 0.2, 0.2),
    oma = c(2.5,2.5,2.5,0),
    # mgp = c(2, 0.5, 0.3), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.02 # Reduce tick length
    # xaxs = "i", yaxs = "i", # Remove plot padding
    # family="serif",
    # ,xpd=NA
) 

      #### b.) Plotting ----
for (j in seq_along(scenarios_sub)){
  # Biomass
  plot(seq_tl[min:max], scenarios_sub[[j]]$`0`$BIOM_MF[min:max], type="n", log = "y", 
       xlab="", ylab="", yaxs="i", xaxs="i"
       ,ylim=c(0.1,100), main=nam[[j]], cex.lab=1.2, yaxt="n")
  axis(2,at=c(0.1,1,10,100))
  
  for (i in seq_along(scenarios_sub[[j]])){
    lines(seq_tl[min:max], scenarios_sub[[j]][[i]][["BIOM_MF"]][min:max], col=cl[i])
  }
  lines(seq_tl[min:max], scenarios_sub[[j]][[1]][["BIOM_MF"]][min:max], col="red", 
        lwd=1.3)
  
  # Catch
  plot(seq_tl[min:max], scenarios_sub[[j]]$`0`$Catches.tot[min:max], type="n", 
       xlab="", ylab="", yaxs="i", xaxs="i"
       ,ylim=c(0.1,60), cex.lab=1.2)
  
  for (i in seq_along(scenarios_sub[[j]])){
    lines(seq_tl[min:max], scenarios_sub[[j]][[i]][["Catches.tot"]][min:max], col=cl[i])
  }
}
mtext('Trophic level', side = 1, outer = TRUE, line = 1)
mtext('Biomass [log]', side = 2, outer = TRUE, las=3, adj = 0.81, line=1)
mtext('Catch', side = 2, outer = TRUE, las=3, adj = 0.22, line=1)

#----
dev.off()

  ## PART II: WE SIMULATE A REALISTIC BALANCED HARVEST SCENARIO WHERE NOT ALL----
#          TROPHIC LEVELS ARE ACCESSIBLE ----
   ### SETTING PLOTTING PARAMETER ----
      #### a.) Color, names and trophic levels  
seq_tl=c(1,seq(2,7,0.1))
min=2                      # trophic level range
max=32

# colours for the plot
cl<-viridis(length(names), direction = -1)         # mE
cl1<-viridis(length(scenarios), direction = -1)    # scenarios

# Plot Numbering
nam<-c(expression("BH"[P/B]), expression("BSH"[])) # headers for the plots

       #### b.) Calculating the inaccessible biomass increase for Figure 7. 

# bumss_perform_sub<-bumss_perform[-2] # here I remove the mE 0.5 from the df_bumss

names(alt_acc_perform)
names(bumss_perform)



# here I need to extract the relative amount of inaccessible biomass from the 
# performance data frame

rel_inaccB_bumss<-c()
rel_inaccB_bumss_P<-c()
rel_inaccB_AH_acc<-c()

for (i in seq_along(bumss_perform)){
  for (j in seq_along(alt_acc_perform)){
    rel_inaccB_bumss[[i]]<-bumss_perform[[i]]$rel.total.B.inacc
    rel_inaccB_bumss_P[[i]]<-bumss_perform_P[[i]]$rel.total.B.inacc
    rel_inaccB_AH_acc[[j]]<-alt_acc_perform[[j]]$rel.total.B.inacc
  }}
sum(df_bumss$`0`$BIOM_MF_acc[-1])
sum(df_alt_acc$`1`$`0`$BIOM_MF_acc[-1])
sum(df_bumss_P$`0`$BIOM_MF_acc[-1])

sum(df_bumss$`0`$BIOM_MF[-1])
sum(df_alt_acc$`1`$`0`$BIOM_MF[-1])
df_alt_acc$`1`$`0`$BIOM_MF/dat_alt_acc$ET_Main$B
df_alt_acc$`1`$`0`$BIOM_MF/df_alt_acc$`10`$`0`$BIOM_MF
df_alt_acc$`1`$`0`$BIOM_MF_acc/dat_alt_acc$ET_Main$B_acc
sum(df_bumss_P$`0`$BIOM_MF[-1])

   ### FIG. 6 : Plotting the fishing mortality, biomass and catch trophic spectra, ----
#     as well as the relative inaccessible biomass ----

tiff("Results/Figures/Integrated implementation/Part II_Realistic balanced harvest/F, C$BTS&AccessB_2.tif", 
     res=300, compression = "lzw", height=150, width=169, units="mm",
     family = "sans", pointsize = 9)

      #### a.) Setting par ----

par(mar = c(0,0,0,0), # Dist' from plot to side of page
    mfcol=c(3,2),
    mgp = c(1, 0.7, 0), # Dist' plot to label
    mai = c(0.2, 0.2, 0.2, 0.2),
    oma = c(2.5,2.5,2.5,0),
    # mgp = c(2, 0.5, 0.3), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.02 # Reduce tick length
    # xaxs = "i", yaxs = "i", # Remove plot padding
    # family="serif",
    # ,xpd=NA
) 

      #### b.) Plotting ----

for (j in seq_along(scenarios_sub_acc)){

  # Biomass
  plot(seq_tl[min:max], scenarios_sub_acc[[j]]$`0`$BIOM_MF[min:max], type="n", 
       log = "y", 
       xlab="", ylab="", main=nam[[j]], cex.main=1.5
       ,ylim=c(0.1,100), yaxs="i", xaxs="i",yaxt="n",
       cex.lab=1.4)
  axis(2, at=c(0.1,1,10,100))
  
  for (i in seq_along(names)){
    lines(seq_tl[min:max], scenarios_sub_acc[[j]][[i]][["BIOM_MF"]][min:max], 
          col=cl[i])
  }
  lines(seq_tl[min:max], scenarios_sub_acc[[j]][[1]][["BIOM_MF"]][min:max], 
        col="red", lwd=1.3)
 
  # Accessible Biomass
  plot(seq_tl[min:max], scenarios_sub_acc[[j]]$`0`$BIOM_MF_acc[min:max], 
       type="n", log = "y", yaxt="n",
       xlab="", ylab="", yaxs="i", xaxs="i"
       ,ylim=c(0.1,11),
       cex.lab=1.4)
  axis(2, at=c(0.1,1,10))
  
  for (i in seq_along(names)){
    lines(seq_tl[min:max], scenarios_sub_acc[[j]][[i]][["BIOM_MF_acc"]][min:max], 
          col=cl[i])
  }
  lines(seq_tl[min:max], scenarios_sub_acc[[j]][[1]][["BIOM_MF_acc"]][min:max], 
        col="red", lwd=1.3)

  # Catch
  plot(seq_tl[min:max], scenarios_sub_acc[[j]]$`0`$Catches.tot[min:max], 
       type="n", xlab="", ylab=""
       ,ylim=c(0.1,4), yaxs="i", xaxs="i",
      cex.lab=1.4)
  
  for (i in seq_along(names)){
    lines(seq_tl[min:max], scenarios_sub_acc[[j]][[i]][["Catches.tot"]][min:max], 
          col=cl[i])
  }
}


mtext('Trophic level', side = 1, outer = TRUE, line = 1)
mtext('Biomass [log]', side = 2, outer = TRUE, las=3, adj = 0.9, line=1)
mtext('Accessible Biomass [log]', side = 2, outer = TRUE, las=3, adj = 0.5, line=1)
mtext('Catch', side = 2, outer = TRUE, las=3, adj = 0.14, line=1)

#----
dev.off()

# WITH 15 ME'S

tiff("Results/Figures/Integrated implementation/Part II_Realistic balanced harvest/F, C$BTS&AccessB_2_mE15.tif", 
     res=300, compression = "lzw", height=150, width=169, units="mm",
     family = "sans", pointsize = 9)

# ----
par(mar = c(0,0,0,0), # Dist' from plot to side of page
    mfcol=c(3,2),
    mgp = c(1, 0.7, 0), # Dist' plot to label
    mai = c(0.2, 0.2, 0.2, 0.2),
    oma = c(2.5,2.5,2.5,0),
    # mgp = c(2, 0.5, 0.3), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.02 # Reduce tick length
    # xaxs = "i", yaxs = "i", # Remove plot padding
    # family="serif",
    # ,xpd=NA
) 
# Biomass
plot(seq_tl[min:max], scenarios_sub_acc$bumss$`0`$BIOM_MF[min:max], type="n", 
     log = "y", 
     xlab="", ylab="", main=expression("BH"[P/B]), cex.main=1.5
     ,ylim=c(0.1,100), yaxs="i", xaxs="i",yaxt="n",
     cex.lab=1.4)
axis(2, at=c(0.1,1,10,100))

for (i in seq_along(names)){
  lines(seq_tl[min:max], scenarios_sub_acc$bumss[[i]][["BIOM_MF"]][min:max], 
        col=cl[i])
}
lines(seq_tl[min:max], scenarios_sub_acc$bumss$`0`[["BIOM_MF"]][min:max], 
      col="red", lwd=1.3)


# Accessible Biomass
plot(seq_tl[min:max], scenarios_sub_acc[[j]]$`0`$BIOM_MF_acc[min:max], 
     type="n", log = "y", yaxt="n",
     xlab="", ylab="", yaxs="i", xaxs="i"
     ,ylim=c(0.1,11),
     cex.lab=1.4)
axis(2, at=c(0.1,1,10))

for (i in seq_along(names)){
  lines(seq_tl[min:max], scenarios_sub_acc[[j]][[i]][["BIOM_MF_acc"]][min:max], 
        col=cl[i])
}
lines(seq_tl[min:max], scenarios_sub_acc[[j]][[1]][["BIOM_MF_acc"]][min:max], 
      col="red", lwd=1.3)


# Catch
plot(seq_tl[min:max], scenarios_sub_acc[[j]]$`0`$Catches.tot[min:max], 
     type="n", xlab="", ylab=""
     ,ylim=c(0.1,4), yaxs="i", xaxs="i",
     cex.lab=1.4)

for (i in seq_along(names)){
  lines(seq_tl[min:max], scenarios_sub_acc[[j]][[i]][["Catches.tot"]][min:max], 
        col=cl[i])
}

# Biomass
namesBSH<-names(df_alt_acc15)
clBSH<-viridis(length(df_alt_acc15), direction = -1)
plot(seq_tl[min:max], df_alt_acc15$`0`$`1`$BIOM_MF[min:max], type="n", 
     log = "y", 
     xlab="", ylab="", main="BSH", cex.main=1.5
     ,ylim=c(0.1,100), yaxs="i", xaxs="i",yaxt="n",
     cex.lab=1.4)
axis(2, at=c(0.1,1,10,100))

for (i in seq_along(namesBSH)){
  lines(seq_tl[min:max], df_alt_acc15[[i]]$`1`[["BIOM_MF"]][min:max], 
        col=clBSH[i])
}
lines(seq_tl[min:max], df_alt_acc15$`0`$`1`[["BIOM_MF"]][min:max], 
      col="red", lwd=1.3)

# Accessible Biomass
plot(seq_tl[min:max], df_alt_acc15$`0`$`1`$BIOM_MF_acc[min:max], type="n", 
     log = "y", 
     xlab="", ylab="", cex.main=1.5
     ,ylim=c(0.1,11), yaxs="i", xaxs="i",yaxt="n",
     cex.lab=1.4)
axis(2, at=c(0.1,1,10))

for (i in seq_along(namesBSH)){
  lines(seq_tl[min:max], df_alt_acc15[[i]]$`1`[["BIOM_MF_acc"]][min:max], 
        col=clBSH[i])
}
lines(seq_tl[min:max], df_alt_acc15$`0`$`1`[["BIOM_MF_acc"]][min:max], 
      col="red", lwd=1.3)

# Catch
plot(seq_tl[min:max], df_alt_acc15$`0`$`1`$Catches.tot[min:max], type="n",
     xlab="", ylab="", cex.main=1.5
     ,ylim=c(0,4), yaxs="i", xaxs="i",
     cex.lab=1.4)
for (i in seq_along(namesBSH)){
  lines(seq_tl[min:max], df_alt_acc15[[i]]$`1`[["Catches.tot"]][min:max], 
        col=clBSH[i])
}

mtext('Trophic level', side = 1, outer = TRUE, line = 1)
mtext('Biomass [log]', side = 2, outer = TRUE, las=3, adj = 0.9, line=1)
mtext('Accessible Biomass [log]', side = 2, outer = TRUE, las=3, adj = 0.5, line=1)
mtext('Catch', side = 2, outer = TRUE, las=3, adj = 0.14, line=1)



#----
dev.off()

   ### FIG. 7 : Plotting the increase in inaccessible biomass with different 
#               midpoints ----
tiff("Results/Figures/Integrated implementation/Part II_Realistic balanced harvest/Relative inaccessible Biomass WITH AH.tif", 
     res=300, compression = "lzw", height=60, width=81, units="mm", 
     family = "sans", pointsize = 9)

      #### a.) Setting par ----

par(mar = c(4, 4, 2, 1), # Dist' from plot to side of page
    mfrow=c(1,1),
    mgp = c(2.7, 0.7, 0), # Dist' plot to label
    # mgp = c(2, 0.5, 0.3), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.02 # Reduce tick length
    # xaxs = "i", yaxs = "i", # Remove plot padding
    # family="serif"
)

      #### b.) Plotting ----

plot(names, rel_inaccB_bumss, type = "n", 
     cex.lab=1,
     xlab="Effort multiplier", 
     ylab=expression("Unexploitable biomass [%]"),
     bg="grey", ylim=c(0.4,0.8), yaxs="i", xaxs="i",
     main=""
)
lines(names, rel_inaccB_bumss, lwd=1.5)
lines(names, rel_inaccB_AH_acc, lwd=1.5, lty=2)


#----
dev.off()


   ## PART III: EXPLORATION OF MULTIPLE DIFFERENT FISHING PATTERNS ----
#           given the assumption of accessible and inaccessible biomass 
   ### SETTING PLOTTING PARAMETER ----
      ## i.) General parameter ----
seq_tl=c(1,seq(2,7,0.1))
min=2                      # trophic level range
max=32
lwd=1.5                    # Widths of the lines plotted in the graphs

# colours for the plot
cl<-viridis(length(names), direction = -1)         # mE
cl1<-viridis(length(scenarios), direction = -1)    # scenarios

# colours for the Figure 9
cl2<-viridis(length(df_fehtla_acc_TLI_sub), direction = -1)
cl3<-viridis(length(1:3), direction = -1)

# Plot Numbering
nam<-c(expression("BH"[P/B]), expression("BSH"[])) # headers for the plots

      ## ii.) Transforming the data frames for the contour plots ----

# adding another column wit the y grid (since the TL50 values are different in 
# both scenarios I need to create another y-grid vector) 
length(perform_fehtla_acc_dat_new_AH)
names(perform_fehtla_acc_dat_new_AH)
perform_fehtla_acc_dat_new_grid_AH<-perform_fehtla_acc_dat_new_AH %>% 
  mutate(xgrid= perform_fehtla_acc_dat_new$Asymptote)
dim(perform_fehtla_acc_dat_new_grid_AH)

perform_fehtla_acc_dat_new_grid<-perform_fehtla_acc_dat_new %>% 
  mutate(xgrid= perform_fehtla_acc_dat_new$Asymptote);dim(perform_fehtla_acc_dat_new_grid)

names(perform_fehtla_acc_dat_new_grid_AH)
names(perform_fehtla_acc_dat_new_grid)
perform_fehtla_acc_dat_new_grid_AH$Asymptote/perform_fehtla_acc_dat_new_grid$xgrid

# We combine the two data frames into one

combined_acc<-merge(perform_fehtla_acc_dat_new_grid, 
                    perform_fehtla_acc_dat_new_grid_AH, 
                    by=c("xgrid","TL50"),
                    suffixes =c("_K", "_AH"));names(combined_acc)
dim(combined_acc)


# now transform the data frame from wide to long
data_long_combined_acc <- gather(combined_acc,measure, 
                                 performance,
                                 Asymptote_K:rel.total.B.inacc_AH, 
                                 factor_key=TRUE);names(data_long_combined_acc)
head(data_long_combined_acc)
# we split the row "measure" again to have the scenario and performance measures in
# two seperate columns

data_long_combined_acc_sep<-separate(data = data_long_combined_acc, 
                                     col = "measure", into = c("measure", "scenario"), 
                                     sep = "_")
data_long_combined_acc_sep
# finally we split the data frame into a list of data frames according to their
# performance measures
dat_final_acc<-split(data_long_combined_acc_sep,
                     f = data_long_combined_acc_sep$measure);names(dat_final_acc)
dim(dat_final_acc$d)
# We attach the columns for the isoline
Bpred_vec<-dat_final_acc$Bpred.ch$performance
y_vec<-dat_final_acc$y$performance
B_vec<-dat_final_acc$B.ch$performance
d_vec<-dat_final_acc$d$performance
Binacc_vec<-dat_final_acc$rel.total.B.inacc$performance


for (i in names(dat_final_acc)){
  dat_final_acc[[i]]$Bpred_vec<-Bpred_vec
  dat_final_acc[[i]][["B_vec"]]<-B_vec
  dat_final_acc[[i]][["y_vec"]]<-y_vec
  dat_final_acc[[i]][["d_vec"]]<-d_vec
  dat_final_acc[[i]][["Binacc_vec"]]<-Binacc_vec
}

# I calculate the minimum catch used for the isoline

perform_fehtla_acc_dat_new_AH
perform_fehtla_acc_dat_new
blub_fehtla_acc_sub<-subset(perform_fehtla_acc_dat_new, B.ch>=0.6 & Bpred.ch>=0.6)
dim(blub_fehtla_acc_sub)
C_max_0.5_pos<-which.max(blub_fehtla_acc_sub$y)
C_max_0.5<-blub_fehtla_acc_sub[C_max_0.5_pos,8]
C_thre<-(C_max_0.5/100)*80;C_thre
max(perform_fehtla_acc_dat_new$y)

blub_fehtla_acc_sub_AH<-subset(perform_fehtla_acc_dat_new_AH, B.ch>=0.6 & Bpred.ch>=0.6)
dim(blub_fehtla_acc_sub_AH)
C_max_0.5_pos_AH<-which.max(blub_fehtla_acc_sub_AH$y)
C_max_0.5_AH<-blub_fehtla_acc_sub_AH[C_max_0.5_pos_AH,8]
C_thre_AH<-(C_max_0.5_AH/100)*80;C_thre_AH
max(perform_fehtla_acc_dat_new_AH$y)

      ## iii.) Extracting the trophic level of first catch of choice ----

# Setting the trophic level of choice
TL_inter<-c(2.5,3.5,4.5)
TLI=2.5
TLII=3.5
TLIII=4.5

# Extract the TL50 of choice for the biomass and catch trophic spectrum plots:
# BH
nam<-c()
for (i in TL_inter){
  nam[[i]]<-grep(i, names(df_fehtla_acc), value = TRUE) 
}
df_fehtla_acc_TLI<-df_fehtla_acc[nam[[2]]]
df_fehtla_acc_TLII<-df_fehtla_acc[nam[[3]]]
df_fehtla_acc_TLIII<-df_fehtla_acc[nam[[4]]]

df_fehtla_acc_TLI_sub<-lapply(df_fehtla_acc_TLI, function(x) x[2])
df_fehtla_acc_TLII_sub<-lapply(df_fehtla_acc_TLII, function(x) x[2])
df_fehtla_acc_TLIII_sub<-lapply(df_fehtla_acc_TLIII, function(x) x[2])

# AH
nam<-c()
for (i in TL_inter){
  nam[[i]]<-grep(i, names(df_fehtla_acc_AH), value = TRUE) 
}
df_fehtla_acc_AH_TLI<-df_fehtla_acc_AH[nam[[2]]]
df_fehtla_acc_AH_TLII<-df_fehtla_acc_AH[nam[[3]]]
df_fehtla_acc_AH_TLIII<-df_fehtla_acc_AH[nam[[4]]]

df_fehtla_acc_AH_TLI_sub<-lapply(df_fehtla_acc_AH_TLI, function(x) x[2])
df_fehtla_acc_AH_TLII_sub<-lapply(df_fehtla_acc_AH_TLII, function(x) x[2])
df_fehtla_acc_AH_TLIII_sub<-lapply(df_fehtla_acc_AH_TLIII, function(x) x[2])


# Extract the TL50 of choice for increase in inaccessible biomass:
# BH~K
namI<-grep(TLI, names(perform_fehtla_acc_sub), value = TRUE) 
perform_fehtla_acc_sub_TLI<-perform_fehtla_acc_sub[namI]

namII<-grep(TLII, names(perform_fehtla_acc_sub), value = TRUE) 
perform_fehtla_acc_sub_TLII<-perform_fehtla_acc_sub[namII]

namIII<-grep(TLIII, names(perform_fehtla_acc_sub), value = TRUE) 
perform_fehtla_acc_sub_TLIII<-perform_fehtla_acc_sub[namIII]

# bind into list
perform_fehtla_acc_sub_TLALL<-list(TLI=perform_fehtla_acc_sub_TLI, 
                                   TLII=perform_fehtla_acc_sub_TLII,
                                   TLIII=perform_fehtla_acc_sub_TLIII)

# AH
namI_K<-grep(TLI, names(perform_fehtla_acc_sub_AH), value = TRUE) 
perform_fehtla_acc_sub_AH_TLI<-perform_fehtla_acc_sub_AH[namI_K]

namII_K<-grep(TLII, names(perform_fehtla_acc_sub_AH), value = TRUE) 
perform_fehtla_acc_sub_AH_TLII<-perform_fehtla_acc_sub_AH[namII_K]

namIII_K<-grep(TLIII, names(perform_fehtla_acc_sub_AH), value = TRUE) 
perform_fehtla_acc_sub_AH_TLIII<-perform_fehtla_acc_sub_AH[namIII_K]

perform_fehtla_acc_sub_AH_TLALL<-list(TLI=perform_fehtla_acc_sub_AH_TLI, 
                                      TLII=perform_fehtla_acc_sub_AH_TLII,
                                      TLIII=perform_fehtla_acc_sub_AH_TLIII)
names(perform_fehtla_acc_sub_AH_TLALL$TLI)

      ## iv.) Calculting the inaccessible biomass increase ----

# extracting the fishing intensity
mul_eff<-K_vec_fehtla_acc
mul_eff_AH<-c(0:10)


# here I need to extract the relative amount of inaccessible biomass from the 
# performance data frame

# BH ~ K
rel_inaccB_bumss_ALL<-list()
for (i in names(perform_fehtla_acc_sub_TLALL)){
  for (j in names(perform_fehtla_acc_sub_TLALL[[i]])){
    rel_inaccB_bumss_ALL[[i]][[j]]<-perform_fehtla_acc_sub_TLALL[[i]][[j]]$rel.total.B.inacc
  }
}

# BH ~ AH
rel_inaccB_alt_acc_ALL<-list()
for (i in names(perform_fehtla_acc_sub_AH_TLALL)){
  for (j in names(perform_fehtla_acc_sub_AH_TLALL[[i]])){
    rel_inaccB_alt_acc_ALL[[i]][[j]]<-perform_fehtla_acc_sub_AH_TLALL[[i]][[j]]$rel.total.B.inacc
  }
}
rel_inaccB_alt_acc_ALL$TLI

      ## v.) Adding a star to the best strategies ----
blub_fehtla_acc_sub<-subset(perform_fehtla_acc_dat_new, B.ch>=0.6 & Bpred.ch>=0.6)
C_max_0.5_pos<-which.max(blub_fehtla_acc_sub$y)
C_max_0.5<-blub_fehtla_acc_sub[C_max_0.5_pos,]
Cpred_max_0.5_pos<-which.max(blub_fehtla_acc_sub$y.pred)
Cpred_max_0.5<-blub_fehtla_acc_sub[Cpred_max_0.5_pos,]

blub_fehtla_acc_sub_AH<-subset(perform_fehtla_acc_dat_new_AH, B.ch>=0.6 & Bpred.ch>=0.6)
C_max_0.5_pos_AH<-which.max(blub_fehtla_acc_sub_AH$y)
C_max_0.5_AH<-blub_fehtla_acc_sub_AH[C_max_0.5_pos_AH,]
Cpred_max_0.5_pos_AH<-which.max(blub_fehtla_acc_sub_AH$y.pred)
Cpred_max_0.5_AH<-blub_fehtla_acc_sub_AH[Cpred_max_0.5_pos_AH,]

pred_best<-blub_fehtla_acc_sub[Cpred_max_0.5_pos,]
w<-subset(dat_final_acc$Bpred.ch, xgrid==0.8 & TL50==3.7 & scenario=="BH[P/B]")
dat_final_acc$Bpred.ch[356,]

tot_best<-blub_fehtla_acc_sub[C_max_0.5_pos,];tot_best
w<-subset(dat_final_acc$Bpred.ch, TL50==1.0)
dat_final_acc$Bpred.ch[124,]

pred_best_AH<-blub_fehtla_acc_sub_AH[Cpred_max_0.5_pos_AH,]
tot_best<-blub_fehtla_acc_sub[C_max_0.5_pos,]
tot_best_AH<-blub_fehtla_acc_sub_AH[C_max_0.5_pos_AH,]

   ### FIG. 8 : Plotting the multiple run ----

tiff("Results/Figures/Integrated implementation/Part III_Multiple scenario runs/Fig.8.tif", 
     res=300, compression = "lzw", height=220, width=169, units="mm", 
     family = "sans", pointsize = 9)

       #### b.) Plotting ----
# now we can run the plots so that the performance measures of each scenario
# share the same legend (facet_wrap)

for (i in seq_along(dat_final_acc)){
  dat_final_acc[[i]]$scenario<- as.factor(dat_final_acc[[i]]$scenario)
levels(dat_final_acc[[i]]$scenario) <- c("BSH", "BH[P/B]")
}

p1_acc<-ggplot(dat_final_acc$d, aes(xgrid, TL50)) + 
  stat_contour(geom = 'polygon', aes(z=performance, fill = ..level..)) +
  geom_tile(aes(fill=performance))+
  viridis::scale_fill_viridis(direction = -1,discrete = F)+
  geom_contour(aes(z = performance), size=0.2) +
  scale_color_continuous(low = "white", high = "white", guide=F) +
  geom_contour(aes(z = Bpred_vec, colour = ..level..), breaks = c(0.6),
               color = "black", cex=1)+
  geom_contour(aes(z = y_vec, colour = ..level..), breaks = c(19), 
               color = "red", cex=1)+
  facet_grid(~scenario)+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_blank())+
  theme(legend.title=element_blank(), strip.background = element_blank(),
        strip.text.x = element_blank())+
  theme(axis.title.x= element_blank(), 
        axis.title.y= element_text(family='sans', size = 12))+
  theme(axis.text = element_text(family='sans', size = 9))+
  ylab(expression(tau[50]))+
  geom_point(data = dat_final_acc$Bpred.ch[356,], shape=8, size=2, color="blue")+
  geom_point(data = dat_final_acc$Bpred.ch[124,], shape=8, size=2, color="black")
  # facet_wrap_custom(~scenario, scales = "free", labeller=labeller(scenario=labs),
  #                   scale_overrides = list(
  #                     scale_override(2, scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
  #                                                          labels = c(0,0.25,0.50,0.75,1.00)))))


p2_acc<-ggplot(dat_final_acc$rel.total.B.inacc, aes(xgrid, TL50)) + 
  stat_contour(geom = 'polygon', aes(z=performance, fill = ..level..)) +
  geom_tile(aes(fill=performance))+
  viridis::scale_fill_viridis(direction = -1,discrete = F)+
  geom_contour(aes(z = performance), size=0.2) +
  scale_color_continuous(low = "white", high = "white", guide=F) +
  geom_contour(aes(z = Bpred_vec, colour = ..level..), breaks = c(0.6), 
               color = "black", cex=1)+
  geom_contour(aes(z = y_vec, colour = ..level..), breaks = c(19), 
               color = "red", cex=1)+
  facet_grid(~scenario)+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_blank())+
  theme(legend.title=element_blank(), strip.background = element_blank(),
        strip.text.x = element_blank())+
  theme(axis.title.x= element_blank(), 
        axis.title.y= element_text(family='sans', size = 12))+
  theme(axis.text = element_text(family='sans', size = 9))+
  ylab(expression(tau[50]))+
  geom_point(data = dat_final_acc$Bpred.ch[356,], shape=8, size=2, color="blue")+
  geom_point(data = dat_final_acc$Bpred.ch[124,], shape=8, size=2, color="black")
  # facet_wrap_custom(~scenario, scales = "free", labeller=labeller(scenario=labs),
  #                   scale_overrides = list(
  #                     scale_override(2, scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
  #                                                          labels = c(0,0.25,0.50,0.75,1.00)))))


p3_acc <- p2_acc %+% dat_final_acc$B.ch
p4_acc <- p2_acc %+% dat_final_acc$Bpred.ch
p5_acc <- p2_acc %+% dat_final_acc$y
p6_acc <- p2_acc %+% dat_final_acc$y.pred


grid.arrange(
  # bottom="Asymptote (Fishing intensity)", 
  bottom=textGrob("Asymptote (Fishing intensity)", gp=gpar(fontsize=11)), ncol=1, nrow=6,
             top = grid::textGrob("",x=0,hjust=0),
                 arrangeGrob(p1_acc, left=textGrob("Disturbance index", gp=gpar(fontsize=11), 
                                               rot = 90)),
             arrangeGrob(p2_acc, left=textGrob("Rel. exploitable B", gp=gpar(fontsize=11), 
                                               rot = 90)),
             arrangeGrob(p3_acc, left=textGrob("B change", gp=gpar(fontsize=11), 
                                               rot = 90)),
             arrangeGrob(p4_acc, left=textGrob("Predator B change", gp=gpar(fontsize=11), 
                                               rot = 90)),
             arrangeGrob(p5_acc, left=textGrob("Total catch", gp=gpar(fontsize=11), 
                                               rot = 90)),
             arrangeGrob(p6_acc, left=textGrob("Predator catch", gp=gpar(fontsize=11), 
                                               rot = 90))
             
             
)
grid.text("BSH", x = unit(0.28, "npc"), 
          y = unit(0.99, "npc"),gp = gpar(fontsize=12))
grid.text(expression("BH"[P/B]), x = unit(0.69, "npc"), 
          y = unit(0.99, "npc"),gp = gpar(fontsize=12))

# 
# # top = grid::textGrob("Title",x=0,hjust=0)
# grid.arrange(arrangeGrob(p,p,p,p,top=textGrob("Sample Title One"),   
#                          ncol=1), arrangeGrob(p,p,p,p,top=textGrob("Sample Title Two"), ncol=1),   
#              ncol = 2)

#-----
dev.off()



# ---------------------- MANUSCRIPT SUPPLEMENTARY FIGURES ----------------------
# 0.) BINDING THE SCENARIOS
# 1.) SUPPLEMENTARY PART I ----
 ## A.) SETTING PLOTTING PARAMETER & DATA TRANSFORMATION ----
   ### i.) Plotting parameter ----
col<-viridis(length(df_bhm12_TE), direction = -1)
mE<-names(df_bhm12_TE$`0.05`)
seq_tl<-c(1,seq(2,7,0.1))

   ### ii.) TOPD : Transforming the data to be used in ggplot2 contour plots ----

perform_topd_alpha_dat_new
names(perform_topd_alpha_dat_new)
dim(perform_topd_alpha_dat_new)

perform_topd_alpha_dat_new_P
names(perform_topd_alpha_dat_new_P)
dim(perform_topd_alpha_dat_new_P)

# now transform the data frame from wide to long
data_long_topd_alpha <- gather(perform_topd_alpha_dat_new,measure, 
                               performance,
                               d:rel.y.forage, 
                               factor_key=TRUE);names(data_long_topd_alpha)
head(data_long_topd_alpha)

data_long_topd_alpha_P <- gather(perform_topd_alpha_dat_new_P,measure, 
                                 performance,
                                 d:rel.y.forage, 
                                 factor_key=TRUE);names(data_long_topd_alpha_P)
head(data_long_topd_alpha_P)

# finally we split the data frame into a list of data frames according to their
# performance measures
dat_final_topd_alpha<-split(data_long_topd_alpha,
                            f = data_long_topd_alpha$measure);names(dat_final_topd_alpha)
dim(dat_final_topd_alpha$d)
names(dat_final_topd_alpha$d)

dat_final_topd_alpha_P<-split(data_long_topd_alpha_P,
                              f = data_long_topd_alpha_P$measure);names(dat_final_topd_alpha_P)
dim(dat_final_topd_alpha_P$d)
names(dat_final_topd_alpha_P$d)
which.max(dat_final_topd_alpha_P$d$performance)
dat_final_topd_alpha_P$d[100,]

   ### iii.) TE : Calculate the ratio between biomass of mE and virgin ecosystem ----
mE_choice<-11

rel_B<-c()
rel_B_P<-c()
for (i in names(df_bhm12_TE)){
  rel_B[[i]]<-df_bhm12_TE[[i]][[mE_choice]]$BIOM_MF/df_bhm12_TE[[i]]$`0`$BIOM_MF
  rel_B_P[[i]]<-df_bhm12_TE_P[[i]][[mE_choice]]$BIOM_MF/df_bhm12_TE_P[[i]]$`0`$BIOM_MF
}

for (i in names(df_bhm12_TE)){
rel_B[[i]]<-df_bhm12_TE[[i]][[1]][["BIOM_MF"]]/
  sum(df_bhm12_TE[[i]][[1]][["BIOM_MF"]])-
  df_bhm12_TE[[i]][[mE_choice]][["BIOM_MF"]]/
    sum(df_bhm12_TE[[i]][[mE_choice]][["BIOM_MF"]])

rel_B_P[[i]]<-df_bhm12_TE_P[[i]][[mE_choice]][["BIOM_MF"]]/
  sum(df_bhm12_TE_P[[i]][[mE_choice]][["BIOM_MF"]])
}


   ### iv.) TE : Calculate the relative predator & forage fish catch and extract ----
      #           the dispersion value for each mE ----
# Calculate the relative predator and low trophic level catch & extract the dispersion

d<-list()
rel.predy<-list()
rel.forage<-list()
y<-list()
ypred<-list()
yforage<-list()
d_P<-list()
rel.predy_P<-list()
rel.forage_P<-list()
y_P<-list()
ypred_P<-list()
yforage_P<-list()
for (i in names(TE_perform)){
  for (j in names(TE_perform$`0.05`)){
    d[[i]][[j]]<-(TE_perform[[i]][[j]]$d)
    y[[i]][[j]]<-(TE_perform[[i]][[j]]$y)
    ypred[[i]][[j]]<-(TE_perform[[i]][[j]]$y.pred)
    yforage[[i]][[j]]<-(TE_perform[[i]][[j]]$y.forage)
    rel.predy[[i]][[j]]<-round(ypred[[i]][[j]]/y[[i]][[j]], digits = 4)
    rel.forage[[i]][[j]]<-round(yforage[[i]][[j]]/y[[i]][[j]], digits = 4)
    
    d_P[[i]][[j]]<-(TE_perform_P[[i]][[j]]$d)
    y_P[[i]][[j]]<-(TE_perform_P[[i]][[j]]$y)
    ypred_P[[i]][[j]]<-(TE_perform_P[[i]][[j]]$y.pred)
    yforage_P[[i]][[j]]<-(TE_perform_P[[i]][[j]]$y.forage)
    rel.predy_P[[i]][[j]]<-round(ypred_P[[i]][[j]]/y_P[[i]][[j]], digits = 4)
    rel.forage_P[[i]][[j]]<-round(yforage_P[[i]][[j]]/y_P[[i]][[j]], digits = 4)
  }}

 ## FIG S1 : Effect of topD and FormD on the dispersion, relative lower and predatory catch ----

tiff("Results/Figures/Integrated implementation/Supplementary/Dispersion_topd&formd.tif", 
     res=300, compression = "lzw", height=120, width=169, units="mm", 
     family = "sans", pointsize = 9)

   ### B) Plotting ----
# now we can run the plots so that the performance measures of each scenario
# share the same legend (facet_wrap)

p1_acc<-ggplot(dat_final_topd_alpha$d, aes(FormD, TopD)) + 
  stat_contour(geom = 'polygon', aes(z=performance, fill = ..level..)) +
  geom_tile(aes(fill=performance))+
  viridis::scale_fill_viridis(direction = -1,discrete = F, begin = 0, end=0.6)+
  geom_contour(aes(z = performance), size=0.2) +
  scale_color_continuous(low = "white", high = "white", guide=F) +
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_blank())+
  theme(legend.title=element_blank(), strip.background = element_blank(),
        strip.text.x = element_blank())+
  theme(axis.title.x= element_text(family='sans', size = 10), 
        axis.title.y= element_text(family='sans', size = 10))+
  theme(axis.text = element_text(family='sans', size = 9))+
  ylab(expression(beta))+
  xlab(expression(gamma))+
  theme(plot.margin=unit(c(0.3,0.1,0.1,0.5), "cm"))


p2_acc<-ggplot(dat_final_topd_alpha$rel.ypred, aes(FormD, TopD)) + 
  stat_contour(geom = 'polygon', aes(z=performance, fill = ..level..)) +
  geom_tile(aes(fill=performance))+
  viridis::scale_fill_viridis(direction = -1,discrete = F, begin = 0, end=0.6)+
  geom_contour(aes(z = performance), size=0.2) +
  scale_color_continuous(low = "white", high = "white", guide=F) +
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_blank())+
  theme(legend.title=element_blank(), strip.background = element_blank(),
        strip.text.x = element_blank())+
  theme(axis.title.x= element_text(family='sans', size = 10), 
        axis.title.y= element_text(family='sans', size = 10))+
  theme(axis.text = element_text(family='sans', size = 9))+
  ylab(expression(beta))+
  xlab(expression(gamma))+
  theme(plot.margin=unit(c(0.3,0,0.1,0.1), "cm"))

p4_acc <- p1_acc %+% dat_final_topd_alpha_P$d
p5_acc <- p2_acc %+% dat_final_topd_alpha_P$rel.ypred


grid.arrange(ncol=2,nrow=2,top = grid::textGrob("",x=0,hjust=0),
             left = grid::textGrob("",just=0,y=0.5), p1_acc, p2_acc, p4_acc, p5_acc)
grid.text("Disturbance index", x = unit(0.23, "npc"), 
          y = unit(0.99, "npc"),gp = gpar(fontsize=10))
grid.text("Relative predatory catch", x = unit(0.71, "npc"), 
          y = unit(0.99, "npc"),gp = gpar(fontsize=10))
grid.text(expression("BH"[P/B]), x = unit(0.01, "npc"), 
          y = unit(0.76, "npc"),gp = gpar(fontsize=10), rot = 90)
grid.text(expression("BH"[P]), x = unit(0.01, "npc"), 
          y = unit(0.27, "npc"),gp = gpar(fontsize=10), rot = 90)


# ----

dev.off()

 ## FIG S2 : Effect of TE on the dispersion, relative lower and predatory catch ----

tiff("Results/Figures/Integrated implementation/Supplementary/TE.tif", 
     res=300, compression = "lzw", height=120, width=169, units="mm", 
     family = "sans", pointsize = 9)

   ### i.) Setting par ----

par(oma = c(2.5,2.5,2.5,0),
    mar = c(2,2,2,2), # Dist' from plot to side of page
    mfrow=c(2,2),
    mgp = c(1.2, 0.7, 0), # Dist' plot to label
    # mgp = c(2, 0.5, 0.3), # Dist' plot to label
    # mai = c(0.2, 0.2, 0.2, 0.2),
    las = 1, # Rotate y-axis text
    tck = -.02 # Reduce tick length
    # xaxs = "i", yaxs = "i", # Remove plot padding
    # family="serif"
) 

   ### B) Plotting ----
# Dispersion plot F~P/B
plot(mE, d$`0.05`, ylim=c(0,0.5), xlab="", 
     main="Disturbance index", type="n", yaxs="i", xaxs="i", ylab="")
for (i in seq_along(df_bhm12_TE)){
  lines(mE, d[[i]], col=col[i])
}
# mtext("(a)", side = 3, outer = TRUE, adj = 0.04, line=-1.65, cex=1.1)

# Relative predatory catch F~P/B
plot(mE, rel.predy$`0.05`, type="n", ylim=c(0,0.1), xlab="", 
     main="Relative predatory catch", yaxs="i", xaxs="i", ylab="")
for (i in seq_along(rel.predy)){
  lines(mE, rel.predy[[i]], col=col[i])
}
# mtext("(b)", side = 3, outer = TRUE, adj = 0.56, line=-1.65, cex=1.1)


# Dispersion plot F~P
plot(mE, d_P$`0.05`, type="n", ylim=c(0,0.3), xlab="", 
     ylab="", yaxs="i", xaxs="i")
for (i in seq_along(df_bhm12_TE_P)){
  lines(mE, d_P[[i]], col=col[i])
}

# mtext('(d)', side = 3, outer=TRUE, adj = 0.56, line=-23, cex=1.1)

# Relative predatory catch F~P
plot(mE, rel.predy_P$`0.05`, type="n", ylim=c(0,0.05), xlab="", 
     ylab="", yaxs="i", xaxs="i")
for (i in seq_along(rel.predy_P)){
  lines(mE, rel.predy_P[[i]], col=col[i])
}
# mtext('(e)', side = 3, outer=TRUE, adj = 0.04, line=-43, cex=1.1)
mtext('Effort Multiplier', side = 1, outer = TRUE, line = 1)
mtext(expression('BH'[P/B]), side = 2, outer = TRUE, las=3, adj = 0.81, line=1)
mtext(expression('BH'[P]), side = 2, outer = TRUE, las=3, adj = 0.22, line=1)



# ----

dev.off()

# 2.) SUPPLEMENTARY PART II ----
 ## A.) SETTING PLOTTING PARAMETER & DATA TRANSFORMATION ----
   ### i.) Setting plotting parameter----

# I set the mE vector
mE<-names(rel.total.B.inacc$`1`)

cl1<-viridis(length(rel.total.B.inacc), direction = -1)

   ### i) Transforming the data for plotting
   ### ii.) Relative inaccessible biomass ----
# Because the relative inaccessible biomass depends on the assumption of the 
# midpoint, I am calculating the increase relative to the baseline (Virgin biomass)
rel.total.B.inacc<-lapply(midpoints_performance, 
                      function(a) lapply(a, function (a) a[["rel.total.B.inacc"]]))

rel.total.B.inacc_alt<-lapply(midpoints_performance_alt, 
                          function(a) lapply(a, function (a) a[["rel.total.B.inacc"]]))

   ### iii.) Relative dispersion ----
# Extract the dispersion values for each midpoint
d_midpoints<-lapply(midpoints_performance, 
                            function(a) sapply(a, function (a) a[["d"]]))
midpoints_performance$`1`$`1`$d
d_midpoints$`1`

d_midpoints_alt<-lapply(midpoints_performance_alt, 
                                function(a) sapply(a, function (a) a[["d"]]))

   ### iv.) Total catch ----
# Extract the total catch for each midpoint
tot_catch_midpoints<-lapply(midpoints_performance, 
                          function(a) sapply(a, function (a) a[["y"]]))
midpoints_performance$`1`$`0`$y

tot_catch_midpoints_alt<-lapply(midpoints_performance_alt, 
                              function(a) sapply(a, function (a) a[["y"]]))
   ### v.) Predator catch ----
# Extract the total catch for each midpoint
predcatch_midpoints<-lapply(midpoints_performance, 
                            function(a) sapply(a, function (a) a[["y.pred"]]))
midpoints_performance$`1`$`0`$y.pred

predcatch_midpoints_alt<-lapply(midpoints_performance_alt, 
                                function(a) sapply(a, function (a) a[["y.pred"]]))


 ## FIG. S3 : Relative inaccessible biomass given different midpoints ----
tiff("Results/Figures/Integrated implementation/Supplementary/Part_II_Different midpoints.tif", 
     res=300, compression = "lzw", height=210, width=169, units="mm", 
     family = "sans", pointsize = 9)
    
    ### a.) Setting par ----
par(mar = c(0, 0, 0, 0), # Dist' from plot to side of page
    mfrow=c(4,2),
    mgp = c(2.7, 0.7, 0), # Dist' plot to label
    mai = c(0.2, 0.4, 0.2, 0.2),
    oma = c(2.5,2.5,2.5,0),
    # mgp = c(2, 0.5, 0.3), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.02 # Reduce tick length
    # xaxs = "i", yaxs = "i", # Remove plot padding
    # family="serif"
) 


    ### b.) Plotting----
# Relative increase of inaccessible biomass
plot(mE, rel.total.B.inacc$`1`, type="n", 
     xlab="", ylab="Relative unexploitable biomass", bg="grey",
     yaxs="i", xaxs="i", cex.main=2, cex.lab=1.5,cex.axis=1.2, ylim=c(0,1))
for (i in seq_along(rel.total.B.inacc)){
  lines(mE, rel.total.B.inacc[[i]],col=cl1[i],lwd=1)
}
# mtext("(a)", side = 3, outer = TRUE, adj = 0.04, line=-1.65, cex=1.1)

plot(mE, rel.total.B.inacc_alt$`1`, type="n", 
     xlab="", ylab="", bg="grey",yaxs="i", xaxs="i",
     main="", cex.main=2, cex.lab=1.5, cex.axis=1.2, ylim=c(0,1))
for (i in seq_along(rel.total.B.inacc_alt)){
  lines(mE, rel.total.B.inacc_alt[[i]],col=cl1[i],lwd=1)
}
# mtext("(b)", side = 3, outer = TRUE, adj = 0.56, line=-1.65, cex=1.1)


# Dispersion
plot(mE, d_midpoints$`1`, type="n", 
     xlab="", ylab="Distubrance index"
     ,ylim=c(0,0.4), yaxs="i", xaxs="i",
     cex.lab=1.5, cex.axis=1.2)
for (i in seq_along(d_midpoints)){
  lines(mE, d_midpoints[[i]], col=cl1[i])
}
# mtext('(c)', side = 3, outer=TRUE, adj = 0.04, line=-23, cex=1.1)


plot(mE, d_midpoints_alt$`1`, type="n", 
     xlab="", ylab=""
     ,ylim=c(0,0.2),yaxs="i", xaxs="i",  cex.lab=1.5, cex.axis=1.2)
for (i in seq_along(d_midpoints_alt)){
  lines(mE, d_midpoints_alt[[i]], col=cl1[i])
}
# mtext('(c)', side = 3, outer=TRUE, adj = 0.04, line=-23, cex=1.1)


# Total catch
plot(mE, tot_catch_midpoints$`1`, type="n",xlab="", ylab="Total catch"
     ,ylim=c(0,220),yaxs="i", xaxs="i",  cex.lab=1.5, cex.axis=1.2)
for (i in seq_along(tot_catch_midpoints)){
  lines(mE, tot_catch_midpoints[[i]], col=cl1[i])
}
# mtext('(e)', side = 3, outer=TRUE, adj = 0.04, line=-43, cex=1.1)

plot(mE, tot_catch_midpoints_alt$`1`, type="n",xlab="", ylab=""
     ,ylim=c(0,220),yaxs="i", xaxs="i", cex.lab=1.5, cex.axis=1.2)
for (i in seq_along(tot_catch_midpoints_alt)){
  lines(mE, tot_catch_midpoints_alt[[i]], col=cl1[i])
}
# mtext('(f)', side = 3, outer=TRUE, adj = 0.56, line=-43, cex=1.1)


# Predatory catch 
plot(mE, predcatch_midpoints$`1`, type="n",xlab="", ylab="Predatory catch"
     ,ylim=c(0,3),yaxs="i", xaxs="i",  cex.lab=1.5, cex.axis=1.2)
for (i in seq_along(predcatch_midpoints)){
  lines(mE, predcatch_midpoints[[i]], col=cl1[i])
}
# mtext('(g)', side = 3, outer=TRUE, adj = 0.04, line=-63, cex=1.1)

plot(mE, predcatch_midpoints_alt$`1`, type="n",xlab="", ylab=""
     ,ylim=c(0,1),yaxs="i", xaxs="i",  cex.lab=1.5, cex.axis=1.2)
for (i in seq_along(predcatch_midpoints_alt)){
  lines(mE, predcatch_midpoints_alt[[i]], col=cl1[i])
}
# mtext('(h)', side = 3, outer=TRUE, adj = 0.56, line=-63, cex=1.1)

mtext('Effort Multiplier', side = 1, outer = TRUE, line = 1, cex=1.1)
mtext(expression('BH'[P/B]), side = 3, outer = TRUE, las=1, adj = 0.22, line=-0.5, 
      cex = 1.4)
mtext(expression('BH'[P]), side = 3, outer = TRUE, las=1, adj = 0.81, line=-0.5, 
      cex=1.4)

#----
dev.off()

# 3.) SUPPLEMENTARY PART III ----
 ## A.) SETTING PLOTTING PARAMETER & DATA TRANSFORMATION ----
   ### i.) Plotting parameter
   ### ii.) Transforming the data to be used in ggplot2 ----
head(perform_fehtla_dat_new)
names(perform_fehtla_dat_new)
perform_fehtla_dat_new$Asymptote
head(perform_fehtla_dat_new_AH)
names(perform_fehtla_dat_new_AH)
perform_fehtla_dat_new_AH$Asymptote

# adding another column with the y grid (since the TL50 values are different in 
# both scenarios I need to create another y-grid vector) 
perform_fehtla_dat_new_grid_AH<-perform_fehtla_dat_new_AH %>% 
  mutate(xgrid= perform_fehtla_dat_new$Asymptote)
dim(perform_fehtla_dat_new_grid_AH)

# Replace the dispersion values with cero
d_new<-rep(0.0004,506)
perform_fehtla_dat_new_grid_AH$d<-d_new


perform_fehtla_dat_new_grid<-perform_fehtla_dat_new %>% 
  mutate(xgrid= perform_fehtla_dat_new$Asymptote);dim(perform_fehtla_dat_new_grid)

names(perform_fehtla_dat_new_grid_AH)
names(perform_fehtla_dat_new_grid)


# We combine the two data frames into one

combined<-merge(perform_fehtla_dat_new_grid, 
                perform_fehtla_dat_new_grid_AH, 
                by=c("xgrid","TL50"),
                suffixes =c("_K", "_AH"));names(combined)
dim(combined)
names(combined)


# now transform the data frame from wide to long
data_long_combined <- gather(combined,measure, 
                             performance,
                             Asymptote_K:rel.total.B.inacc_AH, 
                             factor_key=TRUE);names(data_long_combined)
head(data_long_combined)
names(data_long_combined)

# we split the row "measure" again to have the scenario and performance measures in
# two seperate columns

data_long_combined_sep<-separate(data = data_long_combined, 
                                 col = "measure", into = c("measure", "scenario"), 
                                 sep = "_")
names(data_long_combined_sep)
head(data_long_combined_sep)

# finally we split the data frame into a list of data frames according to their
# performance measures
dat_final<-split(data_long_combined_sep,
                 f = data_long_combined_sep$measure);names(dat_final)
head(dat_final$y.pred)

# We attach the columns for the isoline
Bpred_vec<-dat_final$Bpred.ch$performance
y_vec<-dat_final$y$performance
B_vec<-dat_final$B.ch$performance
d_vec<-dat_final$d$performance
Binacc_vec<-dat_final$rel.total.B.inacc$performance


for (i in names(dat_final)){
  dat_final[[i]]$Bpred_vec<-Bpred_vec
  dat_final[[i]][["B_vec"]]<-B_vec
  dat_final[[i]][["y_vec"]]<-y_vec
  dat_final[[i]][["d_vec"]]<-d_vec
  dat_final[[i]][["Binacc_vec"]]<-Binacc_vec
}

 ## FIG. S4 : Multiruns without accessibilitiy ----
tiff("Results/Figures/Integrated implementation/Supplementary/Multiruns without accessibility.tif", 
     res=300, compression = "lzw", height=220, width=169, units="mm", 
     family = "sans", pointsize = 9)

   ### Running the plot ----
# now we can run the plots so that the performance measures of each scenario
# share the same legend (facet_wrap)
for (i in seq_along(dat_final)){
  dat_final[[i]]$scenario<- as.factor(dat_final[[i]]$scenario)
  levels(dat_final[[i]]$scenario) <- c("BSH", "BH[P/B]")
}

p1<-ggplot(dat_final$d, aes(xgrid, TL50)) + 
  stat_contour(geom = 'polygon', aes(z=performance, fill = ..level..)) +
  geom_tile(aes(fill=performance))+
  viridis::scale_fill_viridis(direction = -1,discrete = F)+
  geom_contour(aes(z = performance), size=0.2) +
  facet_grid(~scenario)+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_blank())+
  theme(legend.title=element_blank(), strip.background = element_blank(),
        strip.text.x = element_blank())+
  theme(axis.title.x= element_blank(), 
        axis.title.y= element_text(family='sans', size = 12))+
  theme(axis.text = element_text(family='sans', size = 9))+
  ylab(expression(tau[50]))


p2<-ggplot(dat_final$Bpred.ch, aes(xgrid, TL50)) + 
  stat_contour(geom = 'polygon', aes(z=performance, fill = ..level..)) +
  geom_tile(aes(fill=performance))+
  viridis::scale_fill_viridis(direction = -1,discrete = F)+
  geom_contour(aes(z = performance), size=0.2) +
  # scale_color_continuous(low = "white", high = "white", guide=F) +
  facet_grid(~scenario)+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_blank())+
  theme(legend.title=element_blank(), strip.background = element_blank(),
        strip.text.x = element_blank())+
  theme(axis.title.x= element_blank(), 
        axis.title.y= element_text(family='sans', size = 14))+
  theme(axis.text = element_text(family='sans', size = 9))+
  ylab(expression(tau[50]))
  

p3 <- p2 %+% dat_final$y
p4 <- p2 %+% dat_final$y.pred


grid.arrange(bottom=textGrob("Asymptote (Fishing intensity)", gp =gpar(fontsize=11)),
             ncol=1, nrow=4,
             top=textGrob("", hjust=0, x=0),
             arrangeGrob(p1, left=textGrob("Disturbance index", gp=gpar(fontsize=11), 
                                               rot = 90)),
             arrangeGrob(p2, left=textGrob("Predator B change", gp=gpar(fontsize=11), 
                                               rot = 90)),
             arrangeGrob(p3, left=textGrob("Total catch", gp=gpar(fontsize=11), 
                                               rot = 90)),
             arrangeGrob(p4, left=textGrob("Predator catch", gp=gpar(fontsize=11), 
                                               rot = 90)))
grid.text("BSH", x = unit(0.28, "npc"), 
          y = unit(0.99, "npc"),gp = gpar(fontsize=12))
grid.text(expression("BH"[P/B]), x = unit(0.69, "npc"), 
          y = unit(0.99, "npc"),gp = gpar(fontsize=12))



#-----

dev.off()


# ---------------------- MANUSCRIPT QUANTITATIVE RESULTS -----------------------
# 1.) PART I ----
 ## A.) Dispersion values ----

bhm12_perform$`0.5`$d
bhm12_perform$`10`$d

bhm12_perform_P$`0.5`$d
bhm12_perform_P$`10`$d

 ## B.) Percentage of catch----

# Total catch 
# differences in catch between F~P and F~K
round(100-100*bhm12_perform_P$`1`$y/bhm12_perform$`1`$y, digits = 1)
round(100-100*bhm12_perform_P$`5`$y/bhm12_perform$`5`$y, digits = 1)
round(100-100*bhm12_perform_P$`10`$y/bhm12_perform$`10`$y, digits = 1)

# differences in catch between F~K and AH
round(100-100*alt_perform$`0.5`$y/bhm12_perform$`0.5`$y, digits = 1)
round(100-100*alt_perform$`5`$y/bhm12_perform$`5`$y, digits = 1)
round(100-100*alt_perform$`10`$y/bhm12_perform$`10`$y, digits = 1)

# Predatory catch 
round(((bhm12_perform$`0.5`$y.pred/bhm12_perform$`0.5`$y)*100), digits = 2)
round(((bhm12_perform$`10`$y.pred/bhm12_perform$`10`$y)*100), digits = 2)

round(((bhm12_perform_P$`1`$y.pred/bhm12_perform_P$`1`$y)*100), digits = 2)
round(((bhm12_perform_P$`10`$y.pred/bhm12_perform_P$`10`$y)*100), digits = 2)

round(((alt_perform$`0.5`$y.pred/alt_perform$`0.5`$y)*100), digits = 2)
round(((alt_perform$`10`$y.pred/alt_perform$`10`$y)*100), digits = 2)

# Forage fish catch
round(((bhm12_perform$`1`$y.forage/bhm12_perform$`1`$y)*100), digits = 2)
round(((bhm12_perform$`10`$y.forage/bhm12_perform$`10`$y)*100), digits = 2)

round(((bhm12_perform_P$`1`$y.forage/bhm12_perform_P$`1`$y)*100), digits = 2)
round(((bhm12_perform_P$`10`$y.forage/bhm12_perform_P$`10`$y)*100), digits = 2)

round(((alt_perform$`0.5`$y.forage/alt_perform$`0.5`$y)*100), digits = 2)
round(((alt_perform$`10`$y.forage/alt_perform$`10`$y)*100), digits = 2)



 ## C.) Difference in exploitation rate between low and high trophic level----
#    under F~P----
df_bhm12_P$`5`$F_loss_acc[2]
df_bhm12_P$`5`$F_loss_acc[32]
round(100*df_bhm12_P$`6`$E_acc[2], digits = 1)
round(100*df_bhm12_P$`6`$E_acc[32], digits = 1)


 ## D.) Sum of fishing mortalities----
sum(df_bhm12$`5`$Fish_mort_acc[-1])
sum(df_bhm12_P$`5`$Fish_mort_acc[-1])
sum(df_alt$`5`$Fish_mort_acc[-1])
 ## E.) Differences of fishing mortality of F~K and F~P for each trophic level----
round((df_bhm12_P$`1`$Fish_mort_acc/df_bhm12$`1`$Fish_mort_acc*100), digits=2)
round((df_bhm12_P$`5`$Fish_mort_acc/df_bhm12$`5`$Fish_mort_acc*100), digits=2)
round((df_bhm12_P$`10`$Fish_mort_acc/df_bhm12$`10`$Fish_mort_acc*100), digits=2)
# 2.) SENSITIVITY ANALYSIS NUMBERS----
 ## A.) Dispersion top-down----
topd_perform<-c()
for (i in seq_along(df_bhm12_topd_P)){
  topd_perform[[i]]<-lapply(df_bhm12_topd_P[[i]], 
                            function (x) perform(data=x, 
                                                 VE=df_bhm12_topd_P$`1`$`0`))
}
names(topd_perform)<-names(df_bhm12_topd_P)

topd_perform$`0`$`0.5`$d
topd_perform$`0`$`10`$d

topd_perform$`1`$`0.5`$d
topd_perform$`1`$`10`$d

 ## B.) Dispersion TE----
TE_perform<-c()
for (i in seq_along(df_bhm12_TE_P)){
  TE_perform[[i]]<-mapply(perform, 
                          data=df_bhm12_TE_P[[i]],VE=df_bhm12_TE_P[[i]][1], 
                          SIMPLIFY = F, USE.NAMES = T)
}
names(TE_perform)<-names(df_bhm12_TE_P)
TE_perform$`0.05`$`0.5`$d
TE_perform$`0.05`$`10`$d

TE_perform$`0.1`$`0.5`$d
TE_perform$`0.1`$`10`$d

TE_perform$`0.2`$`0.5`$d
TE_perform$`0.2`$`10`$d

bhm12_perform_P2$`0.5`$d
bhm12_perform_P2$`10`$d



# 3.) PART II ----
 ## A) Dispersion values----
bumss_perform$`0.5`$d
bumss_perform$`10`$d 

bumss_perform_P$`0.5`$d
bumss_perform_P$`10`$d

alt_acc_perform$`0.5`$d
alt_acc_perform$`10`$d

 ## B) Percentage of catch ----

# Catch compared to a fully accessible system
C_diff0.5<-bumss_perform$`0.5`$y/bhm12_perform$`0.5`$y;C_diff0.5
C_diff10<-bumss_perform$`10`$y/bhm12_perform$`10`$y;C_diff10

C_diff1_alt<-alt_acc_perform$`1`$y/alt_perform$`1`$y;C_diff1_alt
C_diff10_alt<-alt_acc_perform$`10`$y/alt_perform$`10`$y;C_diff10_alt
1-mean(c(C_diff0.5, C_diff10, C_diff1_alt, C_diff10_alt))


# Total catch 
# comparing F~P and F~K
round(100-100*bumss_perform_P$`1`$y/bumss_perform$`1`$y, digits = 1)
round(100-100*bumss_perform_P$`5`$y/bumss_perform$`5`$y, digits = 1)
round(100-100*bumss_perform_P$`10`$y/bumss_perform$`10`$y, digits = 1)

# Comparing AH and F~K
round(100-100*alt_acc_perform$`1`$y/bumss_perform$`1`$y, digits = 1)
round(100-100*alt_acc_perform$`5`$y/bumss_perform$`5`$y, digits = 1)
round(100-100*alt_acc_perform$`10`$y/bumss_perform$`10`$y, digits = 1)

# Predatory catch 
round(((bumss_perform$`0.5`$y.pred/bumss_perform$`0.5`$y)*100), digits = 2)
round(((bumss_perform$`10`$y.pred/bumss_perform$`10`$y)*100), digits = 2)

round(((alt_acc_perform$`0.5`$y.pred/alt_acc_perform$`0.5`$y)*100), digits = 2)
round(((alt_acc_perform$`10`$y.pred/alt_acc_perform$`10`$y)*100), digits = 2)


 ## C) Biomass changes per trophic level ----
df_bumss$`0.5`$BIOM_MF/df_bumss$`0`$BIOM_MF
df_bumss$`10`$BIOM_MF/df_bumss$`0`$BIOM_MF

df_alt_acc$`1`$`1`$BIOM_MF/df_alt_acc$`0`$`0`$BIOM_MF
df_alt_acc$`10`$`1`$BIOM_MF/df_alt_acc$`0`$`1`$BIOM_MF


 ## D) Sum of fishing mortalities ----
sum(df_bhm12$`5`$Fish_mort_acc[-1])
sum(df_bhm12_P$`5`$Fish_mort_acc[-1])
sum(df_alt$`5`$Fish_mort_acc[-1])



#-------------------------------------Remaining (to be erased)

# Belonging to manuscript figures part 1

 ## E) Amount of inaccessible biomass ----
round(bumss_perform$`0`$rel.total.B.inacc, digits = 2)
round(bumss_perform$`10`$rel.total.B.inacc, digits = 2)

round(alt_acc_perform$`0`$rel.total.B.inacc, digits = 2)
round(alt_acc_perform$`10`$rel.total.B.inacc, digits = 2)

# 4.) PART III ----

# Here I will split the data frame first into quantiles of total catch
blub_fehtla_acc<-perform_fehtla_acc_dat_new %>%
  mutate(quan_y = ntile(y, n = 4)) %>% 
  mutate(quan_ypred=ntile(y.pred, n=4)) %>% 
  mutate(quan_d=ntile(d,n=4)) %>% 
  mutate(quan_B=ntile(B.ch, n=4))

blub_fehtla_acc_AH<-perform_fehtla_acc_dat_new_AH %>%
  mutate(quan_y = ntile(y, n = 4)) %>% 
  mutate(quan_ypred=ntile(y.pred, n=4)) %>% 
  mutate(quan_d=ntile(d,n=4)) %>% 
  mutate(quan_B=ntile(B.ch, n=4))

dim(blub_fehtla_acc)
dim(blub_fehtla_acc_AH)


 ## A) Which scenarios is better for d, relative inaccessible biomass and total ----
 #     biomass? AH or BH? ----
blub_fehtla_acc_noncero<-subset(blub_fehtla_acc, Asymptote>0)
blub_fehtla_acc_noncero_AH<-subset(blub_fehtla_acc_AH, Asymptote>0)

   ### i.) Dispersion ----
# BH F~K
d_max<-max(blub_fehtla_acc_noncero$d);d_max
d_min<-min(blub_fehtla_acc_noncero$d);d_min
# AH
d_max_AH<-max(blub_fehtla_acc_AH$d);d_max_AH
d_min_AH<-min(blub_fehtla_acc_AH$d);d_min_AH

   ### ii.) Relative inaccessible biomass ----
# BH F~K
inaccB_max<-max(blub_fehtla_acc_noncero$rel.total.B.inacc);inaccB_max
inaccB_min<-min(blub_fehtla_acc_noncero$rel.total.B.inacc);inaccB_min
# AH
inaccB_max_AH<-max(blub_fehtla_acc_AH$rel.total.B.inacc);inaccB_max_AH
inaccB_min_AH<-min(blub_fehtla_acc_AH$rel.total.B.inacc);inaccB_min_AH

   ### iii.) Total biomass ----
# BH F~K
Btot_max<-max(blub_fehtla_acc_noncero$B.ch);Btot_max
Btot_min<-min(blub_fehtla_acc_noncero$B.ch);Btot_min
# AH
Btot_max_AH<-max(blub_fehtla_acc_AH$B.ch);Btot_max_AH
Btot_min_AH<-min(blub_fehtla_acc_AH$B.ch);Btot_min_AH

# AH has lower maximum and lower minimum values in dispersion 

 ## B) Predatory catch ----
blub_fehtla_acc_sub<-subset(blub_fehtla_acc, B.ch>=0.6 & Bpred.ch>=0.6)
dim(blub_fehtla_acc_sub)
C_maxpred_0.5_pos<-which.max(blub_fehtla_acc_sub$y.pred)
C_maxpred_0.5<-blub_fehtla_acc_sub[C_maxpred_0.5_pos,9]
max(blub_fehtla_acc$y.pred)
check<-subset(blub_fehtla_acc_sub, TL50==1)

# what is the total catch when maximizing predatory catch compared to the maximum total
# catch possible without depleting Biomass
C_max_0.5_pos<-which.max(blub_fehtla_acc_sub$y)
C_max_0.5<-blub_fehtla_acc_sub[C_max_0.5_pos,8]
maxCatmaxCpred<-blub_fehtla_acc_sub[C_maxpred_0.5_pos,8]
1-maxCatmaxCpred/C_max_0.5

# What is the ecological food print when maximizing total catch versus predator C
C_maxpred_0.5<-blub_fehtla_acc_sub[C_maxpred_0.5_pos,]
C_max_0.5<-blub_fehtla_acc_sub[C_max_0.5_pos,]

# Which scenario gives more predator catch the AH or BH?
blub_fehtla_acc_AH_sub<-subset(blub_fehtla_acc_AH, B.ch>=0.5 & Bpred.ch>=0.5)
dim(blub_fehtla_acc_AH_sub)
C_maxpred_0.5_pos_AH<-which.max(blub_fehtla_acc_AH_sub$y.pred)
C_maxpred_0.5_AH<-blub_fehtla_acc_AH_sub[C_maxpred_0.5_pos_AH,9]

1-C_maxpred_0.5_AH/C_maxpred_0.5

 ## C) Maximum total catch ----

# Maximum total catch while limiting the biomass reduction of predators and total
# system
blub_fehtla_acc_sub<-subset(blub_fehtla_acc, B.ch>=0.5 & Bpred.ch>=0.5)
C_max_0.5_pos<-which.max(blub_fehtla_acc_sub$y)
C_max_0.5<-blub_fehtla_acc_sub[C_max_0.5_pos,8]

blub_fehtla_acc_AH_sub<-subset(blub_fehtla_acc_AH, B.ch>=0.5 & Bpred.ch>=0.5)
C_max_0.5_pos_AH<-which.max(blub_fehtla_acc_AH_sub$y)
C_max_0.5_AH<-blub_fehtla_acc_AH_sub[C_max_0.5_pos_AH,8]

1-C_max_0.5_AH/C_max_0.5

 ## D) Balancing between total catch and predatory catch ----

# first I subset to only have fishing strategies that keep biomass above 50 % 

blub_fehtla_acc_sub<-subset(blub_fehtla_acc, B.ch>=0.6 & Bpred.ch>=0.6)
balance<-subset(blub_fehtla_acc_sub, quan_y==3 & quan_ypred==4);dim(balance)

cmax_balance<-max(blub_fehtla_acc_sub$y)
cmax_range<-(cmax_balance/100)*60
cpredmax_balance<-max(blub_fehtla_acc_sub$y.pred)
cpredmax_range<-(cpredmax_balance/100)*60

balance_sub<-subset(blub_fehtla_acc_sub, y>=cmax_range & y.pred>=cpredmax_range);dim(balance_sub)

balance_sub[,9]/max(blub_fehtla_acc_sub$y.pred)
balance_sub[,8]/max(blub_fehtla_acc_sub$y)

blub_fehtla_acc[246,]

# Now, I check what are the ecological performances compared to having the same
# total catch as with balanced harvest
# 247,248,259,250
BH<-subset(blub_fehtla_acc_sub, TL50==1)
blub_fehtla_acc[185,]
blub_fehtla_acc[185,9]/predC_max
blub_fehtla_acc[185,9]/balance_ypred

 ## E) Least impact given certain catch ----
   ### i.) Total catch ----
# Here I will split the data frame first into quantiles of total catch
blub_fehtla_acc<-perform_fehtla_acc_dat_new %>%
  mutate(quan_y = ntile(y, n = 4)) %>% 
  mutate(quan_ypred=ntile(y.pred, n=4)) %>% 
  mutate(quan_d=ntile(d,n=4)) %>% 
  mutate(quan_B=ntile(B.ch, n=4))

blub_fehtla_acc_AH<-perform_fehtla_acc_dat_new_AH %>%
  mutate(quan_y = ntile(y, n = 4)) %>% 
  mutate(quan_ypred=ntile(y.pred, n=4)) %>% 
  mutate(quan_d=ntile(d,n=4)) %>% 
  mutate(quan_B=ntile(B.ch, n=4))

dim(blub_fehtla_acc)
dim(blub_fehtla_acc_AH)

# first I subset to only have fishing strategies (>0) that keep biomass above 50 % 

blub_fehtla_acc_sub<-subset(blub_fehtla_acc, B.ch>=0.6 & Bpred.ch>=0.6 &
                              Asymptote>0)
blub_fehtla_acc_sub_AH<-subset(blub_fehtla_acc_AH, B.ch>=0.6 & Bpred.ch>=0.6 &
                                 Asymptote>0)

# Extracting the balanced harvest strategies
BH<-subset(blub_fehtla_acc_sub, TL50==1)
BH_minY<-min(BH$y)
# yield --> 10.4, 18.3, 24.2

best0<-subset(blub_fehtla_acc_sub, y>=9.9 & y<10.9);dim(best0)
head(best0[order(best0$d),],5) # the lowest dispersion --> TL50=1 (1.1,1.2,1.3,1.4)
head(best0[order(best0$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=1.9 (1.8,1.7,1.6, 1.5)
head(best0[order(best0$B.ch, decreasing = T),],5) # the highest TotB --> TL50=1.9 (1.8...)
head(best0[order(best0$rel.total.B.inacc),],5) # the least change in rel. inacc. 
                                                # B --> TL50=3.2 (2.8, 1.7,1.8, 1.9)
head(best0[order(best0$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B

best1<-subset(blub_fehtla_acc_sub, y>=17.8 & y<18.8);dim(best1)
head(best1[order(best1$d),],5) # the lowest dispersion --> TL50=1 (1.1,1.2,1.3, 1.4)
head(best1[order(best1$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=1.8 (1.7...)
head(best1[order(best1$B.ch, decreasing = T),],5) # the highest TotB --> TL50=1.8 (1.7....)
head(best1[order(best1$rel.total.B.inacc),],5)# the least change in relative inaccessible biomass --> TL50=2.5 (1.8...)
head(best1[order(best1$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B

best2<-subset(blub_fehtla_acc_sub, y>=23.7 & y<24.7);dim(best2)
head(best2[order(best2$d),],5) # the lowest dispersion --> TL50=1 (1.1,1.2,1.3,1.4)
head(best2[order(best2$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=1.8 (1.7...)
head(best2[order(best2$B.ch, decreasing = T),],5) # the highest TotB --> TL50=1.8 (1.7...)
head(best2[order(best2$rel.total.B.inacc),],5) # the least change in relative inaccessible biomass --> TL50=1.8 (1.7...)
head(best1[order(best1$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B


# AH
AH<-subset(blub_fehtla_acc_sub_AHH, TL50==1);AH$y
# yield --> 8, 13.7, 18.1, 21.4, 24, 26, 27.6, 28.8, 29.8, 30.5

best0_AH<-subset(blub_fehtla_acc_sub_AH, y>=7.5 & y<8.5);dim(best0_AH)
head(best0_AH[order(best0_AH$d),],5) # the lowest dispersion --> TL50=1 (1.1, 1.2, 1.3, 1.4)
head(best0_AH[order(best0_AH$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=1.8 (1.5,1.4,1.3,1.2)
head(best0_AH[order(best0_AH$B.ch, decreasing = T),],5) # the highest TotB --> TL50=1.8 (1.5,1.4,1.3,1.2)
head(best0_AH[order(best0_AH$rel.total.B.inacc),],5) 
# the least change in rel. inacc. B. --> TL50=3.2 (2.8,2.9,2.6,2.7)
head(best0_AH[order(best0_AH$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B

best1_AH<-subset(blub_fehtla_acc_sub_AH, y>=17.6 & y<18.6);dim(best1_AH)
head(best1_AH[order(best1_AH$d),],5) # the lowest dispersion --> TL50=1 (1.1, 1.2, 1.3, 1.4)
head(best1_AH[order(best1_AH$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=1.7 (1.4,1.3,1.2,1.1)
head(best1_AH[order(best1_AH$B.ch, decreasing = T),],5) # the highest TotB --> TL50=1.7 (1.4,1.3,1.2,1.1)
head(best1_AH[order(best1_AH$rel.total.B.inacc),],5) 
# the least change in rel. inacc. B. --> TL50=2.8 (2.6,2.4,2.2,1.5)
head(best1_AH[order(best1_AH$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B

best2_AH<-subset(blub_fehtla_acc_sub_AH, y>=23.5 & y<24.5);dim(best2_AH)
head(best2_AH[order(best2_AH$d),],5) # the lowest dispersion --> TL50=1 (1.1, 1.2, 1.3, 1.4)
head(best2_AH[order(best2_AH$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=1.8 (1.4,1.3,1.2,1.1)
head(best2_AH[order(best2_AH$B.ch, decreasing = T),],5) # the highest TotB --> TL50=2.6 (1.4,1.3,1.2,1.1)
head(best2_AH[order(best2_AH$rel.total.B.inacc),],5) 
# the least change in rel. inacc. B. --> TL50=2.6 (2.3,1.5,1.4,1.3)
head(best2_AH[order(best2_AH$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B

best4_AH<-subset(blub_fehtla_acc_sub_AH, y>=28.3 & y<29.3);dim(best4_AH)
head(best4_AH[order(best4_AH$d),],5) # the lowest dispersion --> TL50=1 (1.1, 1.2, 1.3, 1.4)
head(best4_AH[order(best4_AH$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=2.1 (1.5,1.4,1.3,1.2)
head(best4_AH[order(best4_AH$B.ch, decreasing = T),],5) # the highest TotB --> TL50=2.4 (1.5,1.4,1.3,1.2)
head(best4_AH[order(best4_AH$rel.total.B.inacc),],5) 
# the least change in rel. inacc. B. --> TL50=2.4 (1.6,1.5,1.4,1.3)
head(best4_AH[order(best4_AH$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B

best5_AH<-subset(blub_fehtla_acc_sub_AH, y>=30 & y<31);dim(best5_AH)
head(best5_AH[order(best5_AH$d),],5) # the lowest dispersion --> TL50=1 (1.1, 1.2, 1.3, 1.4)
head(best5_AH[order(best5_AH$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=2.2 (1.5,1.4,1.3,1.2)
head(best5_AH[order(best5_AH$B.ch, decreasing = T),],5) # the highest TotB --> TL50=2.2 (2.2,2.3,1.5,1.4)
head(best5_AH[order(best5_AH$rel.total.B.inacc),],5) 
# the least change in rel. inacc. B. --> TL50=2.2 (2.2,1.6,1.5,1.4)
head(best5_AH[order(best5_AH$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B

best6_AH<-subset(blub_fehtla_acc_sub_AH, y>=51.5 & y<52.5);dim(best6_AH)
head(best6_AH[order(best6_AH$d),],5) # the lowest dispersion --> TL50=1 (1.1, 1.2, 1.3, 1.4)
head(best6_AH[order(best6_AH$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=1.6 (1.5,1.4,1.3,1.2)
head(best6_AH[order(best6_AH$B.ch, decreasing = T),],5) # the highest TotB --> TL50=2.2 (1.6,2,1.5,1.4)
head(best6_AH[order(best6_AH$rel.total.B.inacc),],5) # the least change in rel. inacc. B. --> TL50=2.2 (2,1.6,1.5,1.4)
head(best6_AH[order(best6_AH$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B

best7_AH<-subset(blub_fehtla_acc_sub_AH, y>=53.8 & y<54.8);dim(best7_AH)
head(best7_AH[order(best7_AH$d),],5) # the lowest dispersion --> TL50=1 (1.1, 1.2, 1.3, 1.4)
head(best7_AH[order(best7_AH$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=1.7 (1.7,1.5,1.4,1.3)
head(best7_AH[order(best7_AH$B.ch, decreasing = T),],5) # the highest TotB --> TL50=2.2 (1.7,1.6,1.5,1.4)
head(best7_AH[order(best7_AH$rel.total.B.inacc),],5) # the least change in rel. inacc. B. --> TL50=2.2 (1.7,1.6,1.5,1.4)
head(best7_AH[order(best7_AH$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B

best8_AH<-subset(blub_fehtla_acc_sub_AH, y>=55.7 & y<56.7);dim(best8_AH)
head(best8_AH[order(best8_AH$d),],5) # the lowest dispersion --> TL50=1 (1.1, 1.2, 1.3, 1.4)
head(best8_AH[order(best8_AH$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=1.8 (1.7,1.6,1.5,1.4)
head(best8_AH[order(best8_AH$B.ch, decreasing = T),],5) # the highest TotB --> TL50=2.1 (1.8,1.7,1.6,1.5)
head(best8_AH[order(best8_AH$rel.total.B.inacc),],5) # the least change in rel. inacc. B. --> TL50=2.1 (1.8,1.7,1.6,1.5)
head(best8_AH[order(best8_AH$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B

best9_AH<-subset(blub_fehtla_acc_sub_AH, y>=57.2 & y<58.2);dim(best9_AH)
head(best9_AH[order(best9_AH$d),],5) # the lowest dispersion --> TL50=1 (1.1, 1.2, 1.3, 1.4)
head(best9_AH[order(best9_AH$Bpred.ch, decreasing = T),],5) # the highest PredB --> TL50=1.9 (1.8,1.7,1.6,1.5)
head(best9_AH[order(best9_AH$B.ch, decreasing = T),],5) # the highest TotB --> TL50=1.9 (1.8,1.7,1.6,1.5)
head(best9_AH[order(best9_AH$rel.total.B.inacc),],5) # the least change in rel. inacc. B. --> TL50=1.9 (1.8,1.7,1.6,1.5)
head(best9_AH[order(best9_AH$rel.total.B.inacc, decreasing = T),],5) # BH has the highest impact on rel.inacc.B


# I will now broaden the range to see which is the best strategy when maximizing the
# total catch to a range of 90-100% of maximum total catch

maxC<-max(blub_fehtla_acc_sub[,8])
C0.9<-(maxC/100)*90

best0.9<-subset(blub_fehtla_acc_sub, y>=23 & y<25.6);dim(best0.9)
# the lowest dispersion --> TL50=2.3
# the highest PredB --> TL50= 2.3
# the highest TotB --> TL50= 2.3
# the least change in relative inaccessible biomass --> TL50=2.3

maxC_AH<-max(blub_fehtla_acc_sub_AH[,8])
C0.9_AH<-(maxC_AH/100)*90

C0.9_AH<-subset(blub_fehtla_acc_sub_AH, y>=23 & y<25.6);dim(C0.9_AH)


   ### ii.) Which scenario (AH or BH) outperforms when maximizing catch? ----

# Now I compare AH and F~P/B at a given catch

# Which scenario has the highest catch ? 

# first I subset to only have fishing strategies (>0) that keep biomass above 50 % 

blub_fehtla_acc_sub<-subset(blub_fehtla_acc, B.ch>=0.6 & Bpred.ch>=0.6 &
                              Asymptote>0)
blub_fehtla_acc_sub_AH<-subset(blub_fehtla_acc_AH, B.ch>=0.6 & Bpred.ch>=0.6 &
                                 Asymptote>0)

max(blub_fehtla_acc_sub$y)
max(blub_fehtla_acc_sub_AH$y)

AHTL1<-subset(blub_fehtla_acc_sub_AH, TL50==1)
min(AHTL1$y)
max(blub_fehtla_acc_sub$y.pred)
max(blub_fehtla_acc_sub_AH$y.pred)

# Find the same catches (digit of 1) between AH and F~P/B
Lookup = round(blub_fehtla_acc_sub$y, digits = 1)      # Create a vector of catches in F~P/B
sub<-subset(blub_fehtla_acc_sub, y>=7.96);dim(sub)    # Subset for minimum catches of AH for TL50=1
sub_AH<-subset(blub_fehtla_acc_sub_AH, y>=7.96);dim(sub_AH)

# Look-up the catches from F~P/B in AH
Find_AH<-sub_AH[Lookup %in% round(sub_AH$y, digits = 1),];dim(Find_AH)

# The final values of catches that are shard in both strategies
Lookup2 <- round(sub_AH$y, digits = 1)
Lookup2<-unique(Lookup2); length(Lookup2)   # remove duplicates
Lookup2[order(Lookup2)]                     # order them

# Now I can use these values to compare the performance under same catches

comp<-c()
for (i in seq_along(Lookup2)){
  comp[[i]]<-subset(blub_fehtla_acc_sub, y>=Lookup2[i]-0.5 & y<Lookup2[i]+0.5)
}
names(comp)<-Lookup2
comp_sub<-Filter(function(x) dim(x)[1] > 0, comp)
names(comp_sub)

comp_AH<-c()
for (i in seq_along(Lookup2)){
  comp_AH[[i]]<-subset(blub_fehtla_acc_sub_AH, y>=Lookup2[i]-0.5 & y<Lookup2[i]+0.5)
}
names(comp_AH)<-Lookup2
names(comp_sub)
comp_AH$`13.7`
comp_sub$`13.7`

# Comparing the dispersion, relative biomass and inaccessible biomass
min<-list("")
min_AH<-list("")
for (i in names(comp_sub)){
  for (j in names(comp_sub$`13.7`)){
    min[[i]][j]<-min(comp_sub[[i]][j])
    min_AH[[i]][[j]]<-min(comp_AH[[i]][j])
  }
}

min_ord<-min[order(names(min))]
min_AH_ord<-min_AH[order(names(min_AH))]

names(min_ord)
names(min_AH_ord)
min_ord$`10.1`
min_AH_ord$`10.2`
# we subset and only compare 10 catch values

min_sub<-min_ord[c(2,9,16,24,31,38,45,52,59,66,72)]
names(min_sub)
length(min_sub)

min_AH_sub<-min_AH_ord[c(2,9,16,24,31,38,45,52,59,66,72)]
names(min_AH_sub)
length(min_AH_sub)

# 15
min_sub$`8.3`
min_AH_sub$`8.3`
# dispersion --> AH
# pred biomass --> AH
# tot biomass --> AH
# relative inacc B --> BH

# 17.5
min_sub$`13.4`
min_AH_sub$`13.4`
# dispersion --> AH
# pred biomass --> AH
# tot biomass --> AH
# relative inacc B --> BH

# 21.4
min_sub$`16.8`
min_AH_sub$`16.8`
# dispersion --> AH
# pred biomass --> AH
# tot biomass --> AH
# relative inacc B --> BH

# `24.2`
min_sub$`19.1`
min_AH_sub$`19.1`
# dispersion --> AH
# pred biomass --> AH
# tot biomass --> AH
# relative inacc B --> BH

# `26.5`
min_sub$`23.9`
min_AH_sub$`23.9`
# dispersion --> AH
# pred biomass --> AH
# tot biomass --> AH
# relative inacc B --> BH

# `30.4`
min_sub$`30.4`
min_AH_sub$`30.4`
# dispersion --> AH
# pred biomass --> BH
# tot biomass --> AH
# relative inacc B --> BH

# `33.9`
min_sub$`33.9`
min_AH_sub$`33.9`
# dispersion --> AH
# pred biomass --> AH
# tot biomass --> AH
# relative inacc B --> BH

# `37.5`
min_sub$`37.5`
min_AH_sub$`37.5`
# dispersion --> AH
# pred biomass --> AH
# tot biomass --> AH
# relative inacc B --> BH

# `41.8`
min_sub$`41.8`
min_AH_sub$`41.8`
# dispersion --> AH
# pred biomass --> AH
# tot biomass --> AH
# relative inacc B --> BH
  

# `44.8`
min_sub$`44.8`
min_AH_sub$`44.8`
# dispersion --> AH
# pred biomass --> AH
# tot biomass --> AH
# relative inacc B --> BH

# `47.2`
min_sub$`47.2`
min_AH_sub$`47.2`
# dispersion --> AH
# pred biomass --> AH
# tot biomass --> AH
# relative inacc B --> BH


   ### iii.) Predator catch ----
# Here I will split the data frame first into quantiles of total catch
blub_fehtla_acc<-perform_fehtla_acc_dat_new %>%
  mutate(quan_y = ntile(y, n = 4)) %>% 
  mutate(quan_ypred=ntile(y.pred, n=4)) %>% 
  mutate(quan_d=ntile(d,n=4)) %>% 
  mutate(quan_B=ntile(B.ch, n=4))

dim(blub_fehtla_acc)

# first I subset to only have fishing strategies (>0) that keep biomass above 50 % 

blub_fehtla_acc_sub<-subset(blub_fehtla_acc, B.ch>=0.6 & Bpred.ch>=0.6 &
                              Asymptote>0)


# For a given predator catch, does BH perform better or worse? 
# Extracting the balanced harvest strategies
check<-subset(blub_fehtla_acc_sub, TL50==1)

# Balanced harvest predator catch: 1.2; 1.8; 2.2

best_pred0<-subset(blub_fehtla_acc_sub, y.pred>=1.15 & y.pred<1.25);dim(best_pred0)
head(best_pred0[order(best_pred0$d, decreasing = T),],20) 
head(best_pred0[order(best_pred0$Bpred.ch, decreasing = F),],5) # highest impact BH
head(best_pred0[order(best_pred0$B.ch, decreasing = F),],5) # highest impact BH
head(best_pred0[order(best_pred0$rel.total.B.inacc, decreasing = T),],5) # highest impact BH

best_pred1<-subset(blub_fehtla_acc_sub, y.pred>=1.75 & y.pred<1.85);dim(best_pred1);best_pred1
head(best_pred1[order(best_pred1$d, decreasing = T),],4) # BH among the highest impact
head(best_pred1[order(best_pred1$Bpred.ch, decreasing = F),],5) # highest impact BH
head(best_pred1[order(best_pred1$B.ch, decreasing = F),],5) # highest impact BH
head(best_pred1[order(best_pred1$rel.total.B.inacc, decreasing = T),],5) # highest impact BH

best_pred2<-subset(blub_fehtla_acc_sub, y.pred>=2.15 & y.pred<2.25);dim(best_pred2);best_pred2
head(best_pred2[order(best_pred2$d, decreasing = T),],4) # BH among the highest impact
head(best_pred2[order(best_pred2$Bpred.ch, decreasing = F),],5) # highest impact BH
head(best_pred2[order(best_pred2$B.ch, decreasing = F),],5) # highest impact BH
head(best_pred2[order(best_pred2$rel.total.B.inacc, decreasing = T),],5) # highest impact BH

 ## F) Dispersion and catch values when TL50>4 ----

TL50abov4<-subset(blub_fehtla_acc, TL50>=4)
TL50below2<-subset(blub_fehtla_acc, TL50<=2 & Asymptote>0)
max(TL50abov4$d)
min(TL50abov4$d)
dim(TL50abov4)

max(TL50below2$d)
min(TL50below2$d)

TL50abov4_AH<-subset(blub_fehtla_acc_AH, TL50>=4)
max(TL50abov4_AH$d)
min(TL50abov4_AH$d)

# lowest dispersion value for balanced harvest
bh<-subset(blub_fehtla_acc, TL50==1 & Asymptote>0)
bh_d_min<-min(bh$d)

# Which scenarios have a dispersion value below those of balanced harvest minimum?

belowd_bh_d<-subset(blub_fehtla_acc, Asymptote>0 & d<bh_d_min & B.ch>0.6 &
                      Bpred.ch>0.6);dim(belowd_bh_d)
belowd_bh_d_ordered<-belowd_bh_d[order(belowd_bh_d$rel.total.B.inacc),]
belowd_bh_d_ordered_d<-belowd_bh_d[order(belowd_bh_d$d),]
blub_fehtla_acc$Bpred.ch
belowd0.08 %>% 
  filter(Asymptote==0.1)
# TL50 = 1, d = 0.0454, y = 16.8
subset(TL50abov4, d<=0.0454 & Asymptote>0 & y>2)

# Which scenarios have the best relative inaccessible biomass and the best total
# biomass?

besteco<-subset(blub_fehtla_acc, Asymptote>0 & B.ch>0.6 &
                      Bpred.ch>0.6);dim(besteco)
head(besteco[order(besteco$d),],40)      # all are above 4
head(besteco[order(besteco$rel.total.B.inacc,decreasing = F),],40)# all are above 4
head(besteco[order(besteco$B.ch,decreasing = T),],40)# all are above 4
head(besteco[order(besteco$Bpred.ch,decreasing = T),],40)# all are above 4
besteco[which.max(besteco$Bpred.ch),]

# Catch percentage compared with the maximum catch when setting a biomass reduction limit

# I subset to only have fishing strategies that keep biomass above 50 % 

blub_fehtla_acc_sub<-subset(blub_fehtla_acc, B.ch>=0.6 & Bpred.ch>=0.6)
blub_fehtla_acc_AH_sub<-subset(blub_fehtla_acc_AH, B.ch>=0.6 & Bpred.ch>=0.6)

# Maximum catch of a scenarios that does not reduce the biomass below 50 % 
max(blub_fehtla_acc_sub$y)

# Calculating the reduction of catch relative to maximum catch, while preserving 50 % of
# biomass

round((1-max(TL50abov4$y)/max(blub_fehtla_acc_sub$y))*100, digits = 2)
round((1-max(TL50abov4_AH$y)/max(blub_fehtla_acc_AH_sub$y))*100, digits = 2)


max(TL50abov4_AH$y)/max(blub_fehtla_acc_AH_sub$y)
min(TL50abov4_AH$y)






# ----------------- SUPPLEMENTARY ANALYSIS QUANTITATIVE RESULTS ----------------
# SUPPLEMENTARY PART I ----
 # 1.) Effect of topD and FormD on the dispersion, relative lower and predatory catch ----

# Because I run the performance equation on one mE at a time, I do a quick solution
# in renaming the data frames here, while running the code each time with a 
# different mE

# mE=0.5
perform_P_0.5<-perform_topd_alpha_dat_new_P

# mE=1
perform_P_1<-perform_topd_alpha_dat_new_P
perform_1<-perform_topd_alpha_dat_new
# mE=3
perform_P_3<-perform_topd_alpha_dat_new_P
perform_3<-perform_topd_alpha_dat_new

# mE=5
perform_P_5<-perform_topd_alpha_dat_new_P

# mE=7
perform_P_7<-perform_topd_alpha_dat_new_P

# mE=10
perform_P_10<-perform_topd_alpha_dat_new_P

  ## A.) Is the variation in dispersion higher between mE's or topd?----

# I want to know how much the change in the dispersion rate depends on the top-down
# control assumption.

# What if we don't have any top-down control? 
P_0.5<-subset(perform_P_0.5, TopD==0);P_0.5$d
P_1<-subset(perform_P_1, TopD==0);P_1$d
P_3<-subset(perform_P_3, TopD==0);P_3$d
P_5<-subset(perform_P_5, TopD==0);P_5$d
P_7<-subset(perform_P_7, TopD==0);P_7$d
P_10<-subset(perform_P_10, TopD==0);P_10$d
# dispersion:  0.0015; 0.0031; 0.0108; 0.0205; 0.0313; 0.0487

# What if I have the usual top-down control? 
bhm12_perform_P$`0.5`$d
bhm12_perform_P$`1`$d
bhm12_perform_P$`3`$d
bhm12_perform_P$`5`$d
bhm12_perform_P$`7`$d
bhm12_perform_P$`10`$d
# dispersion:  0.0016; 0.0033; 0.012; 0.0233; 0.0365; 0.0583

# What if I have the highest top-down control? 
P_0.5_high<-subset(perform_P_0.5, TopD==0.9 & FormD==0.9);P_0.5_high$d
P_1_high<-subset(perform_P_1, TopD==0.9 & FormD==0.9);P_1_high$d
P_3_high<-subset(perform_P_3, TopD==0.9 & FormD==0.9);P_3_high$d
P_5_high<-subset(perform_P_5, TopD==0.9 & FormD==0.9);P_5_high$d
P_7_high<-subset(perform_P_7, TopD==0.9 & FormD==0.9);P_7_high$d
P_10_high<-subset(perform_P_10, TopD==0.9 & FormD==0.9);P_10_high$d
# dispersion: 0.002; 0.0032; 0.0143; 0.0304; 0.0522; 0.0967

# mE=0.5   -->   0.0015-0.0016-0.002     0.0005
# mE=1     -->   0.0031-0.0033-0.0032    0.0001
# mE=3     -->   0.0108-0.012-0.0143     0.0035
# mE=5     -->   0.0205-0.0233-0.0304    0.0099
# mE=7     -->   0.0313-0.0365-0.0522    0.0209
# mE=10    -->   0.0487-0.0583-0.0967    0.048
# Range from mE=0.5 to mE=10             0.0472;0.0567;0.0947
# Range between mE's is larger than between top-down control

  ## B.) Relative predator and forage fish catch----
# I use an mE of 3 for comparison. 

# What is the range of relative predator catch? 
max(perform_P_3$rel.ypred)     # 0.0053
min(perform_P_3$rel.ypred)     # 0.0052

max(perform_3$rel.ypred)     # 0.022
min(perform_3$rel.ypred)     # 0.021

# So the relative predator catch is remaining low, regardless of the top-down


# F~P
predCmax_P<-which.max(perform_P_3$y.pred); round(max(perform_P_3$y.pred), 
                                                    digits = 4)
perform_P_3[predCmax_P,]
# the highest predator catch is with no top-down control
predCmin_P<-which.min(perform_P_3$y.pred);round(min(perform_P_3$y.pred), 
                                                 digits = 4)
perform_P_3[predCmin_P,]
# the lowest predator catch is with no top-down control

# F~P/B
predCmax<-which.max(perform_3$y.pred); round(max(perform_3$y.pred), 
                                                  digits = 4)
perform_3[predCmax,]
# the highest predator catch is with the highest top-down control
predCmin<-which.min(perform_3$y.pred);round(min(perform_3$y.pred), 
                                                 digits = 4)
perform_3[predCmin,]
# the lowest predator catch is with no top-down control. 

# Forage fish relative catch

# What is the range of relative predator catch? 
max(perform_P_3$rel.y.forage)     # 0.866
min(perform_P_3$rel.y.forage)     # 0.865

max(perform_3$rel.y.forage)     # 0.79
min(perform_3$rel.y.forage)     # 0.78

# So the relative forage fish catch is remaining high, regardless of the top-down

# F~P
predCmax_P<-which.max(perform_P_3$y.forage); round(max(perform_P_3$y.forage), 
                                                 digits = 4)
perform_P_3[predCmax_P,]
# the highest forage fish catch is with no top-down control
predCmin_P<-which.min(perform_P_3$y.forage);round(min(perform_P_3$y.forage), 
                                                digits = 4)
perform_P_3[predCmin_P,]
# the lowest forage fish catch is with no top-down control

# F~P/B
predCmax<-which.max(perform_3$y.forage); round(max(perform_3$y.forage), 
                                             digits = 4)
perform_3[predCmax,]
# the highest forage fish catch is with the highest top-down control
predCmin<-which.min(perform_3$y.forage);round(min(perform_3$y.forage), 
                                            digits = 4)
perform_3[predCmin,]
# the lowest forage fish catch is with no top-down control. 


 # 2.) Effect of different TE values and shapes on the dispersion, relative lower and ----
#     predatory catch----
  ## A.) Dispersion ----
d_TE<-list()
d_TE_P<-list()
for (i in names(TE_perform)){
  for (j in names(TE_perform$`0.05`)){
    d_TE[[i]][j]<-TE_perform[[i]][[j]]$d
    d_TE_P[[i]][j]<-TE_perform_P[[i]][[j]]$d
  }
}
# Range between Effort multipliers
d_TE$`0.05`[21]-d_TE$`0.05`[2] # range is 0.2552
d_TE$`0.1`[21]-d_TE$`0.1`[2]   # range is 0.3228
d_TE$`0.15`[21]-d_TE$`0.15`[2] # range is 0.3937
d_TE$`0.2`[21]-d_TE$`0.2`[2]   # range is 0.4532

d_TE_P$`0.05`[21]-d_TE_P$`0.05`[2] # range is 0.0253
d_TE_P$`0.1`[21]-d_TE_P$`0.1`[2]   # range is 0.0567
d_TE_P$`0.15`[21]-d_TE_P$`0.15`[2] # range is 0.1447
d_TE_P$`0.2`[21]-d_TE_P$`0.2`[2]   # range is 0.2543

# Range between transfer efficiencies
d_TE$`0.2`[2]-d_TE$`0.05`[2] # range is 0.0143
d_TE$`0.2`[11]-d_TE$`0.05`[11]   # range is 0.1374
d_TE$`0.2`[21]-d_TE$`0.05`[21] # range is 0.228

d_TE_P$`0.2`[2]-d_TE_P$`0.05`[2] # range is 0.0091
d_TE_P$`0.2`[11]-d_TE_P$`0.05`[11]   # range is 0.1105
d_TE_P$`0.2`[21]-d_TE_P$`0.05`[21] # range is 0.2381

# the higher the TE the higher the impact on the ecosystem
max(d_TE$`0.05`[[11]])
max(d_TE$`0.1`[[11]])
max(d_TE$`0.15`[[11]])
max(d_TE$`0.2`[[11]])

max(d_TE_P$`0.05`[[11]])
max(d_TE_P$`0.1`[[11]])
max(d_TE_P$`0.15`[[11]])
max(d_TE_P$`0.2`[[11]])

  ## B.) Relative predator & forage fish catch ----
rel.predy.percent<-c()
rel.forage.percent<-c()
rel.predy.percent_P<-c()
rel.forage.percent_P<-c()
for (i in names(rel.predy)){
  rel.predy.percent[[i]]<-100*rel.predy[[i]]
  rel.forage.percent[[i]]<-100*rel.forage[[i]]
  rel.predy.percent_P[[i]]<-100*rel.predy_P[[i]]
  rel.forage.percent_P[[i]]<-100*rel.forage_P[[i]]
}
rel.predy.percent$`0.05`
rel.predy.percent$`0.1`
rel.predy.percent$`0.15`
rel.predy.percent$`0.2`

rel.forage.percent$`0.05`
rel.forage.percent$`0.1`
rel.forage.percent$`0.15`
rel.forage.percent$`0.2`


# SUPPLEMENTARY PART II ----
# 1.) Relative inaccessible biomass given different midpoints----
 ## A) Transforming the data for plotting
   ### i.) Relative inaccessible biomass ----

# Because the relative inaccessible biomass depends on the assumption of the 
# midpoint, I am calculating the increase relative to the baseline (Virgin biomass)
rel.total.B.inacc<-lapply(midpoints_performance, 
                          function(a) lapply(a, function (a) a[["rel.total.B.inacc"]]))

rel.total.B.inacc_alt<-lapply(midpoints_performance_alt, 
                              function(a) lapply(a, function (a) a[["rel.total.B.inacc"]]))

# rel.total.B.inacc_virgin<-lapply(midpoints_performance, 
#                           function(a) a[[1]][["rel.total.B.inacc"]])

# B.inacc_inc<-list()
# for (i in names(rel.total.B.inacc)){
#   for (j in names(rel.total.B.inacc$`1`)){
# B.inacc_inc[[i]][[j]]<-rel.total.B.inacc[[i]][[j]]/rel.total.B.inacc_virgin[[i]]
#   }}

   ### ii.) Relative Biomass----
# I select one effort multiplier (mE=5)
mE_interest<-which(names(rel_B_midpoints$`1`)==5)
rel_B_midpoints_mE_interest<-lapply(rel_B_midpoints, function (a) a[[mE_interest]])
rel_B_midpoints_mE_interest_alt<-lapply(rel_B_midpoints_alt, function (a) a[[mE_interest]])

# I calculate the baseline for an mE of 0 
Ve_rel_B<-rel_B_midpoints$`1`$`0`
Ve_rel_B_alt<-rel_B_midpoints$`1`$`0`
   ### iii.) Total catch ----
# Extract the total catch for each midpoint
tot_catch_midpoints<-lapply(midpoints_performance, 
                            function(a) sapply(a, function (a) a[["y"]]))
midpoints_performance$`1`$`0`$y

tot_catch_midpoints_alt<-lapply(midpoints_performance_alt, 
                                function(a) sapply(a, function (a) a[["y"]]))
   ### iv.) Predator catch ----
# Extract the total catch for each midpoint
predcatch_midpoints<-lapply(midpoints_performance, 
                            function(a) sapply(a, function (a) a[["y.pred"]]))
midpoints_performance$`1`$`0`$y.pred

predcatch_midpoints_alt<-lapply(midpoints_performance_alt, 
                                function(a) sapply(a, function (a) a[["y.pred"]]))


   ### v.) Increase of dispersion with midpoint ----
midpoints_performance$`1`$`1`$d
midpoints_performance$`1.5`$`1`$d
midpoints_performance$`2`$`1`$d
midpoints_performance$`2.5`$`1`$d
midpoints_performance$`3`$`1`$d
midpoints_performance$`3.5`$`1`$d
midpoints_performance$`4`$`1`$d
midpoints_performance$`4.5`$`1`$d
midpoints_performance$`5`$`1`$d

midpoints_performance$`1`$`5`$d
midpoints_performance$`1.5`$`5`$d
midpoints_performance$`2`$`5`$d
midpoints_performance$`2.5`$`5`$d
midpoints_performance$`3`$`5`$d
midpoints_performance$`3.5`$`5`$d
midpoints_performance$`4`$`5`$d
midpoints_performance$`4.5`$`5`$d
midpoints_performance$`5`$`5`$d

midpoints_performance$`1`$`10`$d
midpoints_performance$`1.5`$`10`$d
midpoints_performance$`2`$`10`$d
midpoints_performance$`2.5`$`10`$d
midpoints_performance$`3`$`10`$d
midpoints_performance$`3.5`$`10`$d
midpoints_performance$`4`$`10`$d
midpoints_performance$`4.5`$`10`$d
midpoints_performance$`5`$`10`$d

midpoints_performance_alt$`1`$`1`$d
midpoints_performance_alt$`1.5`$`1`$d
midpoints_performance_alt$`2`$`1`$d
midpoints_performance_alt$`2.5`$`1`$d
midpoints_performance_alt$`3`$`1`$d
midpoints_performance_alt$`3.5`$`1`$d
midpoints_performance_alt$`4`$`1`$d
midpoints_performance_alt$`4.5`$`1`$d
midpoints_performance_alt$`5`$`1`$d

midpoints_performance_alt$`1`$`5`$d
midpoints_performance_alt$`1.5`$`5`$d
midpoints_performance_alt$`2`$`5`$d
midpoints_performance_alt$`2.5`$`5`$d
midpoints_performance_alt$`3`$`5`$d
midpoints_performance_alt$`3.5`$`5`$d
midpoints_performance_alt$`4`$`5`$d
midpoints_performance_alt$`4.5`$`5`$d
midpoints_performance_alt$`5`$`5`$d

midpoints_performance_alt$`1`$`10`$d
midpoints_performance_alt$`1.5`$`10`$d
midpoints_performance_alt$`2`$`10`$d
midpoints_performance_alt$`2.5`$`10`$d
midpoints_performance_alt$`3`$`10`$d
midpoints_performance_alt$`3.5`$`10`$d
midpoints_performance_alt$`4`$`10`$d
midpoints_performance_alt$`4.5`$`10`$d
midpoints_performance_alt$`5`$`10`$d

   ### vi.) Increase in inaccessible biomass ----
midpoints_performance$`1`$`10`$rel.total.B.inacc/midpoints_performance$`1`$`0`$rel.total.B.inacc
midpoints_performance$`1.5`$`10`$rel.total.B.inacc/midpoints_performance$`1.5`$`0`$rel.total.B.inacc
midpoints_performance$`2`$`10`$rel.total.B.inacc/midpoints_performance$`2`$`0`$rel.total.B.inacc
midpoints_performance$`2.5`$`10`$rel.total.B.inacc/midpoints_performance$`2.5`$`0`$rel.total.B.inacc
midpoints_performance$`3`$`10`$rel.total.B.inacc/midpoints_performance$`3`$`0`$rel.total.B.inacc
midpoints_performance$`3.5`$`10`$rel.total.B.inacc/midpoints_performance$`3.5`$`0`$rel.total.B.inacc
midpoints_performance$`4`$`10`$rel.total.B.inacc/midpoints_performance$`4`$`0`$rel.total.B.inacc
midpoints_performance$`4.5`$`10`$rel.total.B.inacc/midpoints_performance$`4.5`$`0`$rel.total.B.inacc
midpoints_performance$`5`$`10`$rel.total.B.inacc/midpoints_performance$`5`$`0`$rel.total.B.inacc

midpoints_performance_alt$`1`$`10`$rel.total.B.inacc/midpoints_performance_alt$`1`$`0`$rel.total.B.inacc
midpoints_performance_alt$`1.5`$`10`$rel.total.B.inacc/midpoints_performance_alt$`1.5`$`0`$rel.total.B.inacc
midpoints_performance_alt$`2`$`10`$rel.total.B.inacc/midpoints_performance_alt$`2`$`0`$rel.total.B.inacc
midpoints_performance_alt$`2.5`$`10`$rel.total.B.inacc/midpoints_performance_alt$`2.5`$`0`$rel.total.B.inacc
midpoints_performance_alt$`3`$`10`$rel.total.B.inacc/midpoints_performance_alt$`3`$`0`$rel.total.B.inacc
midpoints_performance_alt$`3.5`$`10`$rel.total.B.inacc/midpoints_performance_alt$`3.5`$`0`$rel.total.B.inacc
midpoints_performance_alt$`4`$`10`$rel.total.B.inacc/midpoints_performance_alt$`4`$`0`$rel.total.B.inacc
midpoints_performance_alt$`4.5`$`10`$rel.total.B.inacc/midpoints_performance_alt$`4.5`$`0`$rel.total.B.inacc
midpoints_performance_alt$`5`$`10`$rel.total.B.inacc/midpoints_performance_alt$`5`$`0`$rel.total.B.inacc


# SUPPLEMENTARY PART III ----
