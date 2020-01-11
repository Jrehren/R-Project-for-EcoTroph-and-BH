#1.) Emperical model for TL kinetics (Gascuel et al. 2008)
kin <- function(tl,sst=15){
  return(20.19*(tl^(-3.26))*exp(.041*sst))}

#2.) Production equation
geomSeq <- function(start,ratio,len_seq){
  seq<-seq(0,len_seq-1)
  return(start*ratio**(seq))
}

#3.) Biomass equation
biom<-function(prod, kin){
  prod/kin()
}

#4.) Selectivity function of accessible biomass - s-shape
s_selec<-function(x,asymptote,TL50,slope){
  return(asymptote/(1+exp(-slope*(x-TL50+0.05))))
}

#5.) Selectivity function of accessible biomass - bell-shape
bell_selec<-function(x, c, a, b){
  return(1/(1 + (abs((x - c)/a))^(2 * b)))
}        


#2.) Emperical model for TL kinetics accessible (Gascuel et al. 2008)

# For convention, Didier proposed to use the kinetic equation he fit to 
# only fish species, because he says this is most of the part that is accessible,
# in fact the kinetics are similiar for higher trophic levels, but lower for
# lower trophic levels. This fits the assumption that at lower trophic levels
# the accessible part of the biomass is lower. 

kin_acc <- function(tl,sst=15){
  return(2.31*(tl^(-1.72))*exp(0.053*sst))}
