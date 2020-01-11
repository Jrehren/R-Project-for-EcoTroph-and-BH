create.ETdiagnosis.adapted <- function(data, Mul_eff = NULL, fleet.of.interest = NULL, 
                               same.mE = NULL, B.Input=NULL, Beta = NULL, TopD = NULL, 
                               FormD = NULL, TLpred = NULL, fleet_names = "catch.1", Fish_mort, Fish_mort_acc) 
{ 

    ET_Main=data$ET_Main

  
  # 1) Other arguments of the function
  TL_out <- as.numeric(rownames(ET_Main))
  names(TL_out)=1:length(TL_out)
  n.TL=length(TL_out)
  
  #Initialization
  
  # 2) Calculates the number of fleets
  fleet=fleet_names ; n.fleet=length(fleet)
  
  # 3)  Defining what happens to the fishing effort multiplier, when the following arguments 
  #     are not set in the function call. 
  if (is.null(same.mE)){same.mE <- FALSE}                                       
  if(!is.null(fleet.of.interest)){same.mE <- FALSE}
  if (is.null(Mul_eff)){Mul_eff.=list()
  for(i in 1:n.fleet){Mul_eff.[[i]]=c(0, 0.2, 0.4, 0.7, 1, 1.5, 2, 2.5, 3, 4, 5)}
  if(same.mE){Mul_eff.=c(0, 0.2, 0.4, 0.7, 1, 1.5, 2, 2.5, 3, 4, 5)}
  }else{Mul_eff.=list()
  for(i in 1:n.fleet){Mul_eff.[[i]]=Mul_eff}
  if(same.mE){Mul_eff.=Mul_eff}
  print(Mul_eff.)
  }
  
  # 4) Here, it is defined what happens, if the other arguments are not set 
  #    in function call by the user 
  
  if (is.null(B.Input)){B.Input <- FALSE}
  if (is.null(Beta)){Beta <- .2}
  if (is.null(TopD)){TopD <- rep(.4,n.TL)}else{if(length(TopD)==1){TopD=rep(TopD,n.TL)}}
  if (is.null(FormD)){FormD <- rep(.5,n.TL)}else{if(length(FormD)==1){FormD=rep(FormD,n.TL)}}
  if (is.null(TLpred)){TLpred <- 3.5}
  
  
  # 6) Calculates the combination of effort multipliers, either a) using the same 
  #    multipliers for all fleets or b) using only fleet of interest.
  if(!same.mE){
    ff=expand.grid(Mul_eff.)
    if(is.null(fleet.of.interest)){
      for(n in 1:n.fleet){
        if(n==1){FF=ff[,n]}else{FF=paste(FF,'_',ff[,n],sep='')}}
      
    }else{
      colnames(ff)=c('interest','other')
      for(n in 1:n.fleet){
        if(fleet[n]%in%fleet.of.interest){ff.=ff[,'interest']}else{ff.=ff[,'other']}
        if(n==1){FF=ff.}else{FF=paste(FF,'_',ff.,sep='')}}
    }
  }else{
    for(n in 1:n.fleet){
      if(n==1){FF=Mul_eff.}else{FF=paste(FF,'_',Mul_eff.,sep='')}} 
    print(FF)
  }
  
  # 7) Setting the final output list of each fishing effort multiplier
  
  FF=unique(FF)
  comb=as.list(FF)
  names(comb)=FF
  
  # 8) Calculating all TEMPORAL outputs of comb, which are used as input in 
  #    the function mf.diagnosis. mf.diagnosis contains a conditional while loop, 
  #    which is used for an iterative solving of the different functions

  # Fishing mortality and accessible fishing mortality

  for(i in 1:length(comb)){
    comb[[i]]=list()
    comb[[i]][['mf']]=as.numeric(unlist(strsplit(names(comb)[i],'_')))
    names(comb[[i]][['mf']])=paste('mf',1:n.fleet,sep='')
    mf=as.numeric(comb[[i]][['mf']])
    for(j in 1:n.fleet){
      ff=mf[j]*Fish_mort[[j]]
      if(j==1){ff.=ff}else{ff.=ff.+ff}
    }
    comb[[i]][['Fish_mort']]=ff.
    comb[[i]][['Fish_mort_acc']]=ff./ET_Main$Selec
    
    # Kinetics and accessible kinetics
    comb[[i]][['TEMP_Kin']]=ET_Main[,'Kin']-ET_Main[,'Fish_mort']+comb[[i]][['Fish_mort']]
    comb[[i]][['TEMP_Kin_acc']]=ET_Main[,'Kin_acc']-ET_Main[,'Fish_mort_acc']+comb[[i]][['Fish_mort_acc']]
    comb[[i]][['Kin_MF']]=comb[[i]][['TEMP_Kin']]
    comb[[i]][['Kin_MF_acc']]=comb[[i]][['TEMP_Kin_acc']]
    
    # Productivity and accessible productivity
    comb[[i]][['Prod_MF']]=ET_Main[,'P']
    comb[[i]][['Prod_MF_acc']]=ET_Main[,'P_acc']
    
    # Flow and accessible flow
    comb[[i]][['flow']]=ET_Main[,'flow']
    comb[[i]][['flow_acc']]=ET_Main[,'flow_acc']
    
    # Flow and accessible flow
    comb[[i]][['flow_mean']]=ET_Main[,'flow_mean']
    comb[[i]][['flow_mean_acc']]=ET_Main[,'flow_mean_acc']
    
    # Biomass and accessible biomass
    comb[[i]][['BIOM_MF']]=ET_Main[,'B']
    comb[[i]][['BIOM_MF_acc']]=ET_Main[,'B_acc']
  }
  
  
  # 9.) CALCULATING OTHER ARGUMENTS NECESSARY FOR THE FOLLOWING MF.DIAGNOSIS
  tll=names(TL_out[TL_out>=2.8 & TL_out<=3.3])
  range.TLpred=as.numeric(c(tll[1],tll[length(tll)]))-2
  high.tl=abs(TL_out-5.6)
  lim.high.TL=as.numeric(names(high.tl[high.tl==min(high.tl)[1]]))
  tlll=names(TL_out[TL_out>=(TL_out[lim.high.TL]-0.5) & TL_out<=(TL_out[lim.high.TL])])
  range.highTL=as.numeric(c(tlll[1],tlll[length(tlll)]))
  
  
  # 10.) RUN MF.DIAGNOSIS ON EACH MULTIPLIER OF EFFORT, TO SAY ON EACH ELEMENT OF LIST COMB 
  diagn.list=lapply(comb,mf.diagnosis.adapted,ET_Main,TL_out,fleet,n.fleet,B.Input,
                    Beta,TopD,FormD,TLpred,n.TL,range.TLpred,lim.high.TL,range.highTL, Fish_mort, Fish_mort_acc)
  
   # 11.) PRODUCE OUTPUT
  # mf.diagnosis(comb[[10]],ET_Main,TL_out,fleet,n.fleet,Fish_mort_ref,Fish_mort_acc_ref,Beta,TopD,FormD,TLpred)
  names(diagn.list)=names(comb)
  diagn.list[['fleet.of.interest']]=fleet.of.interest
  class(diagn.list)<-"ETdiagnosis"
  return(diagn.list)
}
