create.ETdiagnosis.alt_logsol.integrated <- function(data, Mul_eff = NULL, fleet.of.interest = NULL, 
                               same.mE = NULL, B.Input=NULL, Beta = NULL, TopD = NULL, 
                               FormD = NULL, TLpred = NULL, fleet_names = "catch.1", 
                               Fish_mort, Fish_mort_acc) 
{ 

    ET_Main=data$ET_Main

  
  # 1) Other arguments of the function ----
  TL_out <- as.numeric(rownames(ET_Main))
  names(TL_out)=1:length(TL_out)
  n.TL=length(TL_out)
  
  # 2) Calculating the number of fleets ----
  fleet=fleet_names ; n.fleet=length(fleet)
  
  # 3)  Defining what happens to the arguments of the function, if not ---- 
  #     defined in the function call. ----
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
  
  if (is.null(B.Input)){B.Input <- FALSE}
  if (is.null(Beta)){Beta <- .2}
  if (is.null(TopD)){TopD <- rep(.4,n.TL)}else{if(length(TopD)==1){TopD=rep(TopD,n.TL)}}
  if (is.null(FormD)){FormD <- rep(.5,n.TL)}else{if(length(FormD)==1){FormD=rep(FormD,n.TL)}}
  if (is.null(TLpred)){TLpred <- 3.5}
  
  
  # 4) Calculating the combination of effort multipliers, either a) using the same ----
  #    multipliers for all fleets or b) using only fleet of interest. ----
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
  
  # 6) Setting the final output list of each fishing effort multiplier ----
  
  FF=unique(FF)
  comb=as.list(FF)
  names(comb)=FF
  
  # 7) Calculating all TEMPORAL outputs of comb, which are used as input in ----
  #    the function mf.diagnosis. 

   ## A) Calculating the new fishing mortalities by multiplying the effort ----
   #     multiplier with the reference fishing mortality vector ---- 
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
    comb[[i]][['Fish_mort_acc']]=ff./ET_Main$Selec    # accessible fishing mortality 
    
   ## B) Calculating a starting vector for the Kinetics and accessible kinetics ----
    #    to initiate the iterative solving in the mf.diagnosis function. ----
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
  
  
  # 8.) CALCULATING OTHER ARGUMENTS NECESSARY FOR THE FOLLOWING MF.DIAGNOSIS ----
  tll=names(TL_out[TL_out>=2.8 & TL_out<=3.3])
  range.TLpred=as.numeric(c(tll[1],tll[length(tll)]))-2
  high.tl=abs(TL_out-5.6)
  lim.high.TL=as.numeric(names(high.tl[high.tl==min(high.tl)[1]]))
  tlll=names(TL_out[TL_out>=(TL_out[lim.high.TL]-0.5) & TL_out<=(TL_out[lim.high.TL])])
  range.highTL=as.numeric(c(tlll[1],tlll[length(tlll)]))
  
  # 9.) RUN THE NORMAL NORMAL MF.DIAGNOSIS ON THE mE=0 AND STORE IT SEPERATELY ----
  comb0<-comb[["0"]]
  comb0_res<-mf.diagnosis.logsol.integrated(comb0,ET_Main,TL_out,fleet,n.fleet,B.Input,
                          Beta,TopD,FormD,TLpred,n.TL,range.TLpred,
                          lim.high.TL,range.highTL, Fish_mort, Fish_mort_acc)

  Bvirgin_acc=comb0_res$BIOM_MF_acc
  Bvirgin=comb0_res$BIOM_MF
  
  # 10.) RUN MF.DIAGNOSIS ON EACH MULTIPLIER OF EFFORT, TO SAY ON EACH ELEMENT OF LIST COMB ----
  diagn.list=lapply(comb,mf.diagnosis.alt.logsol.integrated,ET_Main,TL_out,fleet,n.fleet,B.Input,
                    Beta,TopD,FormD,TLpred,n.TL,range.TLpred,lim.high.TL,range.highTL, 
                    Fish_mort, Fish_mort_acc,
                    Bvirgin, Bvirgin_acc)
  
  # 11.) PRODUCE OUTPUT ----
  names(diagn.list)=names(comb)
  diagn.list[['fleet.of.interest']]=fleet.of.interest
  class(diagn.list)<-"ETdiagnosis"
  return(diagn.list)
}
