
F_mort_fixed=function(dat, const=0.1){
  ET_Main=dat$ET_Main
  F_acc_new=c(0,rep(const,51))
  F_new<-F_acc_new*ET_Main$Selec
  Fish_mort=list(catch.1=F_new)
  Fish_mort_acc=list(catch.1=F_acc_new)
  Fish_morts<-list(Fish_mort=Fish_mort, Fish_mort_acc=Fish_mort_acc)
  return(Fish_morts)
  }
