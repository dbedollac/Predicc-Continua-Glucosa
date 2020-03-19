Test.Forecast3_prob<-function(serie,Nobs=288,Track=TRUE,Anderson = FALSE, Parsimony = TRUE, BoxJenkins = TRUE,Horizon=12,InformationCriterion = "AIC",wndow=96, Transform=TRUE){
  if(length(serie)>=((Nobs-1)+wndow)){
    s<-serie
    f_DM<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    level=seq(5,95,5)
    lambdas<-NULL
    
    L_DM5<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM10<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM15<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM20<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM25<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM30<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM35<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM40<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM45<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM50<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM55<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM60<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM65<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM70<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM75<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM80<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM85<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM90<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    L_DM95<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    
    U_DM5<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM10<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM15<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM20<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM25<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM30<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM35<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM40<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM45<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM50<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM55<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM60<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM65<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM70<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM75<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM80<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM85<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM90<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    U_DM95<-matrix(data = NA,ncol = Horizon ,nrow = Nobs)
    
    time<-1
    track<-NULL
    fit<-NULL
    alg_time<-NULL
    Parameters_no<-NULL
    
    while (time<=Nobs) {
      
      x<-first(s,wndow)
      start_time<-Sys.time()
      m<-ARIMA(x,plot = FALSE,Anderson = Anderson, Parsimony = Parsimony, BoxJenkins = BoxJenkins, InformationCriterion = InformationCriterion, transform=Transform)
      lambdas<-c(lambdas,m$BoxCox_lambda)
      mf<-forecast(m$model, h=Horizon,level=level)
      End_time<-Sys.time()
      alg_time<-c(alg_time,as.numeric(difftime(End_time,start_time, units = "secs")))
      Parameters_no<-c(Parameters_no,length(m$model$coef))
      
      if(Track==TRUE){track<-rbind(track,arimaorder(m$model));fit<-c(fit,m$FIT)}
      
      f_DM[time,]<-coredata(mf$mean)
      
      L_DM5[time,]<-coredata(mf$lower[,1])
      L_DM10[time,]<-coredata(mf$lower[,2])
      L_DM15[time,]<-coredata(mf$lower[,3])
      L_DM20[time,]<-coredata(mf$lower[,4])
      L_DM25[time,]<-coredata(mf$lower[,5])
      L_DM30[time,]<-coredata(mf$lower[,6])
      L_DM35[time,]<-coredata(mf$lower[,7])
      L_DM40[time,]<-coredata(mf$lower[,8])
      L_DM45[time,]<-coredata(mf$lower[,9])
      L_DM50[time,]<-coredata(mf$lower[,10])
      L_DM55[time,]<-coredata(mf$lower[,11])
      L_DM60[time,]<-coredata(mf$lower[,12])
      L_DM65[time,]<-coredata(mf$lower[,13])
      L_DM70[time,]<-coredata(mf$lower[,14])
      L_DM75[time,]<-coredata(mf$lower[,15])
      L_DM80[time,]<-coredata(mf$lower[,16])
      L_DM85[time,]<-coredata(mf$lower[,17])
      L_DM90[time,]<-coredata(mf$lower[,18])
      L_DM95[time,]<-coredata(mf$lower[,19])
      
      U_DM5[time,]<-coredata(mf$upper[,1])
      U_DM10[time,]<-coredata(mf$upper[,2])
      U_DM15[time,]<-coredata(mf$upper[,3])
      U_DM20[time,]<-coredata(mf$upper[,4])
      U_DM25[time,]<-coredata(mf$upper[,5])
      U_DM30[time,]<-coredata(mf$upper[,6])
      U_DM35[time,]<-coredata(mf$upper[,7])
      U_DM40[time,]<-coredata(mf$upper[,8])
      U_DM45[time,]<-coredata(mf$upper[,9])
      U_DM50[time,]<-coredata(mf$upper[,10])
      U_DM55[time,]<-coredata(mf$upper[,11])
      U_DM60[time,]<-coredata(mf$upper[,12])
      U_DM65[time,]<-coredata(mf$upper[,13])
      U_DM70[time,]<-coredata(mf$upper[,14])
      U_DM75[time,]<-coredata(mf$upper[,15])
      U_DM80[time,]<-coredata(mf$upper[,16])
      U_DM85[time,]<-coredata(mf$upper[,17])
      U_DM90[time,]<-coredata(mf$upper[,18])
      U_DM95[time,]<-coredata(mf$upper[,19])
      
      s<-s[-1]
      time<-time+1
    }
    
    
    obs<-serie[c(which(index(serie)==last(first(serie,wndow+1),1)%>%index()):length(serie))]
    Estimation<-xts(f_DM,order.by = index(first(obs,Nobs)))
    
    L_DM5<-xts(L_DM5,order.by = index(first(obs,Nobs)))
    L_DM10<-xts(L_DM10,order.by = index(first(obs,Nobs)))
    L_DM15<-xts(L_DM15,order.by = index(first(obs,Nobs)))
    L_DM20<-xts(L_DM20,order.by = index(first(obs,Nobs)))
    L_DM25<-xts(L_DM25,order.by = index(first(obs,Nobs)))
    L_DM30<-xts(L_DM30,order.by = index(first(obs,Nobs)))
    L_DM35<-xts(L_DM35,order.by = index(first(obs,Nobs)))
    L_DM40<-xts(L_DM40,order.by = index(first(obs,Nobs)))
    L_DM45<-xts(L_DM45,order.by = index(first(obs,Nobs)))
    L_DM50<-xts(L_DM50,order.by = index(first(obs,Nobs)))
    L_DM55<-xts(L_DM55,order.by = index(first(obs,Nobs)))
    L_DM60<-xts(L_DM60,order.by = index(first(obs,Nobs)))
    L_DM65<-xts(L_DM65,order.by = index(first(obs,Nobs)))
    L_DM70<-xts(L_DM70,order.by = index(first(obs,Nobs)))
    L_DM75<-xts(L_DM75,order.by = index(first(obs,Nobs)))
    L_DM80<-xts(L_DM80,order.by = index(first(obs,Nobs)))
    L_DM85<-xts(L_DM85,order.by = index(first(obs,Nobs)))
    L_DM90<-xts(L_DM90,order.by = index(first(obs,Nobs)))
    L_DM95<-xts(L_DM95,order.by = index(first(obs,Nobs)))
    
    U_DM5<-xts(U_DM5,order.by = index(first(obs,Nobs)))
    U_DM10<-xts(U_DM10,order.by = index(first(obs,Nobs)))
    U_DM15<-xts(U_DM15,order.by = index(first(obs,Nobs)))
    U_DM20<-xts(U_DM20,order.by = index(first(obs,Nobs)))
    U_DM25<-xts(U_DM25,order.by = index(first(obs,Nobs)))
    U_DM30<-xts(U_DM30,order.by = index(first(obs,Nobs)))
    U_DM35<-xts(U_DM35,order.by = index(first(obs,Nobs)))
    U_DM40<-xts(U_DM40,order.by = index(first(obs,Nobs)))
    U_DM45<-xts(U_DM45,order.by = index(first(obs,Nobs)))
    U_DM50<-xts(U_DM50,order.by = index(first(obs,Nobs)))
    U_DM55<-xts(U_DM55,order.by = index(first(obs,Nobs)))
    U_DM60<-xts(U_DM60,order.by = index(first(obs,Nobs)))
    U_DM65<-xts(U_DM65,order.by = index(first(obs,Nobs)))
    U_DM70<-xts(U_DM70,order.by = index(first(obs,Nobs)))
    U_DM75<-xts(U_DM75,order.by = index(first(obs,Nobs)))
    U_DM80<-xts(U_DM80,order.by = index(first(obs,Nobs)))
    U_DM85<-xts(U_DM85,order.by = index(first(obs,Nobs)))
    U_DM90<-xts(U_DM90,order.by = index(first(obs,Nobs)))
    U_DM95<-xts(U_DM95,order.by = index(first(obs,Nobs)))
    
    for (i in 1:Horizon) {
      Estimation[,i]<-lag(Estimation[,i],k=i-1)
      
      L_DM5[,i]<-lag(L_DM5[,i],k=i-1)
      L_DM10[,i]<-lag(L_DM10[,i],k=i-1)
      L_DM15[,i]<-lag(L_DM15[,i],k=i-1)
      L_DM20[,i]<-lag(L_DM20[,i],k=i-1)
      L_DM25[,i]<-lag(L_DM25[,i],k=i-1)
      L_DM30[,i]<-lag(L_DM30[,i],k=i-1)
      L_DM35[,i]<-lag(L_DM35[,i],k=i-1)
      L_DM40[,i]<-lag(L_DM40[,i],k=i-1)
      L_DM45[,i]<-lag(L_DM45[,i],k=i-1)
      L_DM50[,i]<-lag(L_DM50[,i],k=i-1)
      L_DM55[,i]<-lag(L_DM55[,i],k=i-1)
      L_DM60[,i]<-lag(L_DM60[,i],k=i-1)
      L_DM65[,i]<-lag(L_DM65[,i],k=i-1)
      L_DM70[,i]<-lag(L_DM70[,i],k=i-1)
      L_DM75[,i]<-lag(L_DM75[,i],k=i-1)
      L_DM80[,i]<-lag(L_DM80[,i],k=i-1)
      L_DM85[,i]<-lag(L_DM85[,i],k=i-1)
      L_DM90[,i]<-lag(L_DM90[,i],k=i-1)
      L_DM95[,i]<-lag(L_DM95[,i],k=i-1)
      
      U_DM5[,i]<-lag(U_DM5[,i],k=i-1)
      U_DM10[,i]<-lag(U_DM10[,i],k=i-1)
      U_DM15[,i]<-lag(U_DM15[,i],k=i-1)
      U_DM20[,i]<-lag(U_DM20[,i],k=i-1)
      U_DM25[,i]<-lag(U_DM25[,i],k=i-1)
      U_DM30[,i]<-lag(U_DM30[,i],k=i-1)
      U_DM35[,i]<-lag(U_DM35[,i],k=i-1)
      U_DM40[,i]<-lag(U_DM40[,i],k=i-1)
      U_DM45[,i]<-lag(U_DM45[,i],k=i-1)
      U_DM50[,i]<-lag(U_DM50[,i],k=i-1)
      U_DM55[,i]<-lag(U_DM55[,i],k=i-1)
      U_DM60[,i]<-lag(U_DM60[,i],k=i-1)
      U_DM65[,i]<-lag(U_DM65[,i],k=i-1)
      U_DM70[,i]<-lag(U_DM70[,i],k=i-1)
      U_DM75[,i]<-lag(U_DM75[,i],k=i-1)
      U_DM80[,i]<-lag(U_DM80[,i],k=i-1)
      U_DM85[,i]<-lag(U_DM85[,i],k=i-1)
      U_DM90[,i]<-lag(U_DM90[,i],k=i-1)
      U_DM95[,i]<-lag(U_DM95[,i],k=i-1)
      
    }
    
    L_DM<-list(L5=L_DM5,L10=L_DM10,L15=L_DM15,L20=L_DM20,L25=L_DM25,L30=L_DM30,L35=L_DM35,L40=L_DM40,L45=L_DM45,L50=L_DM50,L55=L_DM55,L60=L_DM60,L65=L_DM65,L70=L_DM70,L75=L_DM75,L80=L_DM80,L85=L_DM85,L90=L_DM90,L95=L_DM95)
    U_DM<-list(U5=U_DM5,U10=U_DM10,U15=U_DM15,U20=U_DM20,U25=U_DM25,U30=U_DM30,U35=U_DM35,U40=U_DM40,U45=U_DM45,U50=U_DM50,U55=U_DM55,U60=U_DM60,U65=U_DM65,U70=U_DM70,U75=U_DM75,U80=U_DM80,U85=U_DM85,U90=U_DM90,U95=U_DM95)
    
    Measures<-NULL
    for(i in 1:Horizon){
      obs<-serie[c((which(index(serie)==last(first(serie,wndow+1),1)%>%index())+i-1):length(serie))]
      l<-accuracy(as.ts(na.omit(Estimation[,i])),as.ts(first(obs,Nobs-(i-1))))
      Measures<-rbind(Measures,l)
    }
    rownames(Measures)<-c(1:Horizon)
    
    return(list(Measures=Measures,track = track,fit=fit,forecast=Estimation, lower=L_DM,upper=U_DM, Time=c(mean(alg_time, na.rm = TRUE),max(alg_time, na.rm = TRUE)),Parameters_Num=c(mean(Parameters_no, na.rm = TRUE),median(Parameters_no, na.rm = TRUE),max(Parameters_no, na.rm = TRUE)),Lambdas=lambdas))
  }else{print("Serie too short")}
}