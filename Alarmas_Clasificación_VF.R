ALARMS.class<-function(Eventos,Alarmas){
  
  Eventos_all<-NULL
  Alarmas_all<-NULL
  
  if(!is.null(Eventos$Fin)){
    
    Eventos_all<-Hypo_ident_Events(Eventos)$Inicio
  }
  
  if(!is.null(Alarmas$Fin)){
    Alarmas_all<-Alarmas$Inicio
    
    False_Alarms<-NULL
    True_Alarms<-NULL
    for (i in 1:length(Alarmas_all)) {
      l<-seq.POSIXt(from=Alarmas_all[i]+60*5,length.out = 12, by = "5 mins")
      if(length(intersect(l,Eventos_all))==0){False_Alarms<-c(False_Alarms,as.character(Alarmas_all[i]))}else{True_Alarms<-c(True_Alarms,as.character(Alarmas_all[i]))}
    }
    if(is.null(False_Alarms)==FALSE){False_Alarms<-as.POSIXct(False_Alarms, format="%Y-%m-%d %H:%M:%OS")}
    if(is.null(True_Alarms)==FALSE){True_Alarms<-as.POSIXct(True_Alarms, format="%Y-%m-%d %H:%M:%OS")}
    
    Eventos_sinAlarma<-NULL
    Reaction_time<-NULL
    if(!is.null(Eventos$Fin)){
      Reaction_time<-xts(x=NULL,order.by = Eventos_all)
      
      for (i in 1:length(Eventos_all)) {
        l<-seq.POSIXt(from=Eventos_all[i]-12*60*5,length.out = 12, by = "5 mins")
        if(length(intersect(l,Alarmas_all))==0){Eventos_sinAlarma<-c(Eventos_sinAlarma,as.character(Eventos_all[i]));Reaction_time[i]<-0}else{Reaction_time[i]<-as.numeric(difftime(Eventos_all[i],l[which(as.numeric(l)==intersect(l,Alarmas_all))[1]],units = "mins"))}
      }
      if(length(Eventos_sinAlarma)>0){Eventos_sinAlarma<-as.POSIXct(Eventos_sinAlarma, format="%Y-%m-%d %H:%M:%OS")}
    }
  }else{True_Alarms<-NULL;False_Alarms<-NULL;Eventos_sinAlarma<-Eventos_all;Reaction_time<-NA}
  
  Sensitivity<-length(True_Alarms)/(length(True_Alarms)+length(Eventos_sinAlarma))
  False_rate<-length(False_Alarms)/length(Alarmas_all)
  Dangerous_rate<-(intersect(Eventos_sinAlarma,Hypo_ident_Events(Eventos)$Inicio)%>%length())/(Hypo_ident_Events(Eventos)$Fin%>%length())
  
  return(list(Sensitivity=Sensitivity,False_Alarm_Rate=False_rate,DangerousRate=Dangerous_rate,TruePositive=True_Alarms,FalsePositive=False_Alarms,Missed_Events=Eventos_sinAlarma,ReactionTime=Reaction_time,HypoObs=length(Eventos_all)))
}