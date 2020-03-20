M_FalseAlarmas_III<-NULL
M_TrueAlarms_III<-NULL
M_MissedEvents_III<-NULL
M_ReactionTime_III<-NULL


M_FalseAlarmas_III_st<-NULL
M_TrueAlarms_III_st<-NULL
M_MissedEvents_III_st<-NULL
M_ReactionTime_III_st<-NULL


#---------------------
S_FalseAlarmas_III<-NULL
S_TrueAlarms_III<-NULL
S_MissedEvents_III<-NULL
S_ReactionTime_III<-NULL

S_FalseAlarmas_III_st<-NULL
S_TrueAlarms_III_st<-NULL
S_MissedEvents_III_st<-NULL
S_ReactionTime_III_st<-NULL
#--------------------
meausures_III<-NULL
meausures_III_st<-NULL


Pax_Features<-read.csv("ITAM/Tesis/Pax_Features.csv")
Pax_Features<-Pax_Features[,-1]
Pax_Features[,2]<-NA
Pax_Features[,3]<-NA
colnames(Pax_Features)<-c("PaxID","M_Hypo","S_Hypo")
PaxID<-Pax_Features[,1]

hours<-12
wndow<-hours*12
error<-NULL


colnames(Sample)<-c(PaxID,PaxID)

for(i in 28:length(PaxID)){
  if(i!=37){
      try<-xts(as.numeric(Sample[,i+length(PaxID)]), order.by = as.POSIXct(Sample[,i],tz="", format="%Y-%m-%d %H:%M:%OS"))
      DatosEmc<-na.omit(try)
  
      obs<-DatosEmc[c((which(index(DatosEmc)==last(first(DatosEmc,wndow+1),1)%>%index())+(6-1)):length(DatosEmc))]
    
    if(is.error(III<-Test.Forecast3_prob(DatosEmc,Parsimony = TRUE,BoxJenkins = TRUE,Anderson = FALSE,Nobs = length(obs), wndow = 12*hours, Transform = TRUE))){error<-c(error,PaxID[i])}else{
      
      meausures_III<-rbind(meausures_III,III$Measures[6,])
      
      aux<-ALARMS_aux_prob(obs,III,Horizon = 6,th=70)
      M_FalseAlarmas_III<-rbind(M_FalseAlarmas_III,aux$A[1,])
      M_TrueAlarms_III<-rbind(M_TrueAlarms_III,aux$A[2,])
      M_MissedEvents_III<-rbind(M_MissedEvents_III,aux$A[3,])
      M_ReactionTime_III<-rbind(M_ReactionTime_III,aux$A[4,])
      M_HypoObs<-aux$HypoglycemiaObs
      
      aux<-ALARMS_aux_prob(obs,III,Horizon = 6,th=54)
      S_FalseAlarmas_III<-rbind(S_FalseAlarmas_III,aux$A[1,])
      S_TrueAlarms_III<-rbind(S_TrueAlarms_III,aux$A[2,])
      S_MissedEvents_III<-rbind(S_MissedEvents_III,aux$A[3,])
      S_ReactionTime_III<-rbind(S_ReactionTime_III,aux$A[4,])
      S_HypoObs<-aux$HypoglycemiaObs
      
      print(c(i,"III"))}

      if(is.error(III_st<-Test.Forecast3_prob(DatosEmc,Parsimony = TRUE,BoxJenkins = TRUE,Anderson = FALSE,Nobs = length(obs), wndow = 12*hours, Transform = FALSE))){error<-c(error,PaxID[i])}else{
        
        meausures_III_st<-rbind(meausures_III_st,III_st$Measures[6,])
        
        aux<-ALARMS_aux_prob(obs,III_st,Horizon = 6,th=70)
        M_FalseAlarmas_III_st<-rbind(M_FalseAlarmas_III_st,aux$A[1,])
        M_TrueAlarms_III_st<-rbind(M_TrueAlarms_III_st,aux$A[2,])
        M_MissedEvents_III_st<-rbind(M_MissedEvents_III_st,aux$A[3,])
        M_ReactionTime_III_st<-rbind(M_ReactionTime_III_st,aux$A[4,])

        aux<-ALARMS_aux_prob(obs,III_st,Horizon = 6,th=54)
        S_FalseAlarmas_III_st<-rbind(S_FalseAlarmas_III_st,aux$A[1,])
        S_TrueAlarms_III_st<-rbind(S_TrueAlarms_III_st,aux$A[2,])
        S_MissedEvents_III_st<-rbind(S_MissedEvents_III_st,aux$A[3,])
        S_ReactionTime_III_st<-rbind(S_ReactionTime_III_st,aux$A[4,])

        print(c(i,"III_st"))}
    
    Pax_Features[i,2]<-M_HypoObs
    Pax_Features[i,3]<-S_HypoObs
    
    print(error)
  }
}
