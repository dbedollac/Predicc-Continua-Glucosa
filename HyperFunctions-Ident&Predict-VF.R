library(lubridate)
#Función para identificar el inicio y fin de las observaciones en Hyperglicemia, con criterio del paper Yang
Hyper_ident<-function(serie,th=240){
  Inicio<-NULL
  Fin<-NULL
  serie[1]<-min(th-1,serie[1])
  for(i in 1:length(serie)-2){
    if(isTRUE(serie[i]<th&&serie[i+1]>=th&&serie[i+2]>=th))
    {inicio<-index(serie)[i+1]}else{inicio<-NULL}
    Inicio<-c(Inicio,as.character(inicio))
  }
  Inicio<-as.POSIXct(Inicio, format="%Y-%m-%d %H:%M:%OS")
  
  if(length(Inicio)>0){
    Fin<-NULL
    for(i in 1:length(Inicio)){
      count<-1
      while(isTRUE(serie[Inicio[i]+count*60*5]>=th)){count<-count+1}
      fin<-Inicio[i]+count*60*5
      Fin<-c(Fin,as.character(fin))
    }
    Fin<-as.POSIXct(Fin, format="%Y-%m-%d %H:%M:%OS")}
  
  return(list(Inicio=Inicio,Fin=Fin))
  
}
#---------------------------------------------------------------------------------------------------------
#Función para identificar el inicio y fin de los eventos hipoglucémicos, se considerará un mismo evento si hay un gap < 10 min entre observaciones
Hyper_ident_Events<-function(x){
  Inicio<-x$Inicio
  Fin<-x$Fin
  InicioE<-NULL
  FinE<-NULL
  if(!is.null(Fin)){
    InicioE<-as.character(Inicio[1])
    if(length(Inicio)>1){
      for(i in 1:(length(Fin)-1)){
        if(difftime(Inicio[i+1],Fin[i], units = "mins") >= 10){
          InicioE<-c(InicioE,as.character(Inicio[i+1]));FinE<-c(FinE,as.character(Fin[i]))}
      }}
    FinE<-FinE<-c(FinE,as.character(last(Fin,1)))
    InicioE<-as.POSIXct(InicioE, format="%Y-%m-%d %H:%M:%OS")
    FinE<-as.POSIXct(FinE, format="%Y-%m-%d %H:%M:%OS")
  }
  aux<-Hyper_predict_aux(InicioE,FinE)
  return(list(Inicio=aux$Inicio,Fin=aux$Fin))
}

#------------------------------------------------------------------------------------------
# Función para dividir las alarmas con máximo 1 hora de duración
Hyper_predict_aux<-function(xi,xf){
  Alarms<-NULL
  Fin<-NULL
  if(length(xi)>0){
    for(i in 1:length(xi)){
      aux<-difftime(xf[i],xi[i], units = "mins")
      if(aux>60){aux<-floor(aux/60)
      j<-0
      while(j<=aux){Alarms<-c(Alarms,as.character(xi[i]+j*60^2+1)); j<-j+1}
      }else{Alarms<-c(Alarms,as.character(xi[i]+1))}
    }
    Alarms<-as.POSIXct(Alarms, format="%Y-%m-%d %H:%M:%OS")
    Alarms<-floor_date(Alarms ,unit = "minute")
    aux1<-as.POSIXct(intersect(as.character(Alarms),as.character(xf)), format="%Y-%m-%d %H:%M:%OS")
    if(length(aux1)>0){
      for(i in 1:length(aux1)){
        auxi<-which(Alarms==aux1[i])
        Alarms<-Alarms[-auxi]}
    }
    
    Fin<-NULL
    for(i in 1:length(Alarms)){
      aux<-seq(Alarms[i],Alarms[i]+60^2,60*5)
      aux1<-intersect(as.character(aux+1),as.character(xf+1))
      if(length(aux1)>0){Fin<-c(Fin,aux1[1])}else{Fin<-c(Fin,as.character(aux[13]))}
    }
    Fin<-as.POSIXct(Fin, format="%Y-%m-%d %H:%M:%OS")
    Fin<-floor_date(Fin ,unit = "minute")
  }
  return(list(Inicio=Alarms,Fin=Fin))
}
#---------------------------------------------------------------------------------------
#Función para identificar el momento de inicio de Alarmas los forecast
Hyper_predict<-function(Estimation,Criterio="p1",Horizon=6,Upper=NULL,Lower=NULL,th=240){
  
  Alarms<-NULL
  Fin<-NULL
  Horizon<-min(dim(Estimation)[2],Horizon)
  
  #Utilizando estimación puntual como paper Yang
  if(Criterio=="p1"){
    fH<-Estimation[,Horizon]
    b<-which(!is.na(fH))[1]-1
    fH[b]<-th-1
    
    for(i in 1:length(fH)-2){
      if(isTRUE(fH[i]<th&&fH[i+1]>=th&&fH[i+2]>=th))
      {Alarm<-index(fH)[i+1]-(60*(Horizon-1)*5)+1}else{Alarm<-NULL}
      Alarms<-c(Alarms,as.character(Alarm))
    }
    
    if(length(Alarms)>0){
      Alarms<-as.POSIXct(Alarms, format="%Y-%m-%d %H:%M:%OS")
      Alarms<-floor_date(Alarms ,unit = "minute")
      
      for(i in 1:length(Alarms)){
        count<-1
        T<-TRUE
        while(T){T<-isTRUE((fH[Alarms[i]+60*(Horizon-1)*5+count*5*60])>=th);count<-count+1}
        fin<-Alarms[i]+(count-2)*60*5+1
        Fin<-c(Fin,as.character(fin))
      }
      Fin<-as.POSIXct(Fin, format="%Y-%m-%d %H:%M:%OS")
      Fin<-floor_date(Fin ,unit = "minute")
    }
  }
  
  #Utilizando estimación de regiones sesgado a reportar Alarma A(H0:BG  <=70mg/dl)
  if(Criterio=="I1"){
    fH<-Estimation[,Horizon]
    U<-Upper[,Horizon]
    L<-Lower[,Horizon]
    
    b<-which(!is.na(L))[1]-1
    U[b]<-th-1
    
    for(i in 1:length(fH)-2){
      T1<-isTRUE(U[i]<th)
      T2<-isTRUE(L[i+1]>=th)
      T3<-isTRUE(L[i+2]>=th)
      
      if(T1*T2*T3){Alarm<-index(fH)[i+1]-(60*(Horizon-1)*5)+1}else{Alarm<-NULL}
      Alarms<-c(Alarms,as.character(Alarm))
    }
    
    if(length(Alarms)>0){
      Alarms<-as.POSIXct(Alarms, format="%Y-%m-%d %H:%M:%OS")
      Alarms<-floor_date(Alarms ,unit = "minute")
      
      for(i in 1:length(Alarms)){
        count<-1
        T<-TRUE
        while(T){T<-isTRUE((L[Alarms[i]+60*(Horizon-1)*5+count*5*60])>=th);count<-count+1}
        fin<-Alarms[i]+(count-2)*60*5+1
        Fin<-c(Fin,as.character(fin))
      }
      Fin<-as.POSIXct(Fin, format="%Y-%m-%d %H:%M:%OS")
      Fin<-floor_date(Fin ,unit = "minute")
    }
  }
  
  return(Hyper_predict_aux(Alarms,Fin))
}