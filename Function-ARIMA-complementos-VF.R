library(forecast)
library(xts)
library(astsa)
library(fUnitRoots)
library(berryFunctions)
library(lubridate)
#-------------------------------------------------------------------------------------------------
#Función para encontrar la lambda para la transformación BoxCox, con base en la metodología de Guerrero
lambdaG<-function(x){
  Nm<-length(x)
  g1<-first(x,Nm/2)
  g2<-last(x,Nm/2)
  
  G<-cbind(coredata(g1),coredata(g2))
  G<-as.matrix(G)
  colnames(G)<-c("G1", "G2")
  
  lm<-seq(-1,1,by = 0.5)
  P<-matrix(rep(0,2*length(lm)),2,length(lm))
  for (i in 1:length(lm)) {
    P[1,i]<-((var(G[,1]))^(1/2))/(mean(G[,1])^(1-i))
    P[2,i]<-((var(G[,2]))^(1/2))/(mean(G[,2])^(1-i))
  }
  colnames(P)<-as.character(lm)
  
  CV<-NULL#Vector de coeficientes de variación
  for (i in 1:length(lm)) {
    CV[i]<-((var(P[,i]))^(1/2))/mean(P[,i])
  }
  P<-rbind(P,t(CV))
  rownames(P)<-c("G1","G2","CV")
  return(list(lambda=lm[which.min(CV)], matriz=P))
}
#-------------------------------------------------------------------------------------------------
#Función para econtrar la "d" necesaria con base en el criterio Anderson (1976, p. 116)
anderson<-function(x){V<-NULL
V[1]<-var(x) #vector de varianzas con diferentes d
for (i in 1:3) {
  V[i+1]<-var(diff(x, differences = i)[-c(1:i)])
}
which(V==min(V, na.rm = TRUE))-1}
#--------------------------------------------------------------------------------------------------
#Función para encontrar los p-1 lags necesarios para una regresión ADF
ADFreg<-function(x,type="nc"){
  p1<-1
  Test<- unitrootTest(x, type = type, lags = p1)
    BoxTestp<-Box.test(Test@test$regression$residuals, lag = 24, type = "Ljung-Box", fitdf = p1+1+isTRUE(type=="c")*1)$p.value
  if(is.na(BoxTestp)){BoxTestp <- 0}
  while(BoxTestp<0.05 && p1<=20){
    p1<-p1+1
    Test<- unitrootTest(x, type = type, lags = p1)
    BoxTestp<-Box.test(Test@test$regression$residuals, lag = 24, type = "Ljung-Box", fitdf = p1+1)$p.value
    if(is.na(BoxTestp)){BoxTestp <- 0}
  }
  return(list(lags=p1,ADFtest=Test,LBtest=BoxTestp))
}
#--------------------------------------------------------------------------------------------------
#Función para evaluar si un modelo cumple con los criterios de estacionaridad e invertibilidad
admisible<-function(x){
  auxi<-which(names(x$coef)=="intercept")
  if(length(auxi)>0){coef<-x$coef[-auxi]}else{coef<-x$coef}
  aux<-NULL
  if(arimaorder(x)[1]>0){auxp<-c(1:arimaorder(x)[1])} else {auxp<-0}
  ar<-as.vector(c(1,-coef[auxp]))
  if(length(auxp)>0){
  ma<-as.vector(c(1,coef[-auxp]))}else{ma<-as.vector(c(1,coef))}
  
  testAR<-which(abs(polyroot(ar))<=1+.Machine$double.eps)
  testMA<-which(abs(polyroot(ma))<=1+.Machine$double.eps)
  
  isTRUE(length(testAR)==0 && length(testMA)==0)
}
#-------------------------------------------------------------------------------------------------
#Función para evaluar si la media de los residuo puede ser considerada 0
meanTEST<-function(x){
  N<-length(x$x)-sum(arimaorder(x)[c(1,2)])
  t<-sum(arimaorder(x)[c(1,2)])+1
  mu<-sum(x$residuals[t:length(x$x)],na.rm = TRUE)/N
  sig<-sqrt(sum((x$residuals[t:length(x$x)]-mu)^2,na.rm = TRUE)/(N-arimaorder(x)[3]))
  test<-sqrt(N)*mu/sig
  isTRUE(abs(test)<2)
}
#-------------------------------------------------------------------------------------------------
#Función para evaluar si los polinomios de retraso tienen algún factor en común
poly_FC<-function(x,epsilon=0.2){
auxi<-which(names(x$coef)=="intercept")
if(length(auxi)==0){
aux<-NULL
if(arimaorder(x)[1]>0){auxp<-c(1:arimaorder(x)[1])} else {auxp<-0}
ar<-as.vector(c(1,-x$coef[auxp]))
if(length(auxp)>0){
ma<-as.vector(c(1,x$coef[-auxp]))}else{ma<-as.vector(c(1,x$coef))}

testAR<-polyroot(ar)
testMA<-polyroot(ma)

testAR_i<-which(abs(Im(testAR))>.Machine$double.eps)
testMA_i<-which(abs(Im(testMA))>.Machine$double.eps)

if(length(testAR_i)>0){aux_AR<-testAR[-testAR_i]}else{aux_AR<-testAR}
if(length(testMA_i)>0){aux_MA<-testMA[-testMA_i]}else{aux_MA<-testMA}

T<-0
if(length(aux_AR)>0 && length(aux_MA)>0){
  for(i in 1:length(aux_AR)){
    for(j in 1:length(aux_MA)){
      T<-c(T,isTRUE(abs(aux_AR[i]-aux_MA[j])<=epsilon)*1)
    }
  }
}
return(list(estabilidad=isTRUE(sum(T)==0),rootsAR=testAR,rootsMA=testMA))}else{
  return(list(estabilidad=TRUE))
}
}
#-------------------------------------------------------------------------------------------------
#Función para evaluar si rk son significativamente distintos de cero (k>q)
Testrk<-function(x,q,k){
  FAC_P<-acf(x, lag = 24, na.action = na.omit, plot = FALSE)
  v<-sqrt((1+2*((q+2)%%(q+1))*sum((FAC_P$acf^2)[2:q+1]))/sum(!is.na(x)))
  return(isTRUE(abs(FAC_P$acf[k+1])>2*v))
}
#Función para evaluar si fik son significativamente distintos de cero
Testfik<-function(x,p){
  FAC_P<-pacf(x, lag = 24, na.action = na.omit, plot = FALSE)
  v<-1/sqrt(sum(!is.na(x)))
  return(isTRUE(abs(FAC_P$acf[p])>2*v))
}
#-------------------------------------------------------------------------------------------------
#Función para fijar el modelo inicial a evaluar inspirado en la intepretación de la FAC y la FACP

BoxJ<-function(x,maxp=7,maxq=7){
  stop<-0
  q<-0
  while(stop==0&&q<=maxq){stop<-(!Testrk(x,q,q+1))*(!Testrk(x,q,q+2));q<-q+1}
  Mq<-max(q-1,0)
  tail_q<-NULL
  for(i in 1:24){tail_q<-c(tail_q,Testrk(x,0,i)*1)}
  tail_q<-sum(tail_q)
  
  stop<-0
  p<-0
  while(stop==0 && p<=maxp){stop<-(!Testfik(x,p+1))*(!Testfik(x,p+2));p<-p+1}
  Mp<-max(p-1,0)
  tail_p<-NULL
  for(i in 1:24){tail_p<-c(tail_p,Testfik(x,i)*1)}
  tail_p<-sum(tail_p)
  
  FAC<-acf(x, lag = 24, na.action = na.omit, plot = FALSE)$acf[-1]
  FACP<-pacf(x, lag = 24, na.action = na.omit, plot = FALSE)$acf
  delta_p<-abs(FACP[-1])-abs(FACP[-24])
  delta_q<-abs(FAC[-1])-abs(FAC[-24])
  
  mp<-which(delta_p<(0.05))[1]
  mp<-mp/(!is.na(mp))
  mq<-which(delta_q<(0.05))[1]
  mq<-mq/(!is.na(mq))
  
  if(tail_q>tail_p){l<-list(p=Mp,q=0)}else{
    if(Mq<mq+mp){l<-list(p=0,q=Mq)}else{
      l<-list(p=mp,q=mq)
    }
  }
  
  return(l)
  
}
#----------------------------------------------------------------------------------------------------------------------
#Función para econtrar el modelo óptimo que cumpla con los supuestos de ruido blanco, siguiendo parsimonía + BoxJenkins
ARIMA_Par_s<-function(x,maxp=7,maxq=3, plot = TRUE, transform = TRUE, Anderson=TRUE, InformationCriterion = "AICc"){
  lambda1<-NULL
  if(transform == TRUE){
    lambda1<-lambdaG(x)$lambda
    x<-BoxCox(x, lambda = lambda1)}
  
  if(Anderson==TRUE){
    d<-anderson(x)}else{d<-0}
  if(d>0){x_d<-diff(x, differences = d); l<-ADFreg(as.ts(x_d),type = "nc");ADFt<-unitrootTest(as.ts(x_d), type = "nc",lags = l$lags)}else{x_d<-x; l<-ADFreg(as.ts(x_d),type = "c");ADFt<-unitrootTest(as.ts(x_d), type = "c",lags = l$lags)}
  while (ADFt@test$p.value[1]>=0.05) {
    d<-d+1
    x_d<-diff(x, differences = d)
    l<-ADFreg(as.ts(x_d),type = "nc")
    ADFt<-unitrootTest(as.ts(x_d), type = "nc",lags = l$lags)
  }
  
  test_pqi<-matrix(data = 0, nrow = maxp+1,ncol = maxq+1)
  pi<-BoxJ(x_d, maxp = maxp, maxq = maxq)$p
  qi<-BoxJ(x_d, maxp = maxp, maxq = maxq)$q
  
  
  fit<-0
  modelsAIC<-matrix(data = NA, nrow = maxp+1,ncol = maxq+1)
  cte<-matrix(data = isTRUE(d==0)*1, nrow = maxp+1, ncol = maxq+1)
  test<-matrix(data = 0, nrow = maxp+1,ncol = maxq+1)
  
  while(fit==0 && sum(test)<length(modelsAIC)){
    
    test_pqi[pi+1,qi+1]<-1
    
    minp_i<-max(pi-1,0, na.rm = TRUE)
    minq_i<-max(qi-1,0, na.rm = TRUE)
    maxp_i<-min(pi+1,maxp, na.rm = TRUE)
    maxq_i<-min(qi+1,maxq, na.rm = TRUE)
    
    p<-minp_i
    q<-minq_i
    
    steps<-length(c(minp_i:maxp_i))*length(c(minq_i:maxq_i))
    count<-0
    track<-NULL
    
    while (fit==0 && count<steps){
      
      if(!is.error(model<-Arima(x, order = c(p,d,q),include.mean = isTRUE(d==0)))&&test[p+1,q+1]==0){
        if(meanTEST(model)==FALSE){if(is.error(model<-Arima(x, order = c(p,d,q),include.mean = TRUE))==TRUE || d==0){fit1<-0}else{
          model<-Arima(x, order = c(p,d,q),include.mean = TRUE)
          if(meanTEST(model)==FALSE){fit1<-0}else{fit1<-1;cte[p+1,q+1]<-1}}}else{fit1<-1}
        
        if(admisible(model)==TRUE && poly_FC(model)$estabilidad==TRUE){fit2<-1}else{fit2<-0}
        
        BoxTestp<-Box.test(model$residuals, lag = 24, type = "Ljung-Box", fitdf = p+q+cte[p+1,q+1])$p.value
        if(is.na(BoxTestp)){BoxTestp <- 0}
        if(BoxTestp<0.05){fit3<-0}else{fit3<-1}
        fit<-fit1*fit2*fit3
        
        modelsAIC[p+1,q+1]<-switch(InformationCriterion,"AIC"=model$aic,"AICc"=model$aicc,"BIC"=model$bic)
        
        
      }
      track<-rbind(track,c(p,d,q))
      
      test[p+1,q+1]<-1
      
      count<-count+1
      if(p==minp_i || q==maxq_i ){if(max(track[,1])==maxp_i){p<-maxp_i;q<-track[max(which(track[,1]==maxp_i)),3]+1}else{p<-max(track[,1])+1;q<-minq_i}}else{p<-p-1;q<-q+1}
    }
    modelsAIC_i<-as.matrix(modelsAIC[c((minp_i+1):(maxp_i+1)),c((minq_i+1):(maxq_i+1))])
    
    if((length(modelsAIC_i)-sum(is.na(modelsAIC_i)))>1){
      pi<-which(modelsAIC==min(modelsAIC_i,na.rm = TRUE), arr.ind = TRUE)[1,1]-1
      qi<-which(modelsAIC==min(modelsAIC_i,na.rm = TRUE), arr.ind = TRUE)[1,2]-1
      
      count_aux<-1
      while(isTRUE(test_pqi[pi+1,qi+1]==1&&count_aux<(length(modelsAIC_i)-sum(is.na(modelsAIC_i))))){
        auxp<-which(modelsAIC_i==min(modelsAIC_i,na.rm = TRUE), arr.ind = TRUE)[1,1]
        auxq<-which(modelsAIC_i==min(modelsAIC_i,na.rm = TRUE), arr.ind = TRUE)[1,2]
        modelsAIC_i[auxp,auxq]<-max(modelsAIC_i, na.rm = TRUE)+1
        pi<-which(modelsAIC==min(modelsAIC_i,na.rm = TRUE), arr.ind = TRUE)[1,1]-1
        qi<-which(modelsAIC==min(modelsAIC_i,na.rm = TRUE), arr.ind = TRUE)[1,2]-1
        count_aux<-count_aux+1
      }
      if(count_aux==(length(modelsAIC_i)-sum(is.na(modelsAIC_i)))){
        pi<-which(test_pqi==0, arr.ind = TRUE)[1,1]-1
        qi<-which(test_pqi==0, arr.ind = TRUE)[1,2]-1
      }
    }else{pi<-which(test_pqi==0, arr.ind = TRUE)[1,1]-1;qi<-which(test_pqi==0, arr.ind = TRUE)[1,2]-1}
    
  }
  
  modelfit<-"The Model fits well"
  if(fit==0){
    modelfit<-"The Model do not fits well"
    p<-which(modelsAIC==min(modelsAIC,na.rm = TRUE), arr.ind = TRUE)[1,1]-1
    q<-which(modelsAIC==min(modelsAIC,na.rm = TRUE), arr.ind = TRUE)[1,2]-1
    model<-Arima(x, order = c(p,d,q), include.mean = (cte[p+1,q+1]==1))
  }
  if(plot == TRUE){checkresiduals(model, lag = 24, test = FALSE)}
  o<-arimaorder(model)
  return(list(BoxCox_lambda=lambda1, FIT=modelfit, model=model, WNCheck = Box.test(model$residuals, lag = 24, type = "Ljung-Box", o[1]+o[3]+cte[o[1]+1,o[3]+1])))
}
#---------------------------------------------------------------------------------------------------------------------------------------
#Función para econtrar el modelo óptimo que cumpla con los supuestos de ruido blanco siguiendo parsimonía
ARIMA_Par<-function(x,maxp=7,maxq=3, plot = TRUE, transform = TRUE, Anderson=TRUE, InformationCriterion = "AICc"){
  lambda1<-NULL
  if(transform == TRUE){
    lambda1<-lambdaG(x)$lambda
    x<-BoxCox(x, lambda = lambda1)}
  
  if(Anderson==TRUE){
    d<-anderson(x)}else{d<-0}
  if(d>0){x_d<-diff(x, differences = d); l<-ADFreg(as.ts(x_d),type = "nc");ADFt<-unitrootTest(as.ts(x_d), type = "nc",lags = l$lags)}else{x_d<-x; l<-ADFreg(as.ts(x_d),type = "c");ADFt<-unitrootTest(as.ts(x_d), type = "c",lags = l$lags)}
  while (ADFt@test$p.value[1]>=0.05) {
    d<-d+1
    x_d<-diff(x, differences = d)
    l<-ADFreg(as.ts(x_d),type = "nc")
    ADFt<-unitrootTest(as.ts(x_d), type = "nc",lags = l$lags)
  }
  
  modelsAIC<-matrix(data = NA, nrow = maxp+1, ncol = maxq+1)
  cte<-matrix(data = isTRUE(d==0)*1, nrow = maxp+1, ncol = maxq+1)
  
  count<-0
  p<-0
  q<-0
  track<-NULL
  fit<-0
  
  while (fit==0 && count<length(modelsAIC)) {
    if(!is.error(model<-Arima(x, order = c(p,d,q),include.mean = isTRUE(d==0)))){
      if(meanTEST(model)==FALSE){if(is.error(model<-Arima(x, order = c(p,d,q),include.mean = TRUE))==TRUE || d==0){fit1<-0}else{
        model<-Arima(x, order = c(p,d,q),include.mean = TRUE)
        if(meanTEST(model)==FALSE){fit1<-0}else{fit1<-1;cte[p+1,q+1]<-1}}}else{fit1<-1}
      
      if(admisible(model)==TRUE && poly_FC(model)$estabilidad==TRUE){fit2<-1}else{fit2<-0}
      
      BoxTestp<-Box.test(model$residuals, lag = 24, type = "Ljung-Box", fitdf = p+q+cte[p+1,q+1])$p.value
      if(is.na(BoxTestp)){BoxTestp <- 0}
      if(BoxTestp<0.05){fit3<-0}else{fit3<-1}
      fit<-fit1*fit2*fit3
      
      modelsAIC[p+1,q+1]<-switch(InformationCriterion,"AIC"=model$aic,"AICc"=model$aicc,"BIC"=model$bic)
      
      track<-rbind(track,arimaorder(model))}
    count<-count+1
    if(p==0 || q==maxq ){if(max(track[,1])==maxp){p<-maxp;q<-track[max(which(track[,1]==maxp)),3]+1}else{p<-max(track[,1])+1;q<-0}}else{p<-p-1;q<-q+1}
  }
  modelfit<-"The Model fits well"
  if(fit==0){
    modelfit<-"The Model do not fits well"
    p<-which(modelsAIC==min(modelsAIC,na.rm = TRUE), arr.ind = TRUE)[1,1]-1
    q<-which(modelsAIC==min(modelsAIC,na.rm = TRUE), arr.ind = TRUE)[1,2]-1
    model<-Arima(x, order = c(p,d,q), include.mean = (cte[p+1,q+1]==1))
  }
  if(plot == TRUE){checkresiduals(model, lag = 24, test = FALSE)}
  o<-arimaorder(model)
  return(list(BoxCox_lambda=lambda1, FIT=modelfit, model=model, WNCheck = Box.test(model$residuals, lag = 24, type = "Ljung-Box", o[1]+o[3]+cte[o[1]+1,o[3]+1])))
}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------