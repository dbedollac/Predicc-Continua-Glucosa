ARIMA<-function(x,Parsimony = TRUE,BoxJenkins = TRUE ,maxp=7,maxq=3, plot = TRUE, transform = TRUE, Anderson=TRUE,InformationCriterion = "AICc"){
  if(Parsimony == TRUE){
    switch (BoxJenkins*1+1,
      sol<-ARIMA_Par(x,maxp=maxp,maxq=maxq, plot = plot, transform = transform, Anderson=Anderson, InformationCriterion = InformationCriterion),
      sol<-ARIMA_Par_s(x,maxp=maxp,maxq=maxq, plot = plot, transform = transform, Anderson=Anderson, InformationCriterion = InformationCriterion)
    )
  }
  else{
    switch (BoxJenkins*1+1,
    sol<-ARIMA_AIC(x,maxp=maxp,maxq=maxq, plot = plot, transform = transform, Anderson=Anderson, InformationCriterion = InformationCriterion),
    sol<-ARIMA_AIC_s(x,maxp=maxp,maxq=maxq, plot = plot, transform = transform, Anderson=Anderson, InformationCriterion = InformationCriterion)
    )
  }
  return(sol)
}
