FalseAlarms<-M_FalseAlarmas_III
TrueAlarms<-M_TrueAlarms_III
MissedEvents<-M_MissedEvents_III
ReactionTime<-M_ReactionTime_III
measures_III<-meausures_III[-37,]

FalseAlarms_st<-M_FalseAlarmas_III_st
TrueAlarms_st<-M_TrueAlarms_III_st
MissedEvents_st<-M_MissedEvents_III_st
ReactionTime_st<-M_ReactionTime_III_st
measures_III_st<-meausures_III_st

colnames(FalseAlarms)<-c("p1",seq(5,95,5))
colnames(TrueAlarms)<-c("p1",seq(5,95,5))
colnames(MissedEvents)<-c("p1",seq(5,95,5))
colnames(ReactionTime)<-c("p1",seq(5,95,5))

colnames(FalseAlarms_st)<-c("p1",seq(5,95,5))
colnames(TrueAlarms_st)<-c("p1",seq(5,95,5))
colnames(MissedEvents_st)<-c("p1",seq(5,95,5))
colnames(ReactionTime_st)<-c("p1",seq(5,95,5))

x<-data.frame(III=measures_III[,"MAE"],III_st=measures_III_st[,"MAE"])
boxplot(x,outline =  FALSE,main="",horizontal = FALSE,cex=0.1,mex=0.05)

x<-data.frame(III=measures_III[,"RMSE"],III_st=measures_III_st[,"RMSE"])
boxplot(x,outline =  FALSE,main="",horizontal = FALSE,cex=0.1,mex=0.05)

median(measures_III[,"MAE"])
median(measures_III_st[,"RMSE"])

Sensitivity<-colSums(TrueAlarms)/(colSums(TrueAlarms)+colSums(MissedEvents))
False_rate<-colSums(FalseAlarms)/(colSums(TrueAlarms)+colSums(FalseAlarms))

Sensitivity_st<-colSums(TrueAlarms_st)/(colSums(TrueAlarms_st)+colSums(MissedEvents_st))
False_rate_st<-colSums(FalseAlarms_st)/(colSums(TrueAlarms_st)+colSums(FalseAlarms_st))

layout(mat = matrix(c(1,3,2,4),2,2))
par(mar=c(4.1, 3.1, 3, 1.1))
plot(x=seq(0,95,5),y=Sensitivity_st, col="blue",type="b", main="(A) Sensibilidad",xlab="",xlim=c(0,95),ylim=c(0,1),cex.main=1, ylab="")
lines(x=seq(0,95,5),y=Sensitivity,col="lightseagreen",type="b")
legend(x="bottomleft",legend = c("III_st","III"),col=c("blue","lightseagreen"),lty = c(1,1),pt.cex = 1,bty="n",y.intersp = 1)

plot(x=seq(0,95,5),y=False_rate_st, col="blue",type="b", main="(B) Porcentaje de falsas alarmas",xlab="",xlim=c(0,95),ylim=c(0,1),cex.main=1,ylab="")
lines(x=seq(0,95,5),y=False_rate,col="lightseagreen",type="b")
legend(x="bottomleft",legend = c("III_st","III"),col=c("blue","lightseagreen"),lty = c(1,1),pt.cex = 1,bty="n",y.intersp = 1)

aux_alarmas<-FalseAlarms+TrueAlarms
PromedioAlarmas_WoA_Pars<-colSums(aux_alarmas,na.rm = TRUE)/119
aux_alarmas<-FalseAlarms_st+TrueAlarms_st
PromedioAlarmas_WoA_Pars_st<-colSums(aux_alarmas,na.rm = TRUE)/119
barplot(rbind(PromedioAlarmas_WoA_Pars,PromedioAlarmas_WoA_Pars_st),main="(C) Alarmas diarias promedio",xlab = "Niveles de predicción",col=c("lightseagreen","blue"),ylim = c(0,8),beside = TRUE,names.arg = seq(0,95,5),cex.main=1,ylab="")
legend(x="topleft",legend = c("III","III_st"),col=c("lightseagreen","blue"),lty = c(1,1),pt.cex = 1,bty="n",y.intersp = 1)

aux_reaction_HA_WoAPars<-ReactionTime*TrueAlarms
aux_reaction_HA_WoAPars<-colSums(aux_reaction_HA_WoAPars, na.rm = TRUE)
reaction_HA_WoAPars<-aux_reaction_HA_WoAPars/colSums(TrueAlarms, na.rm = TRUE)

aux_reaction_HA_WoA_Pars_AICc_30_st<-ReactionTime_st*TrueAlarms_st
aux_reaction_HA_WoA_Pars_AICc_30_st<-colSums(aux_reaction_HA_WoA_Pars_AICc_30_st, na.rm = TRUE)
reaction_HA_WoA_Pars_AICc_30_st<-aux_reaction_HA_WoA_Pars_AICc_30_st/colSums(TrueAlarms_st, na.rm = TRUE)

barplot(rbind(reaction_HA_WoAPars,reaction_HA_WoA_Pars_AICc_30_st),main="(D) Tiempo de reacción promedio",xlab = "Niveles de predicción",col=c("lightseagreen","blue"),ylim = c(0,60),beside = TRUE,names.arg = seq(0,95,5),cex.main=1,ylab="")
legend(x="topleft",legend = c("III","III_st"),col=c("lightseagreen","blue"),lty = c(1,1),pt.cex = 1,bty="n",y.intersp = 1)
