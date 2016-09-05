#calculating unadjusted risk prediction 
library("BayesMendel")

data(brca.fam,BRCApenet.metaDSL.2008,death.othercauses, compriskSurv, CBRCApenet.2012)

temp2<-c(1,3)
proband<-temp2[tt]
fam1n2<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,ttt,"/N2.txt", sep=""), sep="\t",  header=T)	
fam1n2$AgeBreast<-ifelse(fam1n2$AgeBreast>110,110,fam1n2$AgeBreast)

BRCApenet<-BRCApenet.metaDSL.2008
BRCApenet$fFX[,4]<-0
BRCApenet$fFX[,5]<-0
BRCApenet$fFY[,1:5]<-0
BRCApenet$fMY[,1:5]<-0
BRCApenet$fMX[,4]<-0
BRCApenet$fMX[,5]<-0

cr<-compriskSurv
cr[,1:8]<-1

#cleaning N2 families and sorting into two data sets

names(fam1n2)[names(fam1n2)=="Ethnic"] = "ethnic"
fam1n2$ethnic<-paste("Other")
fam1n2Val<-fam1n2[,1:12]
fam1n2Rep<-fam1n2[,c(1:4, 16, 6, 17,8:12) ]
names(fam1n2Rep)[names(fam1n2Rep)=="AgeBreastR"] = "AgeBreast"
names(fam1n2Rep)[names(fam1n2Rep)=="AffectedBreastR"] = "AffectedBreast"


#Running BRCAPro on validated and reported families 
familiesR<-list()
familiesV<-list()
m<-length(unique(fam1n2$FamilyID))
Rep = Val = matrix(0,m,5)
newR<-0
newV<-0
s<-special.index

a1<-(s-1)*10000+1
a2<-s*10000


for (t in a1: a2) {
newR<-fam1n2Rep[which(fam1n2Rep$FamilyID==t),]
newV<-fam1n2Val[which(fam1n2Val$FamilyID==t),]
if (newR$AgeBreast[proband]==110) {newR$AgeBreast[proband]<-105}
fambpR<-as.numeric(try(brcapro(newR, counselee.id=proband, params=brcaparams(penetrance=BRCApenet, comprisk= cr,allef=list(c(1 - 0.00609756097560976, 0.00609756097560976),c(1,0))), print=FALSE)@probs))
Rep[t,1]<-fambpR[1]
Rep[t,2]<-fambpR[2]
Rep[t,3]<-fambpR[3]
Rep[t,4]<-fambpR[4]
Rep[t,5]<-newR$FamilyID[1]
if (newV$AgeBreast[proband]==110) {newV$AgeBreast[proband]<-105}
fambpV<-as.numeric(try(brcapro(newV, counselee.id=proband, params=brcaparams(penetrance=BRCApenet, comprisk= cr, allef=list(c(1 - 0.00609756097560976, 0.00609756097560976),c(1,0))), print=FALSE)@probs))
Val[t,1]<-fambpV[1]
Val[t,2]<-fambpV[2]
Val[t,3]<-fambpV[3]
Val[t,4]<-fambpV[4]
Val[t,5]<-newV$FamilyID[1]
familiesR[[t]]<-newR
familiesV[[t]]<-newV

}

save(Rep, Val, file = paste(getwd(),"/Simulation3/Simulation",tt,ttt,"/resultsUnAdj",special.index,".RData", sep=""))