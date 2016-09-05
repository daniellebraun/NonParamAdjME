#merging results
library(verification)
setwd(getwd())
count<-1
temp<-c('a', 'b','c', 'd', 'e','f','g', 'h','i')
Results<-matrix(0,18,14)
for (tt in 1:2){
	for (ttt in 1:length(temp)){

temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")
			
			Results[count,1]<-tt
			Results[count,2]<-temp[ttt]
			
			#MSE
			l<-length(which(is.na(f3$BRCA1ra)==FALSE))
			MSEc<-sum((f3$BRCA1ra-f3$BRCA1vu)^2, na.rm=T)/l
			MSE<-sum((f3$BRCA1ru-f3$BRCA1vu)^2, na.rm=T)/l
			
			Results[count,3]<-MSE
			Results[count,4]<-MSEc

			Results[count,5]<-sqrt(MSE)*1000
			Results[count,6]<-sqrt(MSEc)*1000
			
			#BRIER
			BSc<-sum((f3$BRCA1ra-f3$BRCA1)^2, na.rm=T)/l
			BS<-sum((f3$BRCA1ru-f3$BRCA1)^2, na.rm=T)/l
			
			Results[count,7]<-BSc
			Results[count,8]<-BS
			
			#OE
			f3<-f3[which(is.na(f3$BRCA1ra)==FALSE),]
			f3<-f3[which(is.na(f3$BRCA1ru)==FALSE),]
			
			OEv<-sum(f3$BRCA1)/sum(f3$BRCA1vu)
			OEr<-sum(f3$BRCA1)/sum(f3$BRCA1ru)
			OEc<-sum(f3$BRCA1)/sum(f3$BRCA1ra)
			
			
			Results[count,9]<-OEv
			Results[count,10]<-OEr
			Results[count,11]<-OEc
			
			Results[count,12]<-roc.area(f3$BRCA1,f3$BRCA1vu)$A
			Results[count,13]<-roc.area(f3$BRCA1,f3$BRCA1ru)$A
			Results[count,14]<-roc.area(f3$BRCA1,f3$BRCA1ra)$A
			
			
			count<-count+1
			}
			
}			


save(Results, file = paste(getwd(),"/Simulation3/ResultsFinal.RData", sep=""))








