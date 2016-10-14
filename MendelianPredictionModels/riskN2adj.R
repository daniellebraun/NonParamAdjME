#calculating adjusted risk prediction for each family 
library("BayesMendel")
set.seed(978543210) 
data(brca.fam,BRCApenet.metaDSL.2008,death.othercauses, compriskSurv, CBRCApenet.2012)

#proband either daughter or mother depending on tt
temp2<-c(1,3)
proband<-temp2[tt]
#obtain families
fam1n2<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,ttt,"/N2.txt", sep=""), sep="\t",  header=T)
fam1n2$AgeBreast<-ifelse(fam1n2$AgeBreast>110,110,fam1n2$AgeBreast)	
#set penetrance for BRCA2 and Ovarian equal to 0. 
BRCApenet<-BRCApenet.metaDSL.2008
BRCApenet$fFX[,4]<-0
BRCApenet$fFX[,5]<-0
BRCApenet$fFY[,1:5]<-0
BRCApenet$fMY[,1:5]<-0
BRCApenet$fMX[,4]<-0
BRCApenet$fMX[,5]<-0
#set competing risk equal to 1
cr<-compriskSurv
cr[,1:8]<-1

#creating fam1n2Val with true family history, and fam1n2Rep with reported error-prone family history
names(fam1n2)[names(fam1n2)=="Ethnic"] = "ethnic"
fam1n2$ethnic<-paste("Other")
fam1n2Val<-fam1n2[,1:12]
fam1n2Rep<-fam1n2[,c(1:4, 16, 6, 17,8:12) ]
names(fam1n2Rep)[names(fam1n2Rep)=="AgeBreastR"] = "AgeBreast"
names(fam1n2Rep)[names(fam1n2Rep)=="AffectedBreastR"] = "AffectedBreast"

#loading measurement error function
source(paste(getwd(),"/Simulation3/Simulation",tt,ttt,"/funSKM.R", sep=""))

m<-length(unique(fam1n2$FamilyID))

#creating matrix to store results
RepFix<-matrix(0,m,5)
#number of samples for integration 
n<-50

s<-special.index
aa1<-(s-1)*500+1
aa2<-s*500


for (t in aa1: aa2) {
	#storing family based on error-prone family history 
newR<-fam1n2Rep[which(fam1n2Rep$FamilyID==t),]


#Number of relatives 
fm<-dim(newR)[1]

#matrix to store reported ep age of diagnosis
tO<-matrix(0,1,fm)
#matrix to store reported ep disease status
dO<-matrix(0,1,fm)
#matrix to store P(H|H*) by age T 1:110 for delta=0, 111:220 for delta=1
pb<-matrix(0,220,fm)
#matrix to store sampled family histories H	
sf<-matrix(0,n,fm)
#matrix to store sampled true disease status	
d1<-matrix(0,n,fm)
#matrix to store sampled true age
t1<-matrix(0,n,fm)
#matrix to store probaility  P(H|H*) for each sampled family history
pr1<-matrix(0,n,fm)

#for non probands and females conduct adjustment
nem<-which(newR$ID!=proband & newR$Gender==0)
for (i in nem) {
#store reported ep age of diagnosis. 
tO[1,i]<-ifelse(newR$AgeBreast[i]==0, 1, ifelse(newR$AgeBreast[i]>110, 110, newR$AgeBreast[i]))
#store reported ep disease status
dO[1,i]<-ifelse(newR$AffectedBreast[i]==2, 1, newR$AffectedBreast[i])

#sample from ME distribution 
a1<-0
a1<-KMden(dO[1,i], tO[1,i]) 
a1[is.na(a1)]<-0
if(sum(a1)!=0) {
	a1<-a1/sum(a1)
#P(H|H*) by age T 1:110 for delta=0, 111:220 for delta=1
pb[,i]<-append(a1[,1], a1[,2])
#sample n family histories H, store as true disease status d1 and true age t1, and pr1, probability of observing P(d1,t1|d0,t0)
sf[,i]<-sample(220,n,replace=TRUE, prob=pb[,i])
	for (j in 1:n){
		if (sf[j,i]<=110){
		d1[j,i]<-0
		t1[j,i]<-sf[j,i]
		pr1[j,i]<-pb[sf[j,i],i]
		}
		if(sf[j,i]>110){
			d1[j,i]<-1
			t1[j,i]<-sf[j,i]-110
			pr1[j,i]<-pb[sf[j,i],i]
		}

	}
}
if (sum(a1)==0){
	d1[,i]<-dO[1,i]
	t1[,1]<-tO[1,i]
}


}
#create matrix to store results
bp<-matrix(0,n,4)
nf<-0	

#calculate risk for each n samples
for (j in 1:n){
	#set family to be error-prone family
	nf<-newR
	#assign family history H for relatives in nem based on samples t1 and d1
	for (i in nem){
	nf$AffectedBreast[i]<-ifelse(newR$AffectedBreast[i]==2 & d1[j,i]==1,2, d1[j,i])
	nf$AgeBreast[i]<-t1[j,i]
	nf$AgeBreastContralateral[i]<-ifelse(nf$AffectedBreast[i]==2,t1[j,i],0)
	}
	if (nf$AgeBreast[proband]==110) {nf$AgeBreast[proband]<-105}
	#calculate risk for each sample, assume AJ allele frequency 
	fambp<-as.numeric(try(brcapro(nf, counselee.id=proband, params=brcaparams(penetrance=BRCApenet, comprisk= cr, allef=list(c(1 - 0.00609756097560976, 0.00609756097560976),c(1,0))), print=FALSE)@probs))
	bp[j,1]<-fambp[1]
	bp[j,2]<-fambp[2]
	bp[j,3]<-fambp[3]
	bp[j,4]<-fambp[4]

	
}

	#average risk across n samples
RepFix[t,1]<-mean(bp[,1], na.rm=TRUE)
RepFix[t,2]<-mean(bp[,2], na.rm=TRUE)
RepFix[t,3]<-mean(bp[,3], na.rm=TRUE)
RepFix[t,4]<-mean(bp[,4], na.rm=TRUE)
RepFix[t,5]<-newR$FamilyID[1]
#familiesR[[t]]<-newR
}

save(RepFix, file = paste(getwd(),"/Simulation3/Simulation",tt,ttt,"/resultsAdj",special.index,".RData", sep=""))


