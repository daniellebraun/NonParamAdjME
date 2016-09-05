#calculating prediction summing model 
s<-special.index


for (a in 1:10){
	for (b in 1:10){
		for (z in s){
		
dat<-read.table(paste(getwd(),"/temp",a,b,"/temp",z,"/N1.txt", sep=""), sep="\t",  header=T)
ppv<-length(which(dat$ercenspfs==1 & dat$censpfs==1))/length(which(dat$ercenspfs==1))
npv<-length(which(dat$ercenspfs==0 & dat$censpfs==0))/length(which(dat$ercenspfs==0))

dat<-read.table(paste(getwd(),"/temp",a,b,"/temp",z,"/N2.txt", sep=""), sep="\t",  header=T)
n<-dim(dat)[1]
#based on gs
GSpred<-c()
GSpred2<-c()
#based on error
EPpred<-c()

load(paste(getwd(),"/temp",a,b,"/temp",z,"/SKM_N3.txt", sep=""))

KMSurv<-function(cenpfs, timepfs){
	pKM<-matrix(0,100,1)
		if (cenpfs==0){
			pKM[,1]<-psurvS[1,,1,(timepfs)]
		}
		if (cenpfs==1){
			pKM[,1]<-psurvS[1,,2,(timepfs)]
			
		}
		
		
	return(pKM)	
}



for (i in 1:n){
	timei<-50
	timeii<-50
GSpred[i]<-KMSurv(dat$censpfs[i], dat$timepfs[i])[timei]
GSpred2[i]<-KMSurv(dat$censpfs[i], dat$timepfs[i])[timeii]
EPpred[i]<-KMSurv(dat$ercenspfs[i], dat$ertimepfs[i])[timei]
}

calc<-matrix(,nrow=25, ncol=2)
Cpred<-c()
#corrected
load(paste(getwd(),"/temp",a,b,"/temp",z,"/SKM_N1.txt", sep=""))

KMden<-function(ercenpfs, ertimepfs){
	pKM<-matrix(0,25,2)
		if (ercenpfs==0){
			pKM[,1]<-psurv[1,,1,(ertimepfs)]*psurv[2,,1,(ertimepfs)]*phaz[2,,1,(ertimepfs)]
		}
		if (ercenpfs==1){
			pKM[,1]<-psurv[1,,2,(ertimepfs)]*psurv[2,,2,(ertimepfs)]*phaz[2,,2,(ertimepfs)]
			
		}
		
	
		if (ercenpfs==0){
			pKM[,2]<-phaz[1,,1,(ertimepfs)]*psurv[1,,1,(ertimepfs)]*psurv[2,,1,(ertimepfs)]
					}
		if (ercenpfs==1){
			pKM[,2]<-phaz[1,,2,(ertimepfs)]*psurv[1,,2,(ertimepfs)]*psurv[2,,2,(ertimepfs)]
			
		}
		
	return(pKM)	
}

#dat$ertimepfs<-ifelse(dat$ertimepfs>=6,6,dat$ertimepfs)
Cpred<-matrix(0, nrow=n*50, ncol=5)
Cpred[,1]<-rep(1:n, each=50)
Cpred[,2]<-rep(1:25, n*2)
Cpred[,3]<-rep(1:2)
for (k in 1:(n*50)){
		timei<-50
		if (sum(KMden(dat$ercenspfs[Cpred[k,1]],dat$ertimepfs[Cpred[k,1]]))!=0) {
		a1<-KMden(dat$ercenspfs[Cpred[k,1]],dat$ertimepfs[Cpred[k,1]])/sum(KMden(dat$ercenspfs[Cpred[k,1]],dat$ertimepfs[Cpred[k,1]]))
		Cpred[k,4]<-KMSurv((Cpred[k,3]-1),(Cpred[k,2]))[timei]*a1[Cpred[k,2],Cpred[k,3]]}	
		
							
		else {
				if (dat$ercenspfs[Cpred[k,1]]==0 & Cpred[k,3]==1)
				{a1<-npv}
				if (dat$ercenspfs[Cpred[k,1]]==0 & Cpred[k,3]==2)
				{a1<-1-npv}
				if (dat$ercenspfs[Cpred[k,1]]==1 & Cpred[k,3]==1)
				{a1<-1-ppv}
				if (dat$ercenspfs[Cpred[k,1]]==1 & Cpred[k,3]==2)
				{a1<-ppv}
				#Cpred[k,4]<-KMSurv((Cpred[k,3]-1),(Cpred[k,2]))[timei]*a1/25}
			Cpred[k,4]<-KMSurv(dat$ercenspfs[Cpred[k,1]],dat$ertimepfs[Cpred[k,1]])[timei]*1/50}

				if (dat$ercenspfs[Cpred[k,1]]==0 & Cpred[k,3]==1)
				{a1<-npv}
				if (dat$ercenspfs[Cpred[k,1]]==0 & Cpred[k,3]==2)
				{a1<-1-npv}
				if (dat$ercenspfs[Cpred[k,1]]==1 & Cpred[k,3]==1)
				{a1<-1-ppv}
				if (dat$ercenspfs[Cpred[k,1]]==1 & Cpred[k,3]==2)
				{a1<-ppv}
				Cpred[k,5]<-KMSurv((Cpred[k,3]-1),(Cpred[k,2]))[timei]*a1/25	
		}		
CP<-c()
CPnt<-c()
CP<-colSums(matrix(Cpred[,4], nrow=50))
CPnt<-colSums(matrix(Cpred[,5], nrow=50))
save(dat, GSpred, GSpred2, EPpred, CP, CPnt, file = paste(getwd(),"/temp",a,b,"/temp",z,"/results_aug2014_final_v4.RData", sep=""))
}}
}

