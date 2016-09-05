#function to calculate measurement error. 
load(paste(getwd(),"/Simulation3/Simulation",tt,ttt,"/SKM_N1.txt", sep=""))


KMden<-function(AffectedBreastR, AgeBreastR){
	pKM<-matrix(0,110,2)
		if (AffectedBreastR==0){
			pKM[,1]<-psurv[1,,1,AgeBreastR]*psurv[2,,1,AgeBreastR]*phaz[2,,1,AgeBreastR]
		}
		if (AffectedBreastR==1){
			pKM[,1]<-psurv[1,,2,AgeBreastR]*psurv[2,,2,AgeBreastR]*phaz[2,,2,AgeBreastR]
			
		}
		
	
		if (AffectedBreastR==0){
			pKM[,2]<-phaz[1,,1,AgeBreastR]*psurv[1,,1,AgeBreastR]*psurv[2,,1,AgeBreastR]
					}
		if (AffectedBreastR==1){
			pKM[,2]<-phaz[1,,2,AgeBreastR]*psurv[1,,2,AgeBreastR]*psurv[2,,2,AgeBreastR]
			
		}
		
	return(pKM)	
}
