load(paste(getwd(),"/SKM_N3.txt", sep=""))

KMSurv<-function(cenpfs, timepfs){
	pKM<-matrix(0,100,1)
		if (cenpfs==0){
			pKM[,1]<-psurvS[1,,1, (timepfs)]
		}
		if (cenpfs==1){
			pKM[,1]<-psurvS[1,,2,(timepfs)]
			
		}
		
		
	return(pKM)	
}
