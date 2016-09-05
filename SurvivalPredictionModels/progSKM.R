load(paste(getwd(),"/SKM_N1.txt", sep=""))

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
