#merging unadjusted results
setwd(getwd())
temp<-c('a', 'b','c', 'd', 'e','f','g', 'h','i')
for (tt in 1:2){
	for (ttt in 1:length(temp)){
RepUnAdj<-matrix(0,50000, 5)
ValUnAdj<-matrix(0,50000, 5)
for(i in 1:5){
    a1<-(i-1)*10000+1
	a2<-i*10000
	load(paste(getwd(), "/Simulation3/Simulation",tt,temp[ttt],"/resultsUnAdj",i,".RData", sep=""))
	RepUnAdj[a1:a2,1]<-Rep[a1:a2,1]
	RepUnAdj[a1:a2,2]<-Rep[a1:a2,2]
	RepUnAdj[a1:a2,3]<-Rep[a1:a2,3]
	RepUnAdj[a1:a2,4]<-Rep[a1:a2,4]
	RepUnAdj[a1:a2,5]<-Rep[a1:a2,5]
	ValUnAdj[a1:a2,1]<-Val[a1:a2,1]
	ValUnAdj[a1:a2,2]<-Val[a1:a2,2]
	ValUnAdj[a1:a2,3]<-Val[a1:a2,3]
	ValUnAdj[a1:a2,4]<-Val[a1:a2,4]
	ValUnAdj[a1:a2,5]<-Val[a1:a2,5]   
}
save(RepUnAdj, ValUnAdj, file = paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
}
}

