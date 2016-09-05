#merging adjusted results
setwd(getwd())
temp<-c('a', 'b','c', 'd', 'e','f','g', 'h','i')
for (tt in 1:2){
	for (ttt in 1:length(temp)){
	RepAdj<-matrix(0,50000, 5)


	for(i in 1:100){
    a1<-(i-1)*500+1
	a2<-i*500
	load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/resultsAdj",i,".RData", sep=""))
	RepAdj[a1:a2,1]<-RepFix[a1:a2,1] 
	RepAdj[a1:a2,2]<-RepFix[a1:a2,2] 
	RepAdj[a1:a2,3]<-RepFix[a1:a2,3] 
	RepAdj[a1:a2,4]<-RepFix[a1:a2,4] 
	RepAdj[a1:a2,5]<-RepFix[a1:a2,5]   
}
save(RepAdj, file = paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

}
}