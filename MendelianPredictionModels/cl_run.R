#running on cluster

setwd(getwd())

temp<-c('a', 'b','c', 'd', 'e','f','g', 'h','i')
for (tt in 1:2){
	for (ttt in 1:length(temp)){
		for(i in 1:100){
     system(paste("bsub -e e.txt -o o.txt -q normal_serial R CMD BATCH ",
	               getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/Run_function",i,".R /dev/null", sep=""));
}}}




