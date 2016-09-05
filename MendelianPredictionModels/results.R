#summarizing results
library("xtable")
aa<-matrix(0, nrow=16, ncol=9)

ResultsA<-matrix(0, 16,14 )
ResultsA<-Results[c(1:8,10:17),]
aa[,1]<-0
aa[,2]<-as.numeric(ResultsA[,5])
aa[,3]<-as.numeric(ResultsA[,6])
aa[,4]<-as.numeric(ResultsA[,9])
aa[,5]<-as.numeric(ResultsA[,10])
aa[,6]<-as.numeric(ResultsA[,11])
aa[,7]<-as.numeric(ResultsA[,12])
aa[,8]<-as.numeric(ResultsA[,13])
aa[,9]<-as.numeric(ResultsA[,14])


xtable(aa, digits=4)