
load(paste(getwd(), "/Results.RData", sep=""))
#MSE
tt1<-Results[,6]-Results[,6]
#sqrt(MSE)*1000
plot(Results[1:10,2], Results[1:10,10], xlim=c(0,10), ylim=c(0,400),type="n",xaxt="n", xlab="", col=rgb(255,0,0,200,maxColorValue=255), ylab=expression(sqrt(MSEP)%*%1000), frame.plot=FALSE, pch=17)

for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	tempp<-Results[a:b,2]+i
	points(tempp, Results[a:b,10], ylim=c(0,400),xlim=c(0,10),col=rgb(255,0,0,200,maxColorValue=255), pch=17)
}

for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	tempp<-Results[a:b,2]+i
	points(tempp, Results[a:b,9], ylim=c(0,400),xlim=c(0,10), pch=16, col=rgb(28,134,238,200,maxColorValue=255))
}

for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	tempp<-Results[a:b,2]+i
	points(tempp, Results[a:b,11], ylim=c(0,400), pch=18, col=rgb(160,32,240,200,maxColorValue=255))
}

for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	tempp<-Results[a:b,2]+i
	points(tempp, tt1[a:b], ylim=c(0,400), pch=15, col=rgb(0,100,0,200,maxColorValue=255))
}


temp<-c(0.1,1, 1.1, 2, 2.1, 3, 3.1, 4, 4.1, 5, 5.1, 6, 6.1, 7, 7.1, 8, 8.1, 9, 9.1, 10)
axis(1,temp, labels=FALSE)
mtext("Sensitivity",1,line=1,at=-1, cex=0.8)

mtext("0.1,..,1",1,line=1,at=0.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=1.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=2.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=3.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=4.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=5.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=6.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=7.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=8.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=9.5, cex=0.8)
mtext("Specificity",1,line=2,at=-1)
mtext("0.1",1,line=2,at=0.5)
mtext("0.2",1,line=2,at=1.5)
mtext("0.3",1,line=2,at=2.5)
mtext("0.4",1,line=2,at=3.5)
mtext("0.5",1,line=2,at=4.5)
mtext("0.6",1,line=2,at=5.5)
mtext("0.7",1,line=2,at=6.5)
mtext("0.8",1,line=2,at=7.5)
mtext("0.9",1,line=2,at=8.5)
mtext("1",1,line=2,at=9.5)
legend("topright", c("Error-Free","Error-Prone", "Full Adjustment", "Time Independent Adjustment"), pch =c(15,16,17,18), col=c(rgb(0,100,0,200,maxColorValue=255),col=rgb(28,134,238,200,maxColorValue=255), rgb(255,0,0,200,maxColorValue=255), rgb(160,32,240,200,maxColorValue=255)), cex=1)
dev.copy2pdf(file="PFS_MSE.pdf")

#ROC
plot(Results[1:10,2], Results[1:10,25], xlim=c(0,10), ylim=c(0,1),type="n",xaxt="n", xlab="", col=rgb(255,0,0,200,maxColorValue=255), ylab="ROC-AUC", frame.plot=FALSE, pch=17)
for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	points(Results[a:b,2]+i, Results[a:b,25], ylim=c(0,1), col=rgb(255,0,0,200,maxColorValue=255), pch=17)
}

for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	points(Results[a:b,2]+i, Results[a:b,24], ylim=c(0,1), col=rgb(28,134,238,200,maxColorValue=255), pch=16)
}

for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	points(Results[a:b,2]+i, Results[a:b,23], ylim=c(0,1), col=rgb(0,100,0,200,maxColorValue=255), pch=15)
}

for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	points(Results[a:b,2]+i, Results[a:b,26], ylim=c(0,1), col=rgb(160,32,240,200,maxColorValue=255), pch=18)
}

temp<-c(0.1,1, 1.1, 2, 2.1, 3, 3.1, 4, 4.1, 5, 5.1, 6, 6.1, 7, 7.1, 8, 8.1, 9, 9.1, 10)
axis(1,temp, labels=FALSE)
mtext("Sensitivity",1,line=1,at=-1, cex=0.8)

mtext("0.1,..,1",1,line=1,at=0.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=1.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=2.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=3.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=4.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=5.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=6.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=7.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=8.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=9.5, cex=0.8)
mtext("Specificity",1,line=2,at=-1)
mtext("0.1",1,line=2,at=0.5)
mtext("0.2",1,line=2,at=1.5)
mtext("0.3",1,line=2,at=2.5)
mtext("0.4",1,line=2,at=3.5)
mtext("0.5",1,line=2,at=4.5)
mtext("0.6",1,line=2,at=5.5)
mtext("0.7",1,line=2,at=6.5)
mtext("0.8",1,line=2,at=7.5)
mtext("0.9",1,line=2,at=8.5)
mtext("1",1,line=2,at=9.5)
legend("topright", c("Error-Free", "Error-Prone", "Full Adjustment", "Time Independent Adjustment"), pch =c(15,16,17, 18), col=c(rgb(0,100,0,200,maxColorValue=255),col=rgb(28,134,238,200,maxColorValue=255),rgb(255,0,0,200,maxColorValue=255), rgb(160,32,240,200,maxColorValue=255)), cex=1)
dev.copy2pdf(file="PFS_ROC.pdf")
#OE
plot(Results[1:10,2], Results[1:10,17], xlim=c(0,10), ylim=c(0.5,1.5),type="n",xaxt="n", xlab="", col=rgb(255,0,0,200,maxColorValue=255), ylab="O/E", frame.plot=FALSE, pch=17)
for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	points(Results[a:b,2]+i, Results[a:b,17], ylim=c(0.5,1.5), col=rgb(255,0,0,200,maxColorValue=255), pch=17)
}

for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	points(Results[a:b,2]+i, Results[a:b,16], ylim=c(0.5,1.5), col=rgb(28,134,238,200,maxColorValue=255), pch=16)
}

for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	points(Results[a:b,2]+i, Results[a:b,15], ylim=c(0.5,1.5), col=rgb(0,100,0,200,maxColorValue=255), pch=15)
}
for (i in 0:9){
	a<-i*10+1
	b<-(i+1)*10
	points(Results[a:b,2]+i, Results[a:b,18], ylim=c(0.5,1.5), col=rgb(160,32,240,200,maxColorValue=255), pch=18)
}
temp<-c(0.1,1, 1.1, 2, 2.1, 3, 3.1, 4, 4.1, 5, 5.1, 6, 6.1, 7, 7.1, 8, 8.1, 9, 9.1, 10)
axis(1,temp, labels=FALSE)
mtext("Sensitivity",1,line=1,at=-1, cex=0.8)

mtext("0.1,..,1",1,line=1,at=0.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=1.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=2.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=3.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=4.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=5.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=6.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=7.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=8.5, cex=0.8)
mtext("0.1,..,1",1,line=1,at=9.5, cex=0.8)
mtext("Specificity",1,line=2,at=-1)
mtext("0.1",1,line=2,at=0.5)
mtext("0.2",1,line=2,at=1.5)
mtext("0.3",1,line=2,at=2.5)
mtext("0.4",1,line=2,at=3.5)
mtext("0.5",1,line=2,at=4.5)
mtext("0.6",1,line=2,at=5.5)
mtext("0.7",1,line=2,at=6.5)
mtext("0.8",1,line=2,at=7.5)
mtext("0.9",1,line=2,at=8.5)
mtext("1",1,line=2,at=9.5)
abline(h=1)
legend("topleft", c("Error-Free", "Error-Prone", "Full Adjustment", "Time Independent Adjustment"), pch =c(15,16,17, 18), col=c(rgb(0,100,0,200,maxColorValue=255),rgb(28,134,238,200,maxColorValue=255),rgb(255,0,0,200,maxColorValue=255), rgb(160,32,240,200,maxColorValue=255)), cex=1)
dev.copy2pdf(file="PFS_OE.pdf")

