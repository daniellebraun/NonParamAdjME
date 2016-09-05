#plots
library(verification)
library(pROC)
setwd(getwd())
count<-1
temp<-c('a', 'b','c', 'd', 'e','f','g', 'h','i')
Results<-matrix(0,18,14)

jpeg('rplot1.jpg', width=5, height=5, units="in", res=1000)
#scenario 1
tt<-1
ttt<-1
temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")

par(mfrow=c(3,2), pty="s", mgp=c(2.5,1,0))
#plot(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="red4", ,cex.lab = 0.8, pch=20)

smoothScatter(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="red4", ,cex.lab = 0.8, pch=20, nrpoints = Inf)
abline(a=0, b=1)

title( main="Scenario 1: Sens=0.954, Spec=0.974", col.main="red4", font.main=4, line=3)
dev.off()
jpeg('rplot2.jpg', width=5, height=5, units="in", res=1000)
tt<-1
ttt<-5
temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")
			
#plot(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="magenta4",cex.lab = 0.8, pch=20)
smoothScatter(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="magenta4",cex.lab = 0.8, pch=20, nrpoints = Inf)

abline(a=0, b=1)
title(main="Scenario 2: Sens=0.649, Spec=0.990", col.main="magenta4", font.main=4, line=3)
dev.off()
jpeg('rplot3.jpg', width=5, height=5, units="in", res=1000)
tt<-1
ttt<-1
temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")

#plot(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="red4",cex.lab = 0.8, pch=20)
smoothScatter(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="red4",cex.lab = 0.8, pch=20, nrpoints = Inf)

abline(a=0, b=1)
dev.off()
jpeg('rplot4.jpg', width=5, height=5, units="in", res=1000)


tt<-1
ttt<-5
temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")
			
#plot(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="magenta4", cex.lab = 0.8, pch=20)
smoothScatter(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="magenta4", cex.lab = 0.8, pch=20, nrpoints = Inf)

abline(a=0, b=1)
dev.off()
jpeg('rplot5.jpg', width=5, height=5, units="in", res=1000)

tt<-1
ttt<-1
temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")

#plot(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="red4", cex.lab = 0.8, pch=20
)
smoothScatter(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="red4", cex.lab = 0.8, pch=20, nrpoints = Inf
)

abline(a=0, b=1)
dev.off()
jpeg('rplot6.jpg', width=5, height=5, units="in", res=1000)


tt<-1
ttt<-5
temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")

#plot(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="magenta4",cex.lab = 0.8, pch=20)
smoothScatter(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="magenta4",cex.lab = 0.8, pch=20, nrpoints = Inf)

abline(a=0, b=1)
dev.off()

smoothScatter(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="magenta4",cex.lab = 0.8, pch=20)
abline(a=0, b=1)
			
#jpeg('rplot.jpg', res=100000)

#dev.off()


# plot(f3$BRCA1vu[f3$BRCA1==1], f3$BRCA1ru[f3$BRCA1==1], xlab="BRCA1 based on True Family History", ylab="BRCA1 based on Error-Prone Family History", main="Carriers", xlim=c(0,1), ylim=c(0,1))
# abline(a=0, b=1)
# plot(f3$BRCA1vu[f3$BRCA1==1], f3$BRCA1ra[f3$BRCA1==1], xlab="BRCA1 based on True Family History", ylab="BRCA1 based on Adjustment", xlim=c(0,1), ylim=c(0,1))
# abline(a=0, b=1)
# plot(f3$BRCA1ru[f3$BRCA1==1], f3$BRCA1ra[f3$BRCA1==1], xlab="BRCA1 based on Error-Prone Family History", ylab="BRCA1 based on Adjustment", xlim=c(0,1), ylim=c(0,1))
# abline(a=0, b=1)
# plot(f3$BRCA1vu[f3$BRCA1==0], f3$BRCA1ru[f3$BRCA1==0], xlab="BRCA1 based on True Family History", ylab="BRCA1 based on Error-Prone Family History", main="Non-Carriers", xlim=c(0,1), ylim=c(0,1))
# abline(a=0, b=1)
# plot(f3$BRCA1vu[f3$BRCA1==0], f3$BRCA1ra[f3$BRCA1==0], xlab="BRCA1 based on True Family History", ylab="BRCA1 based on Adjustment", xlim=c(0,1), ylim=c(0,1))
# abline(a=0, b=1)
# plot(f3$BRCA1ru[f3$BRCA1==0], f3$BRCA1ra[f3$BRCA1==0], xlab="BRCA1 based on Error-Prone Family History", ylab="BRCA1 based on Adjustment", xlim=c(0,1), ylim=c(0,1))
# abline(a=0, b=1)

roc.area(f3$BRCA1,f3$BRCA1vu)$A
roc.area(f3$BRCA1,f3$BRCA1ru)$A
roc.area(f3$BRCA1,f3$BRCA1ra)$A


roc1<-roc(f3$BRCA1,f3$BRCA1vu)
roc2<-roc(f3$BRCA1,f3$BRCA1ru)
roc.test(roc1,roc2)





#scenario 1
tt<-1
ttt<-1
temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")

par(mfrow=c(2,3), pty="s", mgp=c(2,1,0), mar=c(1,3.5,0,0), oma=c(0.5,0.5,0.5,0.5))
#par(mar=c(1,4,4,1))
plot(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="#ff000030", ,cex.lab = 0.8)
abline(a=0, b=1)

plot(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#ff000030",cex.lab = 0.8)
abline(a=0, b=1)


plot(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#ff000030", cex.lab = 0.8
)
abline(a=0, b=1)
#par(mar=c(1,4,4,1))
#title( main="Scenario 1: Sens=0.954, Spec=0.974", col.main="red4", font.main=4, line=3)
tt<-1
ttt<-5
temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")
			
plot(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="#0000ff30",cex.lab = 0.8)
abline(a=0, b=1)
#title(main="Scenario 2: Sens=0.649, Spec=0.990", col.main="magenta4", font.main=4, line=3)
plot(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#0000ff30", cex.lab = 0.8)
abline(a=0, b=1)

plot(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#0000ff30",cex.lab = 0.8)
abline(a=0, b=1)

dev.copy2pdf(file="simprob.pdf")




#ind plots
tt<-1
ttt<-1
temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")

#par(mfrow=c(2,3), pty="s", mgp=c(2,1,0), mar=c(1,3.5,0,0), oma=c(0.5,0.5,0.5,0.5))
#par(mar=c(1,4,4,1))
plot(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="#ff000030", ,cex.lab = 1.3)
abline(a=0, b=1)
dev.copy2pdf(file="simprob1.pdf")
plot(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#ff000030",cex.lab = 1.3)
abline(a=0, b=1)
dev.copy2pdf(file="simprob2.pdf")


plot(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#ff000030", cex.lab = 1.3
)
abline(a=0, b=1)
dev.copy2pdf(file="simprob3.pdf")

#par(mar=c(1,4,4,1))
#title( main="Scenario 1: Sens=0.954, Spec=0.974", col.main="red4", font.main=4, line=3)
tt<-1
ttt<-5
temp2<-c(1,3)
proband<-temp2[tt]

load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsUnAdj",tt,temp[ttt],".RData", sep=""))
load(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/ResultsAdj",tt,temp[ttt],".RData", sep=""))

fam1n1<- read.table(paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/N2.txt", sep=""), sep="\t",  header=T)
fam1n1<-fam1n1[which(fam1n1$ID==proband),]

			#BRCA1
			ru<-cbind(RepUnAdj[,2], RepUnAdj[,5])
			colnames(ru)<-c("BRCA1ru", "FamilyID")
			
			vu<-cbind(ValUnAdj[,2], ValUnAdj[,5])
			colnames(vu)<-c("BRCA1vu", "FamilyID")
			
			ra<-cbind(RepAdj[,2], RepAdj[,5])
			colnames(ra)<-c("BRCA1ra", "FamilyID")
			
			f1<-merge(ru, vu, by="FamilyID")
			f2<-merge(f1, ra, by="FamilyID")
			f3<-merge(f2, fam1n1, by="FamilyID")
			
plot(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="#0000ff30",cex.lab = 1.3)
abline(a=0, b=1)
dev.copy2pdf(file="simprob4.pdf")

#title(main="Scenario 2: Sens=0.649, Spec=0.990", col.main="magenta4", font.main=4, line=3)
plot(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#0000ff30", cex.lab = 1.3)
abline(a=0, b=1)
dev.copy2pdf(file="simprob5.pdf")

plot(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#0000ff30",cex.lab = 1.3)
abline(a=0, b=1)

dev.copy2pdf(file="simprob6.pdf")