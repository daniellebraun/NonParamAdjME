#smoothed plots
setwd(getwd())


#scenario 1
tt<-1
ttt<-1
temp2<-c(1,3)
temp<-c('a', 'b','c', 'd', 'e','f','g', 'h','i')
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
par(mfrow=c(2,3), pty="s", mar=c(0,3.5,0,1), mgp=c(2,1,0))
par(mfrow=c(2,3), pty="s", oma = c(0,0,0,0) + 0.1,
    mar = c(0,2.8,0,1) + 0.1, mgp=c(2,1,0))

#par(mar=c(1,4,4,1))
#plot(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="#ff000030", ,cex.lab = 0.8)
smoothScatter(f3$BRCA1vu, f3$BRCA1ru, nrpoints = 0, colramp = colorRampPalette(c("white", "red","red1", "red2","red3","red4")), nbin=300,xlab="P(BRCA1) using Error-Free FH", ylab="P(BRCA1) using Error-Prone FH", xlim=c(0,1), ylim=c(0,1), col="red4", ,cex.lab = 0.8, cex.axis=0.8)
#axis(1,cex.axis=2)

abline(a=0, b=1)

#plot(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#ff000030",cex.lab = 0.8)
smoothScatter(f3$BRCA1vu, f3$BRCA1ra, nrpoints = 0, colramp = colorRampPalette(c("white", "red","red1", "red2","red3","red4")), nbin=300,xlab="P(BRCA1) using Error-Free FH", ylab="P(BRCA1) using Proposed Adj", xlim=c(0,1), ylim=c(0,1), col="red4", ,cex.lab = 0.8, cex.axis=0.8)
abline(a=0, b=1)


#plot(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#ff000030", cex.lab = 0.8
)

smoothScatter(f3$BRCA1ru, f3$BRCA1ra, nrpoints = 0, colramp = colorRampPalette(c("white", "red","red1", "red2","red3","red4")), nbin=300,xlab="P(BRCA1) using Error-Prone FH", ylab="P(BRCA1) using Proposed Adj", xlim=c(0,1), ylim=c(0,1), col="red4", ,cex.lab = 0.8, cex.axis=0.8)
abline(a=0, b=1)
#abline(a=0, b=1)
#par(mar=c(1,4,4,1))
#title( main="Scenario 1: Sens=0.954, Spec=0.974", col.main="red4", font.main=4, line=3)


smoothScatter(f3$BRCA1vu, f3$BRCA1ru, nrpoints = 0, colramp = colorRampPalette(c("white", "red","red1", "red2","red3","red4")), nbin=300,xlab="P(BRCA1) using Error-Free FH", ylab="P(BRCA1) using Error-Prone FH", xlim=c(0,0.2), ylim=c(0,0.2), col="red4", ,cex.lab = 0.8, cex.axis=0.8)
#axis(1,cex.axis=2)

abline(a=0, b=1)

#plot(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#ff000030",cex.lab = 0.8)
smoothScatter(f3$BRCA1vu, f3$BRCA1ra, nrpoints = 0, colramp = colorRampPalette(c("white", "red","red1", "red2","red3","red4")), nbin=300,xlab="P(BRCA1) using Error-Free FH", ylab="P(BRCA1) using Proposed Adj", xlim=c(0,0.2), ylim=c(0,0.2), col="red4", ,cex.lab = 0.8, cex.axis=0.8)
abline(a=0, b=1)


#plot(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#ff000030", cex.lab = 0.8
)

smoothScatter(f3$BRCA1ru, f3$BRCA1ra, nrpoints = 0, colramp = colorRampPalette(c("white", "red","red1", "red2","red3","red4")), nbin=300,xlab="P(BRCA1) using Error-Prone FH", ylab="P(BRCA1) using Proposed Adj", xlim=c(0,0.2), ylim=c(0,0.2), col="red4", ,cex.lab = 0.8, cex.axis=0.8)
abline(a=0, b=1)
#abline(a=0, b=1)
#par(mar=c(1,4,4,1))
#title( main="Scenario 1: Sens=0.954, Spec=0.974", col.main="red4", font.main=4, line=3)
dev.copy2pdf(file="simprob_smooth_scen1.pdf")


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
			
#plot(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="#0000ff30",cex.lab = 0.8)
smoothScatter(f3$BRCA1vu, f3$BRCA1ru, nrpoints = 0, colramp = colorRampPalette(c("white", "blue","blue1", "blue2","blue3","blue4")), nbin=300,xlab="P(BRCA1) using Error-Free FH", ylab="P(BRCA1) using Error-Prone FH", xlim=c(0,1), ylim=c(0,1), col="red4", ,cex.lab = 0.8, cex.axis=0.8)

abline(a=0, b=1)
#title(main="Scenario 2: Sens=0.649, Spec=0.990", col.main="magenta4", font.main=4, line=3)
#plot(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#0000ff30", cex.lab = 0.8)
smoothScatter(f3$BRCA1vu, f3$BRCA1ra, nrpoints = 0, colramp = colorRampPalette(c("white", "blue","blue", "blue2","blue3","blue4")), nbin=300,xlab="P(BRCA1) using Error-Free FH", ylab="P(BRCA1) using Proposed Adj", xlim=c(0,1), ylim=c(0,1), col="red4", ,cex.lab = 0.8, cex.axis=0.8)

abline(a=0, b=1)

#plot(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#0000ff30",cex.lab = 0.8)
smoothScatter(f3$BRCA1ru, f3$BRCA1ra, nrpoints = 0, colramp = colorRampPalette(c("white", "blue","blue1", "blue2","blue3","blue4")), nbin=300,xlab="P(BRCA1) using Error-Prone FH", ylab="P(BRCA1) using Proposed Adj", xlim=c(0,1), ylim=c(0,1), col="red4", ,cex.lab = 0.8, cex.axis=0.8)

abline(a=0, b=1)

#plot(f3$BRCA1vu, f3$BRCA1ru, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Error-Prone Family History", xlim=c(0,1), ylim=c(0,1), col="#0000ff30",cex.lab = 0.8)
smoothScatter(f3$BRCA1vu, f3$BRCA1ru, nrpoints = 0, colramp = colorRampPalette(c("white", "blue","blue1", "blue2","blue3","blue4")), nbin=300,xlab="P(BRCA1) using Error-Free FH", ylab="P(BRCA1) using Error-Prone FH", xlim=c(0,0.2), ylim=c(0,0.2), col="red4", ,cex.lab = 0.8, cex.axis=0.8)

abline(a=0, b=1)
#title(main="Scenario 2: Sens=0.649, Spec=0.990", col.main="magenta4", font.main=4, line=3)
#plot(f3$BRCA1vu, f3$BRCA1ra, xlab="P(BRCA1) using Error-Free Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#0000ff30", cex.lab = 0.8)
smoothScatter(f3$BRCA1vu, f3$BRCA1ra, nrpoints = 0, colramp = colorRampPalette(c("white", "blue","blue", "blue2","blue3","blue4")), nbin=300,xlab="P(BRCA1) using Error-Free FH", ylab="P(BRCA1) using Proposed Adj", xlim=c(0,0.2), ylim=c(0,0.2), col="red4", ,cex.lab = 0.8, cex.axis=0.8)

abline(a=0, b=1)

#plot(f3$BRCA1ru, f3$BRCA1ra, xlab="P(BRCA1) using Error-Prone Family History", ylab="P(BRCA1) using Proposed Adjustment", xlim=c(0,1), ylim=c(0,1), col="#0000ff30",cex.lab = 0.8)
smoothScatter(f3$BRCA1ru, f3$BRCA1ra, nrpoints = 0, colramp = colorRampPalette(c("white", "blue","blue1", "blue2","blue3","blue4")), nbin=300,xlab="P(BRCA1) using Error-Prone FH", ylab="P(BRCA1) using Proposed Adj", xlim=c(0,0.2), ylim=c(0,0.2), col="red4", ,cex.lab = 0.8, cex.axis=0.8)

abline(a=0, b=1)

dev.copy2pdf(file="simprob_smooth_scen2.pdf")


