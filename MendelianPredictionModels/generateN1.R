library("BayesMendel")
#library("kinship")

#GENERATING N1 Families
data(BRCApenet.metaDSL.2008,death.othercauses, compriskSurv) #load data

## BRCA1 allele frequency in general population 
#x<-0.0005829
x<-0.00609756097560976
## number of families
n=100000

#number of daughters
m1=3

source(paste(getwd(),"/progFamily.R", sep=""))

#generate families for risk prediction
set.seed(1234567) #for 1st set of 100000 families
fam1n1 = fam.gen(x,n,m1)
fam1n1<-fam1n1[which(fam1n1$ID!=4 & fam1n1$ID!=6 & fam1n1$ID!=8),]
write.table(fam1n1, paste(getwd(),"/N1.txt", sep=""), sep="\t")

#generate families to fit measurement error model
set.seed(123456) 
#GENERATING N2a Families
## number of families
n=50000
fam1n1 = fam.gen(x,n,m1)
fam1n1<-fam1n1[which(fam1n1$ID!=4 & fam1n1$ID!=6 & fam1n1$ID!=8),]
write.table(fam1n1, paste(getwd(),"/N2.txt", sep=""), sep="\t")

