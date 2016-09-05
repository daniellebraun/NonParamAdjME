#introducing different types of measurement error
set.seed(9785432) #for n1
fam1n1<-read.table(paste(getwd(),"/Simulation3N1/N1.txt", sep=""), sep="\t",  header=T)
temp2<-c(1,3)
temp<-matrix(0,9, 4)
temp[,1]<-c('a', 'b','c', 'd', 'e','f','g', 'h','i')
temp[,2]<-c(0.954, 0.954,0.954, 0.954, 0.649,0.649,0.649, 0.649,0.954)
#temp[,2]<-as.numeric(temp[,2])
temp[,3]<-c(0.974, 0.974,0.974, 0.974, 0.99,0.99,0.99, 0.99,0.974)
#temp[,3]<-as.numeric(temp[,3])
temp[,4]<-c(5, 3, 1, 'o1', 5, 3, 1, 'o1', 'o2')
for (tt in 1:2){
proband<-temp2[tt]
for (ttt in 1:9){
#Generating "reported" information for fam1n1 "validated" members

#Introducing error in disease diagnosis 
#Introducing error only in breast cancer, only for relatives.


#Number of 1st relatives who have BC

BCA<-length(which((fam1n1$Gender!=1) &(fam1n1$ID!=proband)  & (fam1n1$AffectedBreast==1) ))

#Number of 1st relatives who don't have BC
BCNA<-length(which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) & (fam1n1$AffectedBreast==0)))

#Assigning random number for those who have BC
fam1n1$AffectedBreastA<-0
set.seed(976543)
fam1n1$AffectedBreastA[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) & (fam1n1$AffectedBreast==1))]<-sample(1:BCA,BCA,replace=F)

#Assigning random number for those who don't have BC
fam1n1$AffectedBreastNA<-0
set.seed(976542)
fam1n1$AffectedBreastNA[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) & (fam1n1$AffectedBreast==0))]<-sample(1:BCNA,BCNA,replace=F)

#Sensitivity and specificity from Al's paper for 1st degree relatives
#sens1BC<-0.954
#spec1BC<-0.974
sens1BC<-as.numeric(temp[ttt,2])
spec1BC<-as.numeric(temp[ttt,3])
cut1BCse<-sens1BC*BCA
cut1BCsp<-spec1BC*BCNA

fam1n1$AffectedBreastR<-fam1n1$AffectedBreast

#create reported data for those who have breast cancer
fam1n1$AffectedBreastR[which(fam1n1$AffectedBreastA>=cut1BCse)]<-0

#create reported data for those who do not have breast cancer
fam1n1$AffectedBreastR[which(fam1n1$AffectedBreastNA>=cut1BCsp)]<-1


#Introducing error in age
set.seed(976543)
fam1n1$AgeBreastR<-fam1n1$AgeBreast
if (temp[ttt,4]==1|temp[ttt,4]==3|temp[ttt,4]==5) {
a<-length(which((fam1n1$Gender!=1)&(fam1n1$ID!=proband)))
fam1n1$AgeBreastR[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband))]<-round(fam1n1$AgeBreast[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband))]+rnorm(a,0,as.numeric(temp[ttt,4])))
}
if (temp[ttt,4]=='o1') {
a<-length(which((fam1n1$Gender!=1)&(fam1n1$ID!=proband)))
fam1n1$AgeBreastR[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband))]<-round(fam1n1$AgeBreast[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband))]*rexp(a,1))
}
if(temp[ttt,4]=='o2'){
a<-length(which((fam1n1$Gender!=1)&(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==1)))
fam1n1$AgeBreastR[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) &(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==1))]<-round(fam1n1$AgeBreast[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==1))]+rnorm(a,2.4,1))
a<-length(which((fam1n1$Gender!=1)&(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==0)))
fam1n1$AgeBreastR[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) &(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==0))]<-round(fam1n1$AgeBreast[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==0))]+rnorm(a,5.7,1))
a<-length(which((fam1n1$Gender!=1)&(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==0) & (fam1n1$AffectedBreastR==1)))
fam1n1$AgeBreastR[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) &(fam1n1$AffectedBreast==0) & (fam1n1$AffectedBreastR==1))]<-round(fam1n1$AgeBreast[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==0) & (fam1n1$AffectedBreastR==1))]+rnorm(a,4.2,1))
}
fam1n1$AgeBreastR=ifelse(fam1n1$AgeBreastR>=110,110,fam1n1$AgeBreastR); fam1n1$AgeBreastR=ifelse(fam1n1$AgeBreastR<1,1,fam1n1$AgeBreastR)

write.table(fam1n1, paste(getwd(),"/Simulation3N1/Simulation",tt,temp[ttt,1],"/N1.txt", sep=""), sep="\t")
}}


#set.seed(9785432) #for n1
set.seed(978543) #for n2
fam1n1<-read.table(paste(getwd(),"/N2.txt", sep=""), sep="\t",  header=T)
temp2<-c(1,3)
temp<-matrix(0,9, 4)
temp[,1]<-c('a', 'b','c', 'd', 'e','f','g', 'h','i')
temp[,2]<-c(0.954, 0.954,0.954, 0.954, 0.649,0.649,0.649, 0.649,0.954)
#temp[,2]<-as.numeric(temp[,2])
temp[,3]<-c(0.974, 0.974,0.974, 0.974, 0.99,0.99,0.99, 0.99,0.974)
#temp[,3]<-as.numeric(temp[,3])
temp[,4]<-c(5, 3, 1, 'o1', 5, 3, 1, 'o1', 'o2')
for (tt in 1:2){
proband<-temp2[tt]
for (ttt in 1:9){
#Generating "reported" information for fam1n1 "validated" members

#Introducing error in disease diagnosis 
#Introducing error only in breast cancer, only for relatives.


#Number of 1st relatives who have BC

BCA<-length(which((fam1n1$Gender!=1) &(fam1n1$ID!=proband)  & (fam1n1$AffectedBreast==1) ))

#Number of 1st relatives who don't have BC
BCNA<-length(which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) & (fam1n1$AffectedBreast==0)))

#Assigning random number for those who have BC
fam1n1$AffectedBreastA<-0
set.seed(976543)
fam1n1$AffectedBreastA[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) & (fam1n1$AffectedBreast==1))]<-sample(1:BCA,BCA,replace=F)

#Assigning random number for those who don't have BC
fam1n1$AffectedBreastNA<-0
set.seed(976542)
fam1n1$AffectedBreastNA[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) & (fam1n1$AffectedBreast==0))]<-sample(1:BCNA,BCNA,replace=F)

#Sensitivity and specificity from Al's paper for 1st degree relatives
#sens1BC<-0.954
#spec1BC<-0.974
sens1BC<-as.numeric(temp[ttt,2])
spec1BC<-as.numeric(temp[ttt,3])
cut1BCse<-sens1BC*BCA
cut1BCsp<-spec1BC*BCNA

fam1n1$AffectedBreastR<-fam1n1$AffectedBreast

#create reported data for those who have breast cancer
fam1n1$AffectedBreastR[which(fam1n1$AffectedBreastA>=cut1BCse)]<-0

#create reported data for those who do not have breast cancer
fam1n1$AffectedBreastR[which(fam1n1$AffectedBreastNA>=cut1BCsp)]<-1


#Number of 2nd degree relatives who have BC
#BCA<-length(which((fam1n1$Gender!=1) &(fam1n1$ID!=proband)  & (fam1n1$AffectedBreast==1) & fam1n1$MotherID!=1))

#Number of 2nd degree relatives who don't have BC
#BCNA<-length(which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) & (fam1n1$AffectedBreast==0) & fam1n1$MotherID!=1))

#Assigning random number for those who have BC
#fam1n1$AffectedBreastA<-0
#fam1n1$AffectedBreastA[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) & (fam1n1$AffectedBreast==1) & fam1n1$MotherID!=1)]<-sample(1:BCA,BCA,replace=F)

#Assigning random number for those who don't have BC
#fam1n1$AffectedBreastNA<-0
#fam1n1$AffectedBreastNA[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) & (fam1n1$AffectedBreast==0) & fam1n1$MotherID!=1)]<-sample(1:BCNA,BCNA,replace=F)

#Sensitivity and specificity from Al's paper for 1st degree relatives
#sens1BC<-0.824
#spec1BC<-0.976
#cut1BCse<-sens1BC*BCA
#cut1BCsp<-spec1BC*BCNA

#fam1n1$AffectedBreastR<-fam1n1$AffectedBreast

#create reported data for those who have breast cancer
#fam1n1$AffectedBreastR[which(fam1n1$AffectedBreastA>=cut1BCse)]<-0

#create reported data for those who do not have breast cancer
#fam1n1$AffectedBreastR[which(fam1n1$AffectedBreastNA>=cut1BCsp)]<-1


#Introducing error in age
set.seed(976543)
fam1n1$AgeBreastR<-fam1n1$AgeBreast
if (temp[ttt,4]==1|temp[ttt,4]==3|temp[ttt,4]==5) {
a<-length(which((fam1n1$Gender!=1)&(fam1n1$ID!=proband)))
fam1n1$AgeBreastR[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband))]<-round(fam1n1$AgeBreast[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband))]+rnorm(a,0,as.numeric(temp[ttt,4])))
}
if (temp[ttt,4]=='o1') {
a<-length(which((fam1n1$Gender!=1)&(fam1n1$ID!=proband)))
fam1n1$AgeBreastR[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband))]<-round(fam1n1$AgeBreast[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband))]*rexp(a,1))
}
if(temp[ttt,4]=='o2'){
a<-length(which((fam1n1$Gender!=1)&(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==1)))
fam1n1$AgeBreastR[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) &(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==1))]<-round(fam1n1$AgeBreast[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==1))]+rnorm(a,2.4,1))
a<-length(which((fam1n1$Gender!=1)&(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==0)))
fam1n1$AgeBreastR[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) &(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==0))]<-round(fam1n1$AgeBreast[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==1) & (fam1n1$AffectedBreastR==0))]+rnorm(a,5.7,1))
a<-length(which((fam1n1$Gender!=1)&(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==0) & (fam1n1$AffectedBreastR==1)))
fam1n1$AgeBreastR[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband) &(fam1n1$AffectedBreast==0) & (fam1n1$AffectedBreastR==1))]<-round(fam1n1$AgeBreast[which((fam1n1$Gender!=1) &(fam1n1$ID!=proband)&(fam1n1$AffectedBreast==0) & (fam1n1$AffectedBreastR==1))]+rnorm(a,4.2,1))
}
fam1n1$AgeBreastR=ifelse(fam1n1$AgeBreastR>=110,110,fam1n1$AgeBreastR); fam1n1$AgeBreastR=ifelse(fam1n1$AgeBreastR<1,1,fam1n1$AgeBreastR)

write.table(fam1n1, paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt,1],"/N2.txt", sep=""), sep="\t")
}}