set.seed(4493729)
for (a in 1:10){
	for (b in 1:10){
		for (k in 1:100){
		sens<-a/10
		spec<-b/10
n<-1000

#pfs distribution
ftpfs<-rweibull(n, 1.456, scale = 11.063)

u<-runif(n, 0,1)
ctpfs<- 25*(1-sqrt(1-u))

ftpfs<-ceiling(ftpfs)
ctpfs<-ceiling(ctpfs)

censpfs<-ifelse(ftpfs<=ctpfs,1,0)
timepfs<-ifelse(ftpfs<=ctpfs,ftpfs, ctpfs)
dat<-cbind(timepfs,censpfs, ftpfs, ctpfs)
dat<-as.data.frame(dat)

#measurement error
#Number of patients who have PFS
PFSA<-length(which(dat$censpfs==1))

#Number of patients who don't have PFS
PFSNA<-length(which(dat$censpfs==0))

#Assigning random number for those who have PFS
dat$pfsA<-0
dat$pfsA[which(dat$censpfs==1)]<-sample(1:PFSA,PFSA,replace=F)

#Assigning random number for those who don't have PFS
dat$pfsNA<-0
dat$pfsNA[which(dat$censpfs==0 )]<-sample(1:PFSNA,PFSNA,replace=F)

#Sensitivity and specificity from Korn

cutse<-sens*PFSA
cutsp<-spec*PFSNA

dat$ercenspfs<-dat$censpfs

#create error data for those who have PFS
dat$ercenspfs[which(dat$pfsA>cutse)]<-0

#create error data for those who do not have PFS
dat$ercenspfs[which(dat$pfsNA>cutsp)]<-1

#Introducing error in time
dat$ertimepfs<-dat$timepfs

#those who have non-matching error and true pfs status
dat$ertimepfs[which((dat$ercenspfs==1) &(dat$censpfs==0))]<-dat$ftpfs[which((dat$ercenspfs==1) &(dat$censpfs==0))]

dat$ertimepfs[which((dat$ercenspfs==0) &(dat$censpfs==1))]<-dat$ctpfs[which((dat$ercenspfs==0) &(dat$censpfs==1))]

#those who have matching error and true pfs status that is 1

temp<-which((dat$ercenspfs==1) &(dat$censpfs==1))
temp2<-length(temp)
dat$rand<-NA
dat$rand[temp]<-rbinom(temp2, size=1, prob=0.55)
temp3<-which(dat$rand==0)
l2<-length(temp3)
dat$rl<-NA
dat$rl[temp3]<-rlnorm(l2, meanlog=0, sdlog = log(1.5))
dat$mp<-NA
dat$mp[temp3]<-dat$timepfs[temp3]*dat$rl[temp3]
	
dat$ertimepfs[temp3]<-ifelse(dat$mp[temp3]<=25, dat$mp[temp3], 25)
dat$ertimepfs[which(dat$ertimepfs>25)]<-25
dat$ertimepfs<-ceiling(dat$ertimepfs)
	
#generating OS from PFS
#for those who had PFS, 0.09 were deaths, treat them as progression
dat$timeos<-999
dat$censos<-999

temp1<-which(dat$censpfs==1)
ltemp1<-length(which(dat$censpfs==1))

dat$timeos<-60
dat$censos<-0

#those who had pfs, suppose 50% died, if died assume on average died a year after PFS them OStime of N(12,5)
dat$censos[temp1]<-ifelse(rbinom(ltemp1,1,0.5)==1,1,0)
dat$timeos[which(dat$censos==1)]<-abs(rnorm(length(which(dat$censos==1)), 12, 5))

#for those who had no PFS event, suppose 10% of them died, assume on average died 4 years after PFS, OS=PFS+N(48,5)
dat$censos[which(dat$censpfs==0)]<-ifelse(rbinom(length(which(dat$censpfs==0)),1,0.1)==1,1,0)

dat$timeos[which(dat$censpfs==0 & dat$censos==1)]<-abs(rnorm(length(which(dat$censpfs==0 & dat$censos==1)), 48, 5))

	
	
dat<-round(dat)
dat<-dat[, c(1:4,7:8, 12:13)]
write.table(dat, paste(getwd(),"/temp",a,b,"/temp",k,"/N2.txt", sep=""), sep="\t")
}}}