set.seed(435829)
for (a in 1:10){
	for (b in 1:10){
		for (k in 1:100){
			
		
		sens<-a/10
		spec<-b/10
		
#n1
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
	
	
dat<-round(dat)
dat<-dat[, c(1:4,7:8)]
write.table(dat, paste(getwd(),"/temp",a,b,"/temp",k,"/N1.txt", sep=""), sep="\t")
		
		
	}
}

}