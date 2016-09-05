library("prodlim")
library(survival)
for (a in 1:10){
	for (b in 1:10){
		for (zz in 1:100){
dat<-read.table(paste(getwd(),"/temp",a,b,"/temp",zz,"/N3.txt", sep=""), sep="\t",  header=T)

#fit prediction model for OS based on pfs
dat<-round(dat)
dat$timeosN<-ifelse(dat$timeos>100,100,dat$timeos)
dat$censosN<-ifelse(dat$timeos>100,0,dat$censos)
dat$timeos<-dat$timeosN
dat$censos<-dat$censosN
dat$censosC<-1
dat$censosC[which(dat$censos==1) ]<-0 


dat0<-dat[which(dat$censpfs==0),]
dat1<-dat[which(dat$censpfs==1),]

psurvS<-array(rep(0,100*25*4),dim=c(2,100,2,25))
phazS<-array(rep(0,100*25*4),dim=c(2,100,2,25))

#Calculations for distribution of failure times

#delta*==0
fit1a<-prodlim(Surv(timeos,censos)~timepfs,data=dat0)
covVal<-sort(unique(dat0$timepfs))
c1<-1
ss<-0
for (i in 1:25) {
  if (length(which(covVal==i)>0)) {
  	ss[i]<-fit1a$size.strata[c1]
  	c1<-c1+1
  	}
  	else {
  		ss[i]<-0}
  	}
  		
c2<-1
c3<-1
d<-0
for (i in 1:25){	
	if (ss[i]!=0){
		t1<-fit1a$time[c2:(c2+fit1a$size.strata[c3]-1)]
		surv1<-fit1a$surv[c2:(c2+fit1a$size.strata[c3]-1)]
		haz1<-fit1a$haz[c2:(c2+fit1a$size.strata[c3]-1)]	
		nS<-0
		nH<-0
		
		
		for (m in 1:100) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
		
		NS<-rep(0,100)
		NH<-rep(0,100)
		
		if (length(which(nS==0))==100) {
			tf<-t1[1]
			NS[1:tf]<-1
			NS[(tf+1):100]<-0
			  NH[1:tf]<-0
			  NH[(tf+1):100]<-0		
			}
			else {
			#three different cases
			#first events
			first<-which(nS !=0)[1]
			if (first !=1){
			NS[1:first-1]<-1
			NH[1:first-1]<-0
			}
			
			#last events
			
			last<-which(nS !=0)[length(which(nS!=0))]				
			if(last !=100){
			NS[(last+1):100]<-nS[last]
			NH[(last+1):100]<-nH[last]
			}
			
			#middle events
		
			l=length(which(nS!=0))-1
			if (l>=1){
			for (k in 1:l)
			{
			v1<-which(nS!=0)[k]
			v2<-which(nS!=0)[k+1]
			d<-v2-v1
			if (d==1){
				NS[v1]<-nS[v1]
				NS[v2]<-nS[v2]
				NH[v1]<-nH[v1]
				NH[v2]<-nH[v2]
			}
			
			if (d !=1) {
				beg<-v1+1
				end<-v2-1
			for (j in beg:end){
			NS[j]<-nS[v1]-(nS[v1]-nS[v2])*(j-v1)/d	
			NH[j]<-nH[v1]+(nH[v2]-nH[v1])*(j-v1)/d	
			} 
			NS[v1]<-nS[v1]
			NS[v2]<-nS[v2]
			NH[v1]<-nH[v1]
			NH[v2]<-nH[v2]
			
			}
			}
			}
			if (l==0){
				NS[which(nS!=0)]<-nS[which(nS!=0)]
				NH[which(nS!=0)]<-nH[which(nS!=0)]
			}
	}
		psurvS[1,,1,i]<-NS	
		phazS[1,,1,i]<-NH	

		c2<-c2+fit1a$size.strata[c3]
		c3<-c3+1
		}
}


for (i in 1:25){
	   if( ss[i]==0){
		#t* missing between 1-25
			if (i<=25){			
			
			first<-which(ss !=0)[1]
			if (i < first){
				psurvS[1,,1,i]<-psurvS[1,,1,first]
				phazS[1,,1,i]<-phazS[1,,1,first]
			}
			
			last<-which(ss !=0)[length(which(ss!=0))]				
			if (i > last){
				psurvS[1,,1,i]<-psurvS[1,,1,last]
				phazS[1,,1,i]<-phazS[1,,1,last]
			}

			l=length(which(ss!=0))-1
			for (k in 1:l)
			{
			v1<-which(ss!=0)[k]
			v2<-which(ss!=0)[k+1]
			d<-v2-v1
			if (d !=1) {
				beg<-v1+1
				end<-v2-1
				
			for (j in beg:end){
			psurvS[1,,1,j]<-psurvS[1,,1,v1]-(psurvS[1,,1,v1]-psurvS[1,,1,v1])*(j-v1)/d	
			phazS[1,,1,j]<-phazS[1,,1,v1]+(phazS[1,,1,v2]-phazS[1,,1,v1])*(j-v1)/d	
			} 				
			}
			}}	
}}		

#delta*==1
fit1b<-prodlim(Surv(timeos,censos)~timepfs,data=dat1)
covVal<-sort(unique(dat1$timepfs))
c1<-1
ss<-0
for (i in 1:25) {
  if (length(which(covVal==i)>0)) {
  	ss[i]<-fit1b$size.strata[c1]
  	c1<-c1+1
  	}
  	else {
  		ss[i]<-0}
  	}
  		
c2<-1
c3<-1
d<-0
for (i in 1:25){
	if (ss[i]!=0) {
		t1<-fit1b$time[c2:(c2+fit1b$size.strata[c3]-1)]
		surv1<-fit1b$surv[c2:(c2+fit1b$size.strata[c3]-1)]	
		haz1<-fit1b$haz[c2:(c2+fit1b$size.strata[c3]-1)]	
		nS<-0
		nH<-0
		
		for (m in 1:100) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
		NS<-rep(0,100)
		NH<-rep(0,100)
		
		if (length(which(nS==0))==100) {
			tf<-t1[1]
			NS[1:tf]<-1
			NS[(tf+1):100]<-0
			  NH[1:tf]<-0
			  NH[(tf+1):100]<-0		
			}
		else {			
		
			#three different cases
			#first events
		
			
			first<-which(nS !=0)[1]
			if (first !=1){
			NS[1:first-1]<-1
			NH[1:first-1]<-0
			}
			
			#last events
			last<-which(nS !=0)[length(which(nS!=0))]				
			if(last !=100){
			NS[(last+1):100]<-nS[last]
			NH[(last+1):100]<-nH[last]
			}
			
			#middle events
		
			l=length(which(nS!=0))-1
			if (l>=1){
			for (k in 1:l)
			{
			v1<-which(nS!=0)[k]
			v2<-which(nS!=0)[k+1]
			d<-v2-v1
			if (d==1){
				NS[v1]<-nS[v1]
				NS[v2]<-nS[v2]
				NH[v1]<-nH[v1]
				NH[v2]<-nH[v2]
			}
			
			if (d !=1) {
				beg<-v1+1
				end<-v2-1
			for (j in beg:end){
			NS[j]<-nS[v1]-(nS[v1]-nS[v2])*(j-v1)/d	
			NH[j]<-nH[v1]+(nH[v2]-nH[v1])*(j-v1)/d	
			} 
			NS[v1]<-nS[v1]
			NS[v2]<-nS[v2]
			NH[v1]<-nH[v1]
			NH[v2]<-nH[v2]
			
			}
			}
			}
		if (l==0){
				NS[which(nS!=0)]<-nS[which(nS!=0)]
				NH[which(nS!=0)]<-nH[which(nS!=0)]
			}

}
		psurvS[1,,2,i]<-NS	
		phazS[1,,2,i]<-NH	

		c2<-c2+fit1b$size.strata[c3]
		c3<-c3+1
		}
}			

for (i in 1:25){
		if( ss[i]==0){
			#t* missing between 1-25
			if (i<=25){
				
			first<-which(ss !=0)[1]
			if (i < first){
				psurvS[1,,2,i]<-psurvS[1,,2,first]
				phazS[1,,2,i]<-phazS[1,,2,first]
				}
			
			last<-which(ss !=0)[length(which(ss!=0))]				
			if (i > last){
				psurvS[1,,2,i]<-psurvS[1,,2,last]
				phazS[1,,2,i]<-phazS[1,,2,last]
			}

			l=length(which(ss!=0))-1
			for (k in 1:l)
			{
			v1<-which(ss!=0)[k]
			v2<-which(ss!=0)[k+1]
			d<-v2-v1
			if (d !=1) {
				beg<-v1+1
				end<-v2-1
			for (j in beg:end){
			psurvS[1,,2,j]<-psurvS[1,,2,v1]-(psurvS[1,,2,v1]-psurvS[1,,2,v2])*(j-v1)/d	
			phazS[1,,2,j]<-phazS[1,,2,v1]+(phazS[1,,2,v2]-phazS[1,,2,v1])*(j-v1)/d	
			} 				
			}
			}}		

}	}	

	
#Calculations for distribution of censoring times

#delta*==0
fit1c<-prodlim(Surv(timeos,censosC)~timepfs,data=dat0)
covVal<-sort(unique(dat0$timepfs))
c1<-1
ss<-0
for (i in 1:25) {
  if (length(which(covVal==i)>0)) {
  	ss[i]<-fit1c$size.strata[c1]
  	c1<-c1+1
  	}
  	else {
  		ss[i]<-0}
  	}
  		
c2<-1
c3<-1
d<-0
for (i in 1:25){
	if (ss[i]!=0) {

		t1<-fit1c$time[c2:(c2+fit1c$size.strata[c3]-1)]
		surv1<-fit1c$surv[c2:(c2+fit1c$size.strata[c3]-1)]
		haz1<-fit1c$haz[c2:(c2+fit1c$size.strata[c3]-1)]	
		nS<-0
		nH<-0
		
		
		for (m in 1:100) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
	NS<-rep(0,100)
		NH<-rep(0,100)
		
		if (length(which(nS==0))==100) {
			tf<-t1[1]
			NS[1:tf]<-1
			NS[(tf+1):100]<-0
			  NH[1:tf]<-0
			  NH[(tf+1):100]<-0		
			}
		else {	
			
			#three different cases
			#first events
			first<-which(nS !=0)[1]
			if (first !=1){
			NS[1:first-1]<-1
			NH[1:first-1]<-0
			}
			
			#last events
			last<-which(nS !=0)[length(which(nS!=0))]				
			if (last !=100){
			NS[(last+1):100]<-nS[last]
			NH[(last+1):100]<-nH[last]
			}
			
			#middle events
		
			l=length(which(nS!=0))-1
			
			if (l>=1){
			for (k in 1:l)
			{
			v1<-which(nS!=0)[k]
			v2<-which(nS!=0)[k+1]
			d<-v2-v1
			if (d==1){
				NS[v1]<-nS[v1]
				NS[v2]<-nS[v2]
				NH[v1]<-nH[v1]
				NH[v2]<-nH[v2]
			}
			
			if (d !=1) {
				beg<-v1+1
				end<-v2-1
			for (j in beg:end){
			NS[j]<-nS[v1]-(nS[v1]-nS[v2])*(j-v1)/d	
			NH[j]<-nH[v1]+(nH[v2]-nH[v1])*(j-v1)/d	
			} 
			NS[v1]<-nS[v1]
			NS[v2]<-nS[v2]
			NH[v1]<-nH[v1]
			NH[v2]<-nH[v2]
			
			}
			}
			}
			if (l==0){
				NS[which(nS!=0)]<-nS[which(nS!=0)]
				NH[which(nS!=0)]<-nH[which(nS!=0)]
			}
}
		psurvS[2,,1,i]<-NS	
		phazS[2,,1,i]<-NH	
		c2<-c2+fit1c$size.strata[c3]
		c3<-c3+1
		}
}		

for (i in 1:25){
		if( ss[i]==0){
		
			#t* missing between 1-100
			if (i<=25){
				
			first<-which(ss !=0)[1]
			if (i < first){
				psurvS[2,,1,i]<-psurvS[2,,1,last]
				phazS[2,,1,i]<-phazS[2,,1,last]
			}
			
			last<-which(ss !=0)[length(which(ss!=0))]				
			if (i > last){
				psurvS[2,,1,i]<-psurvS[2,,1,last]
				phazS[2,,1,i]<-phazS[2,,1,last]
			}

			l=length(which(ss!=0))-1
			for (k in 1:l)
			{
			v1<-which(ss!=0)[k]
			v2<-which(ss!=0)[k+1]
			d<-v2-v1
			if (d !=1) {
				beg<-v1+1
				end<-v2-1
			for (j in beg:end){
			psurvS[2,,1,j]<-psurvS[2,,1,v1]-(psurvS[2,,1,v1]-psurvS[2,,1,v2])*(j-v1)/d	
			phazS[2,,1,j]<-phazS[2,,1,v1]+(phazS[2,,1,v2]-phazS[2,,1,v1])*(j-v1)/d	
			} 				
			}
			}}		
	
}}		


#delta*==1
fit1d<-prodlim(Surv(timeos,censosC)~timepfs,data=dat1)
covVal<-sort(unique(dat1$timepfs))
c1<-1
ss<-0
for (i in 1:25) {
  if (length(which(covVal==i)>0)) {
  	ss[i]<-fit1d$size.strata[c1]
  	c1<-c1+1
  	}
  	else {
  		ss[i]<-0}
  	}
  		
c2<-1
c3<-1
d<-0
for (i in 1:25){
	if (ss[i]!=0) {
		
		t1<-fit1d$time[c2:(c2+fit1d$size.strata[c3]-1)]
		surv1<-fit1d$surv[c2:(c2+fit1d$size.strata[c3]-1)]
		haz1<-fit1d$haz[c2:(c2+fit1d$size.strata[c3]-1)]
			
		nS<-0
		nH<-0
		
		
		for (m in 1:100) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
		NS<-rep(0,100)
		NH<-rep(0,100)
		
		if (length(which(nS==0))==100) {
			tf<-t1[1]
			NS[1:tf]<-1
			NS[(tf+1):100]<-0
			  NH[1:tf]<-0
			  NH[(tf+1):100]<-0		
			}
		else {
			#three different cases
			#first events
			first<-which(nS !=0)[1]
			if (first !=1){
			NS[1:first-1]<-1
			NH[1:first-1]<-0
			}
			
			#last events
			last<-which(nS !=0)[length(which(nS!=0))]				
			if (last !=100){
			NS[(last+1):100]<-nS[last]
			NH[(last+1):100]<-nH[last]
			}
			
			#middle events
		
			l=length(which(nS!=0))-1
			if (l>=1){
			for (k in 1:l)
			{
			v1<-which(nS!=0)[k]
			v2<-which(nS!=0)[k+1]
			d<-v2-v1
			if (d==1){
				NS[v1]<-nS[v1]
				NS[v2]<-nS[v2]
				NH[v1]<-nH[v1]
				NH[v2]<-nH[v2]
			}
			
			if (d !=1) {
				beg<-v1+1
				end<-v2-1
			for (j in beg:end){
			NS[j]<-nS[v1]-(nS[v1]-nS[v2])*(j-v1)/d	
			NH[j]<-nH[v1]+(nH[v2]-nH[v1])*(j-v1)/d	
			} 
			NS[v1]<-nS[v1]
			NS[v2]<-nS[v2]
			NH[v1]<-nH[v1]
			NH[v2]<-nH[v2]
			
			}
			}
			}
			if (l==0){
				NS[which(nS!=0)]<-nS[which(nS!=0)]
				NH[which(nS!=0)]<-nH[which(nS!=0)]
			}

		}
		psurvS[2,,2,i]<-NS	
		phazS[2,,2,i]<-NH	
		c2<-c2+fit1d$size.strata[c3]
		c3<-c3+1
		}
}			

for (i in 1:25){
		if( ss[i]==0){
			#t* missing between 10-100
			if (i<=25){
				
			first<-which(ss !=0)[1]
			if (i < first){
				psurvS[2,,2,i]<-psurvS[2,,2,first]
				phazS[2,,2,i]<-phazS[2,,2,first]
			}
			
			last<-which(ss !=0)[length(which(ss!=0))]				
			if (i > last){
				psurvS[2,,2,i]<-psurvS[2,,2,last]
				phazS[2,,2,i]<-phazS[2,,2,last]
			}

			l=length(which(ss!=0))-1
			for (k in 1:l)
			{
			v1<-which(ss!=0)[k]
			v2<-which(ss!=0)[k+1]
			d<-v2-v1
			if (d !=1) {
				beg<-v1+1
				end<-v2-1
			for (j in beg:end){
			psurvS[2,,2,j]<-psurvS[2,,2,v1]-(psurvS[2,,2,v1]-psurvS[2,,2,v2])*(j-v1)/d	
			phazS[2,,2,j]<-phazS[2,,2,v1]+(phazS[2,,2,v2]-phazS[2,,2,v1])*(j-v1)/d	
			} 				
			}
			}}		
	
}	}	




save(phazS, psurvS, file = paste(getwd(),"/temp",a,b,"/temp",zz,"/SKM_N3.txt", sep=""))
}}}

