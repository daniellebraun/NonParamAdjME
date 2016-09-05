library("prodlim")
library(survival)
for (a in 1:10){
	for (b in 1:10){
		for (zz in 1:100){
dat<-read.table(paste(getwd(),"/temp",a,b,"/temp",zz,"/N1.txt", sep=""), sep="\t",  header=T)

#fit measurement error model

#weigh model

dat$censpfsC<-1
dat$censpfsC[which(dat$censpfs==1) ]<-0 

dat0<-dat[which(dat$ercenspfs==0),]
dat1<-dat[which(dat$ercenspfs==1),]

psurv<-array(rep(0,25*25*4),dim=c(2,25,2,25))
phaz<-array(rep(0,25*25*4),dim=c(2,25,2,25))

#Calculations for distribution of failure times

#delta*==0
fit1a<-prodlim(Surv(timepfs,censpfs)~ertimepfs,data=dat0)
covVal<-sort(unique(dat0$ertimepfs))
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
		
		
		for (m in 1:25) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
		
		NS<-rep(0,25)
		NH<-rep(0,25)
		
		if (length(which(nS==0))==25) {
			tf<-t1[1]
			NS[1:tf]<-1
			NS[(tf+1):25]<-0
			  NH[1:tf]<-0
			  NH[(tf+1):25]<-0		
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
			if(last !=25){
			NS[(last+1):25]<-nS[last]
			NH[(last+1):25]<-nH[last]
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
		psurv[1,,1,i]<-NS	
		phaz[1,,1,i]<-NH	

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
				psurv[1,,1,i]<-psurv[1,,1,first]
				phaz[1,,1,i]<-phaz[1,,1,first]
			}
			
			last<-which(ss !=0)[length(which(ss!=0))]				
			if (i > last){
				psurv[1,,1,i]<-psurv[1,,1,last]
				phaz[1,,1,i]<-phaz[1,,1,last]
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
			psurv[1,,1,j]<-psurv[1,,1,v1]-(psurv[1,,1,v1]-psurv[1,,1,v1])*(j-v1)/d	
			phaz[1,,1,j]<-phaz[1,,1,v1]+(phaz[1,,1,v2]-phaz[1,,1,v1])*(j-v1)/d	
			} 				
			}
			}}	
}}		

#delta*==1
fit1b<-prodlim(Surv(timepfs,censpfs)~ertimepfs,data=dat1)
covVal<-sort(unique(dat1$ertimepfs))
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
		
		for (m in 1:25) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
		NS<-rep(0,25)
		NH<-rep(0,25)
		
		if (length(which(nS==0))==25) {
			tf<-t1[1]
			NS[1:tf]<-1
			NS[(tf+1):25]<-0
			  NH[1:tf]<-0
			  NH[(tf+1):25]<-0		
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
			if(last !=25){
			NS[(last+1):25]<-nS[last]
			NH[(last+1):25]<-nH[last]
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
		psurv[1,,2,i]<-NS	
		phaz[1,,2,i]<-NH	

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
				psurv[1,,2,i]<-psurv[1,,2,first]
				phaz[1,,2,i]<-phaz[1,,2,first]
				}
			
			last<-which(ss !=0)[length(which(ss!=0))]				
			if (i > last){
				psurv[1,,2,i]<-psurv[1,,2,last]
				phaz[1,,2,i]<-phaz[1,,2,last]
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
			psurv[1,,2,j]<-psurv[1,,2,v1]-(psurv[1,,2,v1]-psurv[1,,2,v2])*(j-v1)/d	
			phaz[1,,2,j]<-phaz[1,,2,v1]+(phaz[1,,2,v2]-phaz[1,,2,v1])*(j-v1)/d	
			} 				
			}
			}}		

}	}	

	
#Calculations for distribution of censoring times

#delta*==0
fit1c<-prodlim(Surv(timepfs,censpfsC)~ertimepfs,data=dat0)
covVal<-sort(unique(dat0$ertimepfs))
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
		
		
		for (m in 1:25) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
	NS<-rep(0,25)
		NH<-rep(0,25)
		
		if (length(which(nS==0))==25) {
			tf<-t1[1]
			NS[1:tf]<-1
			NS[(tf+1):25]<-0
			  NH[1:tf]<-0
			  NH[(tf+1):25]<-0		
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
			if (last !=25){
			NS[(last+1):25]<-nS[last]
			NH[(last+1):25]<-nH[last]
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
		psurv[2,,1,i]<-NS	
		phaz[2,,1,i]<-NH	
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
				psurv[2,,1,i]<-psurv[2,,1,last]
				phaz[2,,1,i]<-phaz[2,,1,last]
			}
			
			last<-which(ss !=0)[length(which(ss!=0))]				
			if (i > last){
				psurv[2,,1,i]<-psurv[2,,1,last]
				phaz[2,,1,i]<-phaz[2,,1,last]
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
			psurv[2,,1,j]<-psurv[2,,1,v1]-(psurv[2,,1,v1]-psurv[2,,1,v2])*(j-v1)/d	
			phaz[2,,1,j]<-phaz[2,,1,v1]+(phaz[2,,1,v2]-phaz[2,,1,v1])*(j-v1)/d	
			} 				
			}
			}}		
	
}}		


#delta*==1
fit1d<-prodlim(Surv(timepfs,censpfsC)~ertimepfs,data=dat1)
covVal<-sort(unique(dat1$ertimepfs))
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
		
		
		for (m in 1:25) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
		NS<-rep(0,25)
		NH<-rep(0,25)
		
		if (length(which(nS==0))==25) {
			tf<-t1[1]
			NS[1:tf]<-1
			NS[(tf+1):25]<-0
			  NH[1:tf]<-0
			  NH[(tf+1):25]<-0		
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
			if (last !=25){
			NS[(last+1):25]<-nS[last]
			NH[(last+1):25]<-nH[last]
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
		psurv[2,,2,i]<-NS	
		phaz[2,,2,i]<-NH	
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
				psurv[2,,2,i]<-psurv[2,,2,first]
				phaz[2,,2,i]<-phaz[2,,2,first]
			}
			
			last<-which(ss !=0)[length(which(ss!=0))]				
			if (i > last){
				psurv[2,,2,i]<-psurv[2,,2,last]
				phaz[2,,2,i]<-phaz[2,,2,last]
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
			psurv[2,,2,j]<-psurv[2,,2,v1]-(psurv[2,,2,v1]-psurv[2,,2,v2])*(j-v1)/d	
			phaz[2,,2,j]<-phaz[2,,2,v1]+(phaz[2,,2,v2]-phaz[2,,2,v1])*(j-v1)/d	
			} 				
			}
			}}		
	
}	}	

save(phaz, psurv, file = paste(getwd(),"/temp",a,b,"/temp",zz,"/SKM_N1.txt", sep=""))
}}
}
