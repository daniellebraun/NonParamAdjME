#function to fit measurement error model
library("BayesMendel")
library(prodlim)
library("survival")
temp<-c('a', 'b','c', 'd', 'e','f','g', 'h','i')
temp2<-c(1,3)
for (tt in 1:2){
	
	proband<-temp2[tt]
	
for (ttt in 1:length(temp)){
fam1n1<- read.table(paste(getwd(),"/Simulation3N1/Simulation",tt,temp[ttt],"/N1.txt", sep=""), sep="\t",  header=T)

fam1n1<-fam1n1[which(fam1n1$ID!=proband),]

fam1n1<-fam1n1[which(fam1n1$Gender==0),]
fam1n1$AgeBreast<-ifelse(fam1n1$AgeBreast>110,110,fam1n1$AgeBreast)
fam1n1$AffectedBreastCen<-1
fam1n1$AffectedBreastCen[which(fam1n1$AffectedBreast==1) ]<-0 

fam1n10<-fam1n1[which(fam1n1$AffectedBreastR==0),]
fam1n11<-fam1n1[which(fam1n1$AffectedBreastR==1),]

psurv<-array(rep(0,110*110*4),dim=c(2,110,2,110))
phaz<-array(rep(0,110*110*4),dim=c(2,110,2,110))



#Calculations for distribution of failure times

#delta*==0
fit1a<-prodlim(Surv(AgeBreast,AffectedBreast)~AgeBreastR,data=fam1n10)
covVal<-sort(unique(fam1n10$AgeBreastR))
c1<-1
ss<-0
for (i in 1:110) {
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
for (i in 1:110){	
	if (ss[i]!=0){
		t1<-fit1a$time[c2:(c2+fit1a$size.strata[c3]-1)]
		surv1<-fit1a$surv[c2:(c2+fit1a$size.strata[c3]-1)]
		haz1<-fit1a$haz[c2:(c2+fit1a$size.strata[c3]-1)]	
		nS<-0
		nH<-0
		
		
		for (m in 1:110) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
		
		NS<-0
		NH<-0
			#three different cases
			#first events
			#first<-which(nS !=0)[1]
			first<-t1[1]
			if (first !=1){
				if (first<=10)
				{NS[1:(first-1)]<-nS[first]
				NH[1:(first-1)]<-nH[first]}
				if (first >10)
				{
					NS[(first-10):(first-1)]<-nS[first]
					NH[(first-10):(first-1)]<-nH[first]
					if (first!=11){
					NS[1:(first-11)]<-1
					NH[1:(first-11)]<-0
					}
				}	
			}
			
			#last events
			#last<-which(nS !=0)[length(which(nS!=0))]				
			last<-t1[length(t1)]
			if(last !=110){
			#NS[(last+1):110]<-nS[last]
			#NH[(last+1):110]<-nH[last]
				if (last>=100)
				{
					NS[(last+1):110]<-nS[last]
					NH[(last+1):110]<-nH[last]
				}
				if (last <100)
				{
					NS[(last+1):(last+10)]<-nS[last]
					NH[(last+1):(last+10)]<-nH[last]
					NS[(last+11):110]<-0
					NH[(last+11):110]<-0
				}
		
			}
			
			#middle events
		
			#l=length(which(nS!=0))-1
			l=length(t1)-1
			if (l==0) {NS[t1]<-nS[t1]
				NH[t1]<-nH[t1]
				}
			if (l!=0){
			for (k in 1:l)
			{
			#v1<-which(nS!=0)[k]
			#v2<-which(nS!=0)[k+1]
			v1<-t1[k]
			v2<-t1[k+1]
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
			NH[j]<-min(nH[v1], nH[v2])+abs(nH[v2]-nH[v1])*(j-v1)/d	
			} 
			NS[v1]<-nS[v1]
			NS[v2]<-nS[v2]
			NH[v1]<-nH[v1]
			NH[v2]<-nH[v2]
			
			}
			}
}
		psurv[1,,1,i]<-NS	
		phaz[1,,1,i]<-NH	

		c2<-c2+fit1a$size.strata[c3]
		c3<-c3+1
		}
}
			
first<-which(ss !=0)[1]
if (first !=1){	
	if (first<=10)
		{
			psurv[1,,1,(1:(first-1))]<-psurv[1,,1,first]
			phaz[1,,1,(1:(first-1))]<-phaz[1,,1,first]
		}
	if (first >10)
		{
				psurv[1,,1,((first-10):(first-1))]<-psurv[1,,1,first]
				phaz[1,,1,((first-10):(first-1))]<-phaz[1,,1,first]
				if (first !=11){
					psurv[1,,1,(1:(first-11))]<-rep(1,110)
					phaz[1,,1,(1:(first-11))]<-rep(0,110)				
				}	
		}

}			
last<-which(ss !=0)[length(which(ss!=0))]				

if (last!=110){	
	if (last>=100)
		{
			psurv[1,,1,((last+1):110)]<-psurv[1,,1,last]
			phaz[1,,1,((last+1):110)]<-phaz[1,,1,last]
		}
	if (last<100)
	{
		psurv[1,,1,((last+11):110)]<-rep(0,110)
		phaz[1,,1,((last+11):110)]<-rep(0,110)
		psurv[1,,1,((last+1):(last+10))]<-psurv[1,,1,last]
		phaz[1,,1,((last+1):(last+10))]<-phaz[1,,1,last]
	}			

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
			phaz[1,,1,j]<-pmin(phaz[1,,1,v1], phaz[1,,1,v2])+abs(phaz[1,,1,v2]-phaz[1,,1,v1])*(j-v1)/d	
			} 				
			}
			}	


#delta*==1
fit1b<-prodlim(Surv(AgeBreast,AffectedBreast)~AgeBreastR,data=fam1n11)
covVal<-sort(unique(fam1n11$AgeBreastR))
c1<-1
ss<-0
for (i in 1:110) {
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
for (i in 1:110){	
	if (ss[i]!=0){
		t1<-fit1b$time[c2:(c2+fit1b$size.strata[c3]-1)]
		surv1<-fit1b$surv[c2:(c2+fit1b$size.strata[c3]-1)]
		haz1<-fit1b$haz[c2:(c2+fit1b$size.strata[c3]-1)]	
		nS<-0
		nH<-0
		
		
		for (m in 1:110) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
		
		NS<-0
		NH<-0
			#three different cases
			#first events
			#first<-which(nS !=0)[1]
			first<-t1[1]
			if (first !=1){
				if (first<=10)
				{NS[1:(first-1)]<-nS[first]
				NH[1:(first-1)]<-nH[first]}
				if (first >10)
				{
					NS[(first-10):(first-1)]<-nS[first]
					NH[(first-10):(first-1)]<-nH[first]
					if (first!=11){
					NS[1:(first-11)]<-1
					NH[1:(first-11)]<-0
					}
				}	
			}
			
			#last events
			#last<-which(nS !=0)[length(which(nS!=0))]				
			last<-t1[length(t1)]
			if(last !=110){
			#NS[(last+1):110]<-nS[last]
			#NH[(last+1):110]<-nH[last]
				if (last>=100)
				{
					NS[(last+1):110]<-nS[last]
					NH[(last+1):110]<-nH[last]
				}
				if (last <100)
				{
					NS[(last+1):(last+10)]<-nS[last]
					NH[(last+1):(last+10)]<-nH[last]
					NS[(last+11):110]<-0
					NH[(last+11):110]<-0
				}
		
			}
			
			#middle events
		
			#l=length(which(nS!=0))-1
			l=length(t1)-1
				if (l==0) {NS[t1]<-nS[t1]
				NH[t1]<-nH[t1]
				}
			if (l!=0){
			for (k in 1:l)
			{
			#v1<-which(nS!=0)[k]
			#v2<-which(nS!=0)[k+1]
			v1<-t1[k]
			v2<-t1[k+1]
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
			NH[j]<-min(nH[v1], nH[v2])+abs(nH[v2]-nH[v1])*(j-v1)/d	
			} 
			NS[v1]<-nS[v1]
			NS[v2]<-nS[v2]
			NH[v1]<-nH[v1]
			NH[v2]<-nH[v2]
			
			}
			}
}

		psurv[1,,2,i]<-NS	
		phaz[1,,2,i]<-NH	

		c2<-c2+fit1b$size.strata[c3]
		c3<-c3+1
		}
}			


first<-which(ss !=0)[1]
if (first !=1){	
	if (first<=10)
		{
			psurv[1,,2,(1:(first-1))]<-psurv[1,,2,first]
			phaz[1,,2,(1:(first-1))]<-phaz[1,,2,first]
		}
	if (first >10)
		{
				psurv[1,,2,((first-10):(first-1))]<-psurv[1,,2,first]
				phaz[1,,2,((first-10):(first-1))]<-phaz[1,,2,first]
				if (first !=11){
					psurv[1,,2,(1:(first-11))]<-rep(1,110)
					phaz[1,,2,(1:(first-11))]<-rep(0,110)				
				}	
		}

}	

last<-which(ss !=0)[length(which(ss!=0))]
if (last!=110){	
	if (last>=100)
		{
			psurv[1,,2,((last+1):110)]<-psurv[1,,2,last]
			phaz[1,,2,((last+1):110)]<-phaz[1,,2,last]
		}
	if (last<100)
	{
		psurv[1,,2,((last+11):110)]<-rep(0,110)
		phaz[1,,2,((last+11):110)]<-rep(0,110)
		psurv[1,,2,((last+1):(last+10))]<-psurv[1,,2,last]
		phaz[1,,2,((last+1):(last+10))]<-phaz[1,,2,last]
	}			

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
			phaz[1,,2,j]<-pmin(phaz[1,,2,v1], phaz[1,,2,v2])+abs(phaz[1,,2,v2]-phaz[1,,2,v1])*(j-v1)/d	
			} 				
			}
			}		

	
#Calculations for distribution of censoring times

#delta*==0
fit1c<-prodlim(Surv(AgeBreast,AffectedBreastCen)~AgeBreastR,data=fam1n10)
covVal<-sort(unique(fam1n10$AgeBreastR))
c1<-1
ss<-0
for (i in 1:110) {
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
for (i in 1:110){	
	if (ss[i]!=0){
		t1<-fit1c$time[c2:(c2+fit1c$size.strata[c3]-1)]
		surv1<-fit1c$surv[c2:(c2+fit1c$size.strata[c3]-1)]
		haz1<-fit1c$haz[c2:(c2+fit1c$size.strata[c3]-1)]	
		nS<-0
		nH<-0
		
		
		for (m in 1:110) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
		
		NS<-0
		NH<-0
			#three different cases
			#first events
			#first<-which(nS !=0)[1]
			first<-t1[1]
			if (first !=1){
				if (first<=10)
				{NS[1:(first-1)]<-nS[first]
				NH[1:(first-1)]<-nH[first]}
				if (first >10)
				{
					NS[(first-10):(first-1)]<-nS[first]
					NH[(first-10):(first-1)]<-nH[first]
					if (first!=11){
					NS[1:(first-11)]<-1
					NH[1:(first-11)]<-0
					}
				}	
			}
			
			#last events
			#last<-which(nS !=0)[length(which(nS!=0))]				
			last<-t1[length(t1)]
			if(last !=110){
			#NS[(last+1):110]<-nS[last]
			#NH[(last+1):110]<-nH[last]
				if (last>=100)
				{
					NS[(last+1):110]<-nS[last]
					NH[(last+1):110]<-nH[last]
				}
				if (last <100)
				{
					NS[(last+1):(last+10)]<-nS[last]
					NH[(last+1):(last+10)]<-nH[last]
					NS[(last+11):110]<-0
					NH[(last+11):110]<-0
				}
		
			}
			
			#middle events
		
			#l=length(which(nS!=0))-1
			l=length(t1)-1
				if (l==0) {NS[t1]<-nS[t1]
				NH[t1]<-nH[t1]
				}
			if (l!=0){
			for (k in 1:l)
			{
			#v1<-which(nS!=0)[k]
			#v2<-which(nS!=0)[k+1]
			v1<-t1[k]
			v2<-t1[k+1]
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
			NH[j]<-min(nH[v1], nH[v2])+abs(nH[v2]-nH[v1])*(j-v1)/d	
			} 
			NS[v1]<-nS[v1]
			NS[v2]<-nS[v2]
			NH[v1]<-nH[v1]
			NH[v2]<-nH[v2]
			
			}
			}}
		psurv[2,,1,i]<-NS	
		phaz[2,,1,i]<-NH	
		c2<-c2+fit1c$size.strata[c3]
		c3<-c3+1
		}
}		

first<-which(ss !=0)[1]
if (first !=1){	
	if (first<=10)
		{
			psurv[2,,1,(1:(first-1))]<-psurv[2,,1,first]
			phaz[2,,1,(1:(first-1))]<-phaz[2,,1,first]
		}
	if (first >10)
		{
				psurv[2,,1,((first-10):(first-1))]<-psurv[2,,1,first]
				phaz[2,,1,((first-10):(first-1))]<-phaz[2,,1,first]
				if (first !=11){
					psurv[2,,1,(1:(first-11))]<-rep(1,110)
					phaz[2,,1,(1:(first-11))]<-rep(0,110)				
				}	
		}

}	
last<-which(ss !=0)[length(which(ss!=0))]
if (last!=110){	
	if (last>=100)
		{
			psurv[2,,1,((last+1):110)]<-psurv[2,,1,last]
			phaz[2,,1,((last+1):110)]<-phaz[2,,1,last]
		}
	if (last<100)
	{
		psurv[2,,1,((last+11):110)]<-rep(0,110)
		phaz[2,,1,((last+11):110)]<-rep(0,110)
		psurv[2,,1,((last+1):(last+10))]<-psurv[2,,1,last]
		phaz[2,,1,((last+1):(last+10))]<-phaz[2,,1,last]
	}			

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
			phaz[2,,1,j]<-pmin(phaz[2,,1,v1], phaz[2,,1,v2])+abs(phaz[2,,1,v2]-phaz[2,,1,v1])*(j-v1)/d	
			} 				
			}
			}		


#delta*==1
fit1d<-prodlim(Surv(AgeBreast,AffectedBreastCen)~AgeBreastR,data=fam1n11)
covVal<-sort(unique(fam1n11$AgeBreastR))
c1<-1
ss<-0
for (i in 1:110) {
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
for (i in 1:110){	
	if (ss[i]!=0){
		t1<-fit1d$time[c2:(c2+fit1d$size.strata[c3]-1)]
		surv1<-fit1d$surv[c2:(c2+fit1d$size.strata[c3]-1)]
		haz1<-fit1d$haz[c2:(c2+fit1d$size.strata[c3]-1)]	
		nS<-0
		nH<-0
		
		
		for (m in 1:110) {
		if (length(which(t1==m)>0) ) {
			nS[m]<-surv1[which(t1==m)]
			nH[m]<-haz1[which(t1==m)]
		}
		else { 
			nS[m]<-0
			nH[m]<-0
			}
		}
		
		
		NS<-0
		NH<-0
			#three different cases
			#first events
			#first<-which(nS !=0)[1]
			first<-t1[1]
			if (first !=1){
				if (first<=10)
				{NS[1:(first-1)]<-nS[first]
				NH[1:(first-1)]<-nH[first]}
				if (first >10)
				{
					NS[(first-10):(first-1)]<-nS[first]
					NH[(first-10):(first-1)]<-nH[first]
					if (first!=11){
					NS[1:(first-11)]<-1
					NH[1:(first-11)]<-0
					}
				}	
			}
			
			#last events
			#last<-which(nS !=0)[length(which(nS!=0))]				
			last<-t1[length(t1)]
			if(last !=110){
			#NS[(last+1):110]<-nS[last]
			#NH[(last+1):110]<-nH[last]
				if (last>=100)
				{
					NS[(last+1):110]<-nS[last]
					NH[(last+1):110]<-nH[last]
				}
				if (last <100)
				{
					NS[(last+1):(last+10)]<-nS[last]
					NH[(last+1):(last+10)]<-nH[last]
					NS[(last+11):110]<-0
					NH[(last+11):110]<-0
				}
		
			}
			
			#middle events
		
			#l=length(which(nS!=0))-1
			l=length(t1)-1
				if (l==0) {NS[t1]<-nS[t1]
				NH[t1]<-nH[t1]
				}
			if (l!=0){
			for (k in 1:l)
			{
			#v1<-which(nS!=0)[k]
			#v2<-which(nS!=0)[k+1]
			v1<-t1[k]
			v2<-t1[k+1]
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
			NH[j]<-min(nH[v1], nH[v2])+abs(nH[v2]-nH[v1])*(j-v1)/d	
			} 
			NS[v1]<-nS[v1]
			NS[v2]<-nS[v2]
			NH[v1]<-nH[v1]
			NH[v2]<-nH[v2]
			
			}
			}}

		psurv[2,,2,i]<-NS	
		phaz[2,,2,i]<-NH	
		c2<-c2+fit1d$size.strata[c3]
		c3<-c3+1
		}
}			

first<-which(ss !=0)[1]
if (first !=1){	
	if (first<=10)
		{
			psurv[2,,2,(1:(first-1))]<-psurv[2,,2,first]
			phaz[2,,2,(1:(first-1))]<-phaz[2,,2,first]
		}
	if (first >10)
		{
				psurv[2,,2,((first-10):(first-1))]<-psurv[2,,2,first]
				phaz[2,,2,((first-10):(first-1))]<-phaz[2,,2,first]
				if (first !=11){
					psurv[2,,2,(1:(first-11))]<-rep(1,110)
					phaz[2,,2,(1:(first-11))]<-rep(0,110)				
				}	
		}

}	

		
last<-which(ss !=0)[length(which(ss!=0))]				
if (last!=110){	
	if (last>=100)
		{
			psurv[2,,2,((last+1):110)]<-psurv[2,,2,last]
			phaz[2,,2,((last+1):110)]<-phaz[2,,2,last]
		}
	if (last<100)
	{
		psurv[2,,2,((last+11):110)]<-rep(0,110)
		phaz[2,,2,((last+11):110)]<-rep(0,110)
		psurv[2,,2,((last+1):(last+10))]<-psurv[2,,2,last]
		phaz[2,,2,((last+1):(last+10))]<-phaz[2,,2,last]
	}			

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
			phaz[2,,2,j]<-pmin(phaz[2,,2,v1], phaz[2,,2,v2], byrow=T)+abs(phaz[2,,2,v2]-phaz[2,,2,v1])*(j-v1)/d	
			} 				
			}
			}		




save(phaz, psurv, file = paste(getwd(),"/Simulation3/Simulation",tt,temp[ttt],"/SKM_N1.txt", sep=""))
}}