#generate Menedlian families. 
library("BayesMendel")

data(BRCApenet.metaDSL.2008,death.othercauses, compriskSurv)

#calculating distribution of censoring time based on hazards of death from other causes from brcapro.
surDe<-exp(-apply(death.othercauses, 2, cumsum))
dennoF<-surDe[,2]*death.othercauses$femaleBC
dennoM<-surDe[,1]*death.othercauses$maleBC

dennoF<-append(dennoF, 1-sum(dennoF))
dennoM<-append(dennoM, 1-sum(dennoM))

#female BRCA1 penetrance
brca1.pen = BRCApenet.metaDSL.2008$fFX[,1:2]  #no homozygous mutation carriers
brca1.penF<-matrix(,111,2)
brca1.penF[,1]<-append(brca1.pen[,1], 1-sum(brca1.pen[,1]))
brca1.penF[,2]<-append(brca1.pen[,2], 1-sum(brca1.pen[,2]))
brca1.nc= brca1.penF[,1]
brca1.c= brca1.penF[,2]

## male BRCA1 penetrance
brca1.pen.m = BRCApenet.metaDSL.2008$fMX[,1:2]  #no homozygous mutation carriers   
brca1.penM<-matrix(,111,2)
brca1.penM[,1]<-append(brca1.pen.m[,1], 1-sum(brca1.pen.m[,1]))
brca1.penM[,2]<-append(brca1.pen.m[,2], 1-sum(brca1.pen.m[,2]))
brca1.nc.m= brca1.penM[,1]
brca1.c.m= brca1.penM[,2]

fam.gen = function(x,n,m1){
#
# Generating the families

# x: q, allele frequency
# n: target population
# m1: 
# m2: 
##
        q = x

## generating genotypes for parents
##
         geno.f = cbind(rbinom(n,size=1, prob=q), rbinom(n,size=1,prob=q))
         geno.m = cbind(rbinom(n,size=1, prob=q), rbinom(n,size=1,prob=q))

## choosing which allele to pass on for m1 female offspring
##
        g.o = NULL
        for (i in 1:m1){
            # which allele to pass to daughter
            indx.f0 = rbinom(n,size=1,prob=0.5)
            indx.m0 = rbinom(n,size=1,prob=0.5)
            g1 = geno.f * cbind(indx.f0,(1-indx.f0))
            g2 = geno.m * cbind(indx.m0,(1-indx.m0))
            geno.f2=g1+g2
            geno.f2= ifelse(geno.f2>1,1,geno.f2)
            
            # which allele to pass to her husband
            geno.m2=cbind(rbinom(n,size=1, prob=q), rbinom(n,size=1,prob=q))
            

            
            g.o = cbind(g.o,apply(geno.f2,1,sum),apply(geno.m2,1,sum))
        }
        g.o = ifelse(g.o>1,1,g.o)

## parents genotype
        g.f = geno.f[,1]+geno.f[,2]; g.f = ifelse(g.f>1,1,g.f)
        g.m = geno.m[,1]+geno.m[,2]; g.m = ifelse(g.m>1,1,g.m)


## all genotype
        g.all=cbind(g.f, g.m, g.o)
        

tot=2+m1*2
## generating failure time
##

        time=cen=ind=obs=vector(); fail = matrix(0,n,tot)
        
        for ( i in 1:n){  # including parents and offspring
         
			## Sampling censoring/failure time for males. First sample failure time using multinomial, sample censoring based on normal distribution.            
            
            #For father
            wt = brca1.nc.m
            if (g.all[i,2]>=1) wt = brca1.c.m  
            new<-rmultinom(1,1, prob=wt)
            time[2]<-which(new[,1]==1)
            
            # For remaining male relatives 
            for (j in 1:m1){
            relnum<-4+(2)*(j-1)
            wt = brca1.nc.m
            if (g.all[i,relnum]>=1) wt = brca1.c.m
            new<-rmultinom(1,1, prob=wt)
            time[relnum]<-which(new[,1]==1)            
            }            	
            
 
## females
#femtot<-append(2,4:tot)
#femtot<-m1+m1*m2
 			#For mother
 			   wt = brca1.nc
            	if (g.all[i,1]>=1) wt = brca1.c
           	new<-rmultinom(1,1, prob=wt)
            time[1]<-which(new[,1]==1) 
                      		
 			#For m1 daughters
 			
            for (j in 1:m1){
                relnum<-3+(2)*(j-1)
                wt = brca1.nc
            	if (g.all[i,relnum]>=1) wt = brca1.c
            	new<-rmultinom(1,1, prob=wt)
            	time[relnum]<-which(new[,1]==1)
            	}            


			fail[i,] = time

        }
Tot<-n*tot
cen<-round(rnorm(Tot,55,10))
fail<-as.vector(t(fail))
ind = ifelse(((fail <= cen) & (fail!=111)),1,0)
obs = fail
obs[ind==0] = cen[ind==0]

## output the data as the brcapro
          #ID = seq(1, 5*n)
          ID=rep(1:tot,n)
          Gender = rep(c(0,1,rep(c(0,1), times=m1)),times=n)
          
          #FatherID = rep(seq(1,5*n,by=5),each=5)*rep(c(0,0,1,1,1),times=n)
          
		  ftem<-c(0,0,rep(c(2,0), times=m1))
		  FatherID=rep(ftem, times=n)
		
          #MotherID = rep(seq(2,5*n,by=5),each=5)*rep(c(0,0,1,1,1),times=n)
          
          mtem<-c(0,0,rep(c(1,0), times=m1))
          MotherID=rep(mtem, times=n)
		
          AffectedBreast = ind
          AffectedOvary = rep(0,tot*n)
          AgeBreast =obs
          AgeOvary = rep(1,tot*n)
          #AgeOvary=as.vector(t(origcen))
          AgeBreastContralateral = rep(0,tot*n)
       
          BRCA1 = as.vector(t(g.all))
         # FamRisk = rep(omega,each=5)
          
          FamilyID=rep(1:n, each=tot)
          
          Twins= rep(0,tot*n)
  		
  		Ethnic=rep("Other", tot*n)
           dat = data.frame(cbind(ID,Gender, FatherID, MotherID, AffectedBreast,
              AffectedOvary, AgeBreast, AgeOvary, AgeBreastContralateral, Twins, Ethnic, FamilyID,
              BRCA1))
              
            
          return(dat)
}
