# - - - - - - - - - - - - - - -
# Supplementary code for:
# Characterizing the transmission potential of zoonotic infections from minor outbreaks
# Authors: Adam J. Kucharski and W. John Edmunds
# - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - -
# 1 Set up simulation model and inference functions for R0 and susceptibility
# - - - - - - - - - - - - - - -

#  - - - - - - - - - - 
# 1a Define simulation model
#  - - - - - - - - - - 

simulate.data <- function (r0,dispk,mm,spillrates,simruns){
	# Max number of generations modelled
	genrns=500
  #Convert mixing vector into matrix for start in group 1
	mm01=matrix(mm, nrow = 2, ncol = 2, byrow = TRUE)
	# Transpose for start in group 2
	mm02=apply(apply(mm01, 1, rev),1,rev)
  # Store simulation outputs
	simsize=matrix(rep(1, 3*simruns), ncol = 3)
	
# Simulate multiple outbreaks and store distributions
	
	for(nn in 1:simruns){
		
		if(runif(1)<spillrates[1]){seedgp=1}else{seedgp=2} # randomly choose which group spillover occurs in
		if(seedgp==1){mm0=mm01}else{mm0=mm02} # Set appropriate matrix
		
		r00=r0/max(eigen(mm0)$values) # Normalise matrix so dominant eigenvalue=R0
		mm=r00*mm0
		# Define the model
		casemodel <- function (genrns, dispk, mm){ 
			
			cases1=vector(len=genrns) #number each generation
			cases2=vector(len=genrns) #number each generation
			d_ind=1;
			for (d in 1:genrns){
				
				if (d==1){
					cases1[d]=1
					cases2[d]=0
				}else{
					cases1[d]=sum(rnbinom(cases1[d-1],size=dispk,mu=mm[1,1]))+sum(rnbinom(cases2[d-1],size=dispk,mu=mm[2,1]))
					cases2[d]=sum(rnbinom(cases1[d-1],size=dispk,mu=mm[1,2]))+sum(rnbinom(cases2[d-1],size=dispk,mu=mm[2,2]))
				}# number of new cases in that generation
			}
			cbind(cases1,cases2)
		}
    # Simulate model
		simcase=casemodel(genrns, dispk, mm)
		
		if(seedgp==1){simsize[nn,]=c(seedgp,sum(simcase[,1]),sum(simcase[,2]))}else
		{simsize[nn,]=c(seedgp,sum(simcase[,2]),sum(simcase[,1]))}
	}
	
	simsize
}

#  - - - - - - - - - - 
# 1b Define inference model
#  - - - - - - - - - - 

# START Outbreak size distribution - multi-type model

probrn1n2<- function(seedgp,n01,n02,rstar,mm){

# Set up matrices for different spillover groups
	
mm01=matrix(mm, nrow = 2, ncol = 2, byrow = TRUE)
mm02=apply(apply(mm01, 1, rev),1,rev)
	
if(seedgp==1){mm0=mm01}else{mm0=mm02}
if(seedgp==1){n1=n01;n2=n02}else{n2=n01; n1=n02}
	
r00=rstar/max(eigen(mm0)$values)
mm=r00*mm0

# Make sure that sum limits are well defined
if(n2==0|mm[2,1]==0){k1max=1}else{k1max=n1} 

# Store contributions to probability
prob.comb=array(1, dim=c(k1max,n2+1))

for(kk1 in 1:k1max){
	for(kk2 in 1:(n2+1)){
		
		a12=kk2-1 #set cases in group 2 caused by group 1
		if(mm[2,1]==0){a21=0}else{a21=(kk1-1)} #set cases in group 1 caused by group 2
		if (a21>n1-1 | a12>n2){stop('Need k<n')}
		nn=c(n1,n2)
		gg=matrix(c(n1-a21-1,a12,a21,n2-a12), nrow = 2, ncol = 2, byrow = TRUE)
		
# offspring function
		gen.ij<- function(i,j,c1){ 
			t1=nn[i]
			# get product
			prod1=rep(1,c1)
			for(z in 0:(c1-1)){prod1[z+1]<-(t1+z)/(factorial(c1)^(1/(c1-1)))}
			if(c1==0){prod2=1}else{prod2=prod(prod1)}
			# probability t1 infectives generate c1 cases 
			prod2*mm[i,j]^c1*(1+mm[i,j])^(-t1-c1)
		}
		
# offspring probability
		prob<- function(n2t,a21t){
			if (n2t==0){
				gen.ij(1,1,gg[1,1])*gen.ij(1,2,0)/(nn[1])
			}else{
				if (n2t>0 & a21t==0){
					a12*gen.ij(1,1,gg[1,1])*gen.ij(1,2,gg[1,2])*gen.ij(2,1,0)*gen.ij(2,2,gg[2,2])/(nn[1]*nn[2])
				}else{
					a12*gen.ij(1,1,gg[1,1])*gen.ij(1,2,gg[1,2])*gen.ij(2,1,gg[2,1])*gen.ij(2,2,gg[2,2])/(nn[1]*nn[2])
				}
			}
		}
		prob.comb[kk1,kk2]=prob(n2,a21)
	}
}
sum(prob.comb) # total probability n1 cases in group 1 & n2 cases in group 2
}

# END Outbreak size distribution - multi-type model



# START LIKELIHOOD FUNCTION - multi-type model

likelihoodM1<- function(simdata,runs,rstar,mm){
	liks=rep(1,runs)
	for(z in 1:runs){liks[z]=probrn1n2(simdata[z,1],simdata[z,2],simdata[z,3],rstar,mm)}
	sum(log(liks))
}
# END LIKELIHOOD FUNCTION - multi-type model


# START LIKELIHOOD FUNCTION - single-type model
hg.rj<- function(j,rstar){ 
	# get product
	prod1=rep(1,j)
	for(z in 0:(j-2)){prod1[z+1]<-j+z}
	if(j==0){prod2=1}else{prod2=prod(prod1)}
	# probability t1 infectives generate c1 cases
	prod2*rstar^(j-1)*(1+rstar)^(-2*j+1)/factorial(j)
}

likelihoodHG<- function(simdata,runs,rstar,mm){
	
	liks=rep(1,runs)
	for(z in 1:runs){liks[z]=hg.rj(simdata[z,2]+simdata[z,3],rstar)}
	sum(log(liks))
}
# END LIKELIHOOD FUNCTION - single-type model



# Define functions to find point estimate for R0 and S using multi-type model

R0estimateM1v<- function(simdata,mm,rtest,susc){
	simdata=data1
	mm=infmatrix
  # Define boundaries of likelihood calculation
	R0r1=0; R0r2=2; R0range=seq(R0r1,R0r2,by=.001)
	kkr1=0; kkr2=1; kkrange=seq(kkr1,kkr2,by=.001)
	collect1 <- data.frame(matrix(NA, nrow=length(R0range),length(kkrange)))
	for(ii in 1:length(R0range)){
		for(jj in 1:length(kkrange)){
			runs1=length(simdata)/3
			rr=R0range[ii]; v1=kkrange[jj]
			mm2=mm; mm2[2]=mm[2]*v1; mm2[4]=mm[4]*v1
			collect1[ii,jj]=likelihoodM1(simdata,runs1,rr,mm2)
		}
	}
	likm=(collect1==max(collect1))
	r0max=max(seq(1,length(R0range))%*%likm)
	c(R0range[r0max],kkrange[max(likm[r0max,]*seq(1,length(kkrange)))])
}

# Find point estimate for R0 using single-type model

R0estimateHG<- function(simdata,mm){
	meansize=mean(simdata[,2]+simdata[,3])
	1-1/meansize
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2 Simulate data, perform inference of R0 & S, and calculate mean squared & absolute error
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Define 4 different scenarios for R0 and relative susceptibility of udner 20 age group:
susctab=c(0.25,1,0.25,1)
rbasictab=c(0.25,0.25,0.75,0.75)

# Number of outbreak clusters to simulate
ireps=50

# Data frame to store outputs
R0errortable<-data.frame(matrix(NA, nrow=4,ncol=7))
names(R0errortable)=c("R0","vac","relerrorM1","relerrorHG","abserrorM1","abserrorHG","Reff")

for(kk in 1:4){

	susc=susctab[kk]
	rbasic=rbasictab[kk]
	R0table=array(NA, dim=c(ireps,4))
	
	# Define matrix used for inference. In form: from i to i, from i to j, from j to i, from j to j
	infmatrix=c(4.27, 3.04, 1.32, 2.70)
	
	# Define matrix used for simulation
	simmatrix0=c(4.27, 3.04, 1.32, 2.70)
	simmatrix=c(4.27,3.04*susc,1.32,2.70*susc)
              
	 # Define probability of spillover into each group based on UK population distribution
  seeds=c(0.24/(0.24+0.76*susc), 0.76*susc/(0.24+0.76*susc))
	
  # Scale matrix depending on R0
	mm01=matrix(simmatrix0, nrow = 2, ncol = 2, byrow = TRUE)
	r00=rbasic/max(eigen(mm01)$values)
	
	mm02=r00*matrix(simmatrix, nrow = 2, ncol = 2, byrow = TRUE)
	rtest=max(eigen(mm02)$values)
  
  # Simulate data and infer parameters
	
	for(ii in 1:ireps){
		data1=simulate.data(rtest, 1, simmatrix, seeds, 50)
		rtab1=R0estimateM1v(data1,infmatrix,round(rtest, digits = 2),susc)
		rtab2=R0estimateHG(data1,infmatrix)
		
		mm3=infmatrix; mm3[2]=infmatrix[2]*rtab1[2]; mm3[4]=infmatrix[4]*rtab1[2]
		mm3a=matrix(mm3, nrow = 2, ncol = 2, byrow = TRUE)
		mm5=matrix(infmatrix, nrow = 2, ncol = 2, byrow = TRUE)
		
		R0table[ii,1]=rtab1[1]
		R0table[ii,2]=rtab2[1]
		R0table[ii,3]=rtab1[2]
		R0table[ii,4]=rtab1[1]*max(eigen(mm3a)$values)/max(eigen(mm5)$values)
		}
		
	R0table1=na.omit(R0table)
			
	R0errortable[kk,1]=mean(R0table1[,4])
	R0errortable[kk,2]=mean(R0table1[,3])
	R0errortable[kk,3]=sqrt(sum(((R0table1[,1]-rtest)/rtest)^2)/ireps) # relative MSE
	R0errortable[kk,4]=sqrt(sum(((R0table1[,2]-rtest)/rtest)^2)/ireps) # relative MSE
	R0errortable[kk,5]=sum(R0table1[,1]-rtest)/ireps # abs error
	R0errortable[kk,6]=sum(R0table1[,2]-rtest)/ireps # abs error
	R0errortable[kk,7]=rtest
			
	write.csv(R0table1,paste("Vinfer_R0table",kk,".csv",sep=""))
	write.csv(R0errortable,paste("Vinfer_.csv",sep=""))

}

