



desgn = function(v,reco)
	{
   	if(is.numeric(v))
		{
     		va = v
     		mrow = length(va)
      		mcol = reco }
       		X = matrix(data=c(0),nrow=mrow,ncol=mcol)
    		for(i in 1:mrow)
			{
      			ic = va[i]
      			X[i,ic] = 1
			}
     	return(X)  
	}   
Sys.time()
ped=read.table("pat_ped.txt",header=T,na.strings="NA")
dat=read.table("phe_dom.txt",header=T,na.strings="NA")


library("nadiv")
library(MCMCglmm)
library(MASS)
library(psych)

A_m=makeA(ped)
D_m=makeD(ped)
D=D_m$D
A=inverseA(ped)



reco=length(ped$animal)
n=3640
A_in=as.matrix(A$Ainv)
D_in=as.matrix(D_m$Dinv)

X=matrix(1,n,1)
			
			Z=desgn(dat$animal,reco)

			XX = t(X) %*% X
    			XZ = t(X) %*% Z 
    			 ZZ = t(Z) %*% Z
    			Xy = t(X) %*% dat$Phe
   			 Zy = t(Z) %*% dat$Phe
   			 R1 = cbind(XX,XZ,XZ)
			RHS = rbind(t(X),Z,Z)
			 H=cbind(X,Z,Z)		

for ( i in 1:10)
{

	add_var=c(880,910,940,970,1000,1030,850,820,790,1060)
	
	for ( j in 1:10)
	{
		dom_var=c(450,465,480,495,510,525,435,420,405,390)
					
		for ( k in 1:10)
		{
			err_var=c(710,725,740,695,680,665,650,635,620,605)





	
			final_va=numeric(0)
			add_rat=err_var[k]/add_var[i]
			dom_rat=err_var[k]/dom_var[j]
			A_inv=A_in*add_rat
			D_inv=D_in*dom_rat
			
	
					
   			 R2 = cbind(t(XZ),(ZZ+A_inv),ZZ)
   			 R3=  cbind(t(XZ),ZZ,(ZZ+D_inv))
    			 LHS = rbind(R1,R2,R3)
			 sol=solve(LHS,RHS)#inverse(C)*Z'
		
			 bhat=sol%*%dat$Phe
			 S=H%*%sol
			 df=tr(S)
			 breed_t=bhat[2:3641]+bhat[3642:7281]
			 blup=(dat$Phe-breed_t)^2 
			 div=n*(1-(df/n))^2
			 res=sum(blup)/div
			



	final_va=cbind(add_var[i],dom_var[j],err_var[k],res,df)
write.table( final_va,"res_gcv_trace.txt",append=T,col.names=F,row.names=F)

		} #closing for err

	} #closing for domiance

} #closing for additive
Sys.time()


