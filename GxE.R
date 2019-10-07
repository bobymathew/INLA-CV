Sys.time()


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

dat=read.table("barley_rep.txt",header=T,na.strings="NA")
loc=3
sam=80
n=82
gui=loc*n
reco=dim(dat)[1]
library("nadiv")
library(MCMCglmm)
library(MASS)
set.seed(1200)
f=seq(1:sam)
			d=sample(f,sam)
			cv1=d[1:8]
			cv2=d[9:16]
			cv3=d[17:24]
			cv4=d[25:32]
			cv5=d[33:40]
			cv6=d[41:48]
			cv7=d[49:56]
			cv8=d[57:64]
			cv9=d[65:72]
			cv10=d[73:80]
			CV=cbind(cv1,cv2,cv3,cv4,cv5,cv6,cv7,cv8,cv9,cv10)

A_mat=read.table("A_matrix.txt",header=F)
	A_new=as.matrix(A_mat)
	A_in=ginv(A_new)
	D=diag(1,gui,gui)


Z=desgn(dat$animal,n)
Q=desgn(dat$gui,gui)

X1=matrix(1,reco)
			X2=desgn(dat$year,3)
			X=cbind(X1,X2)
			
M=cbind(X,Z,Q)

			
for ( i in 1:10)
{

	add_var=c(9.1,9.6,10.1,10.6,11.1,8.6,8.1,7.6,7.1,6.6)
	
	for ( j in 1:10)
	{
		dom_var=c(3.2,3.5,3.8,4.1,4.4,2.8,2.5,2.2,1.9,4.6)
					
		for ( k in 1:10)
		{
			err_var=c(17,17.5,18,18.5,19,19.5,16.5,16,15.5,15)
				
			add_rat=err_var[k]/add_var[i]	
			dom_rat=err_var[k]/dom_var[j]	
			
			A_inv=A_in*add_rat
			D_inv=D*dom_rat
			
			final=numeric(0)
			
			
			for(l in 1:10)
			
					{
					miss=CV[,l]	
			dat_new=dat;					
			dat_miss=dat[dat$animal %in% miss,];
			dat_new[dat_new$animal %in% miss,]=0;
			new_n= dim(dat_new)[1]
			
			X[721,3]=0
			X[731,2]=0
			
			XX = t(X) %*% X
    			XZ = t(X) %*% Z 
			XQ=t(X) %*% Q 
    			 ZZ = t(Z) %*% Z
			QQ=t(Q) %*% Q
			ZQ=t(Z) %*%Q
    			Xy = t(X) %*% dat_new$TKM
   			 Zy = t(Z) %*% dat_new$TKM
			Qy = t(Q) %*% dat_new$TKM
   			 R1 = cbind(XX,XZ,XQ)
   			 R2 = cbind(t(XZ),(ZZ+A_inv),ZQ)
   			 R3=  cbind(t(XQ),t(ZQ),(QQ+D_inv))
    			LHS = rbind(R1,R2,R3)
    			RHS = rbind(Xy,Zy,Qy)

			C = ginv(LHS)
  			bhat= C %*% RHS
			
			#bhat = solve(LHS,RHS)
			breed_t=M%*%bhat
			blup=sum((dat$TKM[dat_miss$Line]-breed_t[dat_miss$Line])^2)
			final=c(final,blup)

			
			}
			
			} #closing for cv

final_va=cbind(add_var[i],dom_var[j],err_var[k],mean(final))
write.table( final_va,"res_blup_barley.txt",append=T,col.names=F,row.names=F)

		} #closing for err

	} #closing for domiance

} #closing for additive

Sys.time()
			
	
