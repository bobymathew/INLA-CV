MCMCglmm


library(MCMCglmm)
library("nadiv")

ped=read.table("pat_ped.txt",header=T)
 	ped$animal=as.factor(ped$animal)
 	ped$sir=as.factor(ped$sir)
 	ped$dam=as.factor(ped$dam)

 

dat=read.table("phe_dom.txt",header=T,na.string="NA")
	dat$Phe=as.numeric(dat$Phe)
	va=var(dat$Phe)

	dat$animal=as.factor(dat$animal)
	dat$animal.d=dat$animal

prior1.2 <- list(G = list(G1 = list(V =0.5, nu = 0.5),G2=list(V =1, nu = 0.5)), R = list(V =0.5, nu = 0.5))
prior2.2 <- list(G = list(G1 = list(V =1, nu = 1),G2=list(V =1, nu = 1)), R = list(V =1, nu = 1))


A=inverseA(ped)
D=makeD(ped)



mod_mcmc=MCMCglmm(Phe~1,random=~animal+animal.d,ginverse=list(animal=A$Ainv,animal.d=D$Dinv),pr=T,data=dat,nitt=50000,thin=500,burnin=10000,prior=prior2.2)

posterior.mode(mod_mcmc$VCV)

#To geth the breeding values set option 'pr=T'

val=colMeans (mod_mcmc$Sol)
n=length(dat$Phe)

breed=val[2:(n+1)]+val[(n+2):(2*n+1)]

write.table(unname(breed),"breed_mcmc.txt",row.names=F,col.names=F)




INLA

library(MCMCglmm)
library("nadiv")
library("INLA")


dat=read.table("phe_dom.txt",header=T,na.string="NA")
	dat$Phe=as.numeric(dat$Phe)

	ped=read.table("pat_ped.txt",header=T)
 
	ped$animal=as.factor(ped$animal)
	ped$sir=as.factor(ped$sir)
	ped$dam=as.factor(ped$dam)
	dat$animal.d=dat$animal

A=inverseA(ped)
D=makeD(ped)


#Phe=animal+animal.d+error(IndexA.2).


formula = Phe ~  f(animal,model="generic0", Cmatrix=A$Ainv,
constr=TRUE,hyper = list(prec = list(prior="loggamma",param = c(0.5, 0.5), fixed = FALSE))) +f(animal.d,model="generic0", Cmatrix=D$Dinv,
constr=TRUE,hyper = list(prec = list(prior="loggamma",param = c(1, 0.05), fixed = FALSE))) 

model = inla(formula=formula, family="gaussian",data=dat,control.family = list(hyper = list(prec = list(prior="loggamma",param = c(0.5, 0.5), fixed = FALSE))),
only.hyperparam =FALSE,control.compute=list(dic=T))


sigma.animal.a = inla.tmarginal(function(x) 1/x,
model$marginals.hyperpar$"Precision for animal")
e.animal.a=inla.emarginal(function(x) x, sigma.animal.a)


sigma.animal.d = inla.tmarginal(function(x) 1/x,
model$marginals.hyperpar$"Precision for animal.d")
e.animal.d=inla.emarginal(function(x) x, sigma.animal.d)


sigma.animal.e = inla.tmarginal(function(x) 1/x,
model$marginals.hyperpar$"Precision for the Gaussian observations")
e.animal.e=inla.emarginal(function(x) x, sigma.animal.e)

 
#breeding values

breed= model$summary.random$animal.d$mean+ model$summary.random$animal$mean
write.table(unname(breed),"breed.txt",row.names=F,col.names=F)



#ASreml with dominace variance 



library(asreml)
library("nadiv")

ped=read.table("pat_ped.txt",header=T)
	ped$animal=as.factor(ped$animal)

ginvA <- asreml.Ainverse(ped)$ginv
ginvD <- makeD(ped)$listDinv

dat=read.table("dat_dom.txt",header=T)
	dat$Phe=as.numeric(dat$Phe)
	dat$animal=as.factor(dat$animal)
	dat$animal.d=dat$animal




# this code is chcked with asreml command promt and same results 


mod <- asreml(Phe ~ 1, random = ~giv(animal) + giv(animal.d),
ginverse = list(animal = ginvA, animal.d = ginvD), data = dat)
summary(mod)$varcomp


#getting breeding values


#mod$coefficients$random
n=length(dat$Phe)
breed=mod$coefficients$random[1:n]+mod$coefficients$random[(n+1):(2*n)]
write.table(unname(breed),"breed_inla.txt",row.names=F,col.names=F)



# inla with the Z matrix this is the classical mixed model

library(MCMCglmm)
library("nadiv")
library("INLA")

dat=read.table("phe_dom.txt",header=T,na.string="NA")
	dat$Phe=as.numeric(dat$Phe)
	reco=dim(dat)[1]

	ped=read.table("pat_ped.txt",header=T)
 
	ped$animal=as.factor(ped$animal)
	ped$sir=as.factor(ped$sir)
	ped$dam=as.factor(ped$dam)
	dat$animal.d=dat$animal

A=inverseA(ped)
D=makeD(ped)

Z=diag(1,reco,reco)

	
	dat$Phe=as.numeric(dat$Phe)
	animal=dat$animal
	animal_d=dat$animal
 	y=dat$Phe

formula = y ~ 1+f(animal,model="z", Cmatrix=A$Ainv, Z=Z,
constr=TRUE,hyper = list(prec = list(param = c(0.5, 0.5), fixed = FALSE)))+f(animal_d,model="z",Cmatrix=D$Dinv,Z=Z,constr=TRUE,hyper = list(prec = list(param = c(0.5, 0.5), fixed = FALSE)))

model = inla(formula=formula, family="gaussian",data=list(y=y, animal=1:reco,animal_d=1:reco ),control.family = list(hyper = list(prec = list(param = c(0.5, 0.5), fixed = FALSE))),
only.hyperparam =FALSE,control.compute=list(dic=T))


sigma.animal = inla.tmarginal(function(x) 1/x,
model$marginals.hyperpar$"Precision for animal")
e.animal=inla.emarginal(function(x) x, sigma.animal)

sigma.animal_d = inla.tmarginal(function(x) 1/x,
model$marginals.hyperpar$"Precision for animal_d")
e.animal_d=inla.emarginal(function(x) x, sigma.animal_d)


sigma.animal.e = inla.tmarginal(function(x) 1/x,
model$marginals.hyperpar$"Precision for the Gaussian observations")
e.animal.e=inla.emarginal(function(x) x, sigma.animal.e)




	




