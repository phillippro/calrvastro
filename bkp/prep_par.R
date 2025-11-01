iniref <- c()
startvalue <- assign.names(iniref,Np=Np,p=p,q=q)
if(exists('par.opt')){
    if(length(par.opt)==length(iniref)){
        startvalue <- assign.names(par.opt,Np=Np,p=p,q=q)
    }
}
nc <- length(grep('\\dc\\d',names(startvalue)))
na <- length(grep('\\da\\d',names(startvalue)))
nb <- length(grep('\\db\\d',names(startvalue)))
ns <- length(grep('\\ds\\d',names(startvalue)))
if(all(fix.par!='')){
    startvalue <- fix(startvalue,fix.par)
}
###prior hyper parameters
Npar <- length(startvalue)
Sd <- 2.4^2/Npar#hyper par of s prior 
s0 <- 1#hyper par of s prior 
K0 <- 1#hyper par of K prior
Esd <- 0.2#hyper par of the e prior
if(length(trv.all)>100 & median(diff(trv.all))<0.1){
   Esd <- 0.01#hyper par of the e prior
}
###########initial covariance matrix
names(par.max) <- names(par.min) <- names(startvalue)
par.ref <- par.max
######
#############define parameters for naming output files
if(period.par=='P'){
	iniP <- startvalue['per1']
}else if(period.par=='logP'){
	iniP <- exp(startvalue['per1'])
}else{
	iniP <- 1/startvalue['per1']
}
if(noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='PSID' | noise.model=='ARMAPSID'){
    TJAR.name <- paste0(TJAR.sym,Par,TJAR.mode,'AB',TJAR.AB)
}else{
    TJAR.name <- c()
}
#cat('names(startvalue)=',names(startvalue),'\n')
#cat('startvalue=',startvalue,'\n')
#cat('par.min=',par.min,'\n')
#cat('par.max=',par.max,'\n')
if(noise.model=='GP' & gp.single & !gp.type=='sho'){
    cov.red.basic <- rednoise.cov(0.1,l=lmean,tt=trv.all,tol=1e-10,type=gp.type)
}
