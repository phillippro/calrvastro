#########parameters order: Keplerian parameters: {{per},{K},{e},{omega},{Mo}}; trend pars: {a, b}; noise pars (white noise: {s}, Gaussian process red noise: {sigma.red,l}, ARMA red noise:{{phi},alpha,{w},beta}>)
##############################################################################
#####Part I: parameter range
##############################################################################
#cat('par.data[[1]]$q=',par.data[[1]]$q,'\n')
#cat('length(par.data)=',length(par.data),'\n')
#library(Matrix)
source('divide_period.R',local=TRUE)
for(j3 in 1:Nw){
    var <- names(par.data[[ids[j3]]])
    for(k in 1:length(var)){
        assign(var[k],par.data[[ids[j3]]][[var[k]]])
    }
if(!any(var=='cmax')){
    cmax <- cmin <- cini <- c()
}
ind.ok <- which(!is.na(trv))
trv.single <- par.data[[1]]$trv
dur.data <- tmax-tmin
##AK and QP pars
tau.min <- 1
tau.max <- 10*dur.data
logtq.max <- logtau.max <- log(tau.max)
logtq.min <- logtau.min <- log(tau.min)
ta.min <- min(trv[!is.na(trv)])#max(-500,-dur.data)#
ta.max <-  max(trv[!is.na(trv)])#min(-300,dur.data)#
aq.max <- 2*pi
aq.min <- 0
####1.4 other parameter range
#Kmax = 2*max(abs(RV.single-mean(RV.single)))#2129#maximum semi-amplitude of radial velocity; according to Ford & Gregory 07
#amax = Kmax/((tmax-tmin)/time.unit)#unit is m/s/year due to stellar companion etc..
if(!exists('amax|amin|bmax|bmin')){
amax = Kmax/((tmax-tmin)/time.unit)#
#amax <- 1e-3
amin = -amax
bmax = Kmax
bmin = -Kmax
}
smax = Kmax
if(exists('Sindex') & AB){
    ssmax <- Kmax/(max(Sindex)-min(Sindex))#default: 2*max(...)
    ssmin <- 0
}
if(noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='TJ'){
    smin <- 0#-smax
}else{
    smin <- 0
}
lmax = 1e4#100
#lmin = 1#5
lmin = 0.1#5
#lmean <- 393#active-region lifetime of the star
lmean <- 100#active-region lifetime of the star
#lsd <- 30
lsd <- 50
Pgpmax <- 1000
Pgpmin <- 1
Pgp.mean <- 130#rotation period
Pgp.sd <- 5
lp.max <- 1e3
lp.min <- 1e-3
lp.mean <- 0.5#dimensionless; smoothing parameter or maximum number of peaks
lp.sd <- 0.05
#ARMA/TARMA parameters
phi.min <- wmin <- -1
phi.max <- wmax <- 1
if(arma.type=='index'){
    phiS.min <- wS.min <- phiC.min <- wC.min <- -1
    phiS.max <- wS.max <- phiC.max <- wC.max <- 1
}
if(arma.type=='index.exponent'){
    phiS.min <- wS.min <- 0
    phiS.max <- wS.max <- 10
    phiC.max <- wC.max <- 1
    phiC.min <- wC.min <- -1
}
manual <- FALSE
if(!manual){
    beta.up <- log(tmax-tmin)#time span of the data
    beta.low <- log(max(1/24,min(1,min(diff(trv.single)))))#1h or minimum separation
#    beta.low <- log(max(1/24,min(diff(trv.all))))#1h or minimum separation
}else{
    if(grepl('wbin',id)){
        beta.low <- 0
        beta.up <- 4
    }else{
        beta.low <- log(min(diff(trv.single)))#-1
        beta.up <- 1 
    }
}
alpha.max <- beta.max <- beta.up#d; limit the range of beta to avoid multimodal or overfitting
alpha.min <- beta.min <- beta.low#24h
phis.min <- phia.min <- phib.min <- phis.max <- phia.max <- phib.max <- c()
alphas.min <- alphaa.min <- alphab.min <- alphas.max <- alphaa.max <- alphab.max <- c()
if(noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='PSID' | noise.model=='ARMAPSID'){
    phis.max <- min(30,Kmax/(max(Sindex)-min(Sindex)))
    phis.min <- 0#phi.min#-phi.ar.max
#    alpha.ar.min <- 0.005-1e-3
    alphas.min <- alpha.min*time.unit#unit: yr^-1
#    alpha.ar.max <- 0.005+1e-3
    alphas.max <- 1#alpha.max
    if(exists('BIS') & exists('FWHM') & TJAR.AB){
        phia.max <- min(30,Kmax/(max(BIS)-min(BIS)))
        phia.min <- 0
        phib.max <- min(30,Kmax/(max(FWHM)-min(FWHM)))
        phib.min <- 0
        alphaa.max <- alphab.max <- alphas.max
        alphaa.min <- alphab.min <- alphas.min
    }
    if(TJAR.sym=='asym'){
#        psi.ar.max <- phi.max
        psis.max <- min(30,Kmax/(max(Sindex)-min(Sindex)))
        psis.min <- 0#-psi.ar.max
        if(exists('BIS') & exists('FWHM') & TJAR.AB){
            psia.max <- min(30,Kmax/(max(BIS)-min(BIS)))
            psia.min <- 0
            psib.max <- min(30,Kmax/(max(FWHM)-min(FWHM)))
            psib.min <- 0
        }
    }
}
#####################################################
####Part II: initial parameter values
#####################################################
Nkeppar <- 5
if(kep.type=='AK' & prior.type!='e0'){
    Nkeppar <- 7
}else if(kep.type=='AK' & prior.type=='e0'){
    Nkeppar <- 5
}else if(prior.type=='e0'){
    Nkeppar <- 3
}

######2.1 create initial conditions
tmp.kepini <- c()
tmp.gpini <- c()
tmp.arini <- c()
tmp.tjar.ini <- c()
tmp.tarini <- c()
tmp.maini <- c()
tmp.tmaini <- c()
tmp.tjini <- c()
if(Np>0){
    np <- Np
    if(prior.type=='e0'){
        kep1 <- c(per.ini,1, 1)
    }else{
        kep1 <- c(per.ini,runif(1,Kmin,Kmax),0, 1, 1)
    }
    if(grepl('AK',kep.type)){
        kep1 <- c(kep1,c((logtau.max+logtau.min)/2,(ta.min+ta.max)/2))
    }
    if(grepl('QP',kep.type)){
        kep1 <- c(kep1,c((logtq.max+logtq.min)/2,(aq.min+aq.max)/2))
    }
    if(Np>1 & exists('par.opt')){
        tmp.kepini <- c(par.opt[1:((Np-1)*Nkeppar)],kep1)
    }else if(Np>1 & !exists('par.opt')){
        tmp.kepini <- rep(kep1,np)
    }else{
        tmp.kepini <- kep1
    }
}
if(noise.model=='TJ' | noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='ARMATJ'){
    tmp.tjini <- rep(0.1,Ntj)
}
if(noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='PSID'  | noise.model=='ARMAPSID'){
    if(TJAR.sym=='asym'){
        tmp.tjar.ini <- c(rep(c(0.1,0.1),Par),(alphas.min+alphas.max)/2)
        if(exists('BIS') & exists('FWHM') & TJAR.AB){
            tmp.tjar.ini <- c(tmp.tjar.ini,rep(c(0.1,0.1),Par),(alphaa.min+alphaa.max)/2)
            tmp.tjar.ini <- c(tmp.tjar.ini,rep(c(0.1,0.1),Par),(alphab.min+alphab.max)/2)
        }
    }else{
        tmp.tjar.ini <- c(rep(0.1,Par),(alphas.min+alphas.max)/2)
        if(exists('BIS') & exists('FWHM') & TJAR.AB){
            tmp.tjar.ini <- c(tmp.tjar.ini,rep(0.1,Par),(alphaa.min+alphaa.max)/2)
            tmp.tjar.ini <- c(tmp.tjar.ini,rep(0.1,Par),(alphab.min+alphab.max)/2)
        }
    }
}
#####GP or GPR: initial red noise covariance basic matrix
if(noise.model=='GP' | noise.model=='GPR'){
    cov.red.tmp <- rednoise.cov(0.1,l=lmean,tt=trv.single,tol=1e-10,type=gp.type)
#    cov.red.tmp[cov.red.tmp<1e-3] <- 0
#    cov.red.basic <- Matrix(cov.red.tmp,sparse=TRUE)
    cov.red.basic <- cov.red.tmp
    if(any(is.na(cov.red.basic))) cat('NA in cov.red.basic\n')
    if(gp.type=='abs' | gp.type=='sq'){
        tmp.gpini <- c(smin,lmean)
    }
    if(gp.type=='qp'){
        tmp.gpini <- c(smin,lmean,Pgp.mean,lp.mean)
    }
}
if(noise.model=='ARMA' | noise.model=='ARMATJ' | noise.model=='ARMATJAR' | noise.model=='ARMAPSID'){
    if(p>0){
        tmp.arini <- c(rep(0.1,p),(alpha.max+alpha.min)/2)#alpha, beta in [1/1min,1/1year] according Tuomi13
    }
    if(q>0){
        tmp.maini <- c(rep(0.1,q),(beta.max+beta.min)/2)
    }
}
if(noise.model=='TARMA'){
    if(p>0){
        tmp.tarini <- c(rep(0.1,p),rep(0,p),0.5)#alpha, beta in [1/1min,1/1year] according Tuomi13
    }
    if(q>0){
        tmp.tmaini <- c(rep(0.1,q),rep(0,q),0.5)
    }
}
##order of model parameters: kepler, trend, jitter, time-varying jitter (TJ), GP, ARMA, TJAR, TARMA, activity index
if(j3==1){
    iniref <- tmp.kepini
}
#    x <- (trv-min(trv))/time.unit
    x <- (trv-tmin)/time.unit
    polyfit <- lm(RV~poly(x,Npoly,raw=TRUE))
    poly.par <- as.numeric(polyfit$coefficients)
    aini <- poly.par[2:length(poly.par)]
    bini <- poly.par[1]
    if(abs(aini)<10){
        alpha <- 1
    }else{
        alpha <- 0.01
    }
#    cat('aini=',aini,'\n')
#    cat('bini=',bini,'\n')
    amax <- aini+alpha*abs(aini)
    amin <- aini-alpha*abs(aini)
    bmax <- bini+alpha*abs(bini)
    bmin <- bini-alpha*abs(bini)
#aini <- rep(0,Npoly)
#bini <- 0
iniref <- c(iniref,c(aini,bini,0),c(tmp.tjini,tmp.gpini,tmp.arini,tmp.maini,tmp.tjar.ini,tmp.tarini,tmp.tmaini),cini)
if(Nes>1){
    for(i3 in 2:Nes){
        if(any(grepl('^a',ep.par))){
	    iniref <- c(iniref,0)
	}
        if(any(grepl('^b',ep.par))){
	    iniref <- c(iniref,0)
	}
        if(any(grepl('^s',ep.par))){
	    iniref <- c(iniref,0)
	}
        if(any(grepl('ma',ep.par))){
            iniref <- c(iniref,tmp.maini)
        }
        if(any(grepl('^c',ep.par))){
            iniref <- c(iniref,rep(cini[1],length(grep('^c',ep.par))))
        }
        if(any(grepl('d[^c]',ep.par))){
            iniref <- c(iniref,rep(dw.ini.mix[1],length(grep('d[^c]',ep.par))))      
        }
	if(any(grepl('dc',ep.par))){
	    iniref <- c(iniref,rep(dc.ini[1],length(grep('dc',ep.par))))      
	}
    }
}
######2.2 or input initial conditions
####2.3 initial covariance and parameter range
tmp.kep.max <- tmp.kep.min <- tmp.tj.min <- tmp.tj.max <- tmp.gp.max <- tmp.gp.min <- tmp.ar.max <- tmp.tar.max <- tmp.ar.min <- tmp.tar.min <- tmp.ma.max <- tmp.tma.max <- tmp.ma.min <- tmp.tma.min <- tmp.tjar.max <- tmp.tjar.min <-  c()
emin <- 0
emax <- 1
if(Np>0){
    for(jj in 1:Np){
        if(prior.type!='e0'){
            tmp.kep.max <- c(tmp.kep.max,c(per.max,Kmax,emax,2*pi,2*pi))
            tmp.kep.min <- c(tmp.kep.min,c(per.min,Kmin,emin,0,0))
        }else{
            tmp.kep.max <- c(tmp.kep.max,c(per.max,Kmax,2*pi))
            tmp.kep.min <- c(tmp.kep.min,c(per.min,Kmin,0))
        }
        if(grepl('AK',kep.type)){
            tmp.kep.max <- c(tmp.kep.max,c(logtau.max,ta.max))
            tmp.kep.min <- c(tmp.kep.min,c(logtau.min,ta.min))
        }
        if(grepl('QP',kep.type)){
            tmp.kep.max <- c(tmp.kep.max,c(logtq.max,aq.max))
            tmp.kep.min <- c(tmp.kep.min,c(logtq.min,aq.min))
        }
    }
}
if(noise.model=='TJ' | noise.model=='ARMATJ' | noise.model=='TJAR' | noise.model=='ARMATJAR'){
    tmp.tj.max <- rep(ssmax,Ntj)
    tmp.tj.min <- rep(ssmin,Ntj)
}
if(noise.model=='GP' | noise.model=='GPR'){
    if(gp.type=='sq' | gp.type=='abs'){
        tmp.gp.max <- c(smax,lmax)
        tmp.gp.min <- c(smin,lmin)
    }
    if(gp.type=='qp'){
        tmp.gp.max <- c(smax,lmax,Pgpmax,lp.max)
        tmp.gp.min <- c(smin,lmin,Pgpmin,lp.min)
    }
}
if(noise.model=='ARMA' | noise.model=='ARMATJ' | noise.model=='ARMATJAR' | noise.model=='ARMAPSID'){
    if(p>0){
        tmp.ar.max <- c(rep(phi.max,p),alpha.max)
        tmp.ar.min <- c(rep(phi.min,p),alpha.min)
    }
    if(q>0){
        tmp.ma.max <- c(rep(wmax,q),beta.max)
        tmp.ma.min <- c(rep(wmin,q),beta.min)
    }
}
if(noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='PSID' | noise.model=='ARMAPSID'){
    if(Par>0){
        if(TJAR.sym=='asym'){
            tmp.tjar.max <- c(c(rep(c(phis.max,psis.max),Par),alphas.max),c(rep(c(phia.max,psia.max),Par),alphaa.max),c(rep(c(phib.max,psib.max),Par),alphab.max))
            tmp.tjar.min <- c(c(rep(c(phis.min,psis.min),Par),alphas.min),c(rep(c(phia.min,psia.min),Par),alphaa.min),c(rep(c(phib.min,psib.min),Par),alphab.min))
        }else{
            tmp.tjar.max <- c(c(rep(phis.max,Par),alphas.max),c(rep(phia.max,Par),alphaa.max),c(rep(phib.max,Par),alphab.max))
            tmp.tjar.min <- c(c(rep(phis.min,Par),alphas.min),c(rep(phia.min,Par),alphaa.min),c(rep(phib.min,Par),alphab.min))
        }
    }
}
if(noise.model=='TARMA'){
if(p>0){
    tmp.tar.max <- c(rep(phiS.max,p),rep(phiC.max,p),alpha.max)
    tmp.tar.min <- c(rep(phiS.min,p),rep(phiC.min,p),alpha.min)
}
if(q>0){
    tmp.tma.max <- c(rep(wS.max,q),rep(wC.max,q),beta.max)
    tmp.tma.min <- c(rep(wS.min,q),rep(wC.min,q),beta.min)
}
}
if(j3==1){
    par.max <- tmp.kep.max
    par.min <- tmp.kep.min
}
    if(FALSE){
        cat('amin=',amin,'\n')
        cat('amax=',amax,'\n')
        cat('bmin=',bmin,'\n')
        cat('bmax=',bmax,'\n')
    }
par.min <- c(par.min,c(amin,bmin,smin),c(tmp.tj.min,tmp.gp.min,tmp.ar.min,tmp.ma.min,tmp.tjar.min,tmp.tar.min,tmp.tma.min,cmin))
par.max <- c(par.max,c(amax,bmax,smax),c(tmp.tj.max,tmp.gp.max,tmp.ar.max,tmp.ma.max,tmp.tjar.max,tmp.tar.max,tmp.tma.max,cmax))
if(length(nepoch)>1){
    for(i3 in nepoch[2:length(nepoch)]){
        if(any(grepl('^a',ep.par))){
            par.max <- c(par.max,amax)
            par.min <- c(par.min,amin)
        }
        if(any(grepl('^b',ep.par))){
            par.max <- c(par.max,bmax)
            par.min <- c(par.min,bmin)
        }
        if(any(grepl('^s',ep.par))){
            par.max <- c(par.max,smax)
            par.min <- c(par.min,smin)
        }
        if(any(grepl('ma',ep.par))){
            par.max <- c(par.max,tmp.ma.max)
            par.min <- c(par.min,tmp.ma.min)
        }
        if(any(grepl('^c',ep.par))){
            inds <- grep('^c',ep.par)
            if(length(inds)>0){
                index <- as.integer(gsub('^c','',ep.par[inds]))
            }
            par.max <- c(par.max,cmax[index])
            par.min <- c(par.min,cmin[index])
        }
    }
}
}
#########initial conditions
if(!quantify){
    startvalue <- assign.names(iniref,Np=Np,p=p,q=q)
}else{
    startvalue <- assign.names(par.opt,Np=Np,p=p,q=q)
}
nc <- length(grep('\\dc',names(startvalue)))
na <- length(grep('\\da',names(startvalue)))
nb <- length(grep('\\db',names(startvalue)))
ns <- length(grep('\\ds',names(startvalue)))
startvalue <- fix(startvalue,fix.par)
###prior hyper parameters
Npar <- length(startvalue)
Sd <- 2.4^2/Npar#hyper par of s prior 
s0 <- 1#hyper par of s prior 
K0 <- 1#hyper par of K prior
Esd <- 0.02#hyper par of the e prior
###########initial covariance matrix
names(par.max) <- names(par.min) <- names(startvalue)
#par.ref <- (1.1*par.max+par.min)/2
par.ref <- par.max
if(grepl('GP',noise.model)){
    par.ref[grep('\\dl\\d',names(par.ref))] <- lmean
    if(gp.type=='qp'){
        par.ref[grep('\\dlp\\d',names(par.ref))] <- lp.mean
        par.ref[grep('\\dPgp\\d',names(par.ref))] <- Pgp.mean
    }
}
if(Np>0){
    cat('per.ini=',per.ini,'\n')
    cat('Np=',Np,'\n')
    par.ref[1+(1:Np-1)*Nkeppar] <- per.ini#a reference value used to generate initial period values
    if(exists('Kini')){
        cat('Kini=',Kini,'m\n')
        par.ref[2+(1:Np-1)*Nkeppar] <- Kini
    }
}
par.ref <- fix(par.ref,fix.par)
if(Np>1 & !quantify){
    par.ref[1:((Np-1)*Nkeppar)] <- 1e-6*par.ref[1:((Np-1)*Nkeppar)]
    par.ref[((Np-1)*Nkeppar+1):length(par.ref)] <- 1e-2*par.ref[((Np-1)*Nkeppar+1):length(par.ref)]
    if(fixP){
        par.ref[1+(0:(Np-2))*Nkeppar] <- 0#
    }
    cov.start <- abs(diag(par.ref))
}else if(!quantify){
    cov.start <- inicov*abs(diag(par.ref))
}else{
    cov.start <- 1e-6*abs(diag(par.ref))
}
if(exists('cov0')){
    if(nrow(cov0)==Npar){
        cov.start <- cov0
    }
    rm(cov0)
}
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
