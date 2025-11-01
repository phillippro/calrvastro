source('mcmc_func.R')
args <- commandArgs(trailingOnly=TRUE)
fs <- read.table('all_companions.txt',header=TRUE)[,1]
stars <- gsub('_.+','',fs)
showtype <- 'complex'
for(f in fs){
load(f,envir=e0 <- new.env())
Nsig <- e0$Nsig
out <- e0$out
tmin <- e0$tmin
mc <- out$mcmc.opt[[paste0('sig',Nsig)]]

#Npar <- ncol(mc)-2
###modify period
indp <- grep('^per\\d',colnames(mc))
mc[,indp] <- exp(mc[,indp])
c0 <- colnames(mc)
c0[indp] <- paste0('Pd',1:Nsig)
colnames(mc) <- c0
plxs <- out$astrometry[nrow(out$astrometry),'parallax']-mc[,'dplx']
ds <- 1/plxs#kpc

#if(target=='UCAC4569-026385'){
#    out$Mstar <- 2.8
#    out$eMstar <- 0.2
#    mstars <- rnorm(nrow(mc),out$Mstar,out$eMstar)
#}

###calculate mass
mcs <- c()
Tps <- c()
if(any(colnames(mc)=='Mstar')){
    mstars <- mc[,'Mstar']
}else{
    mstars <- rnorm(nrow(mc),out$Mstar,out$eMstar)
}
inds <- which(mstars<0)
if(length(inds)>0)   mstars[inds] <- sample(mstars[-inds],length(inds))

for(j in 1:Nsig){
    Mc <- k2m(mc[,paste0('K',j)],mc[,paste0('Pd',j)],mc[,paste0('e',j)],Ms=mstars,Inc=mc[,paste0('Inc',j)])$ms
    Tp <- M02Tp(M0=mc[,paste0('Mo',j)],T0=tmin,P=mc[,paste0('Pd',j)])
    mcs <- cbind(mcs,Mc)
    Tps <- cbind(Tps,Tp)
}
colnames(Tps) <- paste0('Tp',1:Nsig)
mstars <- t(t(mstars))
colnames(mcs) <- paste0('Mc',1:Nsig)
colnames(mstars) <- 'Mstar'
#cc <- c(paste0('Pd',1:Nsig),paste0('K',1:Nsig),paste0('e',1:Nsig),paste0('Inc',1:Nsig),paste0('omega',1:Nsig),paste0('Omega',1:Nsig),paste0('Mo',1:Nsig))
if(showtype=='simple'){
    dat <- cbind('parallax'=plxs,'distance'=ds,Tps,mc)
}else{
    dat <- cbind(mcs,mstars,'parallax'=plxs,'distance'=ds,Tps,mc)
}
stop()
}