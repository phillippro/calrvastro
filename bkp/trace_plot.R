if(!exists('save.memory')) save.memory <- TRUE
par(mfrow=c(4,4))
if(!exists('mc')) mc <- out[['mcmc.opt']][[paste0('sig',Nsig)]]
#mcmc.opt <- out[['mcmc.opt']][[paste0('sig',Nsig)]]#single chain
Npar <- ncol(mc)-2
ns <- colnames(mc)[1:Npar]
panel <- 0
if(!save.memory){
    lp <- as.numeric(sapply(ind.chain,function(k) out[['mcmc.all']][[paste0('sig',Nsig)]][[k]][,Npar+1]))
}else{
    lp <- as.numeric(mc[,Npar+1])
}
for(j in 1:Npar){
    pars <- mc[,j]
    plot(sample(pars,min(round(length(pars)/2),1e3)),xlab='Iteration',ylab=ns[j],type='l')
#    plot(mcmc.opt[,j],mcmc.opt[,Npar+1],xlab='Iteration',ylab='Posterior')
    if(sd(pars)==0) pars <- rnorm(length(pars),pars[1],1e-10)
    p <- data.distr(x=pars,xlab=ns[j],ylab='Frequency',lp=lp,plotf=TRUE)
    panel <- panel+1
    if(panel==1){
        mtext("Trace Plot and Posterior Distribution", outer=TRUE,  cex=1.5, line=-0.5)
    }
}

panel <- 0
par(mfrow=c(4,4))
###companion mass
for(j in 1:Nsig){
    if(any(dP==j)){
        p <- Prefs[j]+par.opt[paste0('dP',j)]/60/24#day
        ps <- Prefs[j]+mc[,paste0('dP',j)]/60/24#day
    }else{
        p <- exp(par.opt[paste0('per',j)])#day
        ps <- exp(mc[,paste0('per',j)])
    }
    Tp.map <- out$tref-par.opt[paste0('Mo',j)]/(2*pi)*p
    Tps <- out$tref-mc[,paste0('Mo',j)]/(2*pi)*ps-Tp.map
    options(digits=13)
    if(sd(Tps)<1e-10) Tps <- rnorm(nrow(mc),Tps[1],1e-10)
    data.distr(x=Tps,xlab=bquote(T['p'*.(j)]-.(format(Tp.map,7.6))~'[day]'),ylab='Frequency',lp=lp,plotf=TRUE)
    options(digits=10)
    if(any(grepl('^K\\d',colnames(mc)))){
        ks <- mc[,paste0('K',j)]
        if(sd(ks)==0) ks <- rnorm(nrow(mc),ks[1],1e-10)
        if(any(names(par.opt)==paste0('Inc',j))){
            Is <- mc[,paste0('Inc',j)]
            m2s <- k2m(K=ks,P=ps,e=mc[,paste0('e',j)],Ms=mstars,Inc=Is)
            f1=data.distr(x=m2s$ms,xlab=bquote(m['c'*.(j)]*'['*M[sun]*']'),ylab='Frequency',lp=lp,plotf=TRUE)
            f2=data.distr(x=m2s$mj,xlab=bquote(m['c'*.(j)]*'['*M[jup]*']'),ylab='Frequency',lp=lp,plotf=TRUE)
                                        #        f3=data.distr(x=m2s$me,xlab=bquote(m['c'*.(j)]*'['*M[earth]*']'),ylab='Frequency',lp=lp,plotf=TRUE)
        }else{
            Is <- acos(uniform(Nsamp,-1,1))
            m2s <- k2m(K=ks,P=ps,e=mc[,paste0('e',j)],Ms=mstars,Inc=Is)$ms
            f1=data.distr(x=m2s$ms,xlab=bquote(m['c'*.(j)]*'['*M[sun]*']'),ylab='Frequency',lp=lp,plotf=TRUE)
            f2=data.distr(x=m2s$mj,xlab=bquote(m['c'*.(j)]*'['*M[jup]*']'),ylab='Frequency',lp=lp,plotf=TRUE)
            f3=data.distr(x=m2s$me,xlab=bquote(m['c'*.(j)]*'['*M[earth]*']'),ylab='Frequency',lp=lp,plotf=TRUE)
            Is <- rep(pi/2,Nsamp)
            m2s <- k2m(K=ks,P=ps,e=mc[,paste0('e',j)],Ms=mstars,Inc=Is)$ms
            f1=data.distr(x=m2s$ms,xlab=bquote(msini['c'*.(j)]*'['*M[sun]*']'),ylab='Frequency',lp=lp,plotf=TRUE)
            f2=data.distr(x=m2s$mj,xlab=bquote(msini['c'*.(j)]*'['*M[jup]*']'),ylab='Frequency',lp=lp,plotf=TRUE)
                                        #        f3=data.distr(x=m2s$me,xlab=bquote(msini['c'*.(j)]*'['*M[earth]*']'),ylab='Frequency',lp=lp,plotf=TRUE)
        }
        acs <- ((mstars+m2s$ms)*(ps/365.25)^2)^(1/3)#au; this assumes that the semi-major axes are calculated separately
        data.distr(x=acs,xlab='a[au]',ylab='Frequency',lp=lp,plotf=TRUE)
    }
}
if(any(names(par.opt)=='dplx')){
    if(sd(mc[,'dra'])==0){
       mc[,c('dra','ddec','dplx','dpmra','dpmdec')] <- mc[,c('dra','ddec','dplx','dpmra','dpmdec')]+rnorm(nrow(mc),0,1e-6)
    }
    data.distr(x=mc[,'dra'],xlab=paste0('(R.A.- ',out$astrometry[out$iref,'ra'],')*cos(DEC) [mas]'),ylab='Frequency',lp=lp,plotf=TRUE)
    data.distr(x=mc[,'ddec'],xlab=paste0('Decl.- ',out$astrometry[out$iref,'dec'],' [mas]'),ylab='Frequency',lp=lp,plotf=TRUE)
    plxs <- out$astrometry[out$iref,'parallax']-mc[,'dplx']
    data.distr(x=plxs,xlab='parallax [mas]',ylab='Frequency',lp=lp,plotf=TRUE)
    data.distr(x=1/plxs,xlab='Distance [kpc]',ylab='Frequency',lp=lp,plotf=TRUE)
    data.distr(x=out$astrometry[out$iref,'pmra']-mc[,'dpmra'],xlab='pmra [mas/yr]',ylab='Frequency',lp=lp,plotf=TRUE)
    data.distr(x=out$astrometry[out$iref,'pmdec']-mc[,'dpmdec'],xlab='pmdec [mas/yr]',ylab='Frequency',lp=lp,plotf=TRUE)
}
mtext("Posterior Distribution for derived parameters", outer=TRUE,  cex=1.5, line=-0.5)
