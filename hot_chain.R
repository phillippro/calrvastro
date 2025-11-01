tem.up <- 15
#ind <- which(!is.na(tab[,3] & !is.na(tab[,3] & tab[,3]!=0)))
if(Niter0<1e5){
    tem.low <- 1e-3
}else{
    tem.low <- 1e-6
}
#if(any(target==c('CPD-632495','TYC3588-11669-1'))) tem.low <- 1e-12
tem <- tem.low
tmp <- run.metropolis.MCMC(startvalue,cov.start,iterations=1e3,tem=tem,bases=bases,RVonly=RVonly)
acceptance <-  (1-mean(duplicated(tmp$out)))*100
if(acceptance<50) tem.low <- 1e-3*tem.low

tem.min <- tem.low
verbose <- FALSE
#####Find the optimal tempering parameter: tem
Ntem <- round(log(1/tem.min)/log(2))#tem.min*2^n=tem; tem<=1
i1 <- 0
#Ntmp <- Ntmp0 <- 1000
Ntmp <- Ntmp0 <- min(1e3,Niter0/10)
#if(target=='CPD-632495') Ntmp <- Ntmp0 <- 1e4
tems <- c()
for(i0 in 0:Ntem){
    tem <- min(tem.min*2^i0,1)
    tems <- c(tems,tem)
    nburn <- floor(Ntmp/2)
    tmp <- run.metropolis.MCMC(startvalue,cov.start,iterations=Ntmp,tem=tem,bases=bases,RVonly=RVonly)
    out.mcmc <- tmp$out
    Npar <- ncol(out.mcmc)-2
    mcmc.out <- out.mcmc[,1:Npar]
    logpost.out <- out.mcmc[,'logpost']
    loglike.out <- out.mcmc[,'loglike']
    acceptance <-  (1-mean(duplicated(mcmc.out)))*100
#    if(Niter0<1e5){
    if(TRUE){
        cat('Adaptive tempering burning: tem=',format(tem,digit=3),'\n')
        cat('maximum likelihood=',max(loglike.out),'\n')
        cat('acceptance:',acceptance,'\n')
    }
    if(acceptance<(tem.up+5)){
        ##        cat('Finding the optimal tempering with longer chain!\n')
        Ntmp <- Ntmp0*10
        if(i1>0 & acceptance<(tem.up-5)) break()
        i1 <- i1+1
    }else{
        Ntmp <- Ntmp0
    }
    startvalue <- mcmc.out[which.max(loglike.out),]
}
####Run MCMC with the optimally tempered chain
cat('\nRun hot chain with tem =',tem,'\n')
tmp <- run.metropolis.MCMC(startvalue,cov.start,iterations=floor(Niter0/Ncores),tem=tem,bases=bases,RVonly=RVonly)
#mcmc.out <- tmp$out[-(1:floor(Niter0/Ncores/2)),]
mcmc.out <- tmp$out#no burning for hot chains
startvalue <- mcmc.out[which.max(mcmc.out[,'loglike']),1:(ncol(mcmc.out)-2)]
####Initial constraint of a signal with cold chain
cat('\nRun short cold chain to constrain the signal\n')
#save(list=ls(all=TRUE),file=paste0('test',ncore,'.Robj'))
tmp <- run.metropolis.MCMC(startvalue,cov.start,iterations=floor(Niter0/Ncores),tem=1,bases=bases,RVonly=RVonly)
out.mcmc <- tmp$out[-(1:floor(Niter0/Ncores/2)),]
Npar <- ncol(out.mcmc)-2
#mcmc <- out.mcmc[,1:Npar]
logpost <- out.mcmc[,Npar+1]
loglike <- out.mcmc[,Npar+2]
par.hot <- out.mcmc[which.max(loglike),1:Npar]
###optional: update intercept parameters
tmp <- loglikelihood(par.hot,prediction=TRUE)
res <- tmp$res
ll0 <- tmp$loglike
par.hot1 <- par.hot
###update RV parameters
if(out$Nrv>0){
    for(k3 in 1:length(out$ins.rv)){
        par.hot1[paste0('b',k3)] <- par.hot1[paste0('b',k3)]+mean(res[[out$ins.rv[k3]]]$res2)
    }
}
ll1 <- loglikelihood(par.hot1)
#cat('ll0=',ll0,'\n')
#cat('ll1=',ll1,'\n')
if(ll1>ll0) par.hot <- par.hot1

###use solutions from Np=1 runs? no!
if(length(per.prim)>0 & FALSE){
    ind <- which(abs(par.hot[1]-per.prim)<0.1 & ll.prim>max(loglike))
    if(length(ind)>0){
        ind.opt <- ind[which.max(ll.prim[ind])]
        tmp <- out.mcmc$sig1[[ind.opt]]
        ind.max <- which.max(tmp[,ncol(tmp)])
        par.hot <- out.mcmc$sig1[[ind.opt]][ind.max,1:Npar]
    }
}
Ps <- exp(par.hot[grep('per',names(par.hot))])
cat('Hot and short-cold chains:P=',paste(round(Ps,2),collapse=','),'day; lnL=',round(max(loglike),2),'\n')
#unlist(sapply(names(par.opt),function(n) cat(n,'=',par.hot[n],',')))
#startvalue <- par.hot

