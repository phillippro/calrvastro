if(chain.type=='adapt'){
    if(Nbin.per>10){
        tem.up <- 10
    }else if(Nbin.per>1){
        tem.up <- 15
    }else{
        tem.up <- 20
    }
    ind <- which(!is.na(tab[,3] & !is.na(tab[,3] & tab[,3]!=0)))
    tem.min <- min(1e-3,10^floor(log10(min(tab[ind,3]/sd(tab[ind,2])))))
    verbose0 <- verbose
    verbose <- FALSE
#####Find the optimal tempering parameter: tem
#    source('setting_up.R')
    source('adapt_temperature.R',local=TRUE)

#    source('adapt_temperature2.R',local=TRUE)
    verbose <- verbose0
####Run MCMC with the optimally tempered chain
    cat('\nRun chain with tem =',tem,'\n')
    Nsamp <- Niter0/Ncores
    amh <- AMH(round(Nsamp/2),Nsamp)
    acceptance  <- (1-mean(duplicated(mcmc.out)))*100
    out <- amh$out
    out.all <- rbind(out.all,out)
    startvalue <- out[which.max(out[,ncol(out)]),1:Npar]
    cat('Chain size:',Nsamp,'\n')
    cat('acceptance:',acceptance,'\n')

####parameter inference with cold chain
#    data <-data0#find singal signal first
#    source('setting_up.R')
    tem <- 1
    cat('\nRun chain with cold chain\n')
    Nsamp <- Niter0/Ncores
    amh <- AMH(round(Nsamp/2),Nsamp)
    acceptance  <- (1-mean(duplicated(mcmc.out)))*100
    out <- amh$out
    out.all <- rbind(out.all,out)
    cat('Chain size:',Nsamp,'\n')
    cat('acceptance:',acceptance,'\n')
}else if(chain.type=='section'){
    out.section <- foreach(n = 1:Ncores, .combine='rbind') %dopar% {
        nbin.per <- (nbin.per0-1)*Ncores+n
        source('prepare_par.R',local=TRUE)#should always set local=TRUE
        source('mcmc_func.R',local=TRUE)#This is necessary because the range of period should be renewed
        tmp <- run.metropolis.MCMC(startvalue,cov.start,Niter,Mstar=Mstar)
        Nburn <- round(nrow(tmp$chain)/2)
        cbind(tmp$chain[-(1:Nburn),],tmp$post[-(1:Nburn)],tmp$like[-(1:Nburn)])
    }
    out <- out.section
}else if(chain.type=='parallel'){
    cat('\nRun chain with tem =',tem,'\n')
    amh <- AMH(round(Nsamp/2),Nsamp)
    cat('Chain size:',Nsamp,'\n')
    out <- amh$out
}else if(chain.type=='normal'){
    out <- run.metropolis.MCMC(startvalue,cov.start,Niter,Mstar=Mstar)
    Nburn <- round(Niter/2)
    mcmc.out <- out$chain[-(1:Nburn),]
    post.out <- out$post[-(1:Nburn)]
    loglike.out <- out$like[-(1:Nburn)]
}else if(chain.type=='talk'){
    tem <- 1
    likes <- c()
    posts <- c()
    Ps <- c()
    Nstep <- 100
    Ntry <- round(max(min(Niter/Nstep/Ncores,1000),10))
    plot.figure <- FALSE#TRUE
    Nbins <- c()
    Nbin.per <- 1
    nbin.per <- 1
    if(exists('plows')){
        rm(plows)
    }
    source('prepare_par.R',local=TRUE)
    source('mcmc_func.R',local=TRUE)
    mcmc.out <- c()
    post.out <- c()
    loglike.out <- c()
    likes <- c()
    posts <- c()
    cat('Ntry=',Ntry,'\n')
    for(i3 in 1:Ntry){
        tmp <- AMH2(1,Nstep)
        out.amh <- tmp$out
        if(talk.type=='post'){
            ind.max <- which.max(out.amh[,Npar+1])
        }else{
            ind.max <- which.max(out.amh[,Npar+2])
        }
        startvalue <- par.opt <- out.amh[ind.max,1:Npar]
        if(Ntry>100){
            nverbose <- round(Ntry/10)
        }else{
            nverbose <- 10
        }
        if(FALSE){
            cat('i3=',i3,'\n')
            cat('optimal period=',format(exp(startvalue[(0:(Np-1))*Nkeppar+1]),digit=5),'\n')
            cat('max logpost=',format(max(out.amh[,Npar+1]),digit=8),'\n')
            cat('max loglike=',format(max(out.amh[,Npar+2]),digit=8),'\n')
        }
        cov.amh <- tmp$cov
        cov.len <- ncol(cov.amh)
        mcmc.out <- out.amh[,1:Npar]
        post.out <- out.amh[,Npar+1]
        loglike.out <- out.amh[,Npar+2]
        likes <- c(likes,max(loglike.out))
        posts <- c(posts,max(post.out))
        if(i3%%10==0 & i3>1){
            if(max(diff(likes))<0.01){
                break()
            }else{
                likes <- c()
                posts <- c()
            }
        }
    }
#    startvalue <- mcmc.out[which.max(loglike.out),]
    startvalue <- mcmc.out[which.max(post.out),]
    amh <- AMH(round(Nsamp/2),Nsamp)
    out <- amh$out
    cov0 <- amh$cov
    cov0 <- cov0[1:Npar,1:Npar]
}
if(chain.type!='normal'){
    mcmc.out <- out[,1:Npar]
    post.out <- out[,Npar+1]
    loglike.out <- out[,Npar+2]
}
acceptance <- (1-mean(duplicated(mcmc.out)))*100
cat('acceptance percentage:',acceptance,'\n')
