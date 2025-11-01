###This routine is used to conduct adaptive MCMCs to draw posteriors of models for parameter inference

###Helper to lazily load expensive dependencies inside each worker only once
prepare_parallel_worker <- function(){
    if(!exists('.dram_worker_ready', envir = .GlobalEnv, inherits = FALSE) ||
       !isTRUE(get('.dram_worker_ready', envir = .GlobalEnv))){
        source('mcmc_func.R')
        assign('.dram_worker_ready', TRUE, envir = .GlobalEnv)
    }
}

###Make sure the current (master) session also has the dependencies ready
prepare_parallel_worker()

nii <- FALSE
Pmin0 <- Pmin
Pmax0 <- Pmax
#Npar0 <- Npar
logdP <- (log(Pmax0)-log(Pmin0))/Ncores
mcmc.all <- list()
mcmc.opt <- list()
kep.up <- kep.low <- kep.opt <- c()
Nsig <- 0
basis0 <- basis
bases <- c(rep(basis,Ntr),rep(basis2,Nmax-Ntr))
ind.transit0 <- ind.transit
ll.prim <- per.prim <- c()
ll0 <- NA
data.astrometry0 <- out$astrometry
for(np in Nmin:Nmax){
###Parallel signal identification and constraints
    if(!nii){
        fe <- foreach(ncore=1:Ncores,.errorhandling = 'pass',.inorder = FALSE,
                      .export = c('prepare_parallel_worker'))
        backend <- foreach::getDoParName()
        if(backend %in% c('doParallel', 'doParallelSNOW', 'doSNOW')){
            fe <- fe$.options.snow(list(preschedule = FALSE))
        }else if(backend %in% c('doParallelMC', 'doMC')){
            fe <- fe$.options.multicore(list(preschedule = FALSE))
        }
        mcmc <- fe %dopar% {
#        mcmc <- foreach(ncore=1:Ncores,.errorhandling = 'pass') %dopar% {
#        mcmc <- list()
#        for(ncore in 1){
            prepare_parallel_worker()
            replacePar <- FALSE
            if(np>0){
                if(np>Ntr0){
                    Ntr <- 0
                    ind.transit <- 0
                }else{
                    Ntr <- Ntr0
                    ind.transit <- ind.transit0
                }
                Np <- 1
                if(np>Ntr0){
                    basis <- basis2
                }else{
                    basis <- basis0
                }

####hot chain
                cat('\n Find signal',np,'using hot chain!\n')
                                        #            Pmin <- exp(log(Pmin0)+logdP*(ncore-1))
                                        #            Pmax <- exp(log(Pmin0)+logdP*ncore)
                bases0 <- bases
                bases <- basis

###special case
                if(Nmin==Nmax){
                    Np <- Nmin
                    bases <- bases0
                }
###get initial values for free parameters
                source('initial_condition.R',local=TRUE)
####initial MCMC constraints of a signal; output posterior sample and par.opt
                                        #            RVonly <- TRUE
                RVonly <- FALSE
                if(file.exists(fpar)) source('reload_optpar.R',local=TRUE)
                par.hot <- startvalue
                if(!is.null(par.fix)){
                    if(any(par.fix=='Mstar')) par.hot['Mstar'] <- out$Mstar
                }
                if(any(startvalue<par.min | startvalue>par.max)){
                    cat('The following initial values exceed the prior boundary:',names(startvalue)[startvalue<par.min | startvalue>par.max],'\n')
                }
                if(Niter0>=1e3)   source('hot_chain.R',local=TRUE)
                bases <- bases0
            }
####cold chain further constraint
            cat('\n Constrian signal',np,'using cold chain!\n')
            Ntr <- Ntr0
            ind.transit <- ind.transit0
            Np <- np
            Pmin <- Pmin0
            Pmax <- Pmax0
            if(out$Nrv>0){
                for(i in out$ins.rv){
                    out[[i]]$RV[,2] <- out[[i]]$data[,2]
                }
            }
            if(out$Nastro>0) out$astrometry <- data.astrometry0
            source('initial_condition.R',local=TRUE)
            if(Np>1 & Nmin<=1){
                startvalue <- as.numeric(c(kep.opt,par.hot))
                startvalue <- assign.names(startvalue,Np=Np,p=ps,q=qs,n=ns,bases=bases)
            }else{
                startvalue <- par.hot
            }
            tmp <- run.metropolis.MCMC(startvalue,cov.start,iterations=min(1e7,Niter0*max(Nsig,1)),out1=out,tem=1,bases=bases)
###burnin
            if(Niter0>1e5){
                return(tmp$out[-(1:floor(nrow(tmp$out)/2)),])
            }else{
                return(tmp$out)
            }
        }
    }else{
        cat('nii!\n')
        source('initial_condition.R',local=TRUE)
    }
###remove bad parallel chains
    ind <- which(sapply(1:length(mcmc),function(k) is.null(dim(mcmc[[k]]))))
    if(length(ind)>0) mcmc <- mcmc[-ind]

###Save mcmc results
    if(!save.memory){
        mcmc.all[[paste0('sig',np)]] <- mcmc
    }
    llmax <- sapply(1:length(mcmc), function(i) max(mcmc[[i]][,ncol(mcmc[[i]])]))
    cat('Maximum likelihood for all chains:',llmax,'\n')
    mcopt <- mcmc[[which.max(llmax)]]
    mcmc.opt[[paste0('sig',np)]] <- mcopt
###save primary signals
    if(np==1){
        per.prim <- sapply(1:length(mcmc), function(i) mcmc[[i]][which.max(mcmc[[i]][,ncol(mcmc[[i]])-1]),1])
        ll.prim <-llmax
    }
###subtract the best-fit in order to find the next signal
    Npar <- ncol(mcopt)-2
    par.opt <- mcopt[which.max(mcopt[,'loglike']),1:Npar]
    if(np>0) cat('optimal period:',extract.par(par.opt,bases=bases)$P,'days\n')
    if(np>0){
        kep.opt <- par.opt[1:(np*Nkeppar)]
        res <- loglikelihood(par.opt,predict=TRUE)$res
        if(out$Nrv>0){
            for(i in out$ins.rv){
                out[[i]]$RV[,2] <- as.numeric(res[[i]]$res1)
            }
        }
    }

###whether to stop DRAM according to BF criterion
    ll <-  max(mcopt[,ncol(mcopt)])
    if(np>0){
        if(is.na(ll0)) ll0 <- ll-1e2
        lnbf3 <- ll-ll0-1.5*log(length(unlist(res)))
        cat('lnbf3=',lnbf3,'\n')
        if(FALSE){
             break()
        }else{
            Nsig <- Nsig+1
        }
    }
    ll0 <- ll
}
if(Nmin==Nmax) Nsig <- Nmin
###change the data back into raw data
if(out$Nrv>0){
    for(i in out$ins.rv){
        out[[i]]$RV[,2] <- out[[i]]$data[,2]
    }
}
if(any(names(out)=='astrometry')){
    out$astrometry <- data.astrometry0
}
mcopt <- mcmc.opt[[paste0('sig',Nsig)]]
Npar <- ncol(mcopt)-2
par.opt <- mcopt[which.max(mcopt[,'logpost']),1:Npar]
if(Nsig>0){
    Popt <- extract.par(par.opt,out1=out,bases=bases)$P
}else{
    Popt <- NA
}
if(!save.memory){
    out[['mcmc.all']] <- mcmc.all
}
out[['mcmc.opt']] <- mcmc.opt
out$tauT <- RV.kepler(par.opt,out1=out)$tauT
