#NN <- length(out[['mcmc.opt']])
source('convergence_test.R')
#NN <- min(Nsig,Nmax)+1
NN <- Nsig+1

##calculate different types of residuals
res.comb.sigtrend <- res.comb.sig <- res.comb.all <- res.sigtrend <- res.all <- res.sig <- Popts <- par.stat <- list()
for(j in (Nmin+1):NN){
    n <- paste0('sig',j-1)
    if(!save.memory){
        Nchain <- length(out[['mcmc.all']][[n]])
    }else{
        Nchain <- 1
    }
    npar <- ncol(out[['mcmc.opt']][[n]])-2
    par.name <- colnames(out[['mcmc.opt']][[n]])[1:npar]

    mc.opt <- out[['mcmc.opt']][[n]]
    if(sig.type=='likelihood'){
        lp <- lp.opt <- out[['mcmc.opt']][[n]][,npar+2]
        if(!save.memory){
            lp <- sapply(ind.chain,function(k) out[['mcmc.all']][[n]][[k]][,npar+2])
        }
    }else{
        lp <- lp.opt <- out[['mcmc.opt']][[n]][,npar+1]
        if(!save.memory){
            lp <-sapply(ind.chain,function(k) out[['mcmc.all']][[n]][[k]][,npar+1])
        }
    }
    tmp <- c()
    for(i in 1:npar){
        if(!save.memory){
            mm <- as.numeric(sapply(ind.chain,function(k) out[['mcmc.all']][[n]][[k]][,i]))
        }else{
            mm <- out[['mcmc.opt']][[n]][,i]
        }
        tmp <- cbind(tmp,data.distr(mm,lp,plotf=FALSE))
    }
    colnames(tmp) <- par.name
    par.stat[[n]] <- tmp
    if(!any(names(out)=='Nrv')) out$Nrv <- 1
    if(out$Nrv>0){
        t0 <- combine.list(lapply(out$ins.rv,function(i) out[[i]]$RV[,1]))
        index <- sort(t0,index.return=TRUE)$ix
        res <- loglikelihood(par.opt,predict=TRUE)$res
        yy0 <- lapply(out$ins.rv,function(k) res[[k]]$res1)#res1: data-kep
        yy1 <- lapply(out$ins.rv,function(k) res[[k]]$res1-res[[k]]$res2)#res2: data-model; res1-res2=trend+red
        yy2 <- lapply(out$ins.rv,function(k) res[[k]]$res2)#res2: data-model
#        yy3 <- lapply(out$ins.rv,function(k) res[[k]]$res1-res[[k]]$res2+res[[k]]$red)#subtract signal+trend
	yy3 <- list()
	for(k in out$ins.rv){
	    if(length(res[[k]]$red)>0){
	        yy3[[k]] <- res[[k]]$res1-res[[k]]$res2+res[[k]]$red
            }else{
	        yy3[[k]] <- res[[k]]$res1-res[[k]]$res2
	    }
        }
        y0 <- combine.list(yy0)
        y2 <- combine.list(yy2)
        y3 <- combine.list(yy3)
        res.comb.sig[[n]] <- y0[index]
        res.comb.all[[n]] <- y2[index]
        res.comb.sigtrend[[n]] <- y3[index]
        names(yy0) <- names(yy1) <- names(yy2) <- names(yy3) <- out$ins.rv
        res.all[[n]] <- yy2
        res.sigtrend[[n]] <- yy3
        res.sig[[n]] <- yy0
    }

#    if(j>min(1,Nmax-1)){
    if(j>1){
        ns <- colnames(mc.opt)
        ind <- grep('per',ns)
        tmp <- exp(mc.opt[which.max(lp.opt),ind])
        names(tmp) <- gsub('per','Pd',names(tmp))
        Popts[[n]] <- extract.par(par.opt,bases=bases)$P
    }
}
out[['Popt']] <- Popts
out[['par.stat']] <- par.stat
out[['res.sig']] <- res.sig
out[['res.sigtrend']] <- res.sigtrend
out[['res.all']] <- res.all
out[['res.comb.sig']] <- res.comb.sig
out[['res.comb.sigtrend']] <- res.comb.sigtrend
out[['res.comb.all']] <- res.comb.all

