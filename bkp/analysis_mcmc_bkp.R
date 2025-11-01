NN <- length(out[['mcmc.opt']])
res.comb.sig <- res.comb.all <- res.all <- res.sig <- Popt <- par.stat <- list()
t0 <- combine.list(lapply(ins,function(i) out[[i]]$RV[,1]))
index <- sort(t0,index.return=TRUE)$ix
for(j in 1:NN){
    tmp <- out[['mcmc.opt']][[j]]
    npar <- ncol(tmp)-2
    ll <- tmp[,npar+2]
    lp <- tmp[,npar+1]
    mc <- tmp[,1:npar]
    tmp <- sapply(1:npar,function(i) data.distr(mc[,i],lp,plotf=FALSE))
    colnames(tmp) <- colnames(mc)[1:npar]
    par.stat[[paste0('sig',j-1)]] <-tmp
    res.sig[[paste0('sig',j-1)]] <- cal.residual(mc[which.max(lp),])$res.sig
    res.all[[paste0('sig',j-1)]] <- cal.residual(mc[which.max(lp),])$res.all

    yy <- lapply(ins,function(k) res.sig[[paste0('sig',j-1)]][[k]])
    y0 <- combine.list(yy)
    y1 <- combine.list(res.all[[paste0('sig',j-1)]])
    res.comb.sig[[paste0('sig',j-1)]] <- y0[index]
    res.comb.all[[paste0('sig',j-1)]] <- y1[index]
    if(j>1){
        ns <- colnames(mc)
        ind <- grep('per',ns)
        tmp <- exp(mc[which.max(lp),ind])
        names(tmp) <- gsub('per','Pd',names(tmp))
        Popt[[paste0('sig',j-1)]] <- tmp
    }
}
out[['Popt']] <- Popt
out[['par.stat']] <- par.stat
out[['res.sig']] <- res.sig
out[['res.all']] <- res.all
out[['res.comb.sig']] <- res.comb.sig
out[['res.comb.all']] <- res.comb.all
