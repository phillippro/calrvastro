if(!exists('par.data')){
    tmin <- min(trv)
    tmax<- max(trv)
    rmin <- min(RV)
    rmax <- max(RV)
    Ndata <- length(trv)
    small.err <- FALSE
    ids <- id
    noise.models <- noise.model
    par.data <- list()
    par.data <- lapply(1:Nw,function(i) par.data[[ids[i]]] <- list())
    names(par.data) <- ids
    vars <- unique(c('ID','trv','t0','times0','RV','eRV','mom.name','mom','epoch.type','epochs','Nes','Nepoch','nepoch','ep.par','ep.name','noise.model','p','q','data.type','Ntj','gp.type','folder','arma.type','instrument','target','star.name'))
    for(k in 1:length(vars)){
        if(exists(vars[k])){
            par.data[[id]][[vars[k]]] <- eval(parse(text = vars[k]))
####remove variables of previous run
        }
    }
}
ts <- seq(tmin,tmax,by=0.1)
Nsep <- max(round(Niter/10000),1)
fchop <- 0#0.9
Nbins = max(min(round(nrow(mcmc.out)/1e2),1000),10)
Nbreaks = 100#for histograms
#Npeaks = 0
Nc2p <- 0
#binning the parameters and corresponding posteriors and likelihoods to show the posterior distributions
par.bin <- post.bin <- like.bin <- array(data=NA,dim=c(Nbins,ncol(mcmc.out)))
for(i in 1:ncol(mcmc.out)){
if(sd(mcmc.out[,1])!=0){
    tmp <- binning.post2(mcmc.out[,i],post.out,loglike.out,Nbins)
    par.bin[,i] <- tmp$pars
    like.bin[,i] <- tmp$likes
    post.bin[,i] <- tmp$posts
      if(length(tmp$ind.na)>0){
        ind.na <- tmp$ind.na
        #cat('the',ind.na,'th bin of ',colnames(mcmc.out)[i],' has no samples!\n')
        like.bin[ind.na,i] <- min(like.bin[-ind.na,i])
        post.bin[ind.na,i] <- min(post.bin[-ind.na,i])
      }
    }else{
      par.bin[,i] <- mcmc.out[1,i]
      like.bin[,i] <- loglike.out[1]
      post.bin[,i] <- post.out[1]
    }
}
###optimal period
#if(talk.type=='like'){
#ind.opt <- which.max(loglike.out)
#}else{
    ind.opt <- which.max(post.out)
#}
if(length(ind.opt)==0) ind.opt <- 1
par.opt <- mcmc.out[ind.opt,]

ind.opt <- which.max(post.out)
if(length(ind.opt)==0) ind.opt <- 1
par.opt.post <- mcmc.out[ind.opt,]

ind.opt <- which.max(loglike.out)
if(length(ind.opt)==0) ind.opt <- 1
par.opt.like <- mcmc.out[ind.opt,]

per.mul <- par.opt.post[grep('per([[:digit:]]{1})',names(par.opt.post))]
if(period.par=='nv'){
   Popt <- 1/per.mul
}else if(period.par=='logP'){
   Popt <- exp(per.mul)
}else{
   Popt <- per.mul
}
#########report results
cat('optimal period:',Popt,'day\n')
#cat('binning time consumed: ', format((proc.time()-t1)[3],digits=3),'s \n')
ind.per <- c()
if(Np>0){
    ind.per <- grep('per',colnames(mcmc.out))
    per.bin <- as.matrix(par.bin[,ind.per],ncol=length(ind.per))
    if(period.par=='nv'){
        nv.out <- nv <- as.matrix(mcmc.out[,ind.per],ncol=length(ind.per))
        P.out <- P <- 1/nv
        logP.out <- logP <- log(P)
        nv.bin <- per.bin
        P.bin <- 1/nv.bin
        logP.bin <- log(P.bin)
    }else if(period.par=='P'){
        P.out <- P <- as.matrix(mcmc.out[,ind.per],ncol=length(ind.per))
        nv.out <- nv <- 1/P
        logP.out <- logP <- log(P)
        P.bin <- per.bin
        nv.bin <- 1/P.bin
        logP.bin <- log(P.bin)
    }else if(period.par=='logP'){
        logP.out <- logP <- as.matrix(mcmc.out[,ind.per],ncol=length(ind.per))
        P.out <- P <- exp(logP)
        nv.out <- nv <- 1/P
        logP.bin <- per.bin
        P.bin <- exp(logP.bin)
        nv.bin <- 1/P.bin
    }
}
tau.bin <- array(data=NA,dim=c(Nbins,length(nepoch),Nw))
#if(any(grepl('ARMA',noise.models))){
#    qs <- as.integer(gsub('ARMA','',noise.models))
#}else{
#    qs <- 0
#}
for(jj in 1:Nw){
    ind.beta <- grep(paste0('beta',jj),colnames(mcmc.out))
    if(grepl('ARMA',noise.models[jj]) & qs[jj]>0){
        beta.bin <- par.bin[,ind.beta]
        tau.bin[,,jj] <- exp(beta.bin)
    }
}
##########thinning and chopping
index.mcmc <- 1:nrow(mcmc.out)
###thinning
if(Nsep>1){
    ind.out <- Nsep*(1:floor(nrow(mcmc.out)/Nsep))
    mcmc.out1 <- mcmc.out[ind.out,]
    post.out1 <- post.out[ind.out]
    loglike.out1 <- loglike.out[ind.out]
    if(Np>0){
        P.out <- as.matrix(P[ind.out,])
        nv.out <- as.matrix(nv[ind.out,])
        logP.out <- as.matrix(logP[ind.out,])
    }
    index.mcmc <- ind.out
}else{
    mcmc.out1 <- mcmc.out
    post.out1 <- post.out
    loglike.out1 <- loglike.out
    index.mcmc <- 1:length(post.out)
}

Nsubsamp = floor(nrow(mcmc.out1)*(1-fchop))
if(fchop>0){
    len.mcmc <- nrow(mcmc.out1)
    ind.out <- (len.mcmc-Nsubsamp+1):len.mcmc
    mcmc.out1 <- mcmc.out1[ind.out,]
    post.out1 <- mcmc.out1[ind.out]
    loglike.out1 <- loglike.out1[ind.out]
    P.out <- as.matrix(P[ind.out,])
    nv.out <- as.matrix(nv[ind.out,])
    logP.out <- as.matrix(logP[ind.out,])
    index.mcmc <- index.mcmc[ind.out]
}
if(!exists('activity')) activity <- FALSE
#######prepare phase-folded graph
RV.noise <- RV.kepler(pars.kep=par.opt.post,Np.kep = Np,noise.only=TRUE)#RV2= RV-RV.noise
RV.kep <- RV.kepler(pars.kep=par.opt.post)
for(k1 in 1:Nw){
    var <- names(par.data[[ids[k1]]])
    for(k in 1:length(var)){
        assign(var[k],par.data[[ids[k1]]][[var[k]]])
    }
    par.data[[k1]]$RV2 <- par.data[[k1]]$RV-RV.noise[[k1]]
    if(noise.model=='ARMA'){
        tmp <- arma(par.opt.post,rv.kep=RV.kep[[k1]],j3=k1,x=par.data[[k1]]$RV,t=par.data[[k1]]$trv)
        par.data[[k1]]$RV2 <- par.data[[k1]]$RV2-tmp
    }
    if(noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='PSID' | noise.model=='ARMAPSID' & TJAR.mode=='RV'){
        par.can <- par.opt.post[grep(paste0(k1,'$'),names(par.opt.post))]
        phis.opt <- par.can[grep('phis([[:digit:]]{1})',names(par.can))]
        alphas.opt <- par.can['alphas']
        if(exists('BIS') & exists('FWHM') & TJAR.AB){
            phia.opt <- par.can[grep('phia([[:digit:]]{1})',names(par.can))]
            alphaa.opt <- par.can['alphaa']
            phib.opt <- par.can[grep('phib([[:digit:]]{1})',names(par.can))]
            alphab.opt <- par.can['alphab']
        }
        if(TJAR.sym=='asym'){
            psis.opt <- par.can[grep('psis([[:digit:]]{1})',names(par.can))]
           #RV2 <- RV2-tjar(phiars=phis.opt,psiars=psis.opt,alphaar=alphas.opt,par.tj=par.can['ss'])
            par.data[[k1]]$RV2 <- par.data[[k1]]$RV2-tjar(phiars=phis.opt,psi.sab=psis.opt,alpha.sab=alphas.opt,par.tj=1)
            if(exists('BIS') & exists('FWHM') & TJAR.AB){
                psia.opt <- par.can[grep('psia([[:digit:]]{1})',names(par.can))]
                psib.opt <- par.can[grep('psib([[:digit:]]{1})',names(par.can))]
                par.data[[k1]]$RV2 <- par.data[[k1]]$RV2-tjar(phi.sab=phia.opt,psi.sab=psia.opt,alpha.sab=alphaa.opt,par.tj=1)
                par.data[[k1]]$RV2 <- par.data[[k1]]$RV2-tjar(phi.sab=phib.opt,psi.sab=psib.opt,alpha.sab=alphab.opt,par.tj=1)
            }
        }else{
            par.data[[k1]]$RV2 <- par.data[[k1]]$RV2-tjar(phi.sab=phis.opt,psi.sab=phis.opt,alpha.sab=alphas.opt,par.tj=1)
            if(exists('BIS') & exists('FWHM') & TJAR.AB){
                par.data[[k1]]$RV2 <- par.data[[k1]]$RV2-tjar(phi.sab=phia.opt,psi.sab=phia.opt,alpha.sab=alphaa.opt,par.tj=1)
                par.data[[k1]]$RV2 <- par.data[[k1]]$RV2-tjar(phi.sab=phib.opt,psi.sab=phib.opt,alpha.sab=alphab.opt,par.tj=1)
            }
        }
    }
}
ylim <- c(rmin,rmax)
#######residuals
for(j4 in 1:Nw){
    par.data[[j4]]$res <- par.data[[j4]]$RV-RV.kep[[j4]]
    if(noise.model=='ARMA'){
        par.data[[j4]]$res <- par.data[[j4]]$res-arma(par.opt.post,j3=j4,rv.kep=RV.kep[[j4]],x=par.data[[ids[j4]]]$RV,t=par.data[[j4]]$trv)
    }
}
#####################
if(Np>0){
    for(k1 in 1:Nw){
        var <- names(par.data[[ids[k1]]])
        for(k in 1:length(var)){
            assign(var[k],par.data[[ids[k1]]][[var[k]]])
        }
        RV5c <- RV3c <- RV4c <- eRV4c <- t3c <- t4c <- t5c <- c()
        rv.kep <- RV.kepler(pars.kep=par.opt.post,kep.only=TRUE)[[1]]
        for(k3 in 1:Np){
            par.other <- par.target <- par.opt.post
            indk <- grep('K([[:digit:]]{1})',names(par.opt.post))
            if(Np>1){
                par.other[indk[k3]] <- 0
                par.target[indk[-k3]] <- 0
                RV.kep.only <- RV.kepler(pars.kep=par.other,Np.kep = Np,kep.only=TRUE)
                RV3 <- RV2-RV.kep.only[[k1]]
            }else{
                par.other[indk[k3]] <- 0
                RV3 <- RV2
            }
            t3 <- (trv-out$tref)%%Popt[k3]#all time series are combined with respect to the same reference time
####binning the phase folded data
            ind.sort <- sort(t3,index.return=TRUE)$ix
            t3.sort <- t3[ind.sort]
            data.bin <- wtb(t=t3.sort,x=RV3[ind.sort],ex=eRV[ind.sort],dt=Popt[k3]/20)
            t4 <- data.bin[,1]
            RV4 <- data.bin[,2]
            eRV4 <- data.bin[,3]
####model prediction
            t5 <- (ts-out$tref)%%Popt[k3]
            RV5 <- RV.kepler(pars.kep=par.target,tt=ts,kep.only=TRUE)[[1]]
###plot
            ylim <- range(RV4,mean(RV3)+sd(RV3),mean(RV3)-sd(RV3),RV5)

            t3c <- cbind(t3c,t3)
            RV3c <- cbind(RV3c,RV3)
            t4c <- combine.index(t4c,t4)
            RV4c <- combine.index(RV4c,RV4)
            eRV4c <- combine.index(eRV4c,eRV4)
            t5c <- cbind(t5c,t5)
            RV5c <- cbind(RV5c,RV5)
        }
        if(k1==1){
            datac <- list(cbind(t3c,RV3c))
            datas <- list(cbind(t5c,RV5c))
            datab <- list(cbind(t4c,RV4c,eRV4c))
	    datak <- list(rv.kep)
        }else if(k1<=Nw){
            datac <- c(datac,list(cbind(t3c,RV3c)))
            datas <- c(datas,list(cbind(t5c,RV5c)))
            datab <- c(datab,list(cbind(t4c,RV4c,eRV4c)))
            datak <- c(datak,list(rv.kep))
        }
    }
}
