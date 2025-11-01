#plotf <- TRUE
plotf <- FALSE
if(plotf){
pdf('test.pdf')
par(mfrow=c(2,2))
}else{
par(mfrow=c(4,4))
}
panel <- 0
#Pkep <- out$Popt[[paste0('sig',Nsig)]]
Pkep <- par.opt[grep('per',names(par.opt))]
if(!any(names(out$par.stat)==paste0('sig',Nsig)) & Nmax>Nmin){
    par.stat <- matrix.distr(mc[,1:(ncol(mc)-2)],mc[,ncol(mc)],plotf=FALSE)
    out$par.stat[[paste0('sig',Nsig)]] <- par.stat
#    Pkep <- Popt <- extract.par(par.stat['xopt',],Nsig,bases)
}
#par.opt <- out$par.stat[[paste0('sig',Nsig)]]['xopt',]
#Popt <- out$Popt[[paste0('sig',Nsig)]]
Popt <- extract.par(par.opt,bases=bases)$P
#fit.type <- 'multiple'#or multiple
fit.type <- 'multiple'
tcol.percentage <- 80
##calculate residual after subtraction of single signal
if(any(grepl('Omega',names(par.opt)))){
    Npar <- 7
}else if(grepl('apm',names(par.opt))){
    Npar <- 8
}else if(prior.type=='e0'){
    Npar <- 3
}
##combined fit
if(out$Nrv>0){
    tmp <- cal.residual(par.opt,out1=out,bases=bases)
    rvmn <- tmp$res.noise# rv-noise
    rvma <- tmp$res.all#rv-signal-noise
    rvms <- tmp$res.sig#rv-signal
    rv.arma <- tmp$rv.arma#ARMA rv
    rv.trend <- tmp$rv.trend#trend rv
    y0 <- rvmn
    t0 <- combine.list(lapply(out$ins.rv,function(i) out[[i]]$RV[,1]))#only consider primary RVs
    index <- sort(t0,index.return=TRUE)$ix
    rvmn.comb <- y0[index]#
    for(j3 in 1:Nsig){
        psig <- Popt[j3]
        bin.raw <- bin.data <- list()

#        par1 <- rep(1e-10,length(par.opt))
#        par1[(j3-1)*Npar+1:Npar] <- par.opt[(j3-1)*Npar+1:Npar]
#        names(par1) <- names(par.opt)
#        if(any(names(par.opt)=='Mstar')) par1['Mstar'] <- par.opt['Mstar']
        par1 <- par.opt
        if(Nsig>1){
            nps <- (1:Nsig)[-j3]
            par1[paste0('K',nps)] <- 1e-6
        }
        tmp <- cal.residual(par1,bases=bases)# residual
        rv.sig <- tmp$rv.kep
        rv.res <- tmp$res.all
        for(ins in out$ins.rv){
            inds <- out$ind.all$rv[[ins]]
            rv.sig[inds] <- rv.sig[inds]+rvma[inds]
        }
        y0 <- rv.sig
        res0 <- rv.res
        rv.sig.comb <- y0[index]#form a residual combined set
        rv.res.comb <- res0[index]

###calculate residual after subtraction of all signals
                                        #    par.sig <- rep(0,length(par.opt))
        par.sig <- par1
                                        #    names(par.sig) <- names(par.opt)
                                        #    if(!exists('Kmin')) Kmin <- 1e-6
                                        #    if(basis=='natural'){
                                        #        par.sig[(1:Nsig)*5-3] <- Kmin
                                        #    }else{
                                        #        par.sig[(1:Nsig)*5-3] <- log(Kmin)
                                        #    }
                                        #    par.sig[j3*5-3] <- par.opt[j3*5-3]
                                        #    names(par.sig) <- names(par.opt)
        simulate <- RV.kepler(pars.kep=par.sig,tt=out.sim$BJDtdb,out1=out.sim,kep.only=TRUE,bases=bases)
        rv.sim.all <- simulate$rv
        rv.sim <- simulate$rvkep
        t1 <- (out.sim$BJDtdb-out$tref)%%psig
        ix <- sort(t1,index.return=TRUE)$ix
        t2 <- t1[ix]

###single signal fit to individual data set
        for(i in out$ins.rv){
            inds <- out$ind.all$rv[[i]]
            tt <- (out[[i]]$RV[,1]-out$tref)%%psig#day
            rv <- rv.sig[inds]
            erv <- out[[i]]$RV[,3]
            rv2 <- rv.sim[ix]
###add offset rv0_xx
#            rv2 <- ((1+rv2/CMPS)/(1+par.opt[paste0('rv0_',i)]/CKMPS)-1)*CMPS
            ##ylim=range(rv.sig[[i]],rv2)
            ylim <- range(rv.sig)
            plot(tt,rv,xlab=paste0('Orbital Phase [days]'),ylab='RV [m/s]',main=paste0('signal',j3,': P=',round(psig,1),'d; ',i),col='white',xlim=range(tt),ylim=ylim)
            if(exists('Tp')){
                Tp1 <- tmax-(tmax-Tp)%%Popt[1]
                Tp1 <- Tp1+(0:10)*Popt[1]
                Tp2 <- (Tp1-out$tref)%%psig
                abline(v=Tp2,col='blue')
            }
            panel <- panel+1
            if(panel==1){
                mtext("Phase Plot", outer=TRUE,  cex=1.5, line=-0.5)
            }
                                        #        points(t1[ix],rv.sim[ix],col='red',pch='.')
            points(tt,rv.sig[inds],col='grey')
            lines(t2,rv2,col='red')
            ##binnig data
            data.bin <- bin.simple(cbind(tt,rv,erv),10)
            bin.data[[i]] <- data.bin
            raw.bin <- bin.simple(cbind(out[[i]]$RV[,1],rv,erv),10*ceiling((tmax-tmin)/psig))
            bin.raw[[i]] <- raw.bin
            points(data.bin[,1],data.bin[,2],col='black',pch=20)
            try(arrows(data.bin[,1],data.bin[,2]+data.bin[,3],data.bin[,1],data.bin[,2]-data.bin[,3],length=0.05,angle=90,code=3,col='black',pch=20),TRUE)
        }
###single signal fit to combined data
#        ylim <- range(rv.sig.comb,rv2)
        ylim <- range(rv.sig.comb)
        plot((trv.all-out$tref)%%psig,rv.sig.comb,xlab=paste0('Orbital Phase [days]'),ylab='RV [m/s]',col='white',ylim=ylim,main=paste0('RMS=',round(sd(rv.sig.comb),2),'m/s'))
        for(j in 1:length(out$ins.rv)){
            i <- out$ins.rv[j]
            inds <- out$ind.all$rv[[i]]
            points( (out[[i]]$RV[,1]-out$tref)%%psig,rv.sig[inds],col=tcol(colors[j],tcol.percentage))
            data.bin <- bin.data[[i]]
            points(data.bin[,1],data.bin[,2],col=colors[j],pch=20)
            try(arrows(data.bin[,1],data.bin[,2]+data.bin[,3],data.bin[,1],data.bin[,2]-data.bin[,3],length=0.05,angle=90,code=3,col=colors[j]),TRUE)
        }
        lines(t2,rv2,col='red')
        plot((trv.all-out$tref)%%psig,rv.res.comb,xlab=paste0('Orbital Phase [days]'),ylab='RV (data-model) [m/s]',col='white',main=paste0('RMS=',round(sd(rv.res.comb),2),'m/s'))
        for(j in 1:length(out$ins.rv)){
            i <- out$ins.rv[j]
            points( (out[[i]]$RV[,1]-out$tref)%%psig,rv.res[out$ind.all$rv[[i]]],col=tcol(colors[j],tcol.percentage))
        }

                                        #    legend('topright',legend=ins,xpd=NA,inset=c(-0.45,0),col=colors[1:length(ins)],pch=1,bty='n')
        n1 <- ceiling(length(out$ins.rv)/2)
        n2 <- length(out$ins.rv)-n1
        legend('top',legend=out$ins.rv[1:n1],col=colors[1:n1],inset=c(0,-0.2),seg.len = 1,horiz=TRUE,xpd=NA,bty='n',pch=1)
        if(n2>0){
            legend('top',legend=out$ins.rv[(n1+1):length(out$ins.rv)],col=colors[(n1+1):length(out$ins.rv)],seg.len = 1,inset=c(0,-0.1),horiz=TRUE,xpd=NA,bty='n',pch=1)
        }
    }

###combined fit to the whole combined data
#    simulate <- RV.kepler(pars.kep=par.opt,tt=out.sim$BJDtdb,out1=out.sim,kep.only=TRUE,bases=bases)
#    rv.sim <- simulate$rv
#    rv2 <- rv.sim[ix]
    rv0s <- par.opt[grep('rv0',names(par.opt))]
    plot(trv.all,rvmn.comb,xlab=paste0('Time [days]'),ylab='RV [m/s]',col='white',ylim=range(rvmn.comb+max(rv0s),rvmn.comb+min(rv0s),rv.sim.all))
                                        #paste0('signal',j3,': P=',round(psig[j3],1),'d; combined set')
    lines(tsim,rv.sim.all,col='red')
    if(exists('Tp')){
        Tp1 <- tmax-(tmax-Tp)%%Popt[1]-(0:10)*Popt[1]
        abline(v=Tp1,col='blue')
    }

    for(j in 1:length(out$ins.rv)){
        i <- out$ins.rv[j]
        inds <- out$ind.all$rv[[i]]
        y <- ((1+rvmn[inds]/CMPS)*(1+par.opt[paste0('rv0_',i)]/CKMPS)-1)*CMPS
        points(out[[i]]$RV[,1],y,col=tcol(colors[j],tcol.percentage))#
    }
    n1 <- ceiling(length(out$ins.rv)/2)
    n2 <- length(out$ins.rv)-n1
    legend('top',legend=out$ins.rv[1:n1],col=colors[1:n1],inset=c(0,-0.2),seg.len = 1,horiz=TRUE,xpd=NA,bty='n',pch=1)
    if(n2>0){
        legend('top',legend=out$ins.rv[(n1+1):length(out$ins.rv)],col=colors[(n1+1):length(out$ins.rv)],seg.len = 1,inset=c(0,-0.1),horiz=TRUE,xpd=NA,bty='n',pch=1)
    }
}
####plot relativity-induced RVs
if(out$Nrv>0) rvfit <- RV.kepler(par.opt,out1=out)
if(out$relativity & out$Nrv>0){
    ind.rv <- unlist(sapply(out$ins.rv,function(i) which(out$tiall[,'instr']==i & out$tiall[,'type']=='rv')))
    dtt <- out$tiall[ind.rv,1]-trv.all
    res <- out$res.comb.all[[paste0('sig',Nsig)]]
    erv <- out$all[,3]
#    respr <- res+rvfit$rvg[ind.rv]+rvfit$rvs[ind.rv]
#    resps <- res+rvfit$rvs[ind.rv]
#    respg <- res+rvfit$rvg[ind.rv]
#    res <- res-weighted.mean(res)
#    respr <- respr-mean(respr)
#    resps <- resps-mean(resps)
#    respg <- respg-mean(respg)
    loglike1 <- loglikelihood(par.opt)
    out$relativity <- FALSE
    loglike0 <- loglikelihood(par.opt)
    out$relativity <- TRUE
#    loglike0 <- sum(dnorm(respr,0,out$all[,3],log=TRUE))
#    loglike1 <- sum(dnorm(res,0,out$all[,3],log=TRUE))
                                        #loglike0 <- loglikelihood(par.opt,relativity=FALSE)
                                        #loglike1 <- loglikelihood(par.opt,relativity=TRUE)
#    cat('loglike without GR+SR=',loglike0,'; RMS=',sd(respr),'\n')
#    cat('loglike with GR+SR=',loglike1,'; RMS=',sd(res),'\n')
####compare phase of relativisticy RV and Keplerian RV
    plot(trv.all,scale(out$all$RV.all),xlab='BJD',ylab='Scaled RV[m/s]',ylim=c(-10,10))
    lines(tsim,scale(simulate$rvg),col='blue')
    lines(tsim,scale(simulate$rvs),col='green')
    lines(tsim,scale(simulate$rvg+simulate$rvs),col='red')
    lines(tsim,scale(simulate$rv),col='black')
    lines(tsim,-scale(simulate$rv),col='grey')
    legend('top',horiz=TRUE,xpd=NA,legend=c('SR','GR','SR+GR'),col=c('blue','green','red'),bty='n',inset=c(0,-0.3),lty=rep(1,5))
    legend('top',horiz=TRUE,xpd=NA,legend=c('Kepler','-Kepler'),col=c('black','grey'),bty='n',inset=c(0,-0.2),lty=rep(1,5))

    plot(trv.all%%Popt[1],scale(out$all$RV.all),xlab='Phase [day]',ylab='Scaled RV[m/s]',ylim=c(-10,10))
    points(tsim%%Popt[1],scale(simulate$rvg),col='blue',pch='.')
    points(tsim%%Popt[1],scale(simulate$rvs),col='green',pch='.')
    points(tsim%%Popt[1],scale(simulate$rvg+simulate$rvs),col='red',pch='.')
    points(tsim%%Popt[1],scale(simulate$rv),col='black',pch='.')
    points(tsim%%Popt[1],-scale(simulate$rv),col='grey',pch='.')
    legend('top',horiz=TRUE,xpd=NA,legend=c('SR','GR','SR+GR'),col=c('blue','green','red'),bty='n',inset=c(0,-0.3),lty=rep(1,5))
    legend('top',horiz=TRUE,xpd=NA,legend=c('Kepler','-Kepler'),col=c('black','grey'),bty='n',inset=c(0,-0.2),lty=rep(1,5))

####difference between emission time and BJD
    plot(trv.all,(rvfit$tauT[unlist(out$ind.all$rv)]-trv.all)*24*60,xlab='BJD',ylab='tau-BJD [min]',ylim=range(rvfit$tauT[unlist(out$ind.all$rv)]-trv.all,simulate$tauT-tsim)*24*60)
    lines(tsim,(simulate$tauT-tsim)*24*60,col='red')

}
if(length(out$data.binary)>0) source('binary_fit.R')

###planet velocity
if(out$Nrvc>0){
    for(i in 1:Nsig){
        if(any(i==out$Irvc)){
            par1 <- par.opt
            par1[1:(Nsig*Npar)] <- 0
            par1[(i-1)*Npar+1:Npar] <- par.opt[(i-1)*Npar+1:Npar]
            names(par1) <- names(par.opt)
            df <- out$rvc[[paste0('p',i)]]
            sim1 <- RV.kepler(pars.kep=par1,tt=out.sim$BJDtdb,out1=out.sim,kep.only=TRUE,bases=bases)
            rvc.sim <- sim1$rvc[[paste0('p',i)]]
            plot(df[,1],df[,2],xlab='BJD',ylab='RV [km/s]',main=paste('companion',i),ylim=range(df[,2],rvc.sim),xlim=range(tsim,df[,1]))
            arrows(df[,1],df[,2]+df[,3],df[,1],df[,2]-df[,3],length=0.05,angle=90,code=3,col='black',pch=20)
            lines(tsim,rvc.sim,col='red')
        }
    }
}
if(plotf) dev.off()
