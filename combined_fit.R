library(RColorBrewer)
if(!exists('out')){
    load('results/HD74014/HD74014_fix1_relativityFALSE_Niter1000000_Ncores8_ofac2_Nset3_Esd1_transit0_P7044_acc9.2_lnlmax-124.Robj')
#    load('results/HD95089/HD95089_fix1_relativityFALSE_Niter1600000_Ncores8_ofac2_Nset3_Esd1_transit0_P1753d459_acc0.2_lnlmax-391.Robj')
}
colors <- c(brewer.pal(9, "Set1"),brewer.pal(n = 8, name = "Set2"))
source('mcmc_func.R')
###plot
pdf('test.pdf',8,8)
#par(mfrow=c(2,2),mar=c(5,5,1,1))
layout(matrix(1:8, ncol = 2,byrow=TRUE), widths = 1, heights = c(2,4), respect = FALSE)
panel <- 0
#Pkep <- out$Popt[[paste0('sig',Nsig)]]
Pkep <- par.opt[grep('per',names(par.opt))]
if(!any(names(out$par.stat)==paste0('sig',Nsig)) & Nmax>Nmin){
    mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
    par.stat <- matrix.distr(mc[,1:(ncol(mc)-2)],mc[,ncol(mc)],plotf=FALSE)
    out$par.stat[[paste0('sig',Nsig)]] <- par.stat
#    Pkep <- Popt <- extract.par(par.stat['xopt',],Nsig,bases)
}
#par.opt <- out$par.stat[[paste0('sig',Nsig)]]['xopt',]
#Popt <- out$Popt[[paste0('sig',Nsig)]]
Popt <- extract.par(par.opt,bases=bases)$P
#fit.type <- 'multiple'#or multiple
fit.type <- 'multiple'
tcol.percentage <- 50

##combined fit
if(out$Nrv>0){
    tmp <- cal.residual(par.opt,bases=bases)
    rvmn.list <- tmp$res.noise# rv-noise
    rvma.list <- tmp$res.all#rv-signal-noise
    rvms.list <- tmp$res.sig#rv-signal
    rv.arma <- tmp$rv.arma#ARMA rv
    rv.trend <- tmp$rv.trend#trend rv
    y0 <- lapply(ins,function(k) rvmn.list[[k]])#sort
    y0 <- combine.list(y0)#combine list
    t0 <- combine.list(lapply(ins,function(i) out[[i]]$RV[,1]))
    index <- sort(t0,index.return=TRUE)$ix
    rvmn.comb <- y0[index]#
    tsim <- seq(tmin,tmax,length.out=1e4)
    for(j3 in 1:Nsig){
        psig <- Popt[j3]
        bin.raw <- bin.data <- list()

        ##calculate residual after subtraction of single signal
        if(any(grepl('Omega',names(par.opt)))){
            Npar <- 7
        }else{
            Npar <- 5
        }
        if(prior.type=='e0'){
            Npar <- 3
        }
        par1 <- rep(0,length(par.opt))
        par1[(j3-1)*Npar+1:Npar] <- par.opt[(j3-1)*Npar+1:Npar]
        names(par1) <- names(par.opt)
        rv.kep <- cal.residual(par1,bases=bases)$rv.sig# residual
        rv.sig <- rvma.list
                                        #    rv.sig <- rvmn.list
                                        #    rv.sig <- rvma.list+rv.kep
        for(k in 1:length(rv.sig)){
            rv.sig[[k]] <-rv.sig[[k]]+rv.kep[[k]]
        }
        y0 <- sapply(ins,function(k) rv.sig[[k]])#sort
        y0 <- combine.list(y0)#combine list
        rv.sig.comb <- y0[index]#form a residual combined set

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
        rv.sim <- RV.kepler(pars.kep=par.sig,tt=tsim,kep.only=TRUE,bases=bases)$rv
        t1 <- (tsim-tmin)%%psig
        ix <- sort(t1,index.return=TRUE)$ix
        t2 <- t1[ix]
        rv2 <- rv.sim[ix]

###single signal fit to individual data set
        for(i in ins){
            tt <- (out[[i]]$RV[,1]-tmin)%%psig#day
            rv <- rv.sig[[i]]
            erv <- out[[i]]$RV[,3]
                                        #stop()
            if(exists('Tp')){
                Tp1 <- tmax-(tmax-Tp)%%Popt[1]
                Tp1 <- Tp1+(0:10)*Popt[1]
                Tp2 <- (Tp1-tmin)%%psig
                abline(v=Tp2,col='blue')
            }
            panel <- panel+1
                                        #        points(t1[ix],rv.sim[ix],col='red',pch='.')
            ##binnig data
            data.bin <- bin.simple(cbind(tt,rv,erv),10)
            bin.data[[i]] <- data.bin
            raw.bin <- bin.simple(cbind(out[[i]]$RV[,1],rv,erv),10*ceiling((tmax-tmin)/psig))
            bin.raw[[i]] <- raw.bin
        }
###single signal fit to combined data
                                        #paste0('signal',j3,': P=',round(psigpsig[j3],1),'d; combined set')
        for(j in 1:length(ins)){
            i <- ins[j]
            data.bin <- bin.data[[i]]
        }
                                        #    legend('topright',legend=ins,xpd=NA,inset=c(-0.45,0),col=colors[1:length(ins)],pch=1,bty='n')
        n1 <- ceiling(length(ins)/2)
        n2 <- length(ins)-n1
        }

###combined fit to the whole combined data
    sim <- RV.kepler(pars.kep=par.opt,tt=tsim,kep.only=TRUE,bases=bases)
    rv.sim <- sim$rv
    rv2 <- rv.sim[ix]
    plot(trv.all,rvmn.comb,xlab=paste0('Time [days]'),ylab='RV (data-model) [m/s]',col='white',ylim=range(rvmn.comb,rv.sim))
                                        #paste0('signal',j3,': P=',round(psig[j3],1),'d; combined set')
    lines(tsim,rv.sim,col='red')
    if(exists('Tp')){
        Tp1 <- tmax-(tmax-Tp)%%Popt[1]-(0:10)*Popt[1]
        abline(v=Tp1,col='blue')
    }

    for(j in 1:length(ins)){
        i <- ins[j]
        points(out[[i]]$RV[,1],rvmn.list[[i]],col=tcol(colors[j],tcol.percentage))#        raw.bin <- bin.raw[[i]]
                                        #        points(raw.bin[,1],raw.bin[,2],col=colors[j])
                                        #        try(arrows(raw.bin[,1],raw.bin[,2]+raw.bin[,3],raw.bin[,1],raw.bin[,2]-raw.bin[,3],length=0.05,angle=90,code=3,col=colors[j],pch=20),TRUE)
    }
                                        #    legend('topright',legend=ins,xpd=NA,inset=c(-0.45,0),col=colors[1:length(ins)],pch=1,bty='n')
    legend('topright',legend=ins,col=colors[1:length(ins)],inset=c(-0.25,0),seg.len = 1,xpd=NA,bty='n',pch=1)
}


dev.off()
