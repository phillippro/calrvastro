pdf('combined.pdf',14,12)
Nsig <- 2
Pkep <- psig <- Popt <- exp(par.opt['per2'])
size <- 1.2
par(mar=c(5,5,1,5),mfrow=c(2,2),cex.axis=size,cex.lab=size,cex=size)
panel <- 0
t0 <- combine.list(lapply(ins,function(i) out[[i]]$RV[,1]))
index <- sort(t0,index.return=TRUE)$ix
#fit.type <- 'multiple'#or multiple
fit.type <- 'multiple'
tcol.percentage <- 80

##combined fit
Popt <- extract.par(par.opt,Np=Nsig,bases=bases)$P
tmp <- cal.residual(par.opt,bases=bases)
rvmn.list <- tmp$res.noise# rv-noise
rvma.list <- tmp$res.all#rv-signal-noise
rv.arma <- tmp$rv.arma#ARMA rv
rv.trend <- tmp$rv.trend#trend rv
y0 <- lapply(ins,function(k) rvmn.list[[k]])#sort
y0 <- combine.list(y0)#combine list
rvmn.comb <- y0[index]#
tsim <- seq(tmin,tmax,length.out=1e4)
for(j3 in 1:Nsig){
    psig <- Popt[j3]
    bin.raw <- bin.data <- list()

##calculate residual after subtraction of single signal
    par1 <- rep(0,length(par.opt))
    if(prior.type!='e0'){
        par1[(j3-1)*5+1:5] <-par.opt[(j3-1)*5+1:5]
    }else{
        par1[(j3-1)*3+1:3] <-par.opt[(j3-1)*3+1:3]
    }
    names(par1) <- names(par.opt)
    rv.kep <- cal.residual(par1,bases=bases)$rv.sig# residual
    rv.sig <- rvma.list
    for(k in 1:length(rv.sig)){
        rv.sig[[k]] <-rv.sig[[k]]+rv.kep[[k]]
    }
    y0 <- sapply(ins,function(k) rv.sig[[k]])#sort
    y0 <- combine.list(y0)#combine list
    rv.sig.comb <- y0[index]#form a residual combined set

###calculate residual after subtraction of all signals
    par.sig <- par.opt
    if(!exists('Kmin')) Kmin <- 1e-6
    if(basis=='natural'){
        par.sig[(1:Nsig)*5-3] <- Kmin
    }else{
        par.sig[(1:Nsig)*5-3] <- log(Kmin)
    }
    par.sig[j3*5-3] <- par.opt[j3*5-3]
    names(par.sig) <- names(par.opt)
    rv.sim <- RV.kepler(pars.kep=par.sig,tt=tsim,kep.only=TRUE,bases=bases)
    t1 <- (tsim-tmin)%%psig
    ix <- sort(t1,index.return=TRUE)$ix
    t2 <- t1[ix]
    rv2 <- rv.sim[ix]

###single signal fit to individual data set
    for(i in ins){
        tt <- (out[[i]]$RV[,1]-tmin)%%psig#day
        rv <- rv.sig[[i]]
        erv <- out[[i]]$RV[,3]
        plot(tt,rv.sig[[i]],xlab=paste0('Orbital Phase [days]'),ylab='RV (data-model) [m/s]',main=paste0('signal',j3,': P=',round(psig,1),'d; ',i),col='white',xlim=range(t2),ylim=range(rv.sig[[i]],rv2))
        if(exists('Tp')){
            Tp1 <- tmax-(tmax-Tp)%%Popt[1]
            Tp1 <- Tp1+(0:10)*Popt[1]
            Tp2 <- (Tp1-tmin)%%psig
            abline(v=Tp2,col='blue')
        }
    panel <- panel+1
    if(panel==1){
        mtext("Phase Plot", outer=TRUE,  cex=1.5, line=-0.5)
    }
                                        #        points(t1[ix],rv.sim[ix],col='red',pch='.')
        points(tt,rv.sig[[i]],col='grey')
        lines(t2,rv2,col='red')
        ##binnig data
        data.bin <- bin.simple(cbind(tt,rv,erv),50)
        bin.data[[i]] <- data.bin
        raw.bin <- bin.simple(cbind(out[[i]]$RV[,1],rv,erv),10*ceiling((tmax-tmin)/psig))
        bin.raw[[i]] <- raw.bin
        points(data.bin[,1],data.bin[,2],col='black',pch=20)
        try(arrows(data.bin[,1],data.bin[,2]+data.bin[,3],data.bin[,1],data.bin[,2]-data.bin[,3],length=0.05,angle=90,code=3,col='black',pch=20),TRUE)
    }

###single signal fit to combined data
    plot((trv.all-tmin)%%psig,rv.sig.comb,xlab=paste0('Orbital Phase [days]'),ylab='RV (data-model) [m/s]',col='white',ylim=range(mean(rv.sig.comb)-3*sd(rv.sig.comb),mean(rv.sig.comb)+3*sd(rv.sig.comb),rv2))
                                        #paste0('signal',j3,': P=',round(psigpsig[j3],1),'d; combined set')
    ins1 <- ins[ins!='AAT']
    for(j in 1:length(ins1)){
        i <- ins1[j]
        points( (out[[i]]$RV[,1]-tmin)%%psig,rv.sig[[i]],col=tcol(colors[j],tcol.percentage))
        data.bin <- bin.data[[i]]
        points(data.bin[,1],data.bin[,2],col=colors[j],pch=20,cex=2)
        try(arrows(data.bin[,1],data.bin[,2]+data.bin[,3],data.bin[,1],data.bin[,2]-data.bin[,3],length=0.05,angle=90,code=3,col=colors[j]),TRUE)
    }
###add CHIRON
if(FALSE){
#if(TRUE){
#    chiron <- read.table('../data/combined/HD128620v1/HD128620_CHIRON.dat')
#    chiron <- read.table('../data/combined/HD128620v1/HD128620_ES.dat')
    chiron <- read.table('../data/combined/HD128620v1/HD128620_VLCres.dat',header=TRUE)
    y <- lm(chiron[,2]~poly(chiron[,1],3))$residuals
    tc <- (chiron[,1]-tmin)%%psig#day
    ind <- sort(tc,index.return=TRUE)$ix
    cbin <- bin.simple(cbind(tc[ind],y[ind],chiron[ind,3]),20)
    points(tc,y,col=tcol('black',80))
    points(cbin[,1],cbin[,2],col='black',pch=20,cex=2)
    try(arrows(cbin[,1],cbin[,2]+cbin[,3],cbin[,1],cbin[,2]-cbin[,3],length=0.05,angle=90,code=3,col='black'),TRUE)
}

###model prediction
    lines(t2,rv2,col='red')
    legend('topright',legend=ins1,col=colors[1:length(ins1)],inset=c(-0.25,0),seg.len = 1,xpd=NA,bty='n',pch=20)
}

###combined fit to the whole combined data
    rv.sim <- RV.kepler(pars.kep=par.opt,tt=tsim,kep.only=TRUE,bases=bases)
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
        points(out[[i]]$RV[,1],rvmn.list[[i]],col=tcol(colors[j],80))
    }

###legends
    legend('topright',legend=ins,col=colors[1:length(ins)],inset=c(-0.25,0),seg.len = 1,xpd=NA,bty='n',pch=1)
dev.off()
