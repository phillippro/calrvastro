#pdf('test.pdf',8,8)
#par(mfrow=c(2,2))
tmin <- min(trv.all)
tmax <- max(trv.all)
panel <- 0
Popt <- exp(par.opt[(1:Nsig)*5-4])
t0 <- combine.list(lapply(ins,function(i) rv.list[[i]][,1]))
index <- sort(t0,index.return=TRUE)$ix
#fit.type <- 'multiple'#or multiple
fit.type <- 'multiple'
tcol.percentage <- 80
tmp <- cal.residual(par.opt,bases=bases)
rvmn.list <- tmp$res.noise# rv-noise
rvma.list <- tmp$res.all#rv-signal-noise
rv.arma <- cal.residual(par.opt)$rv.arma#ARMA rv
rv.trend <- cal.residual(par.opt)$rv.trend#trend rv
y0 <- lapply(ins,function(k) rvmn.list[[k]])#sort
y0 <- combine.list(y0)#combine list
rvmn.comb <- y0[index]#
tsim <- seq(tmin,tmax,length.out=1e4)
pn <-c('b','c','d','e','f','g','h')
inds <- sort(Popt,index.return=TRUE)$ix
indd <- 1:Nsig
if(target0=='HIP55042') indd <- 2
planets <- c()
for(kk in indd){
    j3  <- inds[kk]
    Nsig0 <- Nsig
    psig <- exp(par.opt[j3*5-4])
    bin.raw <- bin.data <- list()
    t1 <- (tsim-tmin)%%psig
    ix <- sort(t1,index.return=TRUE)$ix
    t2 <- t1[ix]

##calculate residual after subtraction of single signal
    par1 <- rep(0,length(par.opt))
    par1[(j3-1)*5+1:5] <-par.opt[(j3-1)*5+1:5]
    names(par1) <- names(par.opt)
    rv.kep <- cal.residual(par1)$rv.sig# residual
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
    rv.sim <- RV.kepler(pars.kep=par.sig,tt=tsim,kep.only=TRUE,bases=bases)$rv
    rv2 <- rv.sim[ix]

###prepared binned data
    for(i in ins){
        tt <- (rv.list[[i]][,1]-tmin)%%psig#day
        rv <- rv.sig[[i]]
        erv <- rv.list[[i]][,3]
        ##binnig data
        data.bin <- bin.simple(cbind(tt,rv,erv),10)
        bin.data[[i]] <- data.bin
        raw.bin <- bin.simple(cbind(rv.list[[i]][,1],rv,erv),10*ceiling((tmax-tmin)/psig))
        bin.raw[[i]] <- raw.bin
    }

###single signal fit to individual data set
    dx <- 2
    par(mar=c(0,dx,dx,dx))
    np <- np+1
###single signal fit to combined data
    tt1 <- (trv.all-tmin)%%psig
    planet <- paste(target.name,pn[min(length(indd),kk)])
    if(target.name=='GJ 3822' & pn[min(length(indd),kk)]=='b') planet <- 'GJ 3822 rotation'
    if(target.name=='GJ 3822' & pn[min(length(indd),kk)]=='c') planet <- 'GJ 3822 b'
    if(target.name=='GJ 480' & pn[min(length(indd),kk)]=='c') planet <- 'GJ 480 rotation'
    if(target.name=='GJ 382' & pn[min(length(indd),kk)]=='b') planet <- 'GJ 382 rotation harmony'
    if(target.name=='GJ 382' & pn[min(length(indd),kk)]=='c') planet <- 'GJ 382 rotation'
    if(target.name=='GJ 382' & pn[min(length(indd),kk)]=='d') planet <- 'GJ 382 b'

    if(target.name=='GJ 422' & pn[min(length(indd),kk)]=='b') planet <- 'GJ 422 b'
    if(target.name=='GJ 433' & pn[min(length(indd),kk)]=='c') planet <- 'GJ 433 d'
    if(target.name=='GJ 433' & pn[min(length(indd),kk)]=='d') planet <- 'GJ 433 c'
    planets <- c(planets,planet)

    rv.bin.all <- unlist(sapply(names(bin.data),function(i)bin.data[[i]][,2]))
#    ylim=range(rv.sig.comb,rv2)
    ylim=range(rv.bin.all,rv2)
    plot(tt1,rv.sig.comb,xlab='',ylab='RV (data-model) [m/s]',col='white',xaxt='n',main=planet,ylim=ylim,pch=pch.type,cex=pch.size,cex.main=1.5)
    axis(side=1,labels=FALSE,cex.lab=1.5)
    for(j in 1:length(ins)){
        i <- ins[j]
        points( (rv.list[[i]][,1]-tmin)%%psig,rv.sig[[i]],col=tcol(colors[ins.name[j]],tcol.percentage),pch=pch.type,cex=pch.size)
        data.bin <- bin.data[[i]]
        points(data.bin[,1],data.bin[,2],col=colors[ins.name[j]],pch=pch.type,cex=pch.size)
        try(arrows(data.bin[,1],data.bin[,2]+data.bin[,3],data.bin[,1],data.bin[,2]-data.bin[,3],length=0.05,angle=90,code=3,col=colors[ins.name[j]]),TRUE)
    }
    lines(t2,rv2,col='red')
    n1 <- ceiling(length(ins)/2)
    n2 <- length(ins)-n1
#    legend('top',legend=ins[1:n1],col=colors[ins.name[1:n1]],inset=c(0,-0.2),seg.len = 1,horiz=TRUE,xpd=NA,bty='n',pch=1)
    if(n2>0){
#        legend('top',legend=ins[(n1+1):length(ins)],col=colors[ins.name[(n1+1):length(ins)]],seg.len = 1,inset=c(0,-0.1),horiz=TRUE,xpd=NA,bty='n',pch=1)
    }
    if(np%%2==1) mtext(side=2,text='RV [m/s]',xpd=NA,line=2.5)
#    if(np%%2==0) mtext(side=1,text='RV [m/s]',xpd=NA,line=3)

#####residual plot
    ####residual
    par(mar=c(dx,dx,0,dx))
    np <- np+1
    RV.res <- unlist(sapply(ins,function(i)rvma.list[[i]]))
    plot(tt1,RV.res,xlab='',ylab='',col='grey',ylim=range(min(RV.res),max(RV.res)+0.5*(max(RV.res)-min(RV.res))),pch=pch.type,cex=pch.size)
#    axis(4)
    rms <- sqrt(mean(RV.res^2))
    legend('topright',legend=paste0('RMS=',format(round(rms,1),1),' m/s'),bty='n',text.col='black',cex=1.5,inset=c(0,-0.05))
    abline(h=0,lty=2,lwd=1,)
   for(j in 1:length(ins)){
        i <- ins[j]
        points( (rv.list[[i]][,1]-tmin)%%psig,rvma.list[[i]],col=tcol(colors[ins.name[j]],tcol.percentage),pch=pch.type,cex=pch.size)
        data.bin <- bin.data[[i]]
#        points(data.bin[,1],data.bin[,2],col=colors[ins.name[j]],pch=20)
#        try(arrows(data.bin[,1],data.bin[,2]+data.bin[,3],data.bin[,1],data.bin[,2]-data.bin[,3],length=0.05,angle=90,code=3,col=colors[ins.name[j]]),TRUE)
    }
    if(np%%2==0) mtext(side=2,text='O-C [m/s]',xpd=NA,line=2.5)
#    if(np<=8 & np%%2==1) mtext(side=2,text='RV (data-model) [m/s]',xpd=NA,line=3)
#    if(np>=(8*4-7)) mtext(side=4,text='O-C [m/s]',xpd=NA,line=3)
    if(np%%(Nrow*2)==0 | np==Ncol*Nrow*2) mtext(side=1,text='Orbital phase [day]',line=3)
    Nsig <- Nsig0
}
