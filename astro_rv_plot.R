library(magicaxis)
library(graphics)
source('mcmc_func.R')
options(scipen = 0)
Nsig <- length(fs)
Nastro <- 2
#ell.col <- rep('darkgrey',2)
ell.col <- c('orange','darkgreen')
#ell.pch <- c(16,17)
ell.pch <- c(20,20)
cols <- tcol(c('black','blue','green','purple','cyan','brown'),50)
#names(cols) <- c('AAT','APF1','ELODIE','')
#fpdf <- paste0(star0,'_planet1_astro_rv',Nsig,'v',version,'.pdf')
if(plotf){
dir <- paste0('../paper/astrometry/v3/')
fpdf <- paste0(dir,'star',length(stars),'_v',version,'.pdf')
cat('output:\n',fpdf,'\n')
pdf(fpdf,12,3*Nsig)
}
size <- 0.9
show.center <- FALSE
margin <- 1.3
pch.size <- 1
if(Nmc>0){
    line.size <- 1
}else{
    line.size <- 3
}
fit.col <- tcol('red',50)
if(plotf){
    split.screen(c(Nsig, 1))#screen 1 to Nsig
#par(mfrow=c(6,3),oma=c(1,6,1,6),cex=size,cex.lab=size,cex.axis=size,mar=c(margin,2*margin,margin,margin),mgp=c(2,0.5,0),cex.main=1.2*size)
if(Nsig>2){
    oma <- c(4,8,1,6)
}else{
    oma <- c(3,8,0,6)
}
par(oma=oma,cex=size,cex.lab=size,cex.axis=size,mar=c(margin,2*margin,margin,0.5*margin),mgp=c(2,0.5,0),cex.main=1.2*size)
}
#layout(matrix(c(1, 1,2,2), nrow=2, byrow=TRUE))
#Nstar <- length(stars)
Nstar <- Nsig
####RV
for(i1 in 1:Nstar){
    out <- list()
    Mstar <- Mstars[i1]
    out[['astrometry']] <- astrometry.tot[[i1]]
    out$Nm <- 0
    split.screen( figs = c(1,3), screen = i1 )#screen Nsig+(i-1)*5+1...3
    split.screen(rbind(c(0,1,0.3,1),c(0,1,0,0.3)),Nsig+(i1-1)*5+1)#Nsig+(i-1)*5+(4~5)
    if(plotf) par(mar=c(0,0,0,0),cex.axis=1)
    Nset <- Nrv.tot[[i1]]
    trv.all <- tall <- unlist(sapply(1:Nset,function(i) JD.tot[[i1]][[i]]))
    ind <- sort(trv.all,index.return=TRUE)$ix
    tmin <- min(trv.all)
    trv.all <- tall <- tall[ind]%%24e5
    RV.all <- unlist(sapply(1:Nset,function(i) RV.tot[[i1]][[i]]))[ind]
    eRV.all <- unlist(sapply(1:Nset,function(i) eRV.tot[[i1]][[i]]))[ind]
    RV2.all <- unlist(sapply(1:Nset,function(i) RV.tot[[i1]][[i]] + res.tot[[i1]][[i]] - res.sig[[i1]][[i]]))[ind]
    data.astrometry <- astrometry.tot[[i1]]
    cov.astro <- cov.tot[[i1]]
    Nw0 <- Nrv.tot[[i1]]
    tlist <- JD.tot[[i1]]
    erv.list <- eRV.tot[[i1]]
###split the RV plot into two sub-panel
    if(plotf){
    screen(Nsig+(i1-1)*5+4)
    ins <- instr <- ins.tot[[i1]]
    par(mar=c(0,2*margin,margin,margin))
    nqp <- nqp.tot[[i1]]
    for(i3 in 1:Nset){
        t <- JD.tot[[i1]][[i3]]%%24e5
        rv <- RV.tot[[i1]][[i3]] + res.tot[[i1]][[i3]] - res.sig[[i1]][[i3]]
        erv <- eRV.tot[[i1]][[i3]]
        if(i3==1){
            plot(t,rv,xlab=expression('JD-2400000'),ylab='',col=cols[i3],pch=20,xlim=c(min(tall),max(tall)),ylim=range(RV2.all),xaxt='n',yaxt='n',cex=pch.size)
#            magaxis(c(1,2))
            magaxis(2,majorn=2,cex=size)
        }else{
            points(t,rv,pch=20,col=cols[i3],cex=pch.size)
        }
        arrows(t,rv-erv,t,rv+erv,length=0.05,angle=90,code=3,col=cols[i3])
    }
    mtext(side=2,text=expression('RV [m/s]'),line=1.5,xpd=NA,cex=1.2*size)
    ind <- which(!grepl('HARPS',instr))
    if(length(ind)>0) instr[ind] <- toupper(instr[ind])
    if(grepl('HD209100',stars[i1])){
        instr <- c('HARPSpre','LC','VLC','HARPSpost','UVES')
    }
    if(exists('star.name')){
        sn <- star.name
    }else{
        sn <- stars[i1]
    }
    legend('topleft',legend=c(instr,'Hipparcos','Gaia DR2','Best fit'),col=c(cols[1:Nw0],ell.col,fit.col),lty=c(rep(NA,Nw0+2),1),pch=c(rep(20,Nw0),ell.pch,NA),xpd=NA,inset=c(-0.75,0),bty='n',title=sn,lwd=c(rep(NA,Nw0+2),line.size),pt.cex=c(rep(1,Nw0),rep(1,2),NA),y.intersp=1.1)
#    legend('bottomleft',legend=c('Hipparcos','Gaia'),pch=rep(20,2),col=ell.col,xpd=NA,inset=c(-0.7,0),bty='n',)
    }
    tt <- seq(min(tall),max(tall),by=0.1)
    par.opt <- pars <- par1.tot[[i1]]
#    ind <- grep('\\db\\d',rownames(parss))
#    names(pars)[ind] <- paste0(1:length(ind),'b1')
    ind <- grep('per',names(pars))
#    par.opt <- pars[,1]
#    par1.tot[]
    Popt <- exp(par.opt[ind])
    prior.type <- 'mt'
    jitter <- TRUE
    parallax <- FALSE
    if(!exists('bases')) bases <- rep('natural',5)
    astrometry <- 8
    kep.type <- 'pure'
    yr2d <- 365.25
    period.par <- 'logP'
    offset <- TRUE
    time.unit <- 1
    Np <- length(ind)
    ##
    if(plotf){
    lines(tt,RV.kepler(pars.kep=par.opt,tt=tt,kep.only=TRUE,nqp=nqp)[[1]],col=fit.col,lwd=line.size)$rv
    if(Nmc>0){
        for(j in 2:length(pars)){
            lines(tt,RV.kepler(pars.kep=par.opt,tt=tt,kep.only=TRUE,nqp=nqp)$rv[[1]],col=fit.col,lwd=1)
        }
    }
###residual plot
    screen(Nsig+(i1-1)*5+5)
    par(mar=c(margin,2*margin,0,margin))
    rv.planet <- RV.kepler(pars.kep=par.opt,tt=NA,kep.only=TRUE,nqp=nqp)$rv
    drv <- res.tot[[i1]]
#    for(j in 1:Nset){
#        drv[[j]] <- rv.list[[j]]-rv.planet[[j]]
#    }
    ylim <- range(unlist(drv))
    xlim <- range(tall)%%24e5
    for(j in 1:Nset){
        if(j==1){
            plot(tlist[[j]]%%24e5,drv[[j]],xlab='',ylab='',ylim=ylim,xlim=xlim,col=tcol(cols[j],50),pch=20,cex=size,xaxt='n',yaxt='n')
#            axis(1,tcl = 0.5)
            magaxis(1,majorn=2,cex=size)
            magaxis(2,majorn=2,cex=size)
#            magaxis(2,majorn=3,cex=size)
#            axis(2,tck=0.1)
        }else{
            points(tlist[[j]]%%24e5,drv[[j]],col=cols[j],pch=20,cex=1)
        }
#        arrows(tlist[[j]]%%24e5,drv[[j]]-erv.list[[j]],tlist[[j]]%%24e5,drv[[j]]+erv.list[[j]],length=0.05,angle=90,code=3,col=tcol(cols[j],50))
    }
#    par(mgp=c(3,1,0),xpd=NA)
    if(i1==Nsig) mtext(side=1,text='JD-2400000',xpd=NA,line=2,cex=1.2*size)
    mtext(side=2,text=expression('O-C [m/s]'),line=1.5,xpd=NA,cex=1.2*size)
    abline(h=0,lty=2,lwd=2)
    }
####astrometry
    source('reflex_motion.R')
}
if(plotf){
close.screen(all.screens = TRUE)
dev.off()
}
