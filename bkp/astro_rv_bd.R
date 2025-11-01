library(magicaxis)
library(graphics)
source('mcmc_func.R')
options(scipen = 0)
library("plotrix")
source('OrbitFunction.R')
source('timing_function.R')
set.seed(9999)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    star <- as.character(args[1])
}else{
#    star <- 'GJ3222'
    star <- 'HD4113'
}
version <- 1
set.seed(100)
ff <- read.table('companion_files.txt')[,1]
stars <- gsub('\\/.+','',ff)
if(star!=''){
   ind <- which(stars==star)
}else{
#   ind <- 1:length(stars)
   ind <- 1:10
}
if(length(ind)==0){
    stop('No file is found for this target!')	
}
fpdf <- ff[ind]
stars <- stars[ind]
fs <- c()
#stars <- c()
#for(f in fpdf){
for(j in 1:length(stars)){
     f <- fpdf[j]
     s <- stars[j]
     pat <- gsub('.+\\/|\\.pdf','',f)
     pat <- gsub('\\+','\\\\+',pat)
     fs <- c(fs,list.files(paste0('results/',s),pattern=paste0(pat,'.+Robj'),full.name=TRUE)[1])
}
plotf <- TRUE
etaH <- etaG <- 1
fs <- sort(fs)
###modifiy input file names
#Nsig <- length(fs)
planets <- read.csv('../data/code/SigClassify/ranking_complex2.csv')
n0 <- gsub(' ','',as.character(planets[,'Name']))
n1 <- gsub(' ','',as.character(planets[,'ID']))
n2 <- gsub(' ','',as.character(planets[,'StarKnown']))
###collect all data needed for the plot
ins.tot <- list()
res.sig <- res.tot <- list()
eRV.tot <- list()
RV.tot <- list()
nqp.tot <- JD.tot <- list()
Nrv.tot <- cov.tot <- astrometry.tot <- list()
id.tot <- list()
par0.tot <- par1.tot<- list()
for(j in 1:length(fs)){
    star <- stars[j]
    cat('load ',fs[j],'\n')
    load(fs[j],env=e0 <- new.env())
    Nkep <- length(e0$Popt)
    Nsig <- 1#e0$Nsig	
    out <- e0$out
    out$Nm <- 0
    out$Nrvc <- length(out$rvc)
    out$relativity <- FALSE
    out$Nastro <- e0$Nastro
    out$Nrv <- nrow(out$all)
    par1.tot[[j]] <- out$par.stat[[paste0('sig',Nkep)]][1,]
nepoch <- 2
astro.type <- 'astro'
if(any(names(e0)=='Mstar')) Mstars <- e0$Mstar
if(any(names(e0)=='eMstar')) eMstars <- e0$eMstar
if(any(names(out)=='Mstar')){
    Mstars <- Mstar <- out$Mstar
}
if(any(names(out)=='eMstar')){
    eMstars <- eMstar <- out$eMstar
}
if(!exists('eMstar')) eMstars <- eMstar <- 0.1*Mstar
    cov.tot[[j]] <- e0$cov.astro
    ins.tot[[j]] <- ins <- e0$ins
    JD.tot[[j]] <- RV.tot[[j]] <- eRV.tot[[j]] <- res.tot[[j]] <- res.sig[[j]] <- nqp.tot[[j]] <- list()
    Nrv <- length(e0$ins)
    Nrv.tot[[j]] <- Nrv
#    cat('Nrv=',Nrv,'\n')
#    cat('ins=',e0$ins,'\n')
    for(k in 1:Nrv){
        t <- out[[ins[k]]]$RV[,1]
        if(max(t)<24e5) t <- t+24e5
        JD.tot[[j]][[ins[k]]] <- t
        nqp.tot[[j]][[ins[k]]] <- out[[ins[k]]]$noise$nqp
        RV.tot[[j]][[ins[k]]] <- out[[ins[k]]]$RV[,2]
        eRV.tot[[j]][[ins[k]]] <- out[[ins[k]]]$RV[,3]
        res.tot[[j]][[ins[k]]] <- out$res.all[[paste0('sig',Nkep)]][[ins[k]]]
        res.sig[[j]][[ins[k]]] <- out$res.sig[[paste0('sig',Nkep)]][[ins[k]]]
    }
    astrometry.tot[[j]] <- out$astrometry
}
Nastro <- 2
#ell.col <- rep('darkgrey',2)
#ell.col <- c('darkblue','darkgreen')
ell.col <- c('darkblue','darkgreen')
#ell.pch <- c(16,17)
ell.pch <- c(20,20)
cols <- tcol(c('black','blue','green','purple','cyan','brown'),50)
#names(cols) <- c('AAT','APF1','ELODIE','')
#fpdf <- paste0(star0,'_planet1_astro_rv',Nsig,'v',version,'.pdf')
size <- 0.9
show.center <- FALSE
margin <- 1.3
pch.size <- 1
Nmc <- Nsamp <- 0
if(Nmc>0){
    line.size <- 2
}else{
#    line.size <- 3
    line.size <- 1
}
if(plotf){
  fpdf <- paste0('bd_plot/',star,'_N',Nmc,'_fit.pdf')
    cat('output:\n',fpdf,'\n')
    pdf(fpdf,12,3*Nsig)
}
fit.col <- tcol('red',90)
fit.opt <- tcol('red',20)
if(plotf){
    split.screen(c(1, Nsig))#screen 1 to Nsig
#par(mfrow=c(6,3),oma=c(1,6,1,6),cex=size,cex.lab=size,cex.axis=size,mar=c(margin,2*margin,margin,margin),mgp=c(2,0.5,0),cex.main=1.2*size)
if(Nsig>2){
    oma <- c(4,8,1,6)
}else{
    oma <- c(3,8,0,6)
}
par(oma=oma,cex=size,cex.lab=size,cex.axis=size,mar=c(margin,2*margin,margin,0.5*margin),mgp=c(2,0.5,0),cex.main=1.2*size)
}
#layout(matrix(c(1, 1,2,2), nrow=2, byrow=TRUE))
Nstar <- length(stars)
#Nstar <- Nsig
####RV
for(i1 in 1:Nstar){
    j1 <- i1
#    out <- list()
    Mstar <- Mstars[j1]
#    out[['astrometry']] <- astrometry.tot[[j1]]
    split.screen( figs = c(1,3), screen = i1 )#screen Nsig+(i-1)*5+1...3
    split.screen(rbind(c(0,1,0.3,1),c(0,1,0,0.3)),Nsig+(i1-1)*5+1)#Nsig+(i-1)*5+(4~5)
    if(plotf) par(mar=c(0,0,0,0),cex.axis=1)
    Nset <- Nrv.tot[[1]]
    trv.all <- tall <- unlist(sapply(1:Nset,function(i) JD.tot[[1]][[i]]))
    ind <- sort(trv.all,index.return=TRUE)$ix
    tmin <- min(trv.all)
    trv.all <- tall <- tall[ind]%%24e5
    RV.all <- unlist(sapply(1:Nset,function(i) RV.tot[[j1]][[i]]))[ind]
    eRV.all <- unlist(sapply(1:Nset,function(i) eRV.tot[[j1]][[i]]))[ind]
    RV2.all <- unlist(sapply(1:Nset,function(i) RV.tot[[j1]][[i]] + res.tot[[j1]][[i]] - res.sig[[j1]][[i]]))[ind]
    data.astrometry <- astrometry.tot[[j1]]
    cov.astro <- cov.tot[[j1]]
    Nw0 <- Nrv.tot[[j1]]
    tlist <- JD.tot[[j1]]
    erv.list <- eRV.tot[[j1]]
###split the RV plot into two sub-panel
    if(plotf){
    screen(Nsig+(i1-1)*5+4)
    ins <- instr <- ins.tot[[j1]]
    par(mar=c(0,2*margin,margin,margin))
    nqp <- nqp.tot[[j1]]
###simulation for ylim
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
    prior.type <- 'mt'
    par.opt <- pars <- par1.tot[[j1]]
    Popt <- exp(par.opt[ind])
    ind <- grep('per',names(pars))
    t3 <- seq(min(tall),max(tall),length.out=1e4)
    rvs <- RV.kepler(pars.kep=par.opt,tt=t3,kep.only=TRUE,nqp=nqp)[[1]]
    for(i3 in 1:Nset){
        t <- JD.tot[[1]][[i3]]%%24e5-tmin%%24e5
        rv <- RV.tot[[1]][[i3]] + res.tot[[j1]][[i3]] - res.sig[[j1]][[i3]]
        erv <- eRV.tot[[1]][[i3]]

        if(i3==1){
            plot(t,rv,xlab=paste0('JD-',tmin%%24e5+24e5),ylab='',col=cols[i3],pch=20,xlim=c(min(tall%%24e5-tmin%%24e5),max(tall%%24e5-tmin%%24e5)),ylim=range(RV2.all,rvs),xaxt='n',yaxt='n',cex=pch.size)#main='RV fit'
#            magaxis(c(1,2))
#            magaxis(2,majorn=2,cex=size)
            magaxis(2,cex=size)
        }else{
            points(t,rv,pch=20,col=cols[i3],cex=pch.size)
        }
        arrows(t,rv-erv,t,rv+erv,length=0.05,angle=90,code=3,col=cols[i3])
    }
    mtext(side=2,text=expression('RV [m/s]'),line=1.5,xpd=NA,cex=1.2*size)
    ind <- which(!grepl('HARPS',instr))
    if(length(ind)>0) instr[ind] <- toupper(instr[ind])
    if(grepl('HD209100',stars[j1])){
        instr <- c('HARPSpre','LC','VLC','HARPSpost','UVES')
    }
    if(exists('star.name')){
        sn <- star.name
    }else{
        sn <- stars[j1]
    }
    legend('topleft',legend=c(instr,'Hipparcos','Gaia DR2','Best fit','Barycenter'),col=c(cols[1:Nw0],ell.col,fit.opt,'black'),lty=c(rep(NA,Nw0+2),1,NA),pch=c(rep(20,Nw0),ell.pch,NA,3),xpd=NA,inset=c(-0.75,0),bty='n',title=sn,lwd=c(rep(NA,Nw0+2),line.size,NA),pt.cex=c(rep(1,Nw0),rep(1,2),NA,1),y.intersp=1.1)
#    legend('bottomleft',legend=c('Hipparcos','Gaia'),pch=rep(20,2),col=ell.col,xpd=NA,inset=c(-0.7,0),bty='n',)
    }
    tspan <- max(tall)-min(tall)
    tt <- seq(min(tall),max(tall),length.out=1e4)
#            rvs <- RV.kepler(pars.kep=mc[j3,1:Npar],tt=tt,kep.only=TRUE,nqp=nqp)$rv
#    par.opt <- pars <- par1.tot[[j1]]
#    ind <- grep('\\db\\d',rownames(parss))
#    names(pars)[ind] <- paste0(1:length(ind),'b1')
#    par.opt <- pars[,1]
#    par1.tot[]
    ##
    if(plotf){
    if(Nmc>0){
        mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
        inds <- sample(1:nrow(mc),Nsamp)
        Npar <- length(par.opt)
        for(j3 in inds){
#            rvs <- RV.kepler(pars.kep=mc[j3,1:Npar],tt=tt,kep.only=TRUE,nqp=nqp)$rv
            lines(tt%%24e5-tmin%%24e5,rvs,col=fit.col,lwd=1)
#            rvfun <- approxfun(tt,rvs)
#            tsim <- sort(runif(1e4,min(tt),max(tt)))
#            tsim <- seq(min(tt),max(tt),by=0.1)
#            lines(tsim-tmin,rvfun(tsim),col=fit.col,lwd=1)
        }
    }
    lines(tt%%24e5-tmin%%24e5,RV.kepler(pars.kep=par.opt,tt=tt,kep.only=TRUE,nqp=nqp)[[1]],col=fit.opt,lwd=line.size)$rv

####replot the RV data to be on top of the fitting lines
    for(i3 in 1:Nset){
        t <- JD.tot[[j1]][[i3]]%%24e5-tmin%%24e5
        rv <- RV.tot[[j1]][[i3]] + res.tot[[j1]][[i3]] - res.sig[[j1]][[i3]]
        erv <- eRV.tot[[j1]][[i3]]
        points(t,rv,pch=20,col=cols[i3],cex=pch.size)
        arrows(t,rv-erv,t,rv+erv,length=0.05,angle=90,code=3,col=cols[i3])
    }

###residual plot
    screen(Nsig+(i1-1)*5+5)
    par(mar=c(margin,2*margin,0,margin))
    rv.planet <- RV.kepler(pars.kep=par.opt,tt=NA,kep.only=TRUE,nqp=nqp)$rv
    drv <- res.tot[[j1]]
#    for(j in 1:Nset){
#        drv[[j]] <- rv.list[[j]]-rv.planet[[j]]
#    }
    ylim <- range(unlist(drv))
    xlim <- range(tall)%%24e5-tmin%%24e5
    for(j in 1:Nset){
        x <- tlist[[j]]%%24e5-tmin%%24e5
        y <- drv[[j]]
        ey <- erv.list[[j]]
        if(j==1){
            plot(x,y,xlab='',ylab='',ylim=ylim,xlim=xlim,col=tcol(cols[j],50),pch=20,cex=size,xaxt='n',yaxt='n')
#            axis(1,tcl = 0.5)
#            magaxis(1,majorn=2,cex=size)
#            magaxis(2,majorn=2,cex=size)
            magaxis(c(1,2),majorn=3,cex=size)
#            magaxis(2,majorn=3,cex=size)
#            axis(2,tck=0.1)
        }else{
            points(x,y,col=cols[j],pch=20,cex=1)
        }
        arrows(x,y-ey,x,y+ey,length=0.05,angle=90,code=3,col=tcol(cols[j],50))
    }
#    par(mgp=c(3,1,0),xpd=NA)
    if(i1==Nsig) mtext(side=1,text=paste0('JD-',tmin%%24e5+24e5),xpd=NA,line=2,cex=1.2*size)
    mtext(side=2,text=expression('O-C [m/s]'),line=1.5,xpd=NA,cex=1.2*size)
    abline(h=0,lty=2,lwd=2)
    }
####astrometry
    source('reflex_motion_bd.R')
}
if(plotf){
close.screen(all.screens = TRUE)
dev.off()
}


