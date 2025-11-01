library("plotrix")
source('OrbitFunction.R')
source('timing_function.R')
set.seed(9999)
#if(!exists('par0s') | TRUE){
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    star <- as.character(args[1])
}else{
#    star <- 'HD74014'
#    star <- 'HD129191'
#    star <- 'HIP2552'
#    star <- 'HD131664'
#    star <- 'HD16160'
#    star <- 'HD161797'
#    star <- 'HD190360'
#    star <- 'HD190406'
#    star <- 'HD39587'
#    star <- 'HD42581'
    star <- 'HD209100'
#    star <- ''
}
version <- 1
set.seed(100)
ff <- read.table('BD_candidates.txt')[,1]
stars <- gsub('\\/.+','',ff)
ind <- which(stars==star)
if(length(ind)>0){
#fs <- 'HD74014_fix1_relativityFALSE_Niter1000000_Ncores8_ofac2_Nset3_Esd1_transit0_P7044_acc9.2_lnlmax'
    fs <- gsub('.pdf','',gsub('.+\\/','',ff[ind[1]]))
}else{
    stop('Error: No file is found for this target!')
}
plotf <- TRUE
etaH <- etaG <- 1
fs <- sort(fs)
###modifiy input file names
stars <- gsub('_.+','',fs)
ff <- c()
for(k in 1:length(stars)){
    star  <- stars[k]
    path <- paste0('/malibu/ffeng/astro_output/',star)
    if(!file.exists(path)) path <- paste0('results/',star)
    ff <- c(ff,list.files(path=path,pattern=paste0(fs[k],'.+Robj'),full.name=TRUE)[1])
}
fs <- ff
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
    load(fs[j])
    Nkep <- length(Popt)
    out$Nm <- 0
    out$Nrvc <- length(out$rvc)
    out$relativity <- FALSE
    out$Nastro <- Nastro
    out$Nrv <- nrow(out$all)
    par1.tot[[j]] <- out$par.stat[[paste0('sig',Nkep)]][1,]
Nsig <- 1
nepoch <- 2
astro.type <- 'astro'
if(exists('Mstar'))Mstars <- Mstar
if(any(names(out)=='Mstar')){
    Mstars <- Mstar <- out$Mstar
}
if(any(names(out)=='eMstar')){
    eMstars <- eMstar <- out$eMstar
}
if(!exists('eMstar')) eMstars <- eMstar <- 0.1*Mstar
    cov.tot[[j]] <- cov.astro
    ins.tot[[j]] <- ins
    JD.tot[[j]] <- RV.tot[[j]] <- eRV.tot[[j]] <- res.tot[[j]] <- res.sig[[j]] <- nqp.tot[[j]] <- list()
    Nrv <- length(ins)
    Nrv.tot[[j]] <- Nrv
    cat('Nrv=',Nrv,'\n')
    cat('ins=',ins,'\n')
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
library(magicaxis)
library(graphics)
source('mcmc_func.R')
options(scipen = 0)
#Nsig <- length(fs)
Nastro <- 2
#ell.col <- rep('darkgrey',2)
#ell.col <- c('darkblue','darkgreen')
ell.col <- c('darkblue','darkgreen')
#ell.pch <- c(16,17)
ell.pch <- c(20,20)
cols <- tcol(c('black','blue','green','purple','cyan','brown'),50)
#names(cols) <- c('AAT','APF1','ELODIE','')
#fpdf <- paste0(star0,'_planet1_astro_rv',Nsig,'v',version,'.pdf')
if(plotf){
###dir <- paste0('../paper/astrometry/v3/')
    dir <- paste0('results/',target,'/')
###dir <- paste0('/malibu/ffeng/astro_output/',target,'/')
    fpdf <- paste0(dir,target,'_fit.pdf')
    cat('output:\n',fpdf,'\n')
    pdf(fpdf,12,3*Nsig)
}
size <- 0.9
show.center <- FALSE
margin <- 1.3
pch.size <- 1
Nmc <- Nsamp <- 100
if(Nmc>0){
    line.size <- 2
}else{
#    line.size <- 3
    line.size <- 1
}
fit.col <- tcol('red',90)
fit.opt <- tcol('red',20)
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
#    out <- list()
    Mstar <- Mstars[i1]
#    out[['astrometry']] <- astrometry.tot[[i1]]
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
            plot(t,rv,xlab=expression('JD-2400000'),ylab='',col=cols[i3],pch=20,xlim=c(min(tall),max(tall)),ylim=range(RV2.all),xaxt='n',yaxt='n',cex=pch.size)#main='RV fit'
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
    legend('topleft',legend=c(instr,'Hipparcos','Gaia DR2','Best fit','Barycenter'),col=c(cols[1:Nw0],ell.col,fit.opt,'black'),lty=c(rep(NA,Nw0+2),1,NA),pch=c(rep(20,Nw0),ell.pch,NA,3),xpd=NA,inset=c(-0.75,0),bty='n',title=sn,lwd=c(rep(NA,Nw0+2),line.size,NA),pt.cex=c(rep(1,Nw0),rep(1,2),NA,1),y.intersp=1.1)
#    legend('bottomleft',legend=c('Hipparcos','Gaia'),pch=rep(20,2),col=ell.col,xpd=NA,inset=c(-0.7,0),bty='n',)
    }
    tspan <- max(tall)-min(tall)
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
    if(Nmc>0){
        mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
        inds <- sample(1:nrow(mc),Nsamp)
        Npar <- length(par.opt)
        for(j3 in inds){
            rvs <- RV.kepler(pars.kep=mc[j3,1:Npar],tt=tt,kep.only=TRUE,nqp=nqp)$rv
            lines(tt,rvs,col=fit.col,lwd=1)
#            rvfun <- approxfun(tt,rvs)
#            tsim <- sort(runif(1e4,min(tt),max(tt)))
#            tsim <- seq(min(tt),max(tt),by=0.1)
#            lines(tsim,rvfun(tsim),col=fit.col,lwd=1)
        }
    }
    lines(tt,RV.kepler(pars.kep=par.opt,tt=tt,kep.only=TRUE,nqp=nqp)[[1]],col=fit.opt,lwd=line.size)$rv

####replot the RV data to be on top of the fitting lines
    for(i3 in 1:Nset){
        t <- JD.tot[[i1]][[i3]]%%24e5
        rv <- RV.tot[[i1]][[i3]] + res.tot[[i1]][[i3]] - res.sig[[i1]][[i3]]
        erv <- eRV.tot[[i1]][[i3]]
        points(t,rv,pch=20,col=cols[i3],cex=pch.size)
        arrows(t,rv-erv,t,rv+erv,length=0.05,angle=90,code=3,col=cols[i3])
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
        arrows(tlist[[j]]%%24e5,drv[[j]]-erv.list[[j]],tlist[[j]]%%24e5,drv[[j]]+erv.list[[j]],length=0.05,angle=90,code=3,col=tcol(cols[j],50))
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


