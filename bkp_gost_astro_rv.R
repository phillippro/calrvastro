library(magicaxis)
library(graphics)
library(RColorBrewer)
source('mcmc_func.R')
options(scipen = 0)
library(paletteer)
library("plotrix")
source('OrbitFunction.R')
source('timing_function.R')
source('general_function.R')
source('sofa_function.R')
set.seed(9999)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    star <- as.character(args[1])
}else{
#    star <- 'HD42581'
#    star <- 'HD39060'
#    star <- 'HD182488'
#    star <- 'HD215257'
#    star <- c('HD209100','HD22049')
#    star <- 'HD22049'
    star <- 'HD209100'
}
comptype <- 'companion'
Nmc <- Nsamp <- 0
m22 <- read.table('../data/combined/mamajek22.txt',header=TRUE)
ms <- as.numeric(m22[,'Msun'])
mg <- as.numeric(m22[,'M_G'])
ind <- which(!is.na(ms) & !is.na(mg))
mrl.m22 <- approxfun(ms[ind],mg[ind])
mlow.m22 <- min(ms[ind])
mup.m22 <- max(ms[ind])
version <- 1
set.seed(100)
ff <- read.table(paste0(comptype,'_files.txt'))[,1]
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
    source('mcmc_func.R')
    if(!exists('e0')) load(fs[j],env=e0 <- new.env())
#    load(fs[j],env=e0 <- new.env())
    eta0 <- e0$eta0
    tmin <- e0$tmin
    tsim <- e0$tsim
    if(max(tsim)-min(tsim)<max(e0$Popt)){
        tsim <- min(tmin,e0$out$astrometry[,1])+seq(0,max(e0$Popt),length.out=10000)
    }
    Nkep <- length(e0$Popt)
    Nsig <- e0$Nsig
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
    cat('Mstar=',Mstar,'\n')
    if(!exists('eMstar')) eMstars <- eMstar <- 0.1*Mstar
    cov.tot[[j]] <- e0$cov.astro
    ins.tot[[j]] <- ins <- e0$ins
    JD.tot[[j]] <- RV.tot[[j]] <- eRV.tot[[j]] <- res.tot[[j]] <- res.sig[[j]] <- nqp.tot[[j]] <- list()
    Nrv <- length(e0$ins)
    Nrv.tot[[j]] <- Nrv
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
ell.col <- c('grey','orange')
#ell.pch <- c(16,17)
#ell.pch <- c(20,20)
#ell.pch <- 18
ell.pch <- 20
rel.pch <- 15
#cols <- tcol(c('blue','green','purple','cyan','brown','orange','steelblue','yellow',brewer.pal(8,'Accent')),50)
#cols <- tcol(paletteer_c(palette = "viridis::inferno", n = 20),50)
cc <- brewer.pal.info[brewer.pal.info$category=='qual',]
cc1 <- unlist(mapply(brewer.pal,cc$maxcolors,rownames(cc)))
#set.seed(99999)
#cc2 <- tcol(cc1[c(1:11,14:15)],50)
cc2 <- tcol(cc1[c(1:3,16,5:11,14:15)],30)
#cc <- c('black','blue','orange','steelblue','green','darkgreen',rainbow(10))
#cc <- brewer.pal(8,'Accent')
#names(cols) <- c('AAT','APF1','ELODIE','')
#fpdf <- paste0(star0,'_planet1_astro_rv',Nsig,'v',version,'.pdf')
size <- 1
show.center <- FALSE
margin <- 1.3
mm <- c(margin,3*margin,2*margin,0)
pch.size <- 1
if(Nmc>0){
    line.size <- 2
}else{
#    line.size <- 3
    line.size <- 2
}
fit.col <- tcol('red',90)
#fit.opt <- tcol('red',50)
#fit.opt <- c('red','blue','green')
fit.opt <- tcol(c('red','red'),20)
#layout(matrix(c(1, 1,2,2), nrow=2, byrow=TRUE))
Nstar <- length(stars)
#Nstar <- Nsig
####RV
for(i1 in 1:Nstar){
###plotting settings
    j1 <- i1
    imaging <- FALSE
    if(grepl('hg3\\+|astro\\+',fs[i1]) & !grepl('HD39091',fs[i1])) imaging <- TRUE
    if(plotf){
        fpdf <- paste0(comptype,'_plot/',star,'_N',Nmc,'_fit.pdf')
        cat('output:\n',fpdf,'\n')
        if(!imaging){
            pdf(fpdf,12,3*1.2)
        }else{
            pdf(fpdf,14,3*1)
        }
        split.screen(c(1, 1))#screen 1 to Nsig
        if(FALSE){
            oma <- c(4,8,1,6)
        }else{
            oma <- c(3,7,0,7)
        }
        par(oma=oma,cex=size,cex.lab=1.2*size,cex.axis=size,mar=mm,mgp=c(2,0.5,0),cex.main=1.2*size)
    }
                                        #    out <- list()
    Mstar <- Mstars[j1]
                                        #    out[['astrometry']] <- astrometry.tot[[j1]]
    if(imaging){
        split.screen( figs = c(1,4), screen = i1 )
    }else{
        split.screen( figs = c(1,3), screen = i1 )
    }
    split.screen(rbind(c(0,1,0.3,1),c(0,1,0,0.3)),1+(i1-1)*5+1)
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
    pchs <- rep(20,Nw0)
    if(Nw0>10){
        pchs[1:round(Nw0/2)] <- 15
    }
    tlist <- JD.tot[[j1]]
    erv.list <- eRV.tot[[j1]]
###split the RV plot into two sub-panel
    if(plotf){
        if(imaging){
            screen(1+(i1-1)*5+5)
        }else{
            screen(1+(i1-1)*5+4)
        }
        ins <- instr <- ins.tot[[j1]]
        par(mar=c(0,mm[-1]))
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
        ind <- grep('per',names(pars))
        Popt <- exp(par.opt[ind])
        indI <- grep('Inc',names(pars))
        par.opt[indI] <- par.opt[indI]%%pi
        t3 <- seq(min(tall),max(tall),length.out=1e4)
        rvs <- RV.kepler(pars.kep=par.opt,tt=t3,kep.only=TRUE,nqp=nqp)[[1]]
        for(i3 in 1:Nset){
            t <- JD.tot[[1]][[i3]]%%24e5-tmin%%24e5
            rv <- RV.tot[[1]][[i3]] + res.tot[[j1]][[i3]] - res.sig[[j1]][[i3]]
            erv <- eRV.tot[[1]][[i3]]
            if(i3==1){
                plot(t,rv,xlab=paste0('JD-',round(tmin%%24e5+24e5,3)),ylab='',col=cc2[i3],pch=pchs[i3],xlim=c(min(tall%%24e5-tmin%%24e5),max(tall%%24e5-tmin%%24e5)),ylim=range(RV2.all,rvs),xaxt='n',yaxt='n',cex=1.2*pch.size,main='RV sets')#paste0(star,' RV')
                                        #            magaxis(c(1,2))
                                        #            magaxis(2,majorn=2,cex=size)
                magaxis(2,cex=size,majorn=3)
            }else{
                points(t,rv,pch=pchs[i3],col=cc2[i3],cex=pch.size)
            }
            arrows(t,rv-erv,t,rv+erv,length=0.05,angle=90,code=3,col=cc2[i3])
        }
        mtext(side=2,text=expression('RV [m/s]'),line=1.5,xpd=NA,cex=size)
        ind <- which(!grepl('HARPS',instr))
        if(length(ind)>0) instr[ind] <- toupper(instr[ind])
        ind <- grep('UCLES',instr)
        if(length(ind)>0) instr[ind] <- 'AAT'
        ind <- grep('SOPHIE2',instr)
        if(length(ind)>0) instr[ind] <- 'SOPHIEpost'
        ind <- grep('SOPHIE1',instr)
        if(length(ind)>0) instr[ind] <- 'SOPHIEpre'
        ind <- grep('SOPHIEPRE',instr)
        if(length(ind)>0) instr[ind] <- 'SOPHIEpre'
        ind <- grep('SOPHIEPOST',instr)
        if(length(ind)>0) instr[ind] <- 'SOPHIEpost'
        ind <- grep('HARPSPRE',instr)
        if(length(ind)>0) instr[ind] <- 'HARPSpre'
        ind <- grep('HARPSPOST',instr)
        if(length(ind)>0) instr[ind] <- 'HARPSpost'
        ind <- grep('Subaru',instr)
        if(length(ind)>0) instr[ind] <- 'HDS'
        ind <- grep('MCD',instr)
        if(length(ind)>0) instr[ind] <- 'HRS'
        ind <- grep('SUBARU',instr)
        if(length(ind)>0) instr[ind] <- 'HDS'
        ind <- grep('C07|C14|C98',instr)
        if(length(ind)>0) instr[ind] <- gsub('C','COR',instr[ind])
        ind <- grep('CORALIE07|CORALIE14|CORALIE98',instr)
        if(length(ind)>0) instr[ind] <- gsub('CORALIE','COR',instr[ind])
#        if(grepl('HD209100',stars[j1])){
#            instr <- c('HARPSpre','LC','VLC','HARPSpost','UVES')
#        }
        instr[instr=='L20'] <- 'AL19'
        if(exists('star.name')){
            sn <- gsub('GL','GJ',star.name)
        }else{
            sn <- gsub('GL','GJ',stars[j1])
        }
        if(!imaging){
            inset <- c(-0.75,-0.1)
        }else{
            inset <- c(-0.82,-0.1)
        }
#        legend('topleft',legend=c(instr,'Hipparcos','Gaia','Best fit'),col=c(cols[1:Nw0],ell.col,fit.opt,'black'),lty=c(rep(NA,Nw0+2),1,NA),pch=c(rep(20,Nw0),ell.pch,NA,3),xpd=NA,inset=inset,bty='n',title=sn,lwd=c(rep(NA,Nw0+2),line.size,NA),pt.cex=c(rep(1,Nw0),rep(1,2),NA,1),y.intersp=1.1)
        font.size <- 1
        size <- 0.8
        if(!imaging){
#            legend('topleft',legend=c(instr,'Hipparcos IAD','Gaia GOST','Best fit'),col=c(cols[1:Nw0],ell.col,fit.opt[1],'black'),lty=c(rep(NA,Nw0+2),1),pch=c(rep(20,Nw0),rep(ell.pch,2),NA,3),xpd=NA,inset=inset,bty='n',lwd=c(rep(NA,Nw0+2),line.size,NA),pt.cex=c(rep(font.size,Nw0),rep(font.size,2),NA,font.size),y.intersp=1.1,cex=size)#title=sn
            instr.show <- instr
            instr.show[instr.show=='APF'] <- 'APFp'
            instr.show[instr.show=='APFJ'] <- 'APFh'
            instr.show[instr.show=='KECK'] <- 'HIRESp'
            instr.show[instr.show=='KECKJ'] <- 'HIRESh'
            legend('topleft',legend=c(instr.show,'Best fit'),col=c(cc2[1:Nw0],fit.opt[1],'black'),lty=c(rep(NA,Nw0),1),pch=c(pchs,NA),xpd=NA,inset=inset,bty='n',lwd=c(rep(NA,Nw0),line.size),pt.cex=c(rep(font.size,Nw0),NA),y.intersp=1.1,cex=size)
        }else{
            ins0 <- unique(unlist(sapply(1:Nsig,function(n) names(out$rel[[n]]))))
            inss <- toupper(ins0)
            if(sn=='HD42581' & any(inss=='ASTROMETRY')) inss[inss=='ASTROMETRY'] <- 'MB21'
            if(sn=='HD182488' & any(inss=='ASTROMETRY')) inss[inss=='ASTROMETRY'] <- 'BB18'
            if(sn=='HD4113' & any(inss=='SPHERE')) inss[inss=='SPHERE'] <- 'AC18'
            if(sn=='GJ3677' & any(inss=='E10')) inss[inss=='E10'] <- 'DE10'
            if(sn=='HD13724' & any(inss=='R20')) inss[inss=='R20'] <- 'ER20'
            if(sn=='HD7449' & any(inss=='R16')) inss[inss=='R16'] <- 'TR16'
            if(sn=='GJ494' & any(inss=='M19')) inss[inss=='M19'] <- 'AM19'
            if(sn=='HD39060' & any(inss=='G20')) inss[inss=='G20'] <- 'SL21'
            if(sn=='HD39060' & any(inss=='HCI')) inss[inss=='HCI'] <- 'AL20'
            Nins <- length(inss)
            col.rel <- cc[1:Nins]
            names(col.rel) <- ins0
#            Nrel <- length(out$rel)
            Nrel <- 1
            legend('topleft',legend=c(instr,'Hipparcos IAD','Gaia GOST',inss,rep('Best fit',Nrel)),col=c(cc2[1:Nw0],ell.col,col.rel,fit.opt[1:Nrel],'black'),lty=c(rep(NA,Nw0+2+Nins),rep(1,Nrel),NA),pch=c(pchs,rep(ell.pch,2),rep(rel.pch,Nins),rep(NA,Nrel),3),xpd=NA,inset=inset,bty='n',lwd=c(rep(NA,Nw0+2+Nins),rep(line.size,Nrel),NA),pt.cex=c(rep(1,Nw0+2+Nins),rep(NA,Nrel),1),y.intersp=1.1)#title=sn,
        }
    }
    tspan <- max(tall)-min(tall)
    tt <- seq(min(tall),max(tall),length.out=1e4)
                                        #    if(Nmc>0){
    if(FALSE){
        t1 <- min(min(tall),mean(tall)-max(Popt)/2)
        t2 <- max(max(tall),mean(tall)+max(Popt)/2)
        tt <- seq(t1,t2,length.out=1e4)
    }
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
            ind.sample <- inds <- sample(1:nrow(mc),Nsamp)
            Npar <- length(par.opt)
            for(j3 in inds){
                rvs <- RV.kepler(pars.kep=mc[j3,1:Npar],tt=tt,kep.only=TRUE,nqp=nqp)$rv
                lines(tt%%24e5-tmin%%24e5,rvs,col=fit.col,lwd=2)
                                        #            rvfun <- approxfun(tt,rvs)
                                        #            tsim <- sort(runif(1e4,min(tt),max(tt)))
                                        #            tsim <- seq(min(tt),max(tt),by=0.1)
                                        #            lines(tsim-tmin,rvfun(tsim),col=fit.col,lwd=1)
            }
        }
        lines(tt%%24e5-tmin%%24e5,RV.kepler(pars.kep=par.opt,tt=tt,kep.only=TRUE,nqp=nqp)[[1]],col=fit.opt[1],lwd=line.size)$rv

####replot the RV data to be on top of the fitting lines
        if(FALSE){
            for(i3 in 1:Nset){
                t <- JD.tot[[j1]][[i3]]%%24e5-tmin%%24e5
                rv <- RV.tot[[j1]][[i3]] + res.tot[[j1]][[i3]] - res.sig[[j1]][[i3]]
                erv <- eRV.tot[[j1]][[i3]]
                points(t,rv,pch=20,col=cc2[i3],cex=pch.size)
                arrows(t,rv-erv,t,rv+erv,length=0.05,angle=90,code=3,col=cc2[i3])
            }
        }

###residual plot
        if(imaging){
            screen(1+(i1-1)*5+6)
        }else{
            screen(1+(i1-1)*5+5)
        }
        par(mar=c(mm[1:2],0,mm[4]))
        rv.planet <- RV.kepler(pars.kep=par.opt,tt=NULL,kep.only=TRUE,nqp=nqp)$rv
        drv <- res.tot[[j1]]
                                        #    for(j in 1:Nset){
                                        #        drv[[j]] <- rv.list[[j]]-rv.planet[[j]]
                                        #    }
        ylim <- range(unlist(drv))
        xlim <- range(tall)%%24e5-tmin%%24e5
        for(j2 in 1:Nset){
            x <- tlist[[j2]]%%24e5-tmin%%24e5
            y <- drv[[j2]]
            ey <- erv.list[[j2]]
            if(j2==1){
                plot(x,y,xlab='',ylab='',ylim=ylim,xlim=xlim,col=tcol(cc2[j2],50),pch=pchs[j2],cex=size,xaxt='n',yaxt='n')
                                        #            axis(1,tcl = 0.5)
                                        #            magaxis(1,majorn=2,cex=size)
                                        #            magaxis(2,majorn=2,cex=size)
                magaxis(1,majorn=3,cex=size)
                magaxis(2,majorn=2,cex=size)
                                        #            axis(2,tck=0.1)
            }else{
                points(x,y,col=cc2[j2],pch=20,cex=1)
            }
            arrows(x,y-ey,x,y+ey,length=0.05,angle=90,code=3,col=tcol(cc2[j2],50))
        }
                                        #    par(mgp=c(3,1,0),xpd=NA)
        if(i1==1) mtext(side=1,text=paste0('JD-',round(tmin%%24e5+24e5,3)),xpd=NA,line=2,cex=1.2*size)
        mtext(side=2,text=expression('O-C [m/s]'),line=1.5,xpd=NA,cex=1.2*size)
        abline(h=0,lty=2,lwd=2)
    }
####astrometry
    source('paper_ellipse_plot.R')
if(plotf){
close.screen(all.screens = TRUE)
dev.off()
}
}



