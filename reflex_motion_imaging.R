source('mcmc_func.R')
###add function
modifypar <- function(par){
    Nsig <- length(grep('^per',names(par)))
    ps <- exp(par[paste0('per',1:Nsig)])
    ind <- which(ps<1000)
    if(length(ind)>0){
        par[paste0('K',ind)] <- 0
    }
    par
}
#par.opt1 <- modifypar(par.opt)
par.opt1 <- par.opt
###
ind.show <- which.max(par.opt[grep('per',names(par.opt))])
yr2d <- 365.25
##astrometry fit
data.astrometry <- out$astrometry
if(!any(names(out)=='shortP')){
   shortP <- FALSE
}else{
   shortP <- out$shortP
}
tt  <-  RV.kepler(pars.kep=par.opt1,bases=bases,nqp=nqp)$tt
if(!any(names(out)=='plx')) out$plx <- out$astrometry$parallax[2]
astro0 <- astrometry.kepler(par.opt1,bases=bases,tt=tt)$barycenter
tsim <- e0$tsim
Nt <- length(tsim)
planet <- astrometry.kepler(par.opt1,tt=tsim,bases=bases)$planet
dpmra <- data.astrometry[,'pmra']-astro0[,'pmra']
dpmdec <- data.astrometry[,'pmdec']-astro0[,'pmdec']
dpmras <- c()
dpmdecs <- c()
for(j1 in 1:Nastro){
    ell <- error.ellipse(dpmra[j1],dpmdec[j1],cov.astro[4:5,4:5,j],percent=68)
    dpmras <- cbind(dpmras,ell[,1])
    dpmdecs <- cbind(dpmdecs,ell[,2])
}

if(plotf){
screen(1+(i1-1)*5+2)
par(mar=c(margin,2*margin,margin,margin))
plot(dpmra,dpmdec,xlab='',ylab='',xlim=range(planet[,'pmra'],dpmras,dpmra),ylim=range(planet[,'pmdec'],dpmdecs,dpmdec),pch=ell.pch,xaxt='n',yaxt='n',col=ell.col)#main='Proper motion fit'
magaxis(1:2)

###MC lines
Nsamp <- 100
if(Nmc>0){
mc <- out$mcmc.opt[[paste0('sig',Nkep)]]
inds <- sample(1:nrow(mc),Nsamp)
Npar <- length(par.opt)
dra.sim <- ddec.sim <- c()
for(j3 in inds){
    par <- mc[j3,1:Npar]
#    par1 <- modifypar(par)
    par1 <- par
    planet.mc <- astrometry.kepler(par1,tt=tsim,bases=bases)$planet
    lines(planet.mc[,'pmra'],planet.mc[,'pmdec'],col=fit.col)
    dra.sim <- cbind(dra.sim,planet.mc[,'ra'])
    ddec.sim <- cbind(ddec.sim,planet.mc[,'dec'])
}
}
###
if(i1==1) mtext(side=1,text=expression(Delta*mu[alpha]~'[mas/yr]'),xpd=NA,line=2,cex=1.2*size)
mtext(side=2,text=expression(Delta*mu[delta]~'[mas/yr]'),xpd=NA,line=1.5,cex=1.2*size)
if(dpmra[2]>dpmra[1] & dpmdec[2]>dpmdec[1]) pos <- c(4,2)
if(dpmra[2]<dpmra[1] & dpmdec[2]<dpmdec[1]) pos <- c(2,4)
if(dpmra[2]<dpmra[1] & dpmdec[2]>dpmdec[1]) pos <- c(3,1)
if(dpmra[2]>dpmra[1] & dpmdec[2]<dpmdec[1]) pos <- c(1,3)
lines(planet[,'pmra'],planet[,'pmdec'],col=fit.opt[1],lwd=line.size)
lines(c(planet[1,'pmra'],dpmra[1]),c(planet[1,'pmdec'],dpmdec[1]),col=ell.col[1])
lines(c(planet[Nt,'pmra'],dpmra[2]),c(planet[Nt,'pmdec'],dpmdec[2]),col=ell.col[2])
}
for(j1 in 1:Nastro){
    lines(dpmras[,j],dpmdecs[,j],col=ell.col[j1])
}
text(x=dpmra,y=dpmdec,labels=c('H','G'),col=ell.col,pos=pos)
points(0,0,pch='+')

Popt <- exp(par.opt[paste0('per',ind.show)])
Kopt <- par.opt[paste0('K',ind.show)]
eopt <- par.opt[paste0('e',ind.show)]
Iopt <- par.opt[paste0('Inc',ind.show)]%%(pi)
Mopt <- par.opt[paste0('Mo',ind.show)]
Oopt <- par.opt[paste0('Omega',ind.show)]
oopt <- par.opt[paste0('omega',ind.show)]
msini <- k2m(K,P,e,Mstar,Inc=NULL)
mp.opt <- Mp <- msini$mj/sin(Iopt)
Mo <- (Mopt%%(2*pi))*180/pi
inc <- Iopt*180/pi

#####plot position ellipse
dras <- c()
ddecs <- c()
dra <- (data.astrometry[,'ra']-astro0[,'ra'])*3.6e6*cos(data.astrometry[,'dec']/180*pi)#mas
ddec <- (data.astrometry[,'dec']-astro0[,'dec'])*3.6e6#mas
for(j1 in 1:Nastro){
    ell <- error.ellipse(dra[j1],ddec[j1],cov.astro[1:2,1:2,j],percent=68)
    dras <- cbind(dras,ell[,1])
    ddecs <- cbind(ddecs,ell[,2])
}
if(plotf){
screen(1+(i1-1)*5+3)
par(mar=c(margin,2*margin,margin,margin))
xlim <- range(planet[,'ra'],dras)
dxlim <- xlim[2]-xlim[1]
xlim <- rev(c(xlim[1]-0.01*dxlim,xlim[2]+0.01*dxlim))

ylim <- range(planet[,'dec'],ddecs)
dylim <- ylim[2]-ylim[1]
ylim <- c(ylim[1]-0.01*dylim,ylim[2]+0.01*dylim)


plot(dra,ddec,xlab='',ylab='',xlim=xlim,ylim=ylim,pch=ell.pch,xaxt='n',yaxt='n',col=ell.col)#main='Position fit'
magaxis(1:2)
if(Nmc>0){
for(j3 in 1:Nsamp){
    lines(dra.sim[,j3],ddec.sim[,j3],col=fit.col)
}
}
if(i1==1) mtext(side=1,text=expression(Delta*alpha*'* [mas]'),xpd=NA,line=2,cex=1.2*size)
mtext(side=2,text=expression(Delta*delta~'[mas]'),xpd=NA,line=1.5,cex=1.2*size)
if(dra[2]<dra[1] & ddec[2]>ddec[1]) pos <- c(4,2)
if(dra[2]>dra[1] & ddec[2]<ddec[1]) pos <- c(2,4)
if(dra[2]>dra[1] & ddec[2]>ddec[1]) pos <- c(3,1)
if(dra[2]<dra[1] & ddec[2]<ddec[1]) pos <- c(1,3)
lines(planet[,'ra'],planet[,'dec'],col=fit.opt[1],lwd=line.size)
lines(c(planet[1,'ra'],dra[1]),c(planet[1,'dec'],ddec[1]),col=ell.col[1])
lines(c(planet[Nt,'ra'],dra[2]),c(planet[Nt,'dec'],ddec[2]),col=ell.col[2])
for(j2 in 1:Nastro){
    lines(dras[,j],ddecs[,j],col=ell.col[j2])
}
Ndig <- 2
text(x=dra,y=ddec,labels=c('H','G'),col=ell.col,pos=pos)
points(0,0,pch='+')
if(!imaging){
legend('topright',inset=c(-0.5,0),xpd=NA,legend=as.expression(c(bquote(m[c]==.(round(Mp,1))~M[Jup]),bquote(P==.(round(Popt/yr2d,Ndig))~'yr'),bquote(e==.(round(eopt,Ndig))),bquote(I==.(round(inc,Ndig-1))~'deg'),bquote(omega==.(round(oopt*180/pi,Ndig-1))~'deg'),bquote(Omega==.(round(Oopt*180/pi,Ndig-1))~'deg'),bquote(M[0]==.(round(Mo,Ndig-1))~'deg'),bquote(lnJ[hip]==.(round(par.opt['logJ.hip'],Ndig))),bquote(lnJ[gaia]==.(round(par.opt['logJ.gaia'],Ndig))))),bty='n',title='Parameters')
}
}

nc <- 1
#cc <- c('black','blue','orange','steelblue','green','darkgreen',rainbow(10))
target <- star
popts <- exp(par.opt[paste0('per',1:Nsig)])
if(any(names(out)=='rel')){
####simulation
    trel <- out$trel
    ##    Dt <- Popt[j]-(max(trel)-min(trel))
    ##    Dt <- max(0,(max(trel)-min(trel))-Popt)
    Dt <- max(Popt)
#    tsim <- seq(max(24e5,min(trel)-Dt/2),max(max(trel)+Dt/2,sum(time_Yr2jd(2020))+Dt/2),length.out=1e4)
    ##    tsim <- seq(min(trel)-100,max(trel)+100,length.out=1e3)
    ##    tsim <- seq(min(yrs),max(yrs))
    planet <- astrometry.rel(par.opt,tsim=tsim)
    range.raS <- range(unlist(planet$rel$raB))
    range.decS <- range(unlist(planet$rel$decB))
    range.raP <- range(unlist(lapply(names(out$rel),function(n) unlist(lapply(names(out$rel[[n]]),function(i) out$rel[[n]][[i]][,'dra'])))))
    range.decP <- range(unlist(lapply(names(out$rel),function(n) unlist(lapply(names(out$rel[[n]]),function(i) out$rel[[n]][[i]][,'ddec'])))))
####
screen(1+(i1-1)*5+4)
par(mar=c(margin,2*margin,margin,margin))
#pdf('paper_astrometry.pdf',6,6)
#par(mar=c(5,5,1,1))
    for(j1 in 1:Nsig){
        n <- paste0('p',j)
        popt <- exp(par.opt[paste0('per',j1)])
        show.epoch <- FALSE
        if(popt==max(popts)){
            show.epoch <- TRUE
            yrs <- time_Jd2yr(cbind(tsim,0))
#            tts <- seq(max(floor(min(yrs)),1990),min(ceiling(max(yrs)),2030),length.out=5)
            tts <- seq(2020,ceiling((2020+Dt/2/365.25)/10)*10,length.out=3)
            tts <- unique(round(tts))
                                        #    tts <- tts[tts>1990 & tts<2030]
                                        #    tts <- tts[tts>2010 & tts<2030]
                                        #    tts <- c(tts,2022.1)
            ii <- unlist(lapply(tts, function(tt) which.min(abs(yrs-tt))))
            jjs <- trel[ii]
        }
        if(out$Nrel[j1]>0 & length(out$rel[[n]])>0){
            n <- paste0('p',j1)
            inss <- names(out$rel[[n]])
            inss <- inss[inss!='tot']

            tPs <- draPs <- ddecPs <- c()
            for(i in inss){
                draPs <- c(draPs,out$rel[[n]][[i]][,'dra'])
                ddecPs <- c(ddecPs,out$rel[[n]][[i]][,'ddec'])
                tPs <- c(tPs,out$rel[[n]][[i]][,1])
            }
            draS <- planet$rel$raB[[n]]
            ddecS <- planet$rel$decB[[n]]
#            xlim <- 1.2*range(draPs,draS)
#            ylim <- 1.2*range(ddecPs,ddecS)
            xlim <- range(range.raP,range.raS)
            ylim <- range(range.decP,range.decS)
#            xlim <- ylim <- range(xlim,ylim)
            if(nc==1){
                plot(draPs,ddecPs,xlab='',ylab='',xlim=rev(xlim),ylim=ylim,col=NA)
                points(0,0,pch='+',col='black')
            }else{
                points(draPs,ddecPs,xlab='',col=NA)
            }
            for(i in inss){
                crel <- col.rel[i==names(col.rel)]
                cn <- colnames(out$rel[[n]][[i]])
                tP <- out$rel[[n]][[i]][,1]
                draP <- out$rel[[n]][[i]][,'dra']
                ddecP <- out$rel[[n]][[i]][,'ddec']
                edraP <- out$rel[[n]][[i]][,grep('era|edra',cn)]
                eddecP <- out$rel[[n]][[i]][,grep('edec|eddec',cn)]
                points(draP,ddecP,col=crel,cex=0.8,pch=rel.pch)
                inds <- c()
                for(k in 1:length(draP)){
###                inds <- c(inds,which.min(abs(out$rel[[n]]$tot[k,1]-tsim)))
                    inds <- c(inds,which.min(abs(tP[k]-tsim)))
                }
                arrows(draS[inds],ddecS[inds],draP,ddecP,length=0.001,angle=0,code=1,col='grey')
                arrows(draP-edraP,ddecP,draP+edraP,ddecP,length=0.01,angle=90,code=3,col=crel)
                arrows(draP,ddecP-eddecP,draP,ddecP+eddecP,length=0.01,angle=90,code=3,col=crel)
                nc <- nc+1
            }
#            if(j1==1) lines(draS,ddecS,col=fit.opt,lwd=line.size)
            lines(draS,ddecS,col=fit.opt[j1],lwd=line.size)
            if(show.epoch){
                points(draS[ii],ddecS[ii],pch='+',col=tcol('red',50))
                pos <- rep(2,length(ii))
                pos[which(draS[ii]>mean(draS[ii]))] <- 4
                text(draS[ii],ddecS[ii],labels=tts,pos=pos,xpd=NA,col=fit.opt[j1],cex=0.8)
            }

###add additional companion orbit for reference
#            if(length(out$rel)>1){
            if(FALSE){
                kk <- (1:length(out$rel))[-j]
                for(k in kk){
                    ins <- names(out$rel)[k]
                    crel <- col.rel[ins==names(col.rel)]
                    draS <- planet$rel$raB[[k]]
                    ddecS <- planet$rel$decB[[k]]
                    lines(draS,ddecS,col='grey')
                    if(show.epoch){
                    points(draS[ii],ddecS[ii],pch='+',col=crel)
                    pos <- rep(2,length(ii))
                    pos[which(draS[ii]>mean(draS[ii]))] <- 4
                    text(draS[ii],ddecS[ii],labels=tts,pos=pos,xpd=NA,col=cc[k],cex=0.8)
                    }
                }
                if(target=='HD42581' & FALSE){
                    t0 <- sum(time_Yr2jd(2022))
                    t1 <- sum(time_Yr2jd(2024))
                    jj <- which(tsim>t0&tsim<t1)
                    lines(draS[jj],ddecS[jj],col='purple',lwd=4)
                    rho <- sqrt(draS[jj]^2+ddecS[jj]^2)/1000#as
                    cat('rho=',rho,'as\n')
                }
            }
            mtext(side=1,text=expression(Delta*alpha[c]*'* [mas]'),xpd=NA,line=2,cex=1.2*size)
            mtext(side=2,text=expression(Delta*delta[c]~'[mas]'),xpd=NA,line=1.5,cex=1.2*size)
###add legend
            popt <- exp(par.opt[paste0('per',j)])
            K <- par.opt[paste0('K',j)]
            e <- par.opt[paste0('e',j)]
            inc <- par.opt[paste0('Inc',j)]*180/pi
            Mo <- (par.opt[paste0('Mo',j)]%%(2*pi))*180/pi
            omega <- par.opt[paste0('omega',j)]*180/pi
            Omega <- par.opt[paste0('Omega',j)]*180/pi
            msini <- K2msini.full(K,popt,e,Ms=out$Mstar)
            Mp <- msini$mj/sin(par.opt[paste0('Inc',j)])
            legend('topright',inset=c(-0.6,0),xpd=NA,legend=as.expression(c(bquote(m[c]==.(round(Mp,1))~M[Jup]),bquote(P==.(round(Popt/yr2d,Ndig))~'yr'),bquote(e==.(round(eopt,Ndig))),bquote(I==.(round(inc,Ndig-1))~'deg'),bquote(omega==.(round(oopt*180/pi,Ndig-1))~'deg'),bquote(Omega==.(round(Oopt*180/pi,Ndig-1))~'deg'),bquote(M[0]==.(round(Mo,Ndig-1))~'deg'),bquote(lnJ[hip]==.(round(par.opt['logJ.hip'],Ndig))),bquote(lnJ[gaia]==.(round(par.opt['logJ.gaia'],Ndig))))),bty='n',title='Parameters')
#            legend('top',legend=as.expression(c(bquote(P==.(round(popt/yr2d))~'yr'),bquote(M==.(round(Mp,1))~M[Jup]),bquote(e==.(round(e,2))),bquote(I==.(round(inc,1))~'deg'))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.12,x.intersp=0.2)
#            legend('top',legend=as.expression(c(bquote(omega==.(round(omega,1))~'deg'),bquote(Omega==.(round(Omega,1))~'deg'),bquote(M[0]==.(round(Mo,1))~'deg'),bquote(logJ[rel]==.(round(par.opt['logJ.rel'],1))))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.08,x.intersp=0.2)
        }
    }
#dev.off()
}

if(Nmc>0) source('orbit_predict.R')
