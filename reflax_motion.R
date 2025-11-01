Nastro <- nrow(data.astrometry)
tsim0 <- seq(data.astrometry[1,1],data.astrometry[2,1],length.out=1e3)%%2400000
tsim <- tsim0-tmin%%2400000
ra.bary <- dec.bary <- pmra.bary <- pmdec.bary <- ra.ref <- dec.ref <- pmra.ref <- pmdec.ref <- ra.planet <- dec.planet <- pmra.planet <- pmdec.planet <- c()
ra.planet2 <- dec.planet2 <- pmra.planet2 <- pmdec.planet2 <- c()
for(j in 1:ncol(pars)){
    tmp <- astrometry.kepler(pars[,j],tt=tsim)$out
    rab <- tmp[,1]-tmp[,6]/3.6e6
    decb <- tmp[,2]-tmp[,7]/3.6e6
    pmrab <- tmp[,4]-tmp[,8]
    pmdecb <- tmp[,5]-tmp[,9]
    ra.bary <- cbind(ra.bary,rab)#deg
    dec.bary <- cbind(dec.bary,decb)
    pmra.bary <- cbind(pmra.bary,pmrab)
    pmdec.bary <- cbind(pmdec.bary,pmdecb)
    ra.ref <- cbind(ra.ref,rab[c(1,length(rab))])#deg
    dec.ref <- cbind(dec.ref,decb[c(1,length(decb))])
    pmra.ref <- cbind(pmra.ref,pmrab[c(1,length(pmrab))])
    pmdec.ref <- cbind(pmdec.ref,pmdecb[c(1,length(pmdecb))])
    ra.planet <- cbind(ra.planet,tmp[,6])#mas
    dec.planet <- cbind(dec.planet,tmp[,7])
    pmra.planet <- cbind(pmra.planet,tmp[,8])
    pmdec.planet <- cbind(pmdec.planet,tmp[,9])
    dra <- (ra.bary[,j]-ra.bary[,1])*3.6e6
    ddec <- (dec.bary[,j]-dec.bary[,1])*3.6e6
    ra.planet2 <- cbind(ra.planet2,ra.planet[,j]+dra)
    dec.planet2 <- cbind(dec.planet2,dec.planet[,j]+ddec)
    dpmra <- pmra.bary[,j]-pmra.bary[,1]
    dpmdec <- pmdec.bary[,j]-pmdec.bary[,1]
    pmra.planet2 <- cbind(pmra.planet2,pmra.planet[,j]+dpmra)
    pmdec.planet2 <- cbind(pmdec.planet2,pmdec.planet[,j]+dpmdec)
}
psi <- 0.8

pmras <- c()
pmdecs <- c()
jitters <- par.opt[c('jitter.hip','jitter.gaia')]
for(j in 1:Nastro){
#    ell <- error.ellipse(data.astrometry[j,'pmra']-pmra.ref[j,1],data.astrometry[j,'pmdec']-pmdec.ref[j,1],cov.astro[4:5,4:5,j]*(1+jitters[j]),percent=90)
    ell <- error.ellipse(data.astrometry[j,'pmra']-pmra.ref[j,1],data.astrometry[j,'pmdec']-pmdec.ref[j,1],cov.astro[4:5,4:5,j],percent=90)
    pmras <- cbind(pmras,ell[,1])
    pmdecs <- cbind(pmdecs,ell[,2])
}
pmra <- data.astrometry[,'pmra']-pmra.ref[,1]
pmdec <- data.astrometry[,'pmdec']-pmdec.ref[,1]
if(Nmc>0){
    xlim <- range(pmra,pmras,pmra.planet2,0)
    ylim <- range(pmdec,pmdecs,pmdec.planet2,0)
}else{
    xlim <- range(pmra,pmras,pmra.planet2[,1],0)
    ylim <- range(pmdec,pmdecs,pmdec.planet2[,1],0)
}
screen(Nsig+(i1-1)*5+2)
par(mar=c(margin,2*margin,margin,margin))
plot(pmra,pmdec,xlab=expression(mu[alpha]-mu[alpha]^{linear}*' [mas/yr]'),ylab='',xlim=xlim,ylim=ylim,pch=20,xaxt='n',yaxt='n')
magaxis(1:2)
for(j in 1:Nastro){
#    lines(pmras[,j],pmdecs[,j],col=ell.col[j])
    polygon(pmras[,j],pmdecs[,j],col=ell.col[j],border=tcol('white',100))
}
lines(pmra.planet2[,1],pmdec.planet2[,1],col=fit.col,lwd=line.size)
if(Nmc>0){
    for(j in 2:ncol(pars)){
        lines(pmra.planet2[,j],pmdec.planet2[,j],col=fit.col,lwd=line.size)
    }
}
points(pmra,pmdec,col='black',pch=20)
if(show.center){
points(0,0,pch='+',cex=2)
}
msini <- K2msini(par.opt['K1'],Popt,par.opt['e1'],Ms=data.astrometry[1,'mass'])

Mp <- msini$mj/sin(par.opt['Inc1'])
Mo <- (par.opt['Mo1']%%(2*pi))*180/pi
inc <- par.opt['Inc1']*180/pi
#    text(x=psi*pmra,y=psi*pmdec,labels=c('Hipparcos','Gaia'),offset=c(2,1),col=text.col)
if(i1==Nsig) mtext(side=1,text=expression(mu[alpha]-mu[alpha]^{linear}~'[mas/yr]'),xpd=NA,line=2,cex=1.2*size)
mtext(side=2,text=expression(mu[delta]-mu[delta]^{linear}~'[mas/yr]'),xpd=NA,line=1.5,cex=1.2*size)

#####plot position ellipse
screen(Nsig+(i1-1)*5+3)
par(mar=c(margin,2*margin,margin,margin))
if(astrometry>=7){
    ras <- c()
    decs <- c()
    for(j in 1:Nastro){
        ra.mas <- (data.astrometry[j,'ra']-ra.ref[j,1])*3.6e6
        dec.mas <- (data.astrometry[j,'dec']-dec.ref[j,1])*3.6e6
        ell <- error.ellipse(ra.mas,dec.mas,cov.astro[1:2,1:2,j]*(1+jitters[j]),percent=90)
        ras <- cbind(ras,ell[,1])#mas
        decs <- cbind(decs,ell[,2])
    }
##barycentric position
    ra <- (data.astrometry[,'ra']-ra.ref[,1])*3.6e6#mas
    dec <- (data.astrometry[,'dec']-dec.ref[,1])*3.6e6#mas
    if(Nmc>0){
        xlim <- range(ra,ras,ra.planet2,0)
        ylim <- range(dec,decs,dec.planet2,0)
    }else{
        xlim <- range(ra,ras,ra.planet2[,1],0)
        ylim <- range(dec,decs,dec.planet2[,1],0)
    }
    plot(ra,dec,xlab=expression(alpha-alpha^{linear}~'[mas]'),ylab='',pch=20,xaxt='n',yaxt='n',xlim=xlim,ylim=ylim)
if(show.center){
    points(0,0,pch='+',cex=2)
}
    magaxis(1:2)
    for(j in 1:Nastro){
#        lines(ras[,j],decs[,j],col=ell.col[j])
        polygon(ras[,j],decs[,j],col=ell.col[j],border=tcol('white',100))
    }
    lines(ra.planet2[,1],dec.planet2[,1],col=fit.col,lwd=line.size)
    if(Nmc>0){
        for(j in 2:ncol(pars)){
            lines(ra.planet2[,j],dec.planet2[,j],col=fit.col,lwd=line.size)
        }
    }
    points(ra,dec,col='black',pch=20)
    if(i1==Nsig) mtext(side=1,text=expression(alpha-alpha^{linear}~'[mas]'),xpd=NA,line=2,cex=1.2*size)
    mtext(side=2,text=expression(delta-delta^{linear}~'[mas]'),xpd=NA,line=1.5,cex=1.2*size)
    legend('topright',inset=c(-0.5,0),xpd=NA,legend=as.expression(c(bquote(m[p]==.(round(Mp,1))~M[Jup]),bquote(P==.(round(Popt/yr2d,1))~'yr'),bquote(e==.(round(par.opt['e1'],1))),bquote(I==.(round(inc,1))~'deg'),bquote(omega==.(round(par.opt['omega1']*180/pi,1))~'deg'),bquote(Omega==.(round(par.opt['Omega1']*180/pi,1))~'deg'),bquote(M[0]==.(round(Mo,1))~'deg'),bquote(J[hip]==.(round(par.opt['jitter.hip'],1))),bquote(J[gaia]==.(round(par.opt['jitter.gaia'],1))))),bty='n',title='Parameters')
#    text(x=psi*ra,y=psi*dec,labels=c('Hipparcos','Gaia'),offset=c(2,1),col=text.col)
}
