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
Ein  <-  RV.kepler(pars.kep=par.opt1,bases=bases,nqp=nqp)$E
#astro0 <- astrometry.kepler(par.opt1,bases=bases,Ein=Ein)$barycenter
if(!any(names(out)=='plx')) out$plx <- out$astrometry$parallax[2]
astro0 <- astrometry.kepler(par.opt1,bases=bases,Ein=Ein,pa=FALSE)$barycenter
if(star=='HD80869'){
Nsim <- 1e5
}else{
Nsim <- 1e4
}
tsim0 <- seq(data.astrometry[1,1],data.astrometry[2,1],length.out=Nsim)
tsim <- tsim0-tmin
Nt <- length(tsim)
planet <- astrometry.kepler(par.opt1,tt=tsim,bases=bases,pa=FALSE)$planet
dpmra <- data.astrometry[,'pmra']-astro0[,'pmra']
dpmdec <- data.astrometry[,'pmdec']-astro0[,'pmdec']
dpmras <- c()
dpmdecs <- c()
for(j in 1:Nastro){
    ell <- error.ellipse(dpmra[j],dpmdec[j],cov.astro[4:5,4:5,j],percent=68)
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
    planet.mc <- astrometry.kepler(par1,tt=tsim,bases=bases,pa=FALSE)$planet
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
for(j in 1:Nastro){
    lines(dpmras[,j],dpmdecs[,j],col=ell.col[j])
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
msini <- K2msini.full(Kopt,Popt,eopt,Ms=Mstar)
mp.opt <- Mp <- msini$mj/sin(Iopt)
Mo <- (Mopt%%(2*pi))*180/pi
inc <- Iopt*180/pi

#####plot position ellipse
dras <- c()
ddecs <- c()
dra <- (data.astrometry[,'ra']-astro0[,'ra'])*3.6e6*cos(data.astrometry[,'dec']/180*pi)#mas
ddec <- (data.astrometry[,'dec']-astro0[,'dec'])*3.6e6#mas
for(j in 1:Nastro){
    ell <- error.ellipse(dra[j],ddec[j],cov.astro[1:2,1:2,j],percent=68)
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
for(j in 1:Nastro){
    lines(dras[,j],ddecs[,j],col=ell.col[j])
}
Ndig <- 2
text(x=dra,y=ddec,labels=c('H','G'),col=ell.col,pos=pos)
points(0,0,pch='+')
legend('topright',inset=c(-0.5,0),xpd=NA,legend=as.expression(c(bquote(m[c]==.(round(Mp,1))~M[Jup]),bquote(P==.(round(Popt/yr2d,Ndig))~'yr'),bquote(e==.(round(eopt,Ndig))),bquote(I==.(round(inc,Ndig-1))~'deg'),bquote(omega==.(round(oopt*180/pi,Ndig-1))~'deg'),bquote(Omega==.(round(Oopt*180/pi,Ndig-1))~'deg'),bquote(M[0]==.(round(Mo,Ndig-1))~'deg'),bquote(lnJ[hip]==.(round(par.opt['logJ.hip'],Ndig))),bquote(lnJ[gaia]==.(round(par.opt['logJ.gaia'],Ndig))))),bty='n',title='Parameters')
}

if(Nmc>0) source('orbit_predict.R')
