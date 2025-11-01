####################
###plot
####################
if(!exists('out')){
    load('/malibu/ffeng/astro_output/HD74014/HD74014_fix1_relativityFALSE_Niter1000000_Ncores8_ofac2_Nset3_Esd1_transit0_P7044_acc9.2_lnlmax-124.Robj')
}
pdf(pdf.file,8,8)
#par(mar=c(5,5,1,8),cex.lab=size,cex.axis=size)
par(mfrow=c(2,2))
###combined data
plot(trv.all,RV.all,xlab='BJD-2400000',ylab='RV [m/s]',main='',col='white')
for(i in 1:length(ins)){
    points(out[[ins[i]]]$RV[,1],out[[ins[i]]]$RV[,2],col=colors[i])
        try(arrows(out[[ins[i]]]$RV[,1],out[[ins[i]]]$RV[,2]-out[[ins[i]]]$RV[,3],out[[ins[i]]]$RV[,1],out[[ins[i]]]$RV[,2,]+out[[ins[i]]]$RV[,3],length=0.05,angle=90,code=3,col=colors[i]),TRUE)
	}
#try(arrows(trv,RV-eRV,trv,RV+eRV,length=0.05,angle=90,code=3,col='grey'),TRUE)
n1 <- ceiling(length(ins)/2)
n2 <- length(ins)-n1
legend('top',legend=ins[1:n1],col=colors[1:n1],inset=c(0,-0.2),seg.len = 1,horiz=TRUE,xpd=NA,bty='n',pch=1)
if(n2>0){
    legend('top',legend=ins[(n1+1):length(ins)],col=colors[(n1+1):length(ins)],seg.len = 1,inset=c(0,-0.1),horiz=TRUE,xpd=NA,bty='n',pch=1)
}
    
    RV4c <- datab[[k1]][,Np+i3]
    eRV4c <- datab[[k1]][,2*Np+i3]
    t5c <- datas[[k1]][,i3]
    RV5c <- datas[[k1]][,Np+i3]
    cols <- c('black','blue','green','cyan','orange','brown')
    pchc <- 1:Ncol
    Nf <- Nf+1
    ylim <- c(-40,40)
    xlim <- range(trv.all)
#        ylim <- range(RV3c[!is.na(eRV3c)],RV3c[!is.na(eRV3c)],RV4c[!is.na(eRV4c)]-eRV4c[!is.na(eRV4c)],RV4c[!is.na(eRV4c)]+eRV4c[!is.na(eRV4c)],RV5c)
        for(k2 in 1:Nw){
            t <- par.data[[k2]]$trv
            cat('t span:',range(t),'\n')
            rv <- par.data[[k2]]$RV2
            erv<- par.data[[k2]]$eRV
            if(k2==1){
#                rspan <- max(RV)-min(RV)
                plot(t,rv,xlab=paste0('JD-2400000 [day]'),ylab=ylab,pch=pch,cex=pch.size,col=cols[k2],ylim=ylim,xlim=xlim,xaxt='n')#col='black',
                ts <- seq(xlim[1],xlim[2],by=1)
                lines(ts,RV.kepler(pars.kep=par.opt,tt=ts,kep.only=TRUE)[[1]],col="red")
            }else{
                points(t,rv,cex=pch.size,col=cols[k2],pch=pch)
            }
            if(err){
                arrows(t,rv-erv,t,rv-erv,length=0.1,angle=90,code=3,col=cols[k2],cex=size,pch=pch)
            }else{
                arrows(t3c,RV3c-eRV,t3c,RV3c+eRV,length=0.03,angle=90,code=3,col='grey',cex=size,pch=pch)
                points(t3c,RV3c,pch=pch,cex=pch.size)
            }
            points(t5c,RV5c,col="red",pch=pch,cex=pch.size)
                                        #        legend('topright',legend=paste0(i3,'planet'),bty='n')
            if(any(i3==c(1))){
                legend('topleft',legend=paste(format(Popt[i3],digit=3),'d'),bty='n',cex=size)
            }else{
                legend('topleft',legend=paste(format(Popt[i3],digit=3),'d'),bty='n',cex=size)
            }
        }
    legend('topright',inset=c(-0.46,0),legend=ins[c(2:3,1,4,5)],pch=rep(pch,Nw),pt.cex=rep(1,Nw),col=cols[c(2:3,1,4,5)],xpd=NA,cex=size)

    par(mar=c(5,5,0,8))
     for(k2 in 1:Nw){
            t <- par.data[[k2]]$trv
            cat('t span:',range(t),'\n')
            rv <- par.data[[k2]]$res
            erv<- par.data[[k2]]$eRV
            if(k2==1){
#                rspan <- max(RV)-min(RV)
                plot(t,rv,xlab=paste0('JD-2400000 [day]'),ylab='O-C [m/s]',cex=pch.size,col=cols[k2],ylim=ylim,xlim=xlim,yaxt='n',pch=pch)#col='black',
                p <- axis(side=2,at=c(-20,0,20))
                abline(h=0,col='grey',lty=2)
            }else{
                points(t,rv,cex=pch.size,col=cols[k2],pch=pch)
            }
        }
}
source('OrbitFunction.R')
par(mar=c(5,5,1,1),cex.lab=size,cex.axis=size)
layout(matrix(c(1,1,1,1),2,2, byrow=TRUE))
tsim0 <- seq(data.astrometry[1,1],data.astrometry[2,1],by=1)%%2400000
tsim <- tsim0-min(trv.all)%%2400000
tmp <- astrometry.kepler(par.opt,tt=tsim)
pmras <- c()
pmdecs <- c()
for(j in 1:Nastro){
    ell <- error.ellipse(data.astrometry[j,'pmra'],data.astrometry[j,'pmdec'],cov.PPM[,,j],percent=68)
    pmras <- cbind(pmras,ell[,1])
    pmdecs <- cbind(pmdecs,ell[,2])
}
pmra <- data.astrometry[,'pmra']
pmdec <- data.astrometry[,'pmdec']
plot(pmra,pmdec,xlab=expression(mu[alpha]*' [mas/yr]'),ylab=expression(mu[delta]*' [mas/yr]'),xlim=range(tmp[,2],pmras,data.astrometry[,'pmra']),ylim=range(tmp[,3],pmdecs,data.astrometry[,'pmdec']),pch=20)
for(j in 1:Nastro){
    lines(pmras[,j],pmdecs[,j],col='grey')
}
lines(tmp[,2],tmp[,3],col='red')
msini <- K2msini(par.opt['K1'],Popt,par.opt['e1'],Ms=0.762)
Mp <- msini$mj/sin(par.opt['Inc1'])
Mo <- (par.opt['Mo1']%%(2*pi))*180/pi
inc <- par.opt['Inc1']*180/pi
legend('bottomright',legend=as.expression(c(bquote(P==.(round(Popt/yr2d))~'yr'),bquote(M==.(round(Mp,2))~M[Jup]),bquote(e==.(round(par.opt['e1'],2))),bquote(I==.(round(inc,0))~'deg'),bquote(omega==.(round(par.opt['omega1']*180/pi,0))~'deg'),bquote(Omega==.(round(par.opt['Omega1']*180/pi,0))~'deg'),bquote(M[0]==.(round(Mo,0))~'deg'),bquote(jitter/error==.(round(par.opt['jitter.astro'],1))))),bty='n')
text(x=pmra,y=pmdec,labels=c('Hipparcos','Gaia'),offset=c(2,1),pos=c(4,1),col='black')
if(parallax){
    plot(data.astrometry[,1]%%2400000,data.astrometry[,'parallax'],xlab='JD-2400000',ylab=expression(pi*' [mas/yr]'),ylim=range(tmp[,1],data.astrometry[,'parallax']))
    lines(tsim0,tmp[,1],col='red')
}
dev.off()
