library(magicaxis)
options(scipen=0)
etas <- c()
alpha <- 0.6
size <- 3
margin <- 15
mgp=c(6,1,0)
par(cex=size,cex.axis=size,cex.lab=size,cex.main=size,cex.sub=size,mgp=mgp,oma=c(1,1,4,1))
for(j3 in 1:Nsig){
    n <- paste0('sig',j3)
    eta <- c()
    if(any(grepl(n,names(out[['WP']])))){
        ap.set <- names(out[['WP']][[n]])
        for(ap in ap.set){
            xx <- out[['WP']][[n]][[ap]]$xx
            t <- out[['WP']][[n]][[ap]]$t
            y <- out[['WP']][[n]][[ap]]$y
            dy <- out[['WP']][[n]][[ap]]$dy
            yy <- out[['WP']][[n]][[ap]]$yy
            zz <- out[['WP']][[n]][[ap]]$zz
            Nap <- ncol(y)
                                        #        xlim <- c(1,(2*Nap+1)*4)
            xlim <- c(1,Nap)

###MP plot
            zlim <- range(zz)
            layout(matrix(1:4, ncol = 2,byrow=TRUE), widths = 1, heights = c(2,4), respect = FALSE)
            for(j in 1:2){
                sigmaY <- sd(y)
                if(j==1){
                    par(mar = c(0, margin, 1, 0))
                    side <- 2

                }else{
                    par(mar= c(0, 0, 1, margin))
                    side <- 4
                }
                if(FALSE){
                    shifts <- -10+20/(Nap-1)*((1:Nap)-1)
                    plot(t,(y[,1]-mean(y[,1]))/sigmaY+shifts[1],ylab='RV [m/s]',xaxt='n',yaxt='n',cex=size,ylim=range(-12,12),col=colors[1])
                    axis(side=side,at=shifts,labels=paste0(Nap,'AP',1:Nap),las=1)
                    abline(h=shifts[1],col=colors[1],lty=2)
                    for(j in 2:Nap){
                        points(t,(y[,1]-mean(y[,1]))/sigmaY+shifts[j],ylab='RV [m/s]',col=colors[j])
                        abline(h=shifts[j],col=colors[j],lty=2)
                    }
                }else{
                    dy <- sapply(1:Nap,function(i) sd(y[,i]))
                    plot(1:Nap,dy,ylab=expression(sigma[RV]*' [m/s]'),xaxt='n',yaxt='n',cex=size,col=colors[1])
                    if(j==1){
                        axis(side=2)
                    }else{
                                        #                magaxis(side=side,tcl=-0.8,mgp=c(6,2,0))
                                        #                    axis(side=4,tcl=-0.5,mgp=c(4,1,0))
                    }
                }
            }                                       #        par(fig=c(0.5,1,0.7,1),new=TRUE)
#######periodograsm plot
            cols <- rainbow(100,start=alpha)#
            xlab <- 'Aperture'
###show signals; e.g. the signals in the HARPS RVs for HD41248
            for(j in 1:2){
                if(j==1){
                    par(mar = c(margin, margin, 0, 0))
                                        #        par(fig=c(0,0.5,0,0.7),new=TRUE)
                    ylim <- range(log10(yy))
                    yspan <- (ylim[2]-ylim[1])
                    xspan <- (xlim[2]-xlim[1])
                    image(xx,log10(yy),zz,xlab='',ylab='Period [day]',xaxt='n',yaxt='n',col=cols,zlim=zlim,xlim=xlim)
                                        #                options(scipen=-1)
                    magaxis(side=2,unlog=TRUE,tcl=-0.8)
                                        #                p <- axis(side=2,at=labels=log)
                                        #                options(scipen=9)
                    mtext(side=1,outer=TRUE,text=xlab,line=-3,cex=size)
                    abline(h=log10(Pkep[j3]),lty=3,lwd=2,col='darkgrey')
                                        #                dy <- 0.1*(ylim[2]-ylim[1])
                    text(rep(xlim[1]+0.1*xspan,length(Pkep[j3])),log10(Pkep[j3])+0.01*yspan,labels=format(Pkep[j3],digit=3),pos=3,cex=size,col='darkgrey',offset=0,xpd=NA)
                }else{
                                        #        par(fig=c(0.5,1,0,0.7),new=TRUE)
                    par(mar = c(margin, 0, 0, margin))
                    flow <- 0.9
                    fup <- 2
                    inds <- sort(abs(log(Pkep[j3])-log(yy)),index.return=TRUE)$ix
                    inds <- sort(inds[1:max(ceiling(length(yy)/10),10)])
                    ylim <- range(log10(yy[inds]))
                    yspan <- (ylim[2]-ylim[1])
                    xspan <- (xlim[2]-xlim[1])
                    image(xx,log10(yy[inds]),zz[,inds],xlab='',ylab='Period [day]',xaxt='n',yaxt='n',col=cols,zlim=zlim,xlim=xlim)
                    side <- 4
                                        #                p <- axis(side=2,labels=FALSE)
                                        #                p <- axis(side=side,labels=FALSE)
                                        #                axis(side=side)
                    magaxis(side=side,tcl=-0.8,mgp=c(6,2,0),unlog=TRUE)
                                        #                ind <- which((ylim[2]-p)>0.1*yspan & (p-ylim[1])>0.1*yspan)
                                        #                axis(side=side,at=p[ind])
                    abline(h=log10(Pkep[j3]),lty=3,lwd=size,col='darkgrey')
                                        #                text(rep(xlim[1]+0.08*xspan,length(Pkep[j3])),Pkep[j3]+0.01*yspan,labels=format(Pkep[j3],digit=3),pos=3,cex=size,col='darkgrey',offset=0,xpd=NA)
                    text(rep(xlim[2]-0.02*xspan,length(Pkep[j3])),log10(Pkep[j3])+0.03*yspan,labels=format(Pkep[j3],digit=3),pos=2,cex=size,col='darkgrey',offset=0,xpd=NA)
                }
                                        #            p <- axis(side=1,labels=FALSE)
                                        #            ind <- which(p>xlim[1] & p<xlim[2])
                                        #            axis(side=1,at=p[ind],mgp=mgp)
                axis(side=1,mgp=c(4,2,0))
            }
            image.plot(zz,col=cols,legend.only=TRUE,zlim = zlim,legend.mar=7.1)
####
###check whether the signal is consistent with time
            power.sd <- sapply(1:length(yy),function(j3) sd(zz[,j3]))
            ind <- which(power.sd<sd(zz[,which.min(abs(xx-Pkep[j3]))]))
            eta <- c(eta,length(ind)/length(power.sd))
            if(j3==1){
                mtext(paste("Wavelength Periodogram for signal",j3,'for',ap), outer=TRUE,  cex=1.5, line=0.5)
            }
        }
    }else{
         eta <- c(eta,NA)
    }
    etas <- cbind(etas,eta)
}
