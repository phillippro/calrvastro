library(magicaxis)
#options(scipen=20)
options(scipen=0)
if(!exists('MPonly')) MPonly <- FALSE
###Moving periodogram setting up
Nbin0 <- 20
alpha <- 0.6
beta <- 0.05
if(target=='HD42581') alpha <- 0.65
if(target=='HIP80268') alpha <- 0.4
if(target=='GL433') beta <- 0.1
if(target=='HIP21556') beta <- 0.02
if(target=='HIP93206') alpha <- c(0,0.5)
if(target=='HIP117886'){
    alpha <- 0.58
    beta <- 0.02
}
if(target=='HIP22762') alpha <- 0.55
gammas <- c()
size <- 3
margin <- 13
mgp <- c(5,1,0)
np <- c('b','c','d','e','f','g')
for(j3 in 1:Nsig){
    if(MPonly){
        pname <- np[match(Pkep[j3],sort(Pkep))]
        if(target=='GL433'){
           if(j3==2) pname <- 'd'
           if(j3==3) pname <- 'c'
        }
#        fmp <- paste0('results/',target,'/',target,pname,'_MP.pdf')
        fmp <- paste0(dir.in,target,pname,'_MP.pdf')
        cat(fmp,'\n')
        pdf(fmp,16,16)
    }
    par(cex=size,cex.axis=size,cex.lab=size,cex.main=size,cex.sub=size,mgp=mgp,oma=c(1,1,4,1))
    n <- paste0('sig',j3)
    if(any(grepl(n,names(out[['MP']])))){
        xx <- out[['MP']][[n]]$xx
        y <- out[['MP']][[n]]$y
        yy <- out[['MP']][[n]]$yy
        zz <- out[['MP']][[n]]$zz
###MP plot
        zlim <- range(zz)
        xlim <- range(trv.all-out$tref)
        dxlim <- 0.05*(max(xlim)-min(xlim))
        xlim <- c(min(xlim)-dxlim,max(xlim)+dxlim)
        layout(matrix(1:4, ncol = 2,byrow=TRUE), widths = 1, heights = c(2,4), respect = FALSE)
        for(j in 1:2){
            ylab <- 'RV [m/s]'
            if(j==1){
                                        #        par(fig=c(0,0.5,0.7,1),new=TRUE)
                par(mar = c(0, margin, 1, 0))
                plot(trv.all-out$tref,y,ylab=ylab,xaxt='n',yaxt='n',xaxs='i',cex=size,col='white',mgp=c(5,2,0),pch=20,xlim=xlim)
                p <- axis(side=2,labels=FALSE)
                ylim <- range(y)
                dy <- 0.1*diff(ylim)
                ind <- which(p>ylim[1]+dy)
                axis(side=2,at=p[ind])
            }else{
                                        #        par(fig=c(0.5,1,0.7,1),new=TRUE)
                par(mar= c(0, 0, 1, margin))
                plot(trv.all-out$tref,y,ylab='',xaxt='n',yaxt='n',xaxs='i',cex=size,col='white',pch=20,xlim=xlim)
            }
            cc <- c()
            for(kk in 1:length(out$ins.rv)){
                inds <- out[[out$ins.rv[kk]]]$index
                ii <- out$ins.rv[kk]
                if(out$ins.rv[kk]=='HARPS') ii <- 'HARPSpre'
                if(out$ins.rv[kk]=='PFS') ii <- 'PFSpre'
                points(trv.all[inds]-out$tref,y[inds],col=colors[ii],pch=20,cex=2)
                cc <- c(cc,colors[ii])
            }
            t <- trv.all-out$tref
            x <- y
            ex <- eRV.all
            arrows(t,x-ex,t,x+ex,length=0.05,angle=90,code=3,col='grey')
            magaxis(side=3,tcl=-0.5,mgp=c(3,2,0),majorn=3)
            magaxis(side=1,tcl=0.5,mgp=c(3,2,0),labels=FALSE,majorn=3)
        }
        legend('topright',legend=out$ins.rv,xpd=NA,inset=c(-0.42,0),col=cc,pch=20,bty='n',cex=0.8*size)
#######periodograsm plot
        cols <- rainbow(100,start=alpha[min(length(alpha),j3)])#
        xlab <- paste0('BJD-',out$tref)

###show signals; e.g. the signals in the HARPS RVs for HD41248
        for(j in 1:2){
            if(j==1){
                par(mar = c(margin, margin, 0, 0))
                                        #        par(fig=c(0,0.5,0,0.7),new=TRUE)
                ind <- which(diff(yy)<1e-10)
                if(length(ind)>0){
                    yy <- yy[-ind]
                    zz <- zz[,-ind]
                }
                ylim <- range(log10(yy))
                yspan <- (ylim[2]-ylim[1])
                xspan <- (xlim[2]-xlim[1])
                ylab <- 'Period [day]'
                image(xx,log10(yy),zz,xlab='',ylab=ylab,xaxt='n',yaxt='n',col=cols,zlim=zlim,xlim=xlim)
                box()
#                image(xx,log10(yy),t(zz),xlab='',ylab='Period [day]',col=cols,zlim=zlim,xlim=xlim)
                                        #                options(scipen=-1)
                magaxis(side=2,unlog=TRUE,tcl=-0.8)
#                p <- axis(side=2,at=labels=log)
#                options(scipen=9)
                mtext(side=1,outer=TRUE,text=xlab,line=-3,cex=0.8*size)
#                mtext(side=2,outer=TRUE,text=ylab,line=-6,cex=size)
                abline(h=log10(Pkep[j3]),lty=3,lwd=3,col='grey50')
#                dy <- 0.1*(ylim[2]-ylim[1])
                if(min(xx)>500 & length(xx)>5){
                    text(rep(xlim[1]+0.1*xspan,length(Pkep[j3])),log10(Pkep[j3])+0.01*yspan,labels=format(Pkep[j3],digit=3),pos=3,cex=size,col='grey50',offset=0,xpd=NA)
                }else{
                    if(target=='HIP80268' & j3==1){
                        text(rep(xlim[1]+0.08*xspan,length(Pkep[j3])),log10(Pkep[j3])+0.01*yspan,labels=format(Pkep[j3],digit=3),pos=3,cex=size,col='grey50',offset=0,xpd=NA)
                    }else if(target=='HIP93206' & j3==2){
                        text(rep(xlim[1]+0.1*xspan,length(Pkep[j3])),log10(Pkep[j3])+0.01*yspan,labels=format(Pkep[j3],digit=3),pos=3,cex=size,col='grey50',offset=0,xpd=NA)
                    }else{
                        text(rep(xlim[1]-0.05*xspan,length(Pkep[j3])),log10(Pkep[j3])+0.01*yspan,labels=format(Pkep[j3],digit=3),pos=3,cex=size,col='grey50',offset=0,xpd=NA)
                    }
                }
            }else{
                par(mar = c(margin, 0, 0, margin))
                inds <- sort(abs(log(Pkep[j3])-log(yy)),index.return=TRUE)$ix
                inds <- sort(inds[1:ceiling(length(yy)*beta)])
                if(length(inds)<5) inds <- 1:length(yy)
                ylim <- range(log10(yy[inds]))
                yspan <- (ylim[2]-ylim[1])
                xspan <- (xlim[2]-xlim[1])
                image(xx,log10(yy[inds]),zz[,inds],xlab='',ylab='Period [day]',xaxt='n',yaxt='n',col=cols,zlim=zlim,xlim=xlim)
#                magaxis(side=4,tcl=-0.8,mgp=c(6,2,0),unlog=TRUE,logpretty=FALSE,prettybase=1.5,majorn=10)
                magaxis(side=4,tcl=-0.8,mgp=c(6,2,0),unlog=TRUE)
                abline(h=log10(Pkep[j3]),lty=3,lwd=3,col='grey50')
                if(min(xx)>500 & length(xx)>5){
                    text(rep(xlim[2]-0.02*xspan,length(Pkep[j3])),log10(Pkep[j3])+0.03*yspan,labels=format(Pkep[j3],digit=3),pos=2,cex=size,col='grey50',offset=0,xpd=NA)
                }else if(target=='HIP93206' & j3==2){
                    text(rep(xlim[2]+0.1*xspan,length(Pkep[j3])),log10(Pkep[j3])+0.01*yspan,labels=format(Pkep[j3],digit=3),pos=2,cex=size,col='grey50',offset=0,xpd=NA)
                }else{
                    text(rep(xlim[2]+0.1*xspan,length(Pkep[j3])),log10(Pkep[j3])+0.03*yspan,labels=format(Pkep[j3],digit=3),pos=2,cex=size,col='grey50',offset=0,xpd=NA)
                }
            }
#            p <- axis(side=1,labels=FALSE)
#            ind <- which(p>xlim[1] & p<xlim[2])
#            axis(side=1,at=p[ind])
#            axis(side=1,mgp=c(4,2,0))
            magaxis(side=1,tcl=-0.5,mgp=c(4,3,0),majorn=3)
        }
        image.plot(t(zz),col=cols,legend.only=TRUE,zlim = zlim,legend.mar=7.1)
####
###check whether the signal is consistent with time
        power.sd <- sapply(1:length(yy),function(j3) sd(zz[,j3]))
        ind <- which(power.sd<sd(zz[,which.min(abs(xx-Pkep[j3]))]))
        gamma <- length(ind)/length(power.sd)
        if(j3==1 & !MPonly){
            mtext(text=paste("Moving Periodogram for signal",j3), outer=TRUE,  cex=1.5, line=0)
        }
        if(MPonly) dev.off()
    }else{
#        par(mfrow=c(2,2))
        gamma <- NA
    }
    gammas <- c(gammas,gamma)
}
