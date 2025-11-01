library(RColorBrewer)
getLevel <- function(kk,prob=0.95){
    dx <- diff(kk$x[1:2])
    dy <- diff(kk$y[1:2])
    sz <- sort(kk$z)
    c1 <- cumsum(sz) * dx * dy
    approx(c1, sz, xout = 1 - prob)$y
}
sigmas <- c(0.682,0.954, 0.997)

mcmc.opt <- out[['mcmc.opt']][[paste0('sig',Nsig)]]
Npar <- ncol(mcmc.opt)-2
size <- 0.9
par(mfcol=c(Npar,Npar),mar=c(0,0,0,0),oma=c(4,4,4,1),cex.axis=0.8,cex.lab=size,cex=size,mgp=c(4,2,0),las = 2)
inds <- sample(nrow(mcmc.opt),min(1e4,nrow(mcmc.opt)/2))
labs <- nams <- plot.labels.simple(colnames(mcmc.opt[,1:Npar]))
#indP <- grep('per',names(par.opt))
#par.opt[indP] <- exp(par.opt[indP])
for(j in 1:Npar){
#    cat('j=',j,'\n')
    for(k in 1:Npar){
        xlab <- labs[[j]]
        ylab <- labs[[k]]
        x <- mcmc.opt[inds,j]
        y <- mcmc.opt[inds,k]
        if(k>j){
            kk <- try(kde2d(x,y, n=20),TRUE)
            if(class(kk)=='try-error') kk <- kde2d(x,y, n=20,h=c(0.1*(max(x)-min(x)),0.1*(max(y)-min(y))))
            my.cols <- 'black'
            levels <- getLevel(kk,prob=sigmas)
            levels <- levels[!is.na(levels)]
                                        #            levels <- pretty(zlim, n=1)#68,95, and 99% contours
            if(length(levels)>0){
                if(length(inds)<=1e3 & FALSE){
                    plot(x,y,pch=20,cex=0.5,xaxt='n',yaxt='n')
                    contour(kk, drawlabels=FALSE, nlevels=length(levels), levels=levels,col=my.cols,add=TRUE,axes=FALSE)
                }else{
                    contour(kk, drawlabels=FALSE, nlevels=length(levels), levels=levels,col=my.cols,xaxt='n',yaxt='n')
                }
            }
#            legend('topleft',inset=c(-0.2,-0.1),legend=paste0(format(cor(x,y),digit=2)),bty='n',text.col='red',cex=size)
        }else if(k==j){
            hist(x,xlab=labs[[j]],main='',axes=FALSE,breaks=20,col = "gray",border=FALSE)
            tmp <- data.distr(x,plotf=FALSE)
            abline(v=c(tmp[c('xminus.1sig','xplus.1sig')]),lwd=2,lty=1)
#            abline(v=c(tmp['median']),lwd=2,lty=2)
            abline(v=par.opt[j],lwd=2,lty=2)
#            ind <- which.max(post.out)
#            abline(v=mcmc.opt[ind,j],lwd=2,lty=2)
        }else{
            plot(1, type="n", axes=F, xlab="", ylab="",main='')
        }
        if(j==1 & k>=j){
#            axis(side=2,mgp=c(3,1,0),cex=0.5)
            mtext(side=2,text=ylab,line=1,cex=1.3*size,srt = 45)
        }
        if(k==Npar& k>=j){
#            axis(side=1,mgp=c(3,1,0),cex=0.8)
            mtext(side=1,text=xlab,line=1,cex=1.3*size,srt = 45)#, pos = 1, xpd = TRUE)
#            text(seq(1, 10, by=1), par("usr")[3] - 0.2, labels = lablist,
        }
    }
}
mtext("Pair Plot", outer=TRUE,  cex=1.5, line=1)
par(las = 0)
