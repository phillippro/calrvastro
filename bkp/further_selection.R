source('mcmc_func.R')
library(magicaxis)
library(plotrix)

splitpsi <- function(str,ind=1,N=5){
    psis <- as.numeric(unlist(strsplit(str,'d')))
    n <- length(psis)
    out <- c()
    if(n<ind){
        out <- rep(NA,N-ind)
    }else if(n<(N-1)){
        if(ind>1){
            out <- c(psis[-(1:(ind-1))],rep(NA,N-n-1))
        }else{
            out <- c(psis,rep(NA,N-n-1))
        }
    }else{
        if(ind>1){
            out <- psis[-(1:(ind-1))]
        }else{
            out <- psis
        }
    }
    out
}

N <- 5
#tab <- read.table('results/misaligned_starpar_N621.txt',header=TRUE)
tab <- read.csv('results/misaligned_starpar_N627.txt')
nsys <- nrow(tab)
cn <- colnames(tab)
cc <- unlist(lapply(1:4,function(i) lapply((i+1):5,function(j) paste0('psi',i,j))))
psi.map <- psi.lower <- psi.upper <- psi.med <- c()
for(j in 1:nrow(tab)){
    ii.med <- grep('psi.med',cn)
    ii.lower <- grep('psi.lower',cn)
    ii.upper <- grep('psi.upper',cn)
    ii.map <- grep('psi.map',cn)
    jj <- which(!is.na(tab[j,ii.med]))
    tmp1 <- tmp2 <- tmp3 <- tmp4 <- rep(NA,length(cc))
    if(length(jj)==N){
        index <- 1:(N-1)
    }else{
        index <- jj
    }
    for(k in index){
        if(k>1) inds <- sum((N-1):(N-k+1))+(1:(5-k))
        if(k==1) inds <- 1:(N-1)
        tmp1[inds] <- splitpsi(tab[j,ii.med[k]],k,N)
        tmp2[inds] <- splitpsi(tab[j,ii.upper[k]],k,N)
        tmp3[inds] <- splitpsi(tab[j,ii.lower[k]],k,N)
        tmp4[inds] <- splitpsi(tab[j,ii.map[k]],k,N)
    }
    psi.med <- rbind(psi.med,tmp1)
    psi.upper <- rbind(psi.upper,tmp2)
    psi.lower <- rbind(psi.lower,tmp3)
    psi.map <- rbind(psi.map,tmp4)
}
colnames(psi.map) <- colnames(psi.lower) <- colnames(psi.upper) <- colnames(psi.med) <- cc

inds <- list()
thresholds <- c(5,10,20,30,40,50)
#thresholds <- c(5)
#thresholds <- c(30)
#thresholds <- c(50)
tt <- read.table('targets_quality.txt',header=TRUE)
stars <- tab[,1]
for(t in thresholds){
    tmp <- c()
    for(k in 1:nsys){
        jj <- which(tt[,1]==tab[k,1])
        if(any(psi.upper[k,]-psi.lower[k,]<t & tt[jj,2]=='A',na.rm=TRUE)){
            tmp <- c(tmp,k)
        }
    }
    inds[[paste0('dpsi',t)]] <- tmp
    fpdf <- paste0('results/misaligned_ma_dpsi',t,'.pdf')
    ffile <- gsub('pdf','txt',fpdf)
    fid <- gsub('ma','targets',ffile)
    cat(fpdf,'\n')
    cat(ffile,'\n')
    cat(fid,'\n\n')
    pdf(fpdf,8,8)
    Nbreaks <- 10
    size <- 1.2
    par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5,1,0),cex=size,cex.axis=size,cex.lab=size)
    p <- hist(as.numeric(psi.med[tmp,]),xlab=expression(psi*" [deg]"),ylab='Frequency',main='',breaks=Nbreaks)
    i1 <- tmp[which(tab[tmp,'Mstar']>0.9 & tab[tmp,'Mstar']>1.1)]
    i2 <- tmp[which(tab[tmp,'Mstar']<0.9)]
    i3 <- tmp[which(tab[tmp,'Mstar']>1.1)]
    indm <- grep('mp.med',cn)
    inde <- grep('^e.med',cn)
    inda <- grep('^a.med',cn)
    indI <- grep('Inc.med',cn)
    indO <- grep('Omega.med',cn)

    aa <- unlist(lapply(tmp,function(i) max(tab[i,inda],na.rm=TRUE)))
    ac <- 8
    ii1 <- tmp[which(aa<ac)]
    ii2 <- tmp[which(aa>ac)]
#    hist(psi.med[i1,],col='red',add=TRUE,fill=FALSE)
#    hist(psi.med[i2,],col='blue',add=TRUE,fill=FALSE)
#    hist(psi.med[i3,],col='green',add=TRUE,fill=FALSE)
    p1 <- hist(as.numeric(psi.med[ii1,]),plot=FALSE,breaks = p$breaks)
    p2 <- hist(as.numeric(psi.med[ii2,]),plot=FALSE,breaks = p$breaks)
    lines(p1$mids-1,p1$counts,col='red',type='h',lwd=2)
    lines(p2$mids+1,p2$counts,col='blue',type='h',lwd=2)
    legend('top',inset=c(0,-0.1),xpd=NA,bty='n',legend=as.expression(c(bquote(a[max]*'<'*.(ac)~'au'),bquote(a[max]*'>'*.(ac)~'au'))),lty=c(1,1),col=c('red','blue'),horiz=TRUE)

    psi.rad <- psi.med[tmp,1]/180*pi
    plot(unlist(tab[tmp,inda[1]]),cos(psi.rad),xlab='a [au]',ylab='cos(psi)',main='',log='x')
#    hist(cos(psi.med[tmp,]/180*pi),xlab='cos(psi)',ylab='Freq.',main=paste('Threshold: dPsi <',t,'deg'))
    p <- hist(cos(psi.med[tmp,]/180*pi),xlab=expression(cos*psi),ylab='Frequency',main='',breaks=Nbreaks)
    p1 <- hist(as.numeric(cos(psi.med[ii1,]/180*pi)),plot=FALSE,breaks = p$breaks)
    p2 <- hist(as.numeric(cos(psi.med[ii2,]/180*pi)),plot=FALSE,breaks = p$breaks)
    lines(p1$mids-0.01,p1$counts,col='red',type='h',lwd=2)
    lines(p2$mids+0.01,p2$counts,col='blue',type='h',lwd=2)
    legend('top',inset=c(0,-0.1),xpd=NA,bty='n',legend=as.expression(c(bquote(a[max]*'<'*.(ac)~'au'),bquote(a[max]*'>'*.(ac)~'au'))),lty=c(1,1),col=c('red','blue'),horiz=TRUE)

    plot(unlist(tab[tmp,indm[1]]),cos(psi.rad),xlab=expression(m[c]*' ['*M[Jup]*']'),ylab='cos(psi)',main='',log='x')
    plot(unlist(tab[tmp,inde[1]]),cos(psi.rad),xlab='e',ylab='cos(psi)',main='',log='x')
    plot(unlist(tab[tmp,inda]),unlist(tab[tmp,indm]),xlab='a [au]',ylab=expression(m[c]*' ['*M[Jup]*']'),main='',log='xy',pch=20,col='white',cex=0.1,axes=F)
    magaxis(1:2)
    magaxis(3:4,label=FALSE)
    cc <- rainbow(length(tmp))
    for(j0 in 1:length(tmp)){
        j <- tmp[j0]
        ii <- which(!is.na(tab[j,inda]))
#        cat(ii,'\n')
        for(i in ii){
            jj <- ii[ii!=i]
            for(j1 in jj) segments(tab[j,inda[i]],tab[j,indm[i]],x1=tab[j,inda[j1]],y1=tab[j,indm[j1]],col=tcol(cc[j0],50))
            points(tab[j,inda[i]],tab[j,indm[i]],col=cc[j0],pch=20,cex=0.8)
        }
    }

####polar plot
    for(j0 in 1:length(tmp)){
        j <- tmp[j0]
        Iopt <- as.numeric(tab[j,indI]%%(2*pi))*180/pi#deg
        Oopt <- as.numeric(tab[j,indO]%%(2*pi))*180/pi#deg
        r2d <- 180/pi
        if(j0==1)  polar.plot(Iopt,Oopt,main='Polar plot',lwd=3,line.col=4,rp.type='s',point.symbols=20,point.col='white',radial.lim=c(0,max(tab[tmp,indI]*180/pi,na.rm=TRUE)))
        ii <- which(!is.na(tab[j,indI]))
                                        #        cat(ii,'\n')
        for(i in ii){
            jj <- ii[ii!=i]
            i1 <- tab[j,indI[i]]*r2d
            o1 <- tab[j,indO[i]]%%(2*pi)*r2d
            x0 <- i1*cos(o1/r2d)
            y0 <- i1*sin(o1/r2d)
            for(j1 in jj){
                i2 <- tab[j,indI[j1]]*r2d
                o2 <- tab[j,indO[j1]]%%(2*pi)*r2d
                x1 <- i2*cos(o2/r2d)
                y1 <- i2*sin(o2/r2d)
                segments(x0,y0,x1,y1,col=tcol(cc[j0],50))
            }
            points(x0,y0,col=cc[j0],pch=20,cex=0.8)
        }
    }

###polar plot
if(FALSE){
    fpdf <- paste0('results/polar.pdf')
    cat(fpdf,'\n')
    for(j0 in 1:length(tmp)){
        j <- tmp[j0]
        ii <- which(!is.na(Iopt)&!is.na(Oopt))
        polar.plot(Iopt[ii],Oopt[ii],main=stars[k],lwd=3,line.col=4,rp.type='s',point.symbols=20,point.col='steelblue',radial.lim=c(0,max(Iopt[ii])))
    }
    dev.off()
}

    dev.off()
    write.table(gsub('Robj','pdf',tab[tmp,'file']),file=ffile,quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(tab[tmp,1],file=fid,quote=FALSE,row.names=FALSE,col.names=FALSE)
}
