####################
###plot
####################
pdf.file <- paste0('phase_fold_Np',Np,'_',id,'_',noise.model,p,q,'.pdf')
#pdf(pdf.file,15,8)
#pdf(pdf.file,15,8)
cat('Np=',Np,'\n')
Ncol <- Nrow <- floor(sqrt(Np))
cat('Nrow=',Nrow,'\n')
if(Np>Nrow*Nrow & Np<=Nrow*(Nrow+1)){
    Ncol <- Nrow+1
}else if(Np>Nrow*(Nrow+1)){
    Ncol <- Nrow+1
    Nrow <- Nrow+1
}
cat('Nrow=',Nrow,'\n')
cat('Ncol=',Ncol,'\n')
Nrow <- 1
Ncol <- 3
pdf(pdf.file,Ncol*4,Nrow*4)
cat('pdf.file=',pdf.file,'\n')
size <- max(Ncol/2,1.2)
#par(mfrow=c(2,3),mar=c(2.5,2.5,0,0),oma=c(5,5,1,12),cex.lab=2,cex.axis=2,mgp=c(3.5,1,0))
#par(mfrow=c(2,3),mar=c(3,3,0,0),oma=c(5,5,1,19),cex.lab=2.5,cex.axis=2.5,mgp=c(4,1,0))
#par(mfrow=c(Nrow,Ncol),mar=c(3,3,0,0),oma=c(5,5,1,19),cex.lab=2.5,cex.axis=2.5,mgp=c(4,1,0))
par(mfrow=c(Nrow,Ncol),mar=c(4.5,4.5,1,1),cex.axis=size,cex.lab=size)
#####20, 49, 164, 600 signals
data.combine <- FALSE
Nf <- 0
#ind.order <- c(3,1,4,2)
ind.order <- sort(Popt,index.return=TRUE)$ix
#ind.order <- c(4,1,2,3)
#Ncol <- 3
err <- FALSE
#ylab <- 'RV(data-model)[m/s]'
ylab <- 'RV [m/s]'
for(i3 in ind.order){
    t3c <- datac[[k1]][,i3]
    RV3c <- datac[[k1]][,Np+i3]
    t4c <- datab[[k1]][,i3]
    RV4c <- datab[[k1]][,Np+i3]
    eRV4c <- datab[[k1]][,2*Np+i3]
    t5c <- datas[[k1]][,i3]
    RV5c <- datas[[k1]][,Np+i3]
    cols <- c('black','blue','green','orange','cyan')
    pchc <- 1:Ncol
    Nf <- Nf+1
    if(data.combine){
        ylim <- range(RV4c[!is.na(eRV4c)]-eRV4c[!is.na(eRV4c)],RV4c[!is.na(eRV4c)]+eRV4c[!is.na(eRV4c)],RV5c)
        for(k2 in 1:Nw){
            t3c <- datac[[k2]][,i3]
            RV3c <- datac[[k2]][,Np+i3]
            t4c <- datab[[k2]][,i3]
            RV4c <- datab[[k2]][,Np+i3]
            eRV4c <- datab[[k2]][,2*Np+i3]
            t5c <- datas[[k2]][,i3]
            RV5c <- datas[[k2]][,Np+i3]
        if(err){
            cex <- 0.2
            pch <- pchc[k2]
        }else{
            cex <- 1
            pch <- 20
        }
            if(k2==1){
                rspan <- max(RV)-min(RV)
                plot(t3c,RV3c,xlab=paste0('Orbital phase [day]'),ylab=ylab,pch=pch,cex=cex,col=cols[k2],ylim=ylim)#col='black',
            }else{
                points(t3c,RV3c,cex=cex,col=cols[k2])
            }
            if(err){
                points(t4c,RV4c,col=cols[k2],cex=size,pch=20)
                arrows(t4c,RV4c-eRV4c,t4c,RV4c+eRV4c,length=0.1,angle=90,code=3,col=cols[k2],cex=size,pch=20)
            }else{
                arrows(t3c,RV3c-eRV,t3c,RV3c+eRV,length=0.03,angle=90,code=3,col='grey',cex=size,pch=20)
            points(t3c,RV3c,pch=20,cex=cex)
            }
            points(t5c,RV5c,col="red",pch=20,cex=0.2)
                                        #        legend('topright',legend=paste0(i3,'planet'),bty='n')
            if(any(i3==c(1))){
                legend('topright',legend=paste(,format(Popt[i3],digit=3),'d'),bty='n',cex=size)
            }else{
                legend('topright',legend=paste(format(Popt[i3],digit=3),'d'),bty='n',cex=size)
            }
        }
        if(Nf==1){
#            ylab <- 'RV(data-model) [m/s]'
            mtext(ylab,side=2,outer=TRUE,line=1,cex=size)
        }
        if(Nf==round(Ncol+(Ncol+1)/2)){
            mtext(paste0('Orbital phase [day]'),side=1,outer=TRUE,line=1,cex=size)
        }
        if(Nf==Ncol){
            par(xpd=TRUE)
            if(Nw>1){
                legend('topright',inset=c(-0.46,0),legend=c('1AP1','1AP1 binned','Keck','Keck binned','model'),pch=c(20,20,20,20,NA),lty=c(NA,NA,NA,NA,1),pt.cex=c(0.5,1.4,0.5,1.4,NA),lwd=c(rep(NA,4),1),col=c('black','black','blue','blue','red'),xpd=NA,cex=size)
	    }else{
                legend('topright',inset=c(-0.68,0),legend=c('RVs','binned RVs','model'),pch=c(20,20,NA),lty=c(NA,NA,1),pt.cex=c(0.8,1.5,NA),lwd=c(rep(NA,2),1),col=c('black','black','red'),xpd=NA,cex=size)
            }
            par(xpd=FALSE)
        }
    }else{
###k1 is the number of data set
        ylim <- range(RV4c[!is.na(RV4c)]-eRV4c[!is.na(eRV4c)],RV4c[!is.na(RV4c)]+eRV4c[!is.na(eRV4c)],RV5c,RV3c,RV3c+eRV,RV3c-eRV)
        if(err){
            cex <- 0.2
            pch <- pchc[k2]
        }else{
            cex <- 1
            pch <- 20
        }
        plot(t3c,RV3c,xlab=paste0('Orbital phase [day]'),ylab=ylab,pch=pch,ylim=ylim,main='',cex=cex,col=cols[k1])#col='black',
if(err){
        points(t4c,RV4c,cex=size,pch=20,col=cols[k1])
        arrows(t4c,RV4c-eRV4c,t4c,RV4c+eRV4c,length=0.03,angle=90,code=3,col=cols[k1],cex=size,pch=20)
    }else{
            arrows(t3c,RV3c-eRV,t3c,RV3c+eRV,length=0.03,angle=90,code=3,col='grey',cex=size)
            points(t3c,RV3c,pch=20,cex=cex)
}
        points(t5c,RV5c,col="red",pch=20,cex=0.2)
#        legend('topright',legend=paste(format(Popt[i3],digit=3),'d'),bty='n',cex=size)
#        legend('topright',legend=paste0(i3,'planet'),bty='n',cex=size)
#        legend('topright',legend=paste0(i3,'planet'),bty='n',cex=size)
    }
}
dev.off()
