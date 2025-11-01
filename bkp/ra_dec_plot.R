plot(trel,draPs,xlab='BJD',ylab=expression(Delta*alpha*'* [mas]'))
cs <- c()
for(i in inss){
    tP <- out$rel[[n]][[i]][,1]
    draP <- out$rel[[n]][[i]][,'dra']
    ddecP <- out$rel[[n]][[i]][,'ddec']
    edraP <- out$rel[[n]][[i]][,'edra']
    eddecP <- out$rel[[n]][[i]][,'eddec']
    points(tP,draP,col=cc[nc],pch=20,cex=0.8)
    points(0,0,pch='+',col='green',cex=2)
    points(tsim[ii],draS[ii],pch='+',col='red')
    pos <- rep(2,length(ii))
    pos[which(draS[ii]>mean(draS[ii]))] <- 4
    text(tsim[ii],draS[ii],labels=tts,pos=pos,xpd=NA,col='red',cex=0.8)
    inds <- c()
    for(k in 1:length(draP)){
###                inds <- c(inds,which.min(abs(out$rel[[n]]$tot[k,1]-tsim)))
        inds <- c(inds,which.min(abs(tP[k]-tsim)))
    }
#    arrows(draS[inds],ddecS[inds],draP,ddecP,length=0.001,angle=0,code=1,col='grey')
#    arrows(draP-edraP,ddecP,draP+edraP,ddecP,length=0.01,angle=90,code=3,col=cc[nc])
#    arrows(draP,ddecP-eddecP,draP,ddecP+eddecP,length=0.01,angle=90,code=3,col=cc[nc])
    nc <- nc+1
    cs <- c(cs,cc[nc])
}
lab.ins <- inss
                                        #            if(length(inss[inss!='tot'])>1){
if(length(inss)>1){
    legend('topright',legend=lab.ins,col=cs,bty='n',pch=20)
}
###add model prediction
lines(tsim,draS,col='red')

plot(trel,ddecPs,xlab='BJD',ylab=expression(Delta*alpha*'* [mas]'))
cs <- c()
for(i in inss){
    tP <- out$rel[[n]][[i]][,1]
    draP <- out$rel[[n]][[i]][,'dra']
    ddecP <- out$rel[[n]][[i]][,'ddec']
    edraP <- out$rel[[n]][[i]][,'edra']
    eddecP <- out$rel[[n]][[i]][,'eddec']
    points(tP,ddecP,col=cc[nc],pch=20,cex=0.8)
    points(0,0,pch='+',col='green',cex=2)
    points(tsim[ii],ddecS[ii],pch='+',col='red')
    pos <- rep(2,length(ii))
    pos[which(ddecS[ii]>mean(ddecS[ii]))] <- 4
    text(tsim[ii],ddecS[ii],labels=tts,pos=pos,xpd=NA,col='red',cex=0.8)
    inds <- c()
    for(k in 1:length(ddecP)){
###                inds <- c(inds,which.min(abs(out$rel[[n]]$tot[k,1]-tsim)))
        inds <- c(inds,which.min(abs(tP[k]-tsim)))
    }
    nc <- nc+1
    cs <- c(cs,cc[nc])
}
lab.ins <- inss
                                        #            if(length(inss[inss!='tot'])>1){
if(length(inss)>1){
    legend('topright',legend=lab.ins,col=cs,bty='n',pch=20)
}
###add model prediction
lines(tsim,ddecS,col='red')
