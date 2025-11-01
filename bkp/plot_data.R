for(i in out$ins.rv){
    plot(out[[i]]$RV[,1],out[[i]]$RV[,2],xlab='BJD-2400000',ylab='RV [m/s]',main=paste(target,'raw data:',i))
    try(arrows(out[[i]]$RV[,1],out[[i]]$RV[,2]-out[[i]]$RV[,3],out[[i]]$RV[,1],out[[i]]$RV[,2,]+out[[i]]$RV[,3],length=0.05,angle=90,code=3,col='grey'),TRUE)
}
if(length(out$data.binary)>0){
    for(i in out$ins.binary){
        inds <- out$ind.binary[[i]]
        t <- out$data.binary[inds,1]
        rv <- out$data.binary[inds,2]
        erv <- out$data.binary[inds,3]
        plot(t,rv,xlab='BJD-2400000',ylab='RV [m/s]',main=paste(binary,'raw data:',i))
        try(arrows(t,rv-erv,t,rv+erv,length=0.05,angle=90,code=3,col='grey'),TRUE)
    }
}
mtext("Raw RV Data", outer=TRUE,  cex=1.5, line=-0.5)
###combined data
plot(trv.all,RV.all,xlab='BJD-2400000',ylab='RV [m/s]',main='',col='white')
for(i in 1:length(out$ins.rv)){
    points(out[[out$ins.rv[i]]]$RV[,1],out[[out$ins.rv[i]]]$RV[,2],col=colors[i])
    try(arrows(out[[out$ins.rv[i]]]$RV[,1],out[[out$ins.rv[i]]]$RV[,2]-out[[out$ins.rv[i]]]$RV[,3],out[[out$ins.rv[i]]]$RV[,1],out[[out$ins.rv[i]]]$RV[,2,]+out[[out$ins.rv[i]]]$RV[,3],length=0.05,angle=90,code=3,col=colors[i]),TRUE)
}
#try(arrows(trv,RV-eRV,trv,RV+eRV,length=0.05,angle=90,code=3,col='grey'),TRUE)
n1 <- ceiling(length(out$ins.rv)/2)
n2 <- length(out$ins.rv)-n1
legend('top',legend=out$ins.rv[1:n1],col=colors[1:n1],inset=c(0,-0.2),seg.len = 1,horiz=TRUE,xpd=NA,bty='n',pch=1)
if(n2>0){
    legend('top',legend=out$ins.rv[(n1+1):length(out$ins.rv)],col=colors[(n1+1):length(out$ins.rv)],seg.len = 1,inset=c(0,-0.1),horiz=TRUE,xpd=NA,bty='n',pch=1)
}

###combined binary data
if(length(out$data.binary)>0){
    plot(out$data.binary[,1],out$data.binary[,2],xlab='BJD-2400000',ylab='RV [m/s]',main='',col='white')
    for(i in 1:length(out$ins.binary)){
        inds <- out$ind.binary[[out$ins.binary[i]]]
        dd <- out$data.binary[inds,]
        points(dd[,1],dd[,2],col=colors[i])
        try(arrows(dd[,1],dd[,2]-dd[,3],dd[,1],dd[,2]+dd[,3],length=0.05,angle=90,code=3,col=colors[i]),TRUE)
    }
                                        #try(arrows(trv,RV-eRV,trv,RV+eRV,length=0.05,angle=90,code=3,col='grey'),TRUE)
    n1 <- ceiling(length(out$ins.binary)/2)
    n2 <- length(out$ins.binary)-n1
    legend('top',legend=out$ins.binary[1:n1],col=colors[1:n1],inset=c(0,-0.2),seg.len = 1,horiz=TRUE,xpd=NA,bty='n',pch=1)
    if(n2>0){
        legend('top',legend=out$ins.binary[(n1+1):length(out$ins.binary)],col=colors[(n1+1):length(out$ins.binary)],seg.len = 1,inset=c(0,-0.1),horiz=TRUE,xpd=NA,bty='n',pch=1)
    }
}

###plot activity indices
par(mfrow=c(4,4))
panel <- 0
for(i in out$ins.rv){
    act <- out[[i]]$activity
    act.name <- names(act)
    act.name <- act.name[act.name!='window']
    if(length(act.name)>0){
        for(an in act.name){
            x <- act[[an]][,1]
            y <- act[[an]][,2]
            ind <- which(!is.na(x) & !is.na(y))
            if(length(ind)>0){
                plot(x[ind],y[ind],xlab='BJD-2400000',ylab='RV [m/s]',main=paste0(i,';',an))
                panel <- panel+1
                if(panel==1) mtext("Noise Proxy", outer=TRUE,  cex=1.5, line=-0.5)
            }
#            try(arrows(act[[an]][,1],act[[an]][,2]-act[[an]][,3],act[[an]][,1],act[[an]][,2]+act[[an]][,3],length=0.05,angle=90,code=3,col='grey'),TRUE)
        }
    }
}
if(any(grepl('photometry',names(out)))){
    photo <- out$photometry
    for(j in 1:length(photo)){
        plot(photo[[j]][,1],photo[[j]][,2],xlab='Time',ylab='Flux or Magnitude',main=paste0('Photometry:',names(photo)[j]))
###check transit
        if(!all(is.na(Popt))){
            x <- photo[[j]][,2]
            ex <- photo[[j]][,3]
            for(k in 1:Nsig){
                t <- (photo[[j]][,1]-min(photo[[j]][,1]))
                tmp <- wtb(t,x,ex,dt=0.01)
                t <- t%%Popt[k]
                tmp[,1] <- tmp[,1]%%Popt[k]
                plot(t,x,xlab='Phase [day]',ylab='dFlux',pch=1,cex=0.5,main=paste('signal',k,':',round(Popt[k],2),'day'),col=tcol('grey',50),ylim=c(mean(x)-3*sd(x),mean(x)+3*sd(x)))
                points(tmp[,1],tmp[,2],pch=20,cex=0.5,col=tcol('red',20))
            }
        }
    }
}

