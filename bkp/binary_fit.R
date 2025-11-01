par(mfrow=c(4,4))
eys <- ys <- ts <- c()
ysim <- simulate$rvC
tauC <- rv$tauC
resC <- out$data.binary[,2]-rv$rvC
ervC <- out$data.binary[,3]
rvCg <- resC+rv$rvgC
rvCs <- resC+rv$rvsC
rvCgs <- resC+rv$rvsC+rv$rvgC

for(i in out$ins.binary){
    inds <- out$ind.binary[[i]]
    t <- out$data.binary[inds,1]
    y <- out$data.binary[inds,2]-par.opt[[paste0('b_',i)]]
    ey <- out$data.binary[inds,3]
    plot(t%%Popt[1],y,xlab=paste0('Orbital Phase [days]'),ylab='RV (data-model) [m/s]',ylim=range(y,ysim),main=paste0(binary,i))
    ts <- tsim%%Popt[1]
    ii <- sort(ts,index.return=TRUE)$ix
    lines(ts[ii],ysim[ii],col='red')
}

plot(out$data.binary[,1],out$data.binary[,2],col='white',xlab='Time [BJD]',ylab='RV (data-model)',ylim=range(out$data.binary[,2],ysim))
for(j in 1:length(out$ins.binary)){
    i  <-  out$ins.binary[j]
    inds <- out$ind.binary[[i]]
    t <- out$data.binary[inds,1]
    y <- out$data.binary[inds,2]-par.opt[[paste0('b_',i)]]
    ey <- out$data.binary[inds,3]
    data.bin <- bin.simple(cbind(t,y,ey),10)
    points(t,y,col=tcol(colors[j],tcol.percentage))
    try(arrows(data.bin[,1],data.bin[,2]+data.bin[,3],data.bin[,1],data.bin[,2]-data.bin[,3],length=0.05,angle=90,code=3,col=colors[i]),TRUE)
}
lines(tsim,ysim,col='red')
n1 <- ceiling(length(out$ins.binary)/2)
n2 <- length(out$ins.binary)-n1
legend('top',legend=out$ins.binary[1:n1],col=colors[1:n1],inset=c(0,-0.2),seg.len = 1,horiz=TRUE,xpd=NA,bty='n',pch=1)
if(n2>0){
    legend('top',legend=out$ins.binary[(n1+1):length(out$ins.binary)],col=colors[(n1+1):length(out$ins.binary)],seg.len = 1,inset=c(0,-0.1),horiz=TRUE,xpd=NA,bty='n',pch=1)
}
mtext("Binary fit", outer=TRUE,  cex=1.5, line=-0.5)

###relativistic effect
if(relativity){
    for(k in 1:3){
        if(k==1){
            main <- 'GR'
            y <- rvCg
            yc <- simulate$rvgC
        }else if(k==2){
            main <- 'SR'
            y <- rvCs
            yc <- simulate$rvsC
        }else{
            main <- 'SR+GR'
            y <- rvCgs
            yc <- simulate$rvsC+simulate$rvgC
        }

        data.bin <- bin.simple(cbind(tauC,y,ervC),10)
        for(i in 1:2){
            if(i==1){
                ylim <- range(data.bin[,2],yc)
            }else{
                ylim <- range(yc)
            }
            plot(tauC,y,xlab=expression(tau*' [BJD]'),ylab='RV [m/s]',main=main,col='white',ylim=ylim)
            for(j in 1:length(out$ins.binary)){
                i <- out$ins.binary[j]
                ind <- out$ind.binary[[i]]
                t <- tauC[ind]
                erv <- ervC[ind]
                points(t,y[ind],col=tcol(colors[j],tcol.percentage))
                data.bin <- bin.simple(cbind(t,y[ind],erv),10)
                points(data.bin[,1],data.bin[,2],col=colors[j],pch=20)
                try(arrows(data.bin[,1],data.bin[,2]+data.bin[,3],data.bin[,1],data.bin[,2]-data.bin[,3],length=0.05,angle=90,code=3,col=colors[j],pch=20),TRUE)
            }
            lines(tsim,yc,col='red')
        }
    }
}

####difference between emission time and BJD
dt <- (tauC-out$data.binary[,1])*24*60
dt.sim <- (simulate$tauC-tsim)*24*60
plot(tauC,dt,xlab='BJD',ylab='tau-BJD [min]',ylim=range(dt,dt.sim))
lines(tsim,dt.sim,col='red')

####RV variation due to biased astrometry and Keplerian motion of the target star
plot(rv$tauC,rv$rvC,xlab='BJD',ylab='dRV [m/s]',ylim=range(rv$rvC,simulate$rvC),main='Bias due to astrometry+reflex motion')
lines(tsim,simulate$rvC,col='red')
