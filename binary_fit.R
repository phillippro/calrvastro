#plotf <- TRUE
plotf <- FALSE
if(plotf){
pdf('test.pdf')
par(mfrow=c(2,2))
}else{
par(mfrow=c(4,4))
}
eys <- ys <- ts <- ress <- c()
out.sim$indC <- out.sim$indP
simulate.binary <- RV.kepler(pars.kep=par.sig,tt=out.sim$BJDtdb,out1=out.sim,kep.only=TRUE,bases=bases)
out.sim$indC <- c()
rv.sim.all <- simulate.binary$rv
rv.sim <- simulate.binary$rvkep
ysim.kep <- simulate.binary$rvkep
ysim <- simulate.binary$rv

tmp <- c()
reslist <- tlist <- ylist <- eylist <- list()
for(i in out$ins.binary){
#    inds <- out$ind.binary[[i]]
#    t <- out$data.binary[inds,1]
    inds <- out$ind.all$binary[[i]]
    t <- out$BJDtdb[inds]
    ts <- c(ts,t)
    tlist[[i]] <- t
    if(offset){
        b <- par.opt[[paste0('b_',i)]]
    }else{
        b <- 0
    }
#    y <- out$data.binary[inds,2]-b
    inds <- out$ind.all$binary[[i]]
    res <- rvma[inds]
    ress <- c(ress,res)
    reslist[[i]] <- res
    y <- rvfit$rvkep[inds]+res
    ys <- c(ys,y)
    ylist[[i]] <- y
    ey <- out$data.binary[out$data.binary$ins==i,3]
    eys <- c(eys,ey)
    eylist[[i]] <- ey
    tmp <- rbind(tmp,data.frame(t,y,ey,res,inds,i))
    plot(t%%Popt[1],y,xlab=paste0('Orbital Phase [days]'),ylab='RV (data-model) [m/s]',ylim=range(y,ysim.kep),main=paste0(binary,i))
    tsim.fold <- tsim%%Popt[1]
    ii <- sort(tsim.fold,index.return=TRUE)$ix
    lines(tsim.fold[ii],ysim.kep[ii],col='red')
}
colnames(tmp) <- c('t','rvkep.res','erv','res','ind','ins')
###
plot(ts,ys,xlab='JD',ylab='RV [m/s]',col='white',main='combined fit of the binary orbit')
for(j in 1:length(out$ins.binary)){
    i  <-  out$ins.binary[j]
    points(tlist[[i]],ylist[[i]],col=tcol(colors[j],tcol.percentage))
    arrows(tlist[[i]],ylist[[i]]+eylist[[i]],tlist[[i]],ylist[[i]]-eylist[[i]],length=0.05,angle=90,code=3,col=tcol(colors[j],tcol.percentage))
}
lines(out.sim$BJDtdb,ysim.kep,col='red')
##O-C
plot(ts,ress,xlab='JD',ylab='RV [data-model] [m/s]',col='white',main='residuals')
for(j in 1:length(out$ins.binary)){
    i  <-  out$ins.binary[j]
    points(tlist[[i]],reslist[[i]],col=tcol(colors[j],tcol.percentage))
    arrows(tlist[[i]],reslist[[i]]+eylist[[i]],tlist[[i]],reslist[[i]]-eylist[[i]],length=0.05,angle=90,code=3,col=tcol(colors[j],tcol.percentage))
}

####combined fit to raw data
ind.binary <- unlist(out$ind.binary)
plot(out$BJDtdb[ind.binary],out$data.binary[,2],col='white',xlab='Time [BJD]',ylab='RV (data-model)',ylim=range(out$data.binary[,2],ysim))
for(j in 1:length(out$ins.binary)){
    i  <-  out$ins.binary[j]
    inds <- out$ind.binary[[i]]
    t <- out$BJDtdb[inds]
    if(offset){
        b <- par.opt[[paste0('b_',i)]]
    }else{
        b <- 0
    }
#    y <- out$data.binary[inds,2]-b
#    ey <- out$data.binary[inds,3]
    y <- rvfit$rv[inds]+rvma[inds]
    ey <- out$data.binary[out$data.binary$ins==i,3]
    data.bin <- bin.simple(cbind(t,y,ey),10)
#    points(t,y,col=tcol(colors[j],tcol.percentage))
    y <- ((1+rvmn[out$ind.all$binary[[i]]]/CMPS)*(1+par.opt[paste0('rv0b_',i)]/CKMPS)-1)*CMPS
    points(t,y,col=tcol(colors[j],tcol.percentage))#
}
try(arrows(data.bin[,1],data.bin[,2]+data.bin[,3],data.bin[,1],data.bin[,2]-data.bin[,3],length=0.05,angle=90,code=3,col=colors[i]),TRUE)
lines(tsim,ysim,col='red')
n1 <- ceiling(length(out$ins.binary)/2)
n2 <- length(out$ins.binary)-n1
legend('top',legend=out$ins.binary[1:n1],col=colors[1:n1],inset=c(0,-0.2),seg.len = 1,horiz=TRUE,xpd=NA,bty='n',pch=1)
if(n2>0){
    legend('top',legend=out$ins.binary[(n1+1):length(out$ins.binary)],col=colors[(n1+1):length(out$ins.binary)],seg.len = 1,inset=c(0,-0.1),horiz=TRUE,xpd=NA,bty='n',pch=1)
}
mtext("Binary fit", outer=TRUE,  cex=1.5, line=-0.5)
if(plotf) dev.off()
