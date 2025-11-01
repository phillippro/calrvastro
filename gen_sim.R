###simulation
sim <- TRUE
ts <- c(tmin,tmax,out$tiall[,1])
if(length(out$astrometry)>0){
    ts <- c(ts,out$astrometry[,'ref_epoch'])
}
if(length(out$data.ref)>0){
    ts <- c(ts,out$data.ref[,1])
}
if(length(out$data.epoch)>0){
    for(i in out$ins.epoch){
        ts <- c(ts,out$data.epoch[[i]][,1])
    }
}
if(any(names(out)=='rel')){
    tt <- c()
    for(n1 in names(out$rel)){
        for(n2 in names(out$rel[[n1]])){
            tt <- c(tt,out$rel[[n1]][[n2]][[1]])
        }
    }
    ts <- c(ts,tt)
}
if((max(ts)-min(ts))<max(Popt)){
    dT <- max(Popt)-(max(ts)-min(ts))
    ts <- c(ts,min(ts)-dT)
}
if(any(names(out)=='gost')) ts <- c(ts,out$gost[,'BJD'])
tsim <- seq(min(ts),max(ts),length.out=Nsim)
#tsim <- trv.all
#out0 <- out.sim <- out
#out.sim$BJDtcb <- tsim
out0 <- out
#if(barycor){
    if(!is.null(out$tiall)){
        source('barycor.R')
    }else{
        source('barycor_hipgaia.R')
    }
    out.sim <- out
out.sim$indP <- 1:length(tsim)
out.sim$indC <- c()
out <- out0
#}
rm(out0)
#stop()
