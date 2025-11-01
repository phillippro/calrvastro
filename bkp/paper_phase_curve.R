source('mcmc_func.R')
library(magicaxis)
targets <- c('HD210193','HD211369','HD211970','HD218511','HD39855','HIP35173','HD187085','HD208487','HD85683','HIP12961','HIP74995')
ff <- 'results/phase_all.Robj'
if(file.exists(ff)){
    load(ff)
}else{
    out.show <- list()
    out.sim <- list()
}
fout <- 'phase_all.pdf'
cat('output pdf:\n',fout,'\n')
pdf(fout,16,16)
size <- 1.2
gap <-1.5
par(mfrow=c(4,4),mar=c(gap,gap,gap,gap),oma=rep(4,4),cex.lab=size,cex.axis=size,cex=size)
for(kk in 1:length(targets)){
    cat('target:',targets[kk],'\n')
#    plot(e$tt2,e$RV.show,xaxt='n',yaxt='n')
    if(!file.exists(ff)){
        fin <-paste0('results/',targets[kk],'/',targets[kk],'_PFS_white_MA_AR_GP011_Esd0.2_ofac2_xi10_N1e+06.Robj')
        load(fin,envir=e <- new.env())
        trv <-e$trv
        RV <- e$RV
        eRV <- e$eRV
        tsim <- e$tsim
        RV.sim <- e$RV.sim
        out.sim[[kk]] <-cbind(tsim,RV.sim)
        out.show[[kk]] <-cbind(trv,RV,eRV)
    }else{
        trv <- out.show[[kk]][,1]
        RV <- out.show[[kk]][,2]
        eRV <- out.show[[kk]][,3]
        tsim <-out.sim[[kk]][,1]
        RV.sim <-out.sim[[kk]][,2]
    }
    plot(trv-min(trv),RV,xaxt='n',yaxt='n',pch=20,cex=0.5,main=targets[kk])
    magaxis(side=1,cex=1.0)
    magaxis(side=3,label=FALSE)
    magaxis(side=c(2,4))
    try(arrows(trv-min(trv),RV+eRV,trv-min(trv),RV-eRV,length=0.05,angle=90,code=3,col='grey'),TRUE)
    points(tsim-min(trv),RV.sim,xaxt='n',yaxt='n',col='red',pch='.')
}
mtext(side=1,text='Orbital Phase [day]',outer=TRUE,cex=1.8,line=1.8)
mtext(side=c(2,4),text='K [m/s]',outer=TRUE,cex=2,line=2)
dev.off()
save(out.show,out.sim,file='phase_all.Robj')
