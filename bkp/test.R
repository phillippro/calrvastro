pdf('test.pdf')
                                        #        inds <- 1:length(tsim)
i <- 'FORS2'
inds <- which(tsim>min(out$data.epoch$FORS2[,1]) & tsim<max(out$data.epoch$FORS2[,1]))
imax <- which.max(out$astrometry[,'ref_epoch'])
                                        #        imax <- which.min(out$astrometry[,'ref_epoch'])
#dt <- (tsim[inds]-out$astrometry[imax,'ref_epoch'])/365.25#yr
dt <- (tsim[inds]-out$data.epoch$FORS2[1,1])/365.25#yr
pmra <- out$astrometry[imax,'pmra']#-32
pmdec <- out$astrometry[imax,'pmdec']#-35
if(TRUE){
    dra <- plx.ra[inds]*out$astrometry[imax,'parallax']+pmra*dt
    ddec <- plx.dec[inds]*out$astrometry[imax,'parallax']+pmdec*dt
    dra0 <- out$data.epoch[[i]][,'dra']
    ddec0 <- out$data.epoch[[i]][,'ddec']
#    ddec0 <- (out$data.epoch[[i]][,'dec']-out$astrometry[imax,'dec'])*3.6e6
    plot(dra,ddec,xlab='dRA* [mas]',ylab='dDEC [mas]',type='l',xlim=rev(range(dra,dra0)),ylim=range(ddec,ddec0),col='red')
    points(dra0,ddec0)
}else{
    dt <- (tsim[inds]-out$astrometry[imax,'ref_epoch'])#day
    obs1 <- obs.lin.prop(unlist(out$astrometry[imax,]),dt)
    ra <- obs1[,'ra']+plx.ra[inds]*out$astrometry[imax,'parallax']/cos(obs1[,'dec']/180*pi)/3.6e6
    dec <- obs1[,'dec']+plx.dec[inds]*out$astrometry[imax,'parallax']/3.6e6
    ra0 <- out$data.epoch[[i]][,'ra']
    dec0 <- out$data.epoch[[i]][,'dec']
    plot(ra,dec,xlab='dRA* [mas]',ylab='dDEC [mas]',type='l',xlim=rev(range(ra,ra0)),ylim=range(dec,dec0),col='red')
    points(ra0,dec0)
}
dev.off()
