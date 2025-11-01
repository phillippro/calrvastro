plotf <- FALSE
#plotf <- TRUE
if(plotf){
cat('test.pdf\n')
pdf('test.pdf',8,8)
par(mfrow=c(2,2))
}
rvp  <-  RV.kepler(pars.kep=par.opt,bases=bases)
rv.sim  <-  RV.kepler(pars.kep=par.opt,tt=tsim,bases=bases)
nt <- grep('^t',names(par.opt))
ng <- grep('^g',names(par.opt))
nh <- grep('^h',names(par.opt))
for(n1 in names(out$timing)){
    for(n2 in names(out$timing[[n1]])){
        if(grepl('pulsation|eb',n2)){
             bjd <- out$timing[[n1]][[n2]][,1]
#             dt.yr.sim <- (tsim-out$timing[[n1]][[n2]][1,1])/365.25
#             dt.yr <- (out$timing[[n1]][[n2]][,1]-out$timing[[n1]][[n2]][1,1])/365.25
             dt.yr.sim <- (tsim-out$tref)/365.25
             dt.yr <- (out$timing[[n1]][[n2]][,1]-out$tref)/365.25
             dt <- out$timing[[n1]][[n2]][,2]
             dt.trend <- par.opt[nt]+par.opt[ng]*dt.yr+par.opt[nh]*dt.yr^2
             dt.sig <- rvp$tt-rvp$tauT
             dt.res <- dt-dt.trend-dt.sig
             et <- out$timing[[n1]][[n2]][,3]
             n3 <- gsub('pulsation|eb','',n2)
#             tsim1 <- seq(min(out$timing[[n1]][[n2]][,1]),max(out$timing[[n1]][[n2]][,1]),length.out=1e3)
             if(grepl('pulsation',n2)){
                 dt.sim <- (rv.sim$tt-rv.sim$tauT)*24*3600+par.opt[paste0('t',n3)]+par.opt[paste0('g',n3)]*dt.yr.sim+par.opt[paste0('h',n3)]*dt.yr.sim^2
                 dt.trend.sim <- par.opt[grep('^t',names(par.opt))]+par.opt[grep('^g',names(par.opt))]*dt.yr.sim+par.opt[grep('^h',names(par.opt))]*dt.yr.sim^2
                 dt.sig.sim <- dt.sim-dt.trend.sim
                 ylab <- 'Delay [second]'
             }else if(grepl('eb',n2)){
                 dt.trend.sim <- (par.opt[nt]+par.opt[ng]*dt.yr.sim+par.opt[nh]*dt.yr.sim^2)*3600*24
                 dt.sig.sim <- (rv.sim$tt-rv.sim$tauT)*3600*24
                 dt.sim <- dt.trend.sim+dt.sig.sim
                 dt <- dt*3600*24
                 et <- et*3600*24
                 dt.sig <- dt.sig*3600*24
                 dt.res <- dt.res*3600*24
                 dt.trend <- dt.trend*3600*24
                 ylab <- 'Delay [second]'
             }
             if(TRUE){
                 bin <- wtb(bjd,dt.sig+dt.res,et,dt=1)
                 plot(bin[,1],bin[,2],xlab='BJD',ylab=ylab,main=paste('planet',gsub('p','',n1),';',n2),pch=20,col='white',ylim=range(bin[,2]+bin[,3],bin[,2]-bin[,3]))
                 arrows(bin[,1],bin[,2]+bin[,3],bin[,1],bin[,2]-bin[,3],length=0.05,angle=90,code=3,col=tcol('black',50))
                 points(bin[,1],bin[,2],pch=20,col=tcol('black',50))
                 lines(tsim,dt.sig.sim,col='red')
             }else{
                 plot(bjd,dt.sig+dt.res,xlab='BJD',ylab=ylab,main=paste('planet',gsub('p','',n1),';',n2),pch=20,col='white',ylim=range(dt.sig+et,dt.sig-et))
                 arrows(bjd,dt.sig+et,bjd,dt.sig-et,length=0.05,angle=90,code=3,col=tcol('black',50))
                 points(bjd,dt.sig,pch=20,col=tcol('black',50))
                 lines(tsim,dt.sig.sim,col='red')
             }
             ###phase-folded
             d2s <- 3600*24
             if(grep('pulsation|eb',timetype)){
                 dt.hat <- (rvp$tt-rvp$tauT)*d2s+par.opt[nt]+par.opt[ng]*dt.yr+par.opt[nh]*dt.yr^2#second
             }else{
                 dt.hat <- (rvp$tt-rvp$tauT)+par.opt[nt]+par.opt[ng]*dt.yr+par.opt[nh]*dt.yr^2
             }
             for(jj in 1:Nsig){
                 pp <- par.opt
                 if(length(Nsig)>1 & any(grepl('^K\\d',names(pp)))) pp[,paste('K',(1:Nsig)[-jj])] <- 1e-6
                 if(length(Nsig)>1 & any(grepl('^arc\\d',names(pp)))) pp[,paste('arc',(1:Nsig)[-jj])] <- 1e-6
                 rvpj  <-  RV.kepler(pars.kep=pp,bases=bases)
                 rv.simj  <-  RV.kepler(pars.kep=pp,tt=tsim,bases=bases)
                 if(grep('pulsation|eb',timetype)){
                     dtj.sig <- (out$timing[[n1]][[n2]][,'dt']-dt.hat)+(rvpj$tt-rvpj$tauT)*d2s
                     dtj.sim <- (rv.simj$tt-rv.simj$tauT)*d2s
                 }else{
                     dtj.sig <- ((out$timing[[n1]][[n2]][,'dt']-dt.hat)+rvpj$tt-rvpj$tauT)*3600*24
                     dtj.sim <- (rv.simj$tt-rv.simj$tauT)*3600*24#s
                 }
                 tt <- bjd%%Popt[jj]
                 ind <- sort(tt,index.return=TRUE)$ix
                 bin <- wtb(tt[ind],dtj.sig[ind],et[ind],dt=1)
                 plot(bin[,1],bin[,2],xlab='BJD',ylab='O-C [s]')
                 arrows(bin[,1],bin[,2]+bin[,3],bin[,1],bin[,2]-bin[,3],length=0.05,angle=90,code=3,col='grey')
                 points(bin[,1],bin[,2],pch=1)
                 inds <- sort(tsim%%Popt[jj],index.return=TRUE)$ix
                 lines((tsim%%Popt[jj])[inds],dtj.sim[inds],col='red')
             }
        }else{
#            index <- out$ind.all$timing[[paste0(n1,'_',n2)]]
#            lte  <-  (rvp$tt[index]-rvp$tauT[index])*24*60#min
            lte.sim <- (rv.sim$tt-rv.sim$tauT)*24*60#min
            ip <- as.integer(gsub('p','',n1))
            if(any(dP==ip)){
                p <- Prefs[ip]+par.opt[paste0('dP',ip)]/60/24#day
            }else{
                p <- exp(par.opt[paste0('per',ip)])#day
            }
            dt <- out$timing[[n1]][[n2]][,1]-out$ttv0
            if(grepl('occult',n2)) dt <- dt+p/2
            ttv <- (dt-round(dt/p)*p)*24*60
            ettv <- out$timing[[n1]][[n2]][,2]*24*60
            bjd <- out$timing[[n1]][[n2]][,1]
            oc <- ttv-par.opt['t0']
            plot(bjd,oc,xlab='BJD',ylab='TTV [min]',main=paste('transiting planet',gsub('p','',n1),';',n2),pch=20,col='white',ylim=range(oc+ettv,oc-ettv))
            arrows(bjd,oc+ettv,bjd,oc-ettv,length=0.05,angle=90,code=3,col=tcol('black',50))
            points(bjd,oc,pch=20,col=tcol('black',50))
            lines(tsim,lte.sim,col='red')
        }
    }
}

if(plotf){
dev.off()
}
