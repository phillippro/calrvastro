##library("viridis")
#plotf <- TRUE
plotf <- FALSE
if(plotf){
pdf('test.pdf')
par(mfrow=c(2,2))
}else{
par(mfrow=c(4,4))
}
yr2d <- 365.25
ng <- 'G'
nast <- length(out$astro.index)
nastro <- nrow(out$astrometry)
poss <- (1:4)[1:nastro]
#ell.col <- c('orange','darkgreen')
#ell.col <- rainbow(nastro)
ell.col <- terrain.colors(nastro)
if(nastro>2) ng <- out$cats
cc0 <- cc <- c('black','blue','orange','steelblue','green','darkgreen',rainbow(10))
#cc <- c('black','blue','cyan','steelblue','orange','darkgreen')
##astrometry fit
imax <- which.max(Popt)
###add binary motion
nc <- 1
if(any(names(par.opt)=='Mstar')){
    mstar <- par.opt['Mstar']
}else{
    mstar <- out$Mstar
}

tmp <- astrometry.kepler(par.opt,bases=bases)
astro0 <- tmp$barycenter
if(any(names(out)=='rel')){
    source('timing_function.R')
    source('general_function.R')
    source('sofa_function.R')

####simulation
    trel <- out$trel
#    Dt <- Popt[j]-(max(trel)-min(trel))
#    Dt <- max(0,(max(trel)-min(trel))-Popt)
    Dt <- Popt[imax]
#    tsim <- seq(max(24e5,min(trel)-Dt/2),min(25e5,max(trel)+Dt/2),length.out=1e3)

    yrs <- time_Jd2yr(cbind(tsim,0))
#    tts <- seq(max(floor(min(yrs)),1990),min(ceiling(max(yrs)),2030),length.out=10)
    tts <- seq(max(floor(min(yrs)),1990),min(ceiling(max(yrs)),2030),length.out=10)
    tts <- unique(round(tts))
#    tts <- tts[tts>1990 & tts<2030]
#    tts <- tts[tts>2010 & tts<2030]
#    tts <- c(tts,2022.1)
    ii <- unlist(lapply(tts, function(tt) which.min(abs(yrs-tt))))
    jjs <- trel[ii]
    popt <- exp(par.opt['per1'])
    planet <- astrometry.rel(par.opt,tt=out.sim$BJDtdb,out1=out.sim)
    range.raS <- range(unlist(planet$rel$dra))
    range.decS <- range(unlist(planet$rel$ddec))
    for(n in names(out$rel)){
        for(i in names(out$rel[[n]])){
            if(out$typerel[[n]][[i]]=='sp'){
                rho <- out$rel[[n]][[i]][,'sep']
                erho <- out$rel[[n]][[i]][,'esep']
                pa <- out$rel[[n]][[i]][,'pa']
                epa <- out$rel[[n]][[i]][,'epa']
                rd <- ps2rd(rho,erho,pa,epa)
                out$rel[[n]][[i]][,'dra'] <- rd[,'dra']
                out$rel[[n]][[i]][,'era'] <- rd[,'edra']
                out$rel[[n]][[i]][,'ddec'] <- rd[,'ddec']
                out$rel[[n]][[i]][,'edec'] <- rd[,'eddec']
            }
        }
    }
    range.raP <- range(unlist(lapply(names(out$rel),function(n) unlist(lapply(names(out$rel[[n]]),function(i) out$rel[[n]][[i]][,'dra'])))))
    range.decP <- range(unlist(lapply(names(out$rel),function(n) unlist(lapply(names(out$rel[[n]]),function(i) out$rel[[n]][[i]][,'ddec'])))))
####
#pdf('paper_astrometry.pdf',6,6)
#par(mar=c(5,5,1,1))
    for(j in 1:Nsig){
        n <- paste0('p',j)
        if(out$Nrel[j]>0 & length(out$rel[[n]])>0){
            n <- paste0('p',j)
            inss <- names(out$rel[[n]])
            inss <- inss[inss!='tot']

            tPs <- draPs <- ddecPs <- c()

            for(i in inss){
                draPs <- c(draPs,out$rel[[n]][[i]][,'dra'])
                ddecPs <- c(ddecPs,out$rel[[n]][[i]][,'ddec'])
                tPs <- c(tPs,out$rel[[n]][[i]][,1])
            }
            draS <- planet$rel$dra[[n]]
            ddecS <-planet$rel$ddec[[n]]
#            xlim <- 1.2*range(draPs,draS)
#            ylim <- 1.2*range(ddecPs,ddecS)
            xlim <- range(range.raP,range.raS)
            ylim <- range(range.decP,range.decS)
#            xlim <- ylim <- range(xlim,ylim)
            plot(draPs,ddecPs,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),xlim=rev(xlim),ylim=ylim,col='white')
            cs <- c()
            for(i in inss){
                nc <- nc+1
                tP <- out$rel[[n]][[i]][,1]
                draP <- out$rel[[n]][[i]][,'dra']
                ddecP <- out$rel[[n]][[i]][,'ddec']
                edraP <- out$rel[[n]][[i]][,'era']
                eddecP <- out$rel[[n]][[i]][,'edec']
                points(draP,ddecP,col=cc0[nc],pch=20,cex=0.8)
                points(0,0,pch='+',col='green',cex=2)
                points(draS[ii],ddecS[ii],pch='+',col='red')
                if(target=='HD42581'){
                    t0 <- sum(time_Yr2jd(2022))
                    t1 <- sum(time_Yr2jd(2024))
                    jj <- which(tsim>t0&tsim<t1)
                    lines(draS[jj],ddecS[jj],col='purple',lwd=4)
                    rho <- sqrt(draS[jj]^2+ddecS[jj]^2)/1000#as
                    cat('rho=',rho,'as\n')
                }
                pos <- rep(2,length(ii))
                pos[which(draS[ii]>mean(draS[ii]))] <- 4
                text(draS[ii],ddecS[ii],labels=tts,pos=pos,xpd=NA,col='red',cex=0.8)
                inds <- c()
                for(k in 1:length(draP)){
###                inds <- c(inds,which.min(abs(out$rel[[n]]$tot[k,1]-tsim)))
                    inds <- c(inds,which.min(abs(tP[k]-tsim)))
                }
#                arrows(draS[inds],ddecS[inds],draP,ddecP,length=0.001,angle=0,code=1,col='grey')
#                arrows(draP-edraP,ddecP,draP+edraP,ddecP,length=0.01,angle=90,code=3,col=cc0[nc])
#		cat(i,';edraP=',edraP,';eddecP=',eddecP,'\n')
                arrows(draP,ddecP-eddecP,draP,ddecP+eddecP,length=0.01,angle=90,code=3,col=cc0[nc])
                cs <- c(cs,cc0[nc])
            }
            lab.ins <- inss
#            if(length(inss[inss!='tot'])>1){
            if(length(inss)>1){
                legend('topright',legend=lab.ins,col=cs,bty='n',pch=20)
            }
###add model prediction
            lines(draS,ddecS,col='red')

###add additional companion orbit for reference
            if(length(out$rel)>1){
                kk <- (1:length(out$rel))[-j]
                for(k in kk){
                    draS <- planet$rel$dra[[k]]
                    ddecS <- planet$rel$ddec[[k]]
                    lines(draS,ddecS,col='grey')
                    points(draS[ii],ddecS[ii],pch='+',col=cc[k])
                    pos <- rep(2,length(ii))
                    pos[which(draS[ii]>mean(draS[ii]))] <- 4
                    text(draS[ii],ddecS[ii],labels=tts,pos=pos,xpd=NA,col=cc[k],cex=0.8)
                }
                if(target=='HD42581'){
                    t0 <- sum(time_Yr2jd(2022))
                    t1 <- sum(time_Yr2jd(2024))
                    jj <- which(tsim>t0&tsim<t1)
                    lines(draS[jj],ddecS[jj],col='purple',lwd=4)
                    rho <- sqrt(draS[jj]^2+ddecS[jj]^2)/1000#as
                    cat('rho=',rho,'as\n')
                }
            }

###add legend
            popt <- exp(par.opt[paste0('per',j)])
            K <- par.opt[paste0('K',j)]
            e <- par.opt[paste0('e',j)]
            inc <- par.opt[paste0('Inc',j)]*180/pi
            Mo <- (par.opt[paste0('Mo',j)]%%(2*pi))*180/pi
            omega <- par.opt[paste0('omega',j)]*180/pi
            Omega <- par.opt[paste0('Omega',j)]*180/pi
            Mp <- k2m(K,popt,e,Ms=mstar,par.opt[paste0('Inc',j)])$mj
            Mo <- (par.opt[paste0('Mo',imax)]%%(2*pi))*180/pi
            legend('top',legend=as.expression(c(bquote(P==.(round(popt/yr2d))~'yr'),bquote(M==.(round(Mp,1))~M[Jup]),bquote(e==.(round(e,2))),bquote(I==.(round(inc,1))~'deg'))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.12,x.intersp=0.2)
            legend('top',legend=as.expression(c(bquote(omega==.(round(omega,1))~'deg'),bquote(Omega==.(round(Omega,1))~'deg'),bquote(M[0]==.(round(Mo,1))~'deg'),bquote(logJ[rel]==.(round(par.opt['logJ.rel'],1))))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.08,x.intersp=0.2)
        }
    }
#dev.off()
}

####add fit to epoch data
if(length(out$data.epoch)>0){
    astro <- astrometry.kepler(pars.kep=par.opt,bases=bases)
    astro.sim <- astrometry.kepler(par.opt,tt=out.sim$BJDtdb,out1=out.sim)
    epoch <- astro$epoch
    epoch.sim <- astro.sim$epoch
    for(j3 in 1:length(out$ins.epoch)){
        i  <- out$ins.epoch[j3]
        reflex.sim <- astrometry.epoch(par.opt,tt=out.sim$BJDtcb,out1=out.sim,comp=out.sim$comp.epoch[j3])$epoch
        pred.res <- loglikelihood(par.opt,prediction=TRUE)
        par.opt1 <- par.opt
        par.opt1[grepl('^K\\d',names(par.opt1))] <- 1e-6
###The following approach may not be good to show the residual relatied to long period companion because the stellar reflex motion is subtraced from the Gaia epoch and then propagate the Gaia epoch to Hipparcos epoch which will result in huge offset if without reflex motion
                                        #        pred.sig <- loglikelihood(par.opt1,prediction=TRUE)
#        tepoch <- out$data.epoch[[i]][,1]
        tepoch <- out$BJDtdb[out$ind.all$epoch1[[i]]]
        i2 <- which(tsim>min(tepoch) & tsim<max(tepoch))
        if(i=='hip1'){
            absci <- out$data.epoch[[i]][,'dv']
            eabsci <- out$data.epoch[[i]][,'ev']
            cpsi <- out$data.epoch[[i]][,'pv_pra']
            spsi <- out$data.epoch[[i]][,'pv_pdec']
            dra0 <- out$data.epoch[[i]][,'dv']*out$data.epoch[[i]][,'pv_pra']
            era <- era0 <- out$data.epoch[[i]][,'ev']*out$data.epoch[[i]][,'pv_pra']
            ddec0 <- out$data.epoch[[i]][,'dv']*out$data.epoch[[i]][,'pv_pdec']
            edec <- edec0 <- out$data.epoch[[i]][,'ev']*out$data.epoch[[i]][,'pv_pdec']
            dra.res <- pred.res$res[[paste0('epoch_',i)]]*cpsi
            dra.dec <- pred.res$res[[paste0('epoch_',i)]]*spsi
            dra.obs <- dra <- dra.res+epoch[[i]][,'dra']
            ddec.obs  <- ddec <- dra.dec+epoch[[i]][,'ddec']
        }else if(i=='hip2'){
            cpsi <- out$data.epoch[[i]][,'CPSI']
            spsi <- out$data.epoch[[i]][,'SPSI']
            absci <- out$data.epoch[[i]][,'RES']
            eabsci <- out$data.epoch[[i]][,'SRES']
            dra0 <- absci*cpsi
            era <- era0 <- eabsci*cpsi
            ddec0 <- absci*spsi
            edec <- edec0 <- eabsci*spsi
            dra.res <- pred.res$res[[paste0('epoch_',i)]]*cpsi
            ddec.res <- pred.res$res[[paste0('epoch_',i)]]*spsi
            dra.obs <- dra <- dra.res+epoch[[i]][,'dra']
            ddec.obs  <- ddec <- ddec.res+epoch[[i]][,'ddec']
        }else{
#            nepoch <- length(out$ind.all$epoch[[i]])
            dra.res <- pred.res$res[[paste0('epoch_',i)]][,1]
            ddec.res <- pred.res$res[[paste0('epoch_',i)]][,2]
            dra.obs <- out$data.epoch[[i]][,'dra']
            ddec.obs <- out$data.epoch[[i]][,'ddec']
            dra <- dra.res+epoch[[i]][,'dra']
            ddec <- ddec.res+epoch[[i]][,'ddec']
            era <- out$data.epoch[[i]][,'era']
            edec <- out$data.epoch[[i]][,'edec']
        }
        dra.sim <- epoch.sim[[i]][,'dra']
        ddec.sim <- epoch.sim[[i]][,'ddec']
        dra.pred <- epoch[[i]][,'dra']
        ddec.pred <- epoch[[i]][,'ddec']

###full fit
        plot(dra.obs,ddec.obs,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),xlim=rev(range(dra,dra.sim[i2])),ylim=range(ddec,ddec.sim[i2]),main=paste(i,'full'),pch=20)
        points(0,0,pch='+')
        ncolor <- 20
        cc <- paletteer_c(palette = "viridis::inferno", n = ncolor)
        jj <- as.factor( as.numeric( cut((tepoch-tepoch[1])%%Popt[imax], ncolor)))
        if(grepl('hip',i)){
            segments(dra.obs-era,ddec.obs-edec,dra.obs+era,ddec.obs+edec,col=tcol('grey',50))
        }else{
            arrows(dra.obs-era,ddec.obs,dra.obs+era,ddec.obs,length=0.001,angle=0,code=1,col='grey')
            arrows(dra.obs,ddec.obs-edec,dra.obs,ddec.obs+edec,length=0.001,angle=0,code=1,col='grey')
        }
        ii <- sapply(tepoch,function(t) which.min(abs(tsim-t)))
        segments(dra.sim[ii],ddec.sim[ii],dra.obs,ddec.obs,col=tcol(cc[jj],50))
        lines(dra.sim,ddec.sim,col='red')
        points(dra.sim[ii],ddec.sim[ii],col=cc[jj],pch=20,cex=0.5)
###projected orbit
#        tobs <- out$data.epoch[[i]][,1]
        tobs <- out$BJDtdb[out$ind.all$epoch1[[i]]]
        plot(tobs,dra.obs,xlab='BJD',ylab=expression(Delta*alpha*'* [mas]'))
        arrows(tobs,dra.obs-era,tobs,dra.obs+era,length=0.001,angle=0,code=1,col='grey')
        lines(tsim,dra.sim,col='red')

        plot(tobs,ddec.obs,xlab='BJD',ylab=expression(Delta*delta))
        arrows(tobs,ddec.obs-edec,tobs,ddec.obs+edec,length=0.001,angle=0,code=1,col='grey')
        lines(tsim,ddec.sim,col='red')

        Mp <- k2m(par.opt[paste0('K',imax)],Popt[imax],par.opt[paste0('e',imax)],Ms=mstar,Inc=par.opt[paste0('Inc',imax)])$mj
        Mo <- (par.opt[paste0('Mo',imax)]%%(2*pi))*180/pi
        inc <- par.opt[paste0('Inc',imax)]*180/pi
        legend('top',legend=as.expression(c(bquote(P==.(round(Popt[imax]/yr2d))~'yr'),bquote(M==.(round(Mp,1))~M[Jup]),bquote(e==.(round(par.opt['e1'],2))),bquote(I==.(round(inc,1))~'deg'),bquote(omega==.(round(par.opt['omega1']*180/pi,1))~'deg'))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.12,x.intersp=0.15,cex=0.9)
        legend('top',legend=as.expression(c(bquote(omega==.(round(par.opt['omega1']*180/pi,1))~'deg'),bquote(Omega==.(round(par.opt['Omega1']*180/pi,1))~'deg'),bquote(M[0]==.(round(Mo,1))~'deg'),bquote(logJ[hip]==.(round(par.opt[grep('J_hip',names(par.opt))],2))),bquote(logJ[gaia]==.(round(par.opt['logJ_gaia'],2))))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.08,x.intersp=0.15,cex=0.9)

####individual plot
        if(grepl('hip',i)){
            dy <- 0
            plot(tepoch,absci+dy,xlab='Time [BJD]',ylab='Abscissa [mas]',ylim=range(absci+dy,pred.res$res[[paste0('epoch_',i)]]+dy))
            arrows(tepoch,absci-eabsci+dy,tepoch,absci+eabsci+dy,length=0.001,angle=0,code=1,col='grey')
                                        #        res.absci <- absci-(epoch[,'dra']+par.opt['dra'])*cpsi-(epoch[,'ddec']+par.opt['ddec'])*spsi
            points(tepoch,pred.res$res[[paste0('epoch_',i)]]+dy,col='red')
            plot(tepoch,dra,xlab='Time [BJD]',ylab=expression(Delta*alpha*'* [mas]'))
            arrows(tepoch,dra-era,tepoch,dra+era,length=0.001,angle=0,code=1,col='grey')
            lines(tsim,epoch.sim[[i]][,'dra'],col='red')
            points(tepoch,epoch.sim[[i]][ii,'dra'],col=cc[jj],pch=20,cex=0.8)

            plot(tepoch,ddec,xlab='Time [BJD]',ylab=expression(Delta*delta*' [mas]'))
            arrows(tepoch,ddec-edec,tepoch,ddec+edec,length=0.001,angle=0,code=1,col='grey')
            lines(tsim,epoch.sim[[i]][,'ddec'],col='red')
            points(tepoch,epoch.sim[[i]][ii,'ddec'],col=cc[jj],pch=20,cex=0.8)
        }else{
            reflex <- astrometry.epoch(par.opt,tt=tepoch,comp=out$comp.epoch[j3])$epoch
            dra.lin.sim <- dra.sim-reflex.sim[,'dra']
            ddec.lin.sim <- ddec.sim-reflex.sim[,'ddec']
            dra.lin <- dra-reflex[,'dra']
            ddec.lin <- ddec-reflex[,'ddec']
#out$astrometry[imax,'pmdec']*dt.sim+out$astrometry[out$iref,'parallax']*out$pf[,'dec']+par.opt[paste0('ddec_',i)]
#out$astrometry[imax,'pmra']*dt.sim+out$astrometry[out$iref,'parallax']*out$pf[,'ra']+par.opt[paste0('dra_',i)]
            plot(dra.lin,ddec.lin,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),xlim=rev(range(dra.lin,dra.lin.sim[i2])),ylim=range(ddec.lin,ddec.lin.sim[i2]),main=paste(i,'linear'),pch=20)
            arrows(dra.lin-era,ddec.lin,dra.lin+era,ddec.lin,length=0.001,angle=0,code=1,col='grey')
            arrows(dra.lin,ddec.lin-edec,dra.lin,ddec.lin+edec,length=0.001,angle=0,code=1,col='grey')
            points(0,0,pch='+')
            ncolor <- 20
            cc <- paletteer_c(palette = "viridis::inferno", n = ncolor)
            jj <- as.factor( as.numeric( cut((tepoch-tepoch[1])%%Popt[imax], ncolor)))
            if(grepl('hip',i)){
                segments(dra.lin-era,ddec.lin-edec,dra.lin+era,ddec.lin+edec,col=tcol('grey',50))
            }else{
                arrows(dra.lin-era,ddec.lin,dra.lin+era,ddec.lin,length=0.001,angle=0,code=1,col='grey')
                arrows(dra.lin,ddec.lin-edec,dra.lin,ddec.lin+edec,length=0.001,angle=0,code=1,col='grey')
            }
            ii <- sapply(tepoch,function(t) which.min(abs(tsim-t)))
            segments(dra.lin.sim[ii],ddec.lin.sim[ii],dra.lin,ddec.lin,col=tcol(cc[jj],50))
            lines(dra.lin.sim,ddec.lin.sim,col='red')
            points(dra.lin.sim[ii],ddec.lin.sim[ii],col=cc[jj],pch=20,cex=0.5)

###
            dt <- (tepoch-tepoch[1])/365.25
            if(grepl('hip',i)){
                dra0 <- (astro$barycenter[,'ra']-out$astrometry[out$iref,'ra'])*3.6e6*cos(out$astrometry[out$iref,'dec']/180*pi)+out$pf[out$ind.all$epoch[[i]],'ra']*astro$barycenter[,'parallax']#already considered perspective acceleration
                ddec0 <- (astro$barycenter[,'dec']-out$astrometry[out$iref,'dec'])*3.6e6+out$pf[out$ind.all$epoch[[i]],'dec']*astro$barycenter[,'parallax']
#                dra0 <- astro$barycenter[out$iref,'pmra']*dt+astro$barycenter[out$iref,'parallax']*out$pf[out$ind.epoch[[i]],'ra']+par.opt[paste0('dra_',i)]
#                ddec0 <- out$astrometry[out$iref,'pmdec']*dt+out$astrometry[out$iref,'parallax']*out$pf[out$ind.epoch[[i]],'dec']+par.opt[paste0('ddec_',i)]
                dra.nl <- dra-dra0
                ddec.nl <- ddec-ddec0
            }else{
###parallax projected onto R.A. and Decl.
                pfra.sim <- plx.opt*out.sim$pf[,1]
                pfdec.sim <- plx.opt*out.sim$pf[,2]
#                tobs <- out$tiall[out$ind.all$epoch[[i]],1]
                tobs <- out$BJDtdb[out$ind.all$epoch1[[i]]]
                pfra.obs <- plx.opt*out$pf[out$ind.all$epoch[[i]],1]+dra.res
                pfdec.obs <- plx.opt*out$pf[out$ind.all$epoch[[i]],2]+ddec.res
                plot(out.sim$BJDtdb,pfra.sim,xlab='BJD',ylab=expression(Delta*alpha*'* [mas]'),xlim=range(tobs),ylim=range(pfra.obs,pfra.sim),main=paste(i,'parallax'),type='l',col='red')
                points(tobs,pfra.obs)
                arrows(tobs,pfra.obs-era,tobs,pfra.obs+era,length=0.001,angle=0,code=1,col='grey')
                plot(out.sim$BJDtdb,pfdec.sim,xlab='BJD',ylab=expression(Delta*delta*' [mas]'),xlim=range(tobs),ylim=range(pfdec.obs,pfdec.sim),main=paste(i,'parallax'),type='l',col='red')
                points(tobs,pfdec.obs)
                arrows(tobs,pfdec.obs-era,tobs,pfdec.obs+era,length=0.001,angle=0,code=1,col='grey')
                dra.nl <- dra.res+reflex[,'dra']
                ddec.nl <- ddec.res+reflex[,'ddec']
            }
            dra.nl.sim <- reflex.sim[,'dra']
            ddec.nl.sim <- reflex.sim[,'ddec']

            plot(dra.nl,ddec.nl,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),xlim=rev(range(dra.nl,dra.nl.sim)),ylim=range(ddec.nl,ddec.nl.sim),main=paste(i,'orbital'),pch=20)
            arrows(dra.nl-era,ddec.nl,dra.nl+era,ddec.nl,length=0.001,angle=0,code=1,col='grey')
            arrows(dra.nl,ddec.nl-edec,dra.nl,ddec.nl+edec,length=0.001,angle=0,code=1,col='grey')
            points(0,0,pch='+')
            ncolor <- 20
            cc <- paletteer_c(palette = "viridis::inferno", n = ncolor)
            jj <- as.factor( as.numeric( cut((tepoch-tepoch[1])%%Popt[imax], ncolor)))
#            segments(dra.nl-era,ddec.nl-edec,dra.nl+era,ddec.nl+edec,col=tcol('grey',50))
            ii <- sapply(out$data.epoch[[i]][,1],function(t) which.min(abs(out.sim$BJDtdb-t)))
#            segments(dra.sim[ii],ddec.sim[ii],dra,ddec,col=tcol(cc[jj],50))
            lines(dra.nl.sim,ddec.nl.sim,col='red')
            points(dra.nl.sim[ii],ddec.nl.sim[ii],col=cc[jj],pch=20,cex=0.5)
        }
    }
}
###add GOST fit
if(length(out$gost)>0){
    reflex.gost <- astrometry.epoch(par.opt,tt=out$gost[,'BJD'],bases=bases)$epoch
    bary.gost <- astrometry.bary(par.opt,tt=out$gost[,'BJD'],bases=bases)
    kep <- astrometry.kepler(par.opt,bases=bases)
###observed position and proper motion relative to the predicted barycentric position and proper motion
    dra.obs <- ((out$astrometry[out$astro.index,'ra']-kep$barycenter[out$astro.index,'ra']))*cos(out$astrometry[out$astro.index,'dec']/180*pi)*3.6e6
    ddec.obs <- ((out$astrometry[out$astro.index,'dec']-kep$barycenter[out$astro.index,'dec']))*3.6e6
    dpmra.obs <- out$astrometry[out$astro.index,'pmra']-kep$barycenter[out$astro.index,'pmra']
    dpmdec.obs <- out$astrometry[out$astro.index,'pmdec']-kep$barycenter[out$astro.index,'pmdec']
    dplx.obs <- out$astrometry[out$astro.index,'parallax']-kep$barycenter[out$astro.index,'parallax']
###the position and proper motion reconstructed by a linear fit to gost-based observations generated by accounting for reflex motion
    dra.model <- dra.obs-kep$cats[,1]
    ddec.model <- ddec.obs-kep$cats[,2]
    dplx.model <- dplx.obs-kep$cats[,3]
    dpmra.model <- dpmra.obs-kep$cats[,4]
    dpmdec.model <- dpmdec.obs-kep$cats[,5]

    dra.obs.range <- range(sapply(1:length(dra.obs), function(k) dra.obs[k]+dpmra.obs[k]*c(min(out$dtg[[k]]),max(out$dtg[[k]]))))
    ddec.obs.range <- range(sapply(1:length(dra.obs), function(k) ddec.obs[k]+dpmdec.obs[k]*c(min(out$dtg[[k]]),max(out$dtg[[k]]))))
    dra.model.range <- range(sapply(1:length(dra.model), function(k) dra.model[k]+dpmra.model[k]*c(min(out$dtg[[k]]),max(out$dtg[[k]]))))
    ddec.model.range <- range(sapply(1:length(dra.model), function(k) ddec.model[k]+dpmdec.model[k]*c(min(out$dtg[[k]]),max(out$dtg[[k]]))))
    dra.range <- range(dra.obs.range,dra.model.range)
    ddec.range <- range(ddec.obs.range,ddec.model.range)

####model fit to catalog data visualization
##combined
    dt <- 0.1
    epochs <- range(gdr1.epoch,gdr2.epoch,gdr3.epoch)
    reflex.sim <- astrometry.epoch(par.opt,tt=out.sim$BJDtdb,out1=out.sim)$epoch
    index <- which(out.sim$BJDtdb>=(min(out$astrometry[,'ref_epoch'])-100) & out.sim$BJDtdb<=(max(epochs)+100))
    dtsim <- (out.sim$BJDtdb-out$astrometry[out$iref,'ref_epoch'])/365.25#yr
    pfra.sim <- out.sim$pf[index,'ra']
    pfdec.sim <- out.sim$pf[index,'dec']
    dtsim <- dtsim[index]
    reflex.sim <- reflex.sim[index,]
    bary.sim <- astrometry.bary(par.opt,out1=out.sim,bases=bases,tt=out.sim$BJDtdb[index])
#out.sim$BJDtdb
    rd <- obs.lin.prop(out$astrometry[out$iref,],t=dtsim*365.25,PA=TRUE)
    dras <- (rd[,1]-out$astrometry[out$iref,'ra'])*cos(out$astrometry[out$iref,'dec']/180*pi)*3.6e6+out$astrometry[out$iref,'parallax']*pfra.sim
    ddecs <- (rd[,2]-out$astrometry[out$iref,'dec'])*3.6e6+out$astrometry[out$iref,'parallax']*pfdec.sim
    dra.linear.sim <- (bary.sim[,'ra']-out$astrometry[out$iref,'ra'])*cos(out$astrometry[out$iref,'dec']*pi/180)*3.6e6+bary.sim[,'parallax']*pfra.sim
    ddec.linear.sim <- (bary.sim[,'dec']-out$astrometry[out$iref,'dec'])*3.6e6+bary.sim[,'parallax']*pfdec.sim
    dra.all.sim <- dra.linear.sim+reflex.sim$dra
    ddec.all.sim <- ddec.linear.sim+reflex.sim$ddec
    tgost <- (out$gost[,'BJD']-out$astrometry[out$iref,'ref_epoch'])/365.25#yr
    theta <- out$gost[,'psi']-pi/2
    pfra.gost <- out$gost[,'parf']*cos(theta)+out$gost[,'parx']*sin(theta)
    pfdec.gost <- out$gost[,'parf']*sin(-theta)+out$gost[,'parx']*cos(theta)
    dra.linear.gost <- (bary.gost[,'ra']-out$astrometry[out$iref,'ra'])*cos(out$astrometry[out$iref,'dec']*pi/180)*3.6e6+bary.gost[,'parallax']*pfra.gost
    ddec.linear.gost <- (bary.gost[,'dec']-out$astrometry[out$iref,'dec'])*3.6e6+bary.gost[,'parallax']*pfdec.gost
    dra.all.gost <- dra.linear.gost+reflex.gost[,'dra']
    ddec.all.gost <- ddec.linear.gost+reflex.gost[,'ddec']

    dra1 <- out$astro.gost[,'dra']#-par.opt['dra']
    ddec1 <- out$astro.gost[,'ddec']#-par.opt['ddec']
    plx1 <- out$astro.gost[,'parallax']#-par.opt['dplx']
    pmra1 <- out$astro.gost[,'pmra']#-par.opt['dpmra']
    pmdec1 <- out$astro.gost[,'pmdec']#-par.opt['dpmdec']
    astro.gost <- cbind(dra1,ddec1,plx1,pmra1,pmdec1)
    colnames(astro.gost) <- colnames(out$astro.gost)
    ind.valid <- which(kep$cats[,3]!=0)

    dtobs <- (out$astrometry[out$astro.index,'ref_epoch']-out$astrometry[out$iref,'ref_epoch'])/365.25#yr
    era.obs <- out$astrometry[out$astro.index,'ra_error']
    edec.obs <- out$astrometry[out$astro.index,'dec_error']
    gost.fit <- astro.gost-kep$cats

    if(binary!='' & length(out[[binary]]$dtg)>0){
        mb <- k2m(par.opt['K1'],exp(par.opt['per1']),par.opt['e1'],par.opt['Mstar'],Inc=par.opt['Inc1'])$ms
        mratio <- mb/par.opt['Mstar']
        dra.obs.binary <- ((out[[binary]]$astrometry[out[[binary]]$astro.index,'ra']-kep$barycenter[out[[binary]]$astro.index,'ra']))*cos(out[[binary]]$astrometry[out[[binary]]$astro.index,'dec']/180*pi)*3.6e6
        ddec.obs.binary <- ((out[[binary]]$astrometry[out[[binary]]$astro.index,'dec']-kep$barycenter[out[[binary]]$astro.index,'dec']))*3.6e6
        dpmra.obs.binary <- out[[binary]]$astrometry[out[[binary]]$astro.index,'pmra']-kep$barycenter[out[[binary]]$astro.index,'pmra']
        dpmdec.obs.binary <- out[[binary]]$astrometry[out[[binary]]$astro.index,'pmdec']-kep$barycenter[out[[binary]]$astro.index,'pmdec']
        dplx.obs.binary <- out[[binary]]$astrometry[out[[binary]]$astro.index,'parallax']-kep$barycenter[out[[binary]]$astro.index,'parallax']
###the position and proper motion reconstructed by a linear fit to gost-based observations generated by accounting for reflex motion
        dra.model.binary <- dra.obs.binary-kep$binarycats[,1]
        ddec.model.binary <- ddec.obs.binary-kep$binarycats[,2]
        dplx.model.binary <- dplx.obs.binary-kep$binarycats[,3]
        dpmra.model.binary <- dpmra.obs.binary-kep$binarycats[,4]
        dpmdec.model.binary <- dpmdec.obs.binary-kep$binarycats[,5]

        dra.obs.range.binary <- range(sapply(1:length(dra.obs.binary), function(k) dra.obs.binary[k]+dpmra.obs.binary[k]*c(min(out[[binary]]$dtg[[k]]),max(out[[binary]]$dtg[[k]]))))
        ddec.obs.range.binary <- range(sapply(1:length(ddec.obs.binary), function(k) ddec.obs.binary[k]+dpmdec.obs.binary[k]*c(min(out[[binary]]$dtg[[k]]),max(out[[binary]]$dtg[[k]]))))
        dra.model.range.binary <- range(sapply(1:length(dra.model.binary), function(k) dra.model.binary[k]+dpmra.model.binary[k]*c(min(out[[binary]]$dtg[[k]]),max(out[[binary]]$dtg[[k]]))))
        ddec.model.range.binary <- range(sapply(1:length(dra.model.binary), function(k) ddec.model.binary[k]+dpmdec.model.binary[k]*c(min(out[[binary]]$dtg[[k]]),max(out[[binary]]$dtg[[k]]))))
        dra.range.binary <- range(dra.obs.range.binary,dra.model.range.binary)
        ddec.range.binary <- range(ddec.obs.range.binary,ddec.model.range.binary)

        par.opt1 <- par.opt
        if(Nsig>1){
            ind.rm <- 8:(Nsig*7)
            par.opt1 <- par.opt[-ind.rm]
        }
        reflex.gost.binary <- -astrometry.epoch(par.opt1,out1=out,tt=out[[binary]]$gost[,'BJD'],bases=bases)$epoch/mratio
        index.binary <- which(out.sim$BJDtdb>=(min(out[[binary]]$astrometry[,'ref_epoch'])-100) & out.sim$BJDtdb<=(max(epochs)+100))
        rd.binary <- obs.lin.prop(out[[binary]]$astrometry[out[[binary]]$iref,],out.sim$BJDtdb[index.binary]-out[[binary]]$astrometry[out[[binary]]$iref,'ref_epoch'],PA=TRUE)
        dras.binary <- (rd.binary[,1]-out[[binary]]$astrometry[out[[binary]]$iref,'ra'])*cos(out[[binary]]$astrometry[out[[binary]]$iref,'dec']/180*pi)*3.6e6+out[[binary]]$astrometry[out[[binary]]$iref,'parallax']*out$pf[index.binary,'ra']
        ddecs.binary <- (rd.binary[,2]-out[[binary]]$astrometry[out[[binary]]$iref,'dec'])*3.6e6+out[[binary]]$astrometry[out[[binary]]$iref,'parallax']*out$pf[index.binary,'dec']
        reflex.sim.binary <- -reflex.sim/mratio
        dra.linear.sim.binary <- (bary.sim[,'ra']-out[[binary]]$astrometry[out[[binary]]$iref,'ra'])*cos(out[[binary]]$astrometry[out[[binary]]$iref,'dec']*pi/180)*3.6e6+bary.sim[,'parallax']*pfra.sim
        ddec.linear.sim.binary <- (bary.sim[,'dec']-out[[binary]]$astrometry[out[[binary]]$iref,'dec'])*3.6e6+bary.sim[,'parallax']*pfdec.sim
        dra.all.sim.binary <- dra.linear.sim.binary+reflex.sim.binary$dra[index.binary]
        ddec.all.sim.binary <- ddec.linear.sim.binary+reflex.sim.binary$ddec[index.binary]
        tgost.binary <- (out[[binary]]$gost[,'BJD']-out[[binary]]$astrometry[out[[binary]]$iref,'ref_epoch'])/365.25#yr
        theta <- out[[binary]]$gost[,'psi']-pi/2
        pfra.gost.binary <- out[[binary]]$gost[,'parf']*cos(theta)+out[[binary]]$gost[,'parx']*sin(theta)
        pfdec.gost.binary <- out[[binary]]$gost[,'parf']*sin(-theta)+out[[binary]]$gost[,'parx']*cos(theta)

        dra.linear.sim.binary <- (bary.sim[,'ra']-out[[binary]]$astrometry[out[[binary]]$iref,'ra'])*cos(out[[binary]]$astrometry[out[[binary]]$iref,'dec']*pi/180)*3.6e6+bary.sim[,'parallax']*pfra.sim
        ddec.linear.sim.binary <- (bary.sim[,'dec']-out[[binary]]$astrometry[out[[binary]]$iref,'dec'])*3.6e6+bary.sim[,'parallax']*pfdec.sim
        dra.linear.gost.binary <- (bary.gost[,'ra']-out[[binary]]$astrometry[out[[binary]]$iref,'ra'])*cos(out[[binary]]$astrometry[out[[binary]]$iref,'dec']*pi/180)*3.6e6+bary.gost[,'parallax']*pfra.gost.binary
        ddec.linear.gost.binary <- (bary.gost[,'dec']-out[[binary]]$astrometry[out[[binary]]$iref,'dec'])*3.6e6+bary.gost[,'parallax']*pfdec.gost.binary
        dra.all.gost.binary <- dra.linear.gost.binary+reflex.gost.binary[,'dra']
        ddec.all.gost.binary <- ddec.linear.gost.binary+reflex.gost.binary[,'ddec']
        dra1.binary <- out[[binary]]$astro.gost[,'dra']#-par.opt['dra']
        ddec1.binary <- out[[binary]]$astro.gost[,'ddec']#-par.opt['ddec']
        plx1.binary <- out[[binary]]$astro.gost[,'parallax']#-par.opt['dplx']
        pmra1.binary <- out[[binary]]$astro.gost[,'pmra']#-par.opt['dpmra']
        pmdec1.binary <- out[[binary]]$astro.gost[,'pmdec']#-par.opt['dpmdec']
        astro.gost.binary <- cbind(dra1.binary,ddec1.binary,plx1.binary,pmra1.binary,pmdec1.binary)
        colnames(astro.gost.binary) <- colnames(out[[binary]]$astro.gost)
        ind.valid.binary <- which(kep$binarycats[,3]!=0)

        dtobs.binary <- (out[[binary]]$astrometry[out[[binary]]$astro.index,'ref_epoch']-out[[binary]]$astrometry[out[[binary]]$iref,'ref_epoch'])/365.25#yr
        era.obs.binary <- out[[binary]]$astrometry[out[[binary]]$astro.index,'ra_error']
        edec.obs.binary <- out[[binary]]$astrometry[out[[binary]]$astro.index,'dec_error']
        gost.fit.binary <- astro.gost.binary-kep$binarycats
    }

#if(FALSE){
if(TRUE){
    ra.lab <- paste0('(R.A.- ',out$astrometry[out$iref,'ra'],')*cos(DEC) [mas]')
    dec.lab <- paste0('Decl.- ',out$astrometry[out$iref,'dec'],' [mas]')
    tlab <- 'Time - J2016 [yr]'
    pos <- 1
    ccs <- c('red','blue','black','grey','green','darkgrey','pink','steelblue')

    for(i3 in 1:2){
        if(i3==2 & min(dtobs)<min(epochs)){
            zoomin <- TRUE
            dtlim <- range((epochs-out$astrometry[out$iref,'ref_epoch'])/365.25)
            ind.sim <- which(dtsim>=dtlim[1])
            ind.gost <- which(dtobs>=dtlim[1])
            dra.lim <- range(ddec.all.sim[ind.sim],astro.gost[ind.gost,'dra']-era.obs[ind.gost],astro.gost[ind.gost,'dra']+era.obs[ind.gost],gost.fit[ind.gost,'dra']-gost.fit[ind.gost,'pmra']*dt,gost.fit[ind.gost,'dra']+gost.fit[ind.gost,'pmra']*dt,dra.linear.sim)
            ddec.lim <- range(ddec.all.sim[ind.sim],astro.gost[ind.gost,'ddec']-edec.obs[ind.gost],astro.gost[ind.gost,'ddec']+edec.obs[ind.gost],gost.fit[ind.gost,'ddec']-edec.obs[ind.gost],gost.fit[ind.gost,'ddec']-gost.fit[ind.gost,'pmdec']*dt,gost.fit[ind.gost,'ddec']+gost.fit[ind.gost,'pmdec']*dt,ddec.linear.sim)
        }else{
            zoomin <- FALSE
            dra.lim <- range(dra.all.sim,astro.gost[,'dra']-era.obs,astro.gost[,'dra']+era.obs,gost.fit[,'dra']-gost.fit[,'pmra']*dt,gost.fit[,'dra']+gost.fit[,'pmra']*dt,dra.linear.sim)
            ddec.lim <- range(ddec.all.sim,astro.gost[,'ddec']-edec.obs,astro.gost[,'ddec']+edec.obs,gost.fit[,'ddec']-edec.obs,gost.fit[,'ddec']-gost.fit[,'pmdec']*dt,gost.fit[,'ddec']+gost.fit[,'pmdec']*dt,ddec.linear.sim)
            dtlim <- range(dtsim)
        }
        plot(dra.all.sim,ddec.all.sim,xlab=ra.lab,ylab=dec.lab,type='l',col=ccs[1],xlim=rev(dra.lim),ylim=ddec.lim)#barycenter+reflex motion
        lines(dras,ddecs,col=ccs[8])#steelblue;
        if(zoomin) legend('topright',bty='n',legend='Zoom into GDRs')
        points(dra.all.gost,ddec.all.gost,col=ccs[4],pch=20)#grey; gost emulation
        lines(dra.linear.sim,ddec.linear.sim,col=ccs[2])#blue; barycentric motion
        points(dra1,ddec1,pch=20,col=ccs[3])#black;
        arrows(dra1,ddec1-edec.obs,dra1,ddec1+edec.obs,length=0.001,angle=0,code=1,col=ccs[3])
        arrows(dra1-era.obs,ddec1,dra1+era.obs,ddec1,length=0.001,angle=0,code=1,col=ccs[3])
        dastros <- list()
        nsamp <- 10
        for(kk in 1:nsamp){
            dastro <- t(sapply(1:nrow(astro.gost), function(k3) mvtnorm::rmvnorm(1,mean=astro.gost[k3,], sigma=out$cov.astro[1:5,1:5,out$astro.index[k3]])))
            colnames(dastro) <- c('dra','ddec','plx','pmra','pmdec')
            dastros[[kk]] <- dastro
            dastro <- dastro[ind.valid,,drop=FALSE]
                                        #        points(dastro[,'dra'],dastro[,'ddec'],pch=20,col=tcol(ccs[6],50))
            segments(dastro[,'dra']-dastro[,'pmra']*dt,dastro[,'ddec']-dastro[,'pmdec']*dt,dastro[,'dra']+dastro[,'pmra']*dt,dastro[,'ddec']+dastro[,'pmdec']*dt,col=tcol(ccs[6],50))
        }
        segments(dra1[ind.valid]-pmra1[ind.valid]*dt,ddec1[ind.valid]-pmdec1[ind.valid]*dt,dra1[ind.valid]+pmra1[ind.valid]*dt,ddec1[ind.valid]+pmdec1[ind.valid]*dt,col=ccs[3])
        points(gost.fit[,'dra'],gost.fit[,'ddec'],pch=20,col=ccs[7])
        segments(gost.fit[ind.valid,'dra']-gost.fit[ind.valid,'pmra']*dt,gost.fit[ind.valid,'ddec']-gost.fit[ind.valid,'pmdec']*dt,gost.fit[ind.valid,'dra']+gost.fit[ind.valid,'pmra']*dt,gost.fit[ind.valid,'ddec']+gost.fit[ind.valid,'pmdec']*dt,col=ccs[7])
        text(astro.gost[,'dra'],astro.gost[,'ddec'],labels=out$cats,pos=pos,xpd=NA,col=ccs[5],cex=0.8)
        legend('top',horiz=TRUE,xpd=NA,legend=c('Cobmined motion','Barycentric motion','GDR3 solution'),col=ccs[c(1:2,8)],lty=c(1,1,1),pch=c(NA,NA,NA),inset=c(0,-0.3),bty='n')
        legend('top',horiz=TRUE,xpd=NA,legend=c('Catalog Data','GOST Epoch'),col=ccs[3:4],lty=c(1,NA),pch=c(20,20),inset=c(0,-0.2),bty='n')
        legend('top',horiz=TRUE,xpd=NA,legend=c('Catalog Sample','GOST Fit'),col=ccs[6:7],lty=c(1,1),pch=c(NA,20),inset=c(0,-0.1),bty='n')

        plot(dtsim,dra.all.sim,xlab=tlab,ylab=ra.lab,type='l',col=ccs[1],ylim=dra.lim,xlim=dtlim)
        lines(dtsim,dras,col=ccs[8])
        if(zoomin) legend('topright',bty='n',legend='Zoom into GDRs')
        points(tgost,dra.all.gost,col=ccs[4],pch=20)
        lines(dtsim,dra.linear.sim,col=ccs[2])
        points(dtobs,astro.gost[,'dra'],pch=20,col=ccs[3])
        arrows(dtobs,astro.gost[,'dra']-era.obs,dtobs,astro.gost[,'dra']+era.obs,length=0.001,angle=0,code=1,col=ccs[3])
        segments(dtobs[ind.valid]-dt,dra1[ind.valid]-pmra1[ind.valid]*dt,dtobs[ind.valid]+dt,dra1[ind.valid]+pmra1[ind.valid]*dt,col=ccs[3])
        for(kk in 1:nsamp){
            dastro <- dastros[[kk]][ind.valid,,drop=FALSE]
            segments(dtobs[ind.valid]-dt,dastro[,'dra']-dastro[,'pmra']*dt,dtobs[ind.valid]+dt,dastro[,'dra']+dastro[,'pmra']*dt,col=tcol(ccs[6],50))
        }
        points(dtobs,gost.fit[,'dra'],pch=20,col=ccs[7])
        segments(dtobs[ind.valid]-dt,gost.fit[ind.valid,'dra']-gost.fit[ind.valid,'pmra']*dt,dtobs[ind.valid]+dt,gost.fit[ind.valid,'dra']+gost.fit[ind.valid,'pmra']*dt,col=ccs[7])
        text(dtobs,astro.gost[,'dra'],labels=out$cats,pos=pos,xpd=NA,col=ccs[5],cex=0.8)
        plot(dtsim,ddec.all.sim,xlab=tlab,ylab=dec.lab,type='l',col=ccs[1],ylim=ddec.lim,xlim=dtlim)
        lines(dtsim,ddecs,col=ccs[8])
        if(zoomin) legend('topright',bty='n',legend='Zoom into GDRs')
        lines(dtsim,ddec.linear.sim,col=ccs[2])
        points(tgost,ddec.all.gost,col=ccs[4],pch=20)
        points(dtobs,astro.gost[,'ddec'],pch=20,col=ccs[3])
        arrows(dtobs,astro.gost[,'ddec']-edec.obs,dtobs,astro.gost[,'ddec']+edec.obs,length=0.001,angle=0,code=1,col=ccs[3])
        segments(dtobs[ind.valid]-dt,ddec1[ind.valid]-pmdec1[ind.valid]*dt,dtobs[ind.valid]+dt,ddec1[ind.valid]+pmdec1[ind.valid]*dt,col=ccs[3])
        for(kk in 1:nsamp){
            dastro <- dastros[[kk]][ind.valid,,drop=FALSE]
            segments(dtobs[ind.valid]-dt,dastro[,'ddec']-dastro[,'pmdec']*dt,dtobs[ind.valid]+dt,dastro[,'ddec']+dastro[,'pmdec']*dt,col=tcol(ccs[6],50))
        }
        points(dtobs,gost.fit[,'ddec'],pch=20,col=ccs[7])
        segments(dtobs[ind.valid]-dt,gost.fit[ind.valid,'ddec']-gost.fit[ind.valid,'pmdec']*dt,dtobs[ind.valid]+dt,gost.fit[ind.valid,'ddec']+gost.fit[ind.valid,'pmdec']*dt,col=ccs[7])
        text(dtobs,astro.gost[,'ddec'],labels=out$cats,pos=pos,xpd=NA,col=ccs[5],cex=0.8)
    }
##reflex motion fit
#    reflex.sim <- astrometry.epoch(par.opt,tt=out$gost[,'BJD'])$epoch
    plot(reflex.sim[,'dra'],reflex.sim[,'ddec'],xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),type='l',col='grey',xlim=rev(range(reflex.sim[,'dra'],dra.range)),ylim=range(reflex.sim[,'ddec'],ddec.range))
#    points(0,0,pch='+')
    jj <- as.factor(((out$gost[,1]-out$gost[1,1])%%max(Popt))/max(Popt))
    cc <- paletteer_c(palette = "viridis::inferno", n = length(jj))
    points(reflex.gost[,'dra'],reflex.gost[,'ddec'],col=cc[jj],pch=20)
#    source('goodness_fit.R')
    if(length(out$cats)==0) out$cats <- c('GDR2','GDR3')
    phases <- seq(0,1,by=0.2)
    cl <- paletteer_c(palette = "viridis::inferno", n = length(phases))
    legend('top',legend=phases,pch=20,col=cl[as.factor(phases)],horiz=TRUE,bty='n',inset=c(0,-0.25),xpd=NA,title='Phase')
    ref.epochs <- c()
    cols <- c('darkgreen','blue','red','green','cyan')
#    dra.model.comb <- c()
    for(k in 1:length(out$astro.index)){
        ttsim <- seq(min(out$dtg[[k]]),max(out$dtg[[k]]),by=0.1)#year
        dras.obs <- dra.obs[k]+dpmra.obs[k]*ttsim
        ddecs.obs <- ddec.obs[k]+dpmdec.obs[k]*ttsim
        if(!any(out$cats[k]==c('GDR1','TYC'))) lines(dras.obs,ddecs.obs,col=cols[k],lty=2)
        points(dra.obs[k],ddec.obs[k],col=cols[k],pch=1)

        dras.model <- dra.model[k]+dpmra.model[k]*ttsim
        ddecs.model <- ddec.model[k]+dpmdec.model[k]*ttsim
        if(!any(out$cats[k]==c('GDR1','TYC'))) lines(dras.model,ddecs.model,col=cols[k],lty=1)
        points(dra.model[k],ddec.model[k],col=cols[k],pch=20)
    }
    legend('top',legend=out$cats,col=cols,lty=1,cex=0.8,bty='n',inset=c(0,-0.1),xpd=NA,horiz=TRUE)
#####individual fit
    trefs <- out$astrometry[out$astro.index,'ref_epoch']
    obs <- cbind(dra.obs,ddec.obs,dplx.obs,dpmra.obs,dpmdec.obs)
    eobs <- out$astrometry[out$astro.index,c('ra_error','dec_error','parallax_error','pmra_error','pmdec_error')]
    model <- cbind(dra.model,ddec.model,dplx.model,dpmra.model,dpmdec.model)
    ylabs <- c(expression(Delta*alpha*'* [mas]'),expression(Delta*delta*' [mas]'),'parallax [mas]',expression(Delta*mu[alpha]*' [mas/yr]'),expression(Delta*mu[delta]*' [mas/yr]'))
    ind.rm <- c(out$ityc,out$igdr1)
    for(j in 1:ncol(obs)){
        inds <- 1:nrow(obs)
        if(j>2 & length(ind.rm)>0){
            inds <- inds[-ind.rm]
        }
        plot(trefs[inds],obs[inds,j],xlab='BJD',ylab=as.expression(ylabs[j]),ylim=range(obs[inds,j],model[inds,j]),col=cc[inds])
        arrows(trefs[inds],obs[inds,j]-eobs[inds,j],trefs[inds],obs[inds,j]+eobs[inds,j],length=0.001,angle=0,code=1,col='grey')
        points(trefs[inds],model[inds,j],col='red')
    }
}
###################
####binary plot
    if(binary!='' & length(out[[binary]]$dtg)>0){
####proper motion + parallax + reflex motion
        for(i3 in 1:2){
            if(i3==2 & min(dtobs.binary)<min(epochs)){
                zoomin <- TRUE
                dtlim.binary <- range((epochs-out[[binary]]$astrometry[out[[binary]]$iref,'ref_epoch'])/365.25)
                ind.sim.binary <- which(dtsim>=dtlim.binary[1])
                ind.gost.binary <- which(dtobs.binary>=dtlim.binary[1])
                dra.lim.binary <- range(ddec.all.sim.binary[ind.sim.binary],astro.gost.binary[ind.gost,'dra']-era.obs.binary[ind.gost],astro.gost.binary[ind.gost,'dra']+era.obs.binary[ind.gost],gost.fit.binary[ind.gost,'dra']-gost.fit.binary[ind.gost,'pmra']*dt,gost.fit.binary[ind.gost,'dra']+gost.fit.binary[ind.gost,'pmra']*dt)
                ddec.lim.binary <- range(ddec.all.sim[ind.sim.binary],astro.gost.binary[ind.gost,'ddec']-edec.obs.binary[ind.gost],astro.gost.binary[ind.gost,'ddec']+edec.obs.binary[ind.gost],gost.fit.binary[ind.gost,'ddec']-edec.obs.binary[ind.gost],gost.fit.binary[ind.gost,'ddec']-gost.fit.binary[ind.gost,'pmdec']*dt,gost.fit.binary[ind.gost,'ddec']+gost.fit.binary[ind.gost,'pmdec']*dt)
            }else{
                zoomin <- FALSE
                dra.lim.binary <- range(dra.all.sim,astro.gost.binary[,'dra']-era.obs,astro.gost.binary[,'dra']+era.obs,gost.fit.binary[,'dra']-gost.fit.binary[,'pmra']*dt,gost.fit.binary[,'dra']+gost.fit.binary[,'pmra']*dt)
                ddec.lim.binary <- range(ddec.all.sim,astro.gost.binary[,'ddec']-edec.obs,astro.gost.binary[,'ddec']+edec.obs,gost.fit.binary[,'ddec']-edec.obs,gost.fit.binary[,'ddec']-gost.fit.binary[,'pmdec']*dt,gost.fit.binary[,'ddec']+gost.fit.binary[,'pmdec']*dt)
                dtlim.binary <- range(dtsim)
            }
                                        #        dra.lim.binary <- rev(-dra.lim/mratio)
                                        #        ddec.lim.binary <- sort(-ddec.lim/mratio)
            plot(dra.all.sim.binary,ddec.all.sim.binary,xlab=ra.lab,ylab=dec.lab,type='l',col=ccs[1],xlim=rev(dra.lim.binary),ylim=ddec.lim.binary,main=binary)
            lines(dras.binary,ddecs.binary,col=ccs[8])
            if(zoomin) legend('topright',bty='n',legend='Zoom into GDRs')
            points(dra.all.gost.binary,ddec.all.gost.binary,col=ccs[4],pch=20)
            lines(dra.linear.sim.binary,ddec.linear.sim.binary,col=ccs[2])
            points(dra1.binary,ddec1.binary,pch=20,col=ccs[3])
            arrows(dra1.binary,ddec1.binary-edec.obs.binary,dra1,ddec1+edec.obs,length=0.001,angle=0,code=1,col=ccs[3])
            arrows(dra1.binary-era.obs.binary,ddec1.binary,dra1.binary+era.obs.binary,ddec1.binary,length=0.001,angle=0,code=1,col=ccs[3])
            dastros.binary <- list()
            nsamp <- 10
            for(kk in 1:nsamp){
                dastro <- t(sapply(1:nrow(astro.gost.binary), function(k3) mvtnorm::rmvnorm(1,mean=astro.gost.binary[k3,], sigma=out[[binary]]$cov.astro[1:5,1:5,out[[binary]]$astro.index[k3]])))
                colnames(dastro) <- c('dra','ddec','plx','pmra','pmdec')
                dastros.binary[[kk]] <- dastro
                dastro <- dastro[ind.valid.binary,,drop=FALSE]
                                        #        points(dastro[,'dra'],dastro[,'ddec'],pch=20,col=tcol(ccs[6],50))
                segments(dastro[,'dra']-dastro[,'pmra']*dt,dastro[,'ddec']-dastro[,'pmdec']*dt,dastro[,'dra']+dastro[,'pmra']*dt,dastro[,'ddec']+dastro[,'pmdec']*dt,col=tcol(ccs[6],50))
            }
            segments(dra1.binary[ind.valid.binary]-pmra1.binary[ind.valid.binary]*dt,ddec1.binary[ind.valid.binary]-pmdec1.binary[ind.valid.binary]*dt,dra1.binary[ind.valid.binary]+pmra1.binary[ind.valid.binary]*dt,ddec1[ind.valid.binary]+pmdec1[ind.valid.binary]*dt,col=ccs[3])
            points(gost.fit.binary[,'dra'],gost.fit.binary[,'ddec'],pch=20,col=ccs[7])
            segments(gost.fit.binary[ind.valid.binary,'dra']-gost.fit.binary[ind.valid.binary,'pmra']*dt,gost.fit.binary[ind.valid.binary,'ddec']-gost.fit.binary[ind.valid.binary,'pmdec']*dt,gost.fit.binary[ind.valid.binary,'dra']+gost.fit.binary[ind.valid.binary,'pmra']*dt,gost.fit.binary[ind.valid.binary,'ddec']+gost.fit.binary[ind.valid.binary,'pmdec']*dt,col=ccs[7])
            text(astro.gost.binary[,'dra'],astro.gost.binary[,'ddec'],labels=out[[binary]]$cats,pos=pos,xpd=NA,col=ccs[5],cex=0.8)
            legend('top',horiz=TRUE,xpd=NA,legend=c('Cobmined motion','Barycentric motion','GDR3 solution'),col=ccs[c(1:2,8)],lty=c(1,1,1),pch=c(NA,NA,NA),inset=c(0,-0.3),bty='n')
            legend('top',horiz=TRUE,xpd=NA,legend=c('Catalog Data','GOST Epoch'),col=ccs[3:4],lty=c(1,NA),pch=c(20,20),inset=c(0,-0.2),bty='n')
            legend('top',horiz=TRUE,xpd=NA,legend=c('Catalog Sample','GOST Fit'),col=ccs[6:7],lty=c(1,1),pch=c(NA,20),inset=c(0,-0.1),bty='n')

            plot(dtsim,dra.all.sim.binary,xlab=tlab,ylab=ra.lab,type='l',col=ccs[1],ylim=dra.lim.binary,xlim=dtlim,main=binary)
            lines(dtsim,dras.binary,col=ccs[8])
            if(zoomin) legend('topright',bty='n',legend='Zoom into GDRs')
            points(tgost.binary,dra.all.gost.binary,col=ccs[4],pch=20)
            lines(dtsim,dra.linear.sim.binary,col=ccs[2])
            points(dtobs.binary,astro.gost.binary[,'dra'],pch=20,col=ccs[3])
            arrows(dtobs.binary,astro.gost.binary[,'dra']-era.obs.binary,dtobs.binary,astro.gost.binary[,'dra']+era.obs.binary,length=0.001,angle=0,code=1,col=ccs[3])
            segments(dtobs.binary[ind.valid.binary]-dt,dra1.binary[ind.valid]-pmra1.binary[ind.valid.binary]*dt,dtobs.binary[ind.valid]+dt,dra1.binary[ind.valid]+pmra1.binary[ind.valid]*dt,col=ccs[3])
            for(kk in 1:nsamp){
                dastro.binary <- dastros.binary[[kk]][ind.valid.binary,,drop=FALSE]
                segments(dtobs.binary[ind.valid.binary]-dt,dastro.binary[,'dra']-dastro.binary[,'pmra']*dt,dtobs.binary[ind.valid.binary]+dt,dastro.binary[,'dra']+dastro.binary[,'pmra']*dt,col=tcol(ccs[6],50))
            }
            points(dtobs.binary,gost.fit.binary[,'dra'],pch=20,col=ccs[7])
            segments(dtobs.binary[ind.valid.binary]-dt,gost.fit.binary[ind.valid.binary,'dra']-gost.fit.binary[ind.valid.binary,'pmra']*dt,dtobs.binary[ind.valid.binary]+dt,gost.fit.binary[ind.valid.binary,'dra']+gost.fit.binary[ind.valid.binary,'pmra']*dt,col=ccs[7])
            text(dtobs.binary,astro.gost.binary[,'dra'],labels=out$cats,pos=pos,xpd=NA,col=ccs[5],cex=0.8)

            plot(dtsim,ddec.all.sim.binary,xlab=tlab,ylab=dec.lab,type='l',col=ccs[1],ylim=ddec.lim,xlim=dtlim)
            lines(dtsim,ddecs.binary,col=ccs[8])
            if(zoomin) legend('topright',bty='n',legend='Zoom into GDRs')
            lines(dtsim,ddec.linear.sim.binary,col=ccs[2])
            points(tgost.binary,ddec.all.gost.binary,col=ccs[4],pch=20)
            points(dtobs.binary,astro.gost.binary[,'ddec'],pch=20,col=ccs[3])
            arrows(dtobs.binary,astro.gost.binary[,'ddec']-edec.obs.binary,dtobs.binary,astro.gost.binary[,'ddec']+edec.obs.binary,length=0.001,angle=0,code=1,col=ccs[3])
            segments(dtobs.binary[ind.valid.binary]-dt,ddec1.binary[ind.valid.binary]-pmdec1.binary[ind.valid.binary]*dt,dtobs.binary[ind.valid.binary]+dt,ddec1.binary[ind.valid.binary]+pmdec1.binary[ind.valid.binary]*dt,col=ccs[3])
            for(kk in 1:nsamp){
                dastro.binary <- dastros.binary[[kk]][ind.valid.binary,,drop=FALSE]
                segments(dtobs.binary[ind.valid.binary]-dt,dastro.binary[,'ddec']-dastro.binary[,'pmdec']*dt,dtobs.binary[ind.valid]+dt,dastro.binary[,'ddec']+dastro.binary[,'pmdec']*dt,col=tcol(ccs[6],50))
            }
            points(dtobs.binary,gost.fit.binary[,'ddec'],pch=20,col=ccs[7])
            segments(dtobs.binary[ind.valid.binary]-dt,gost.fit.binary[ind.valid.binary,'ddec']-gost.fit.binary[ind.valid.binary,'pmdec']*dt,dtobs.binary[ind.valid.binary]+dt,gost.fit.binary[ind.valid.binary,'ddec']+gost.fit.binary[ind.valid.binary,'pmdec']*dt,col=ccs[7])
            text(dtobs.binary,astro.gost.binary[,'ddec'],labels=out[[binary]]$cats,pos=pos,xpd=NA,col=ccs[5],cex=0.8)
        }

####reflex motion fit
        plot(reflex.sim.binary[,'dra'],reflex.sim.binary[,'ddec'],xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),type='l',col='grey',xlim=rev(range(reflex.sim.binary[,'dra'],dra.range.binary)),ylim=range(reflex.sim.binary[,'ddec'],ddec.range.binary),main=binary)
        points(0,0,pch='+')
        jj <- as.factor(((out[[binary]]$gost[,1]-out[[binary]]$gost[1,1])%%max(Popt))/max(Popt))
        cc <- paletteer_c(palette = "viridis::inferno", n = length(jj))
        points(reflex.gost.binary[,'dra'],reflex.gost.binary[,'ddec'],col=cc[jj],pch=20)
                                        #    source('goodness_fit.R')
        if(length(out[[binary]]$cats)==0) out[[binary]]$cats <- c('GDR2','GDR3')
        phases <- seq(0,1,by=0.2)
        cl <- paletteer_c(palette = "viridis::inferno", n = length(phases))
        legend('top',legend=phases,pch=20,col=cl[as.factor(phases)],horiz=TRUE,bty='n',inset=c(0,-0.25),xpd=NA,title='Phase')
        ref.epochs <- c()
        cols <- c('darkgreen','blue','red','green','cyan')
                                        #    dra.model.comb <- c()
        if(length(out[[binary]]$astro.index)>0){
            for(k in 1:length(out[[binary]]$astro.index)){
                ttsim <- seq(min(out[[binary]]$dtg[[k]]),max(out[[binary]]$dtg[[k]]),by=0.1)#year
                dras.obs.binary <- dra.obs.binary[k]+dpmra.obs.binary[k]*ttsim
                ddecs.obs.binary <- ddec.obs.binary[k]+dpmdec.obs.binary[k]*ttsim
                if(!any(out[[binary]]$cats[k]==c('GDR1','TYC'))) lines(dras.obs.binary,ddecs.obs.binary,col=cols[k],lty=2)
                points(dra.obs.binary[k],ddec.obs.binary[k],col=cols[k],pch=1)

                dras.model.binary <- dra.model.binary[k]+dpmra.model.binary[k]*ttsim
                ddecs.model.binary <- ddec.model.binary[k]+dpmdec.model.binary[k]*ttsim
                if(!any(out[[binary]]$cats[k]==c('GDR1','TYC'))) lines(dras.model.binary,ddecs.model.binary,col=cols[k],lty=1)
                points(dra.model.binary[k],ddec.model.binary[k],col=cols[k],pch=20)
            }
            legend('top',legend=out[[binary]]$cats,col=cols,lty=1,cex=0.8,bty='n',inset=c(0,-0.1),xpd=NA,horiz=TRUE)
#####individual fit
            trefs.binary <- out[[binary]]$astrometry[out[[binary]]$astro.index,'ref_epoch']
            obs.binary <- cbind(dra.obs.binary,ddec.obs.binary,dplx.obs.binary,dpmra.obs.binary,dpmdec.obs.binary)
            eobs.binary <- out[[binary]]$astrometry[out[[binary]]$astro.index,c('ra_error','dec_error','parallax_error','pmra_error','pmdec_error')]
            model.binary <- cbind(dra.model.binary,ddec.model.binary,dplx.model.binary,dpmra.model.binary,dpmdec.model.binary)
            ylabs <- c(expression(Delta*alpha*'* [mas]'),expression(Delta*delta*' [mas]'),'parallax [mas]',expression(Delta*mu[alpha]*' [mas/yr]'),expression(Delta*mu[delta]*' [mas/yr]'))
            ind.rm <- c(out[[binary]]$ityc,out[[binary]]$igdr1)
            for(j in 1:ncol(obs.binary)){
                inds <- 1:nrow(obs.binary)
                if(j>2 & length(ind.rm)>0){
                    inds <- inds[-ind.rm]
                }
                plot(trefs.binary[inds],obs.binary[inds,j],xlab='BJD',ylab=as.expression(ylabs[j]),ylim=range(obs.binary[inds,j],model.binary[inds,j]),col=cc[inds],main=binary)
                arrows(trefs.binary[inds],obs.binary[inds,j]-eobs.binary[inds,j],trefs.binary[inds],obs.binary[inds,j]+eobs.binary[inds,j],length=0.001,angle=0,code=1,col='grey')
                points(trefs.binary[inds],model.binary[inds,j],col='red')
            }
        }
    }
}
if(ruweDR>1){
    kep <- astrometry.kepler(par.opt)
    ts <- uwes <- ruwes <- sds <- c()
    if(ruweDR==2 | ruweDR==23){
        ruwes <- c(ruwes,ruwe2)
        dabs2 <- kep$dabs$GDR2
        uwes <- c(uwes,calc.uwe(dabs2,sfov2,Nfov2,Nccd2,Npar=5))
        sds <- c(sds,0.25)
        ts <- c(ts,astrometry[igdr2,'ref_epoch'])
    }
    if(ruweDR==3 | ruweDR==23){
        ruwes <- c(ruwes,ruwe3)
        dabs3 <- kep$dabs$GDR3
        uwes <- c(uwes,calc.uwe(dabs3,sfov3,Nfov3,Nccd3,Npar=5))
        ts <- c(ts,astrometry[igdr3,'ref_epoch'])
        sds <- c(sds,0.14)
    }
    plot(ts,ruwes,xlab='Epoch',ylab='RUWE(black) or UWE(red)',log='y',ylim=range(ruwes,exp(log(ruwes)-sds),uwes,exp(log(ruwes)+sds)))
    arrows(ts,exp(log(ruwes)-sds),ts,exp(log(ruwes)+sds),length=0.001,angle=0,code=1,col='grey')
    points(ts,uwes,col='red')
}
if(plotf){
    dev.off()
}
