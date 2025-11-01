#library("viridis")
plotf <- TRUE
#plotf <- FALSE
if(plotf){
cat('test.pdf\n')
pdf('test.pdf',8,8)
par(mfrow=c(2,2))
}
yr2d <- 365.25
ng <- 'G'
nast <- length(out$astro.index)
nastro <- nrow(out$astrometry)
poss <- (1:4)[1:nastro]
#ell.col <- c('orange','darkgreen')
#ell.col <- rainbow(nastro)
ell.col <- terrain.colors(nastro)
if(nastro>2) ng <- out$gdrs
cc <- c('black','blue','orange','steelblue','green','darkgreen',rainbow(10))
#cc <- c('black','blue','cyan','steelblue','orange','darkgreen')
##astrometry fit
imax <- which.max(Popt)
if(out$Nastro>0 | length(out$data.epoch)>0){
    data.astrometry <- out$astrometry
    tmp <- astrometry.kepler(par.opt,bases=bases)
    astro0 <- tmp$barycenter
    rel <- tmp$rel
    Nt <- length(tsim)
#    rel.sim <- astrometry.rel(pars.kep=par.opt,tt=tsim,bases=bases,eta=par.opt['eta'])
    planet <- astrometry.ref(pars.kep=par.opt,tt=tsim,bases=bases,Pmin=1000)$ref
    planet.dra <- planet$dra
    planet.ddec <- planet$ddec
    planet.dpmra <- planet$dpmra
    planet.dpmdec <- planet$dpmdec
#    planet.dra <- rowSums(planet$dra)
#    planet.ddec <- rowSums(planet$ddec)
#    planet.dpmra <- rowSums(planet$dpmra)
#    planet.dpmdec <- rowSums(planet$dpmdec)
    if(grepl('hgca',astro.type)){
        dpmra <- data.astrometry[,'pmra']-astro0[c(1,3),'pmra']
        dpmdec <- data.astrometry[,'pmdec']-astro0[c(2,4),'pmdec']
    }else{
        dpmra <- data.astrometry[,'pmra']-astro0[,'pmra']
        dpmdec <- data.astrometry[,'pmdec']-astro0[,'pmdec']
    }
    dpmras <- c()
    dpmdecs <- c()
    for(j in 1:Nastro){
        ell <- error.ellipse(dpmra[j],dpmdec[j],cov.astro[4:5,4:5,j],percent=68)
        dpmras <- cbind(dpmras,ell[,1])
        dpmdecs <- cbind(dpmdecs,ell[,2])
    }
    if(Popt[imax]>1000){
        plot(dpmra,dpmdec,xlab=expression(Delta*mu[alpha]*' [mas/yr]'),ylab=expression(Delta*mu[delta]*' [mas/yr]'),xlim=range(planet.dpmra,dpmras,dpmra),ylim=range(planet.dpmdec,dpmdecs,dpmdec),pch=20)
        for(j in 1:nastro){
            lines(dpmras[,j],dpmdecs[,j],col='grey')
        }
        lines(planet.dpmra,planet.dpmdec,col='red')
        for(k in 1:nastro){
            i1 <- which.min(abs(tsim-out$astrometry[k,'ref_epoch']))
            lines(c(planet.dpmra[i1],dpmra[k]),c(planet.dpmdec[i1],dpmdec[k]),col=ell.col[k])
        }
        points(0,0,pch='+')
                                        #    if(target=='HD42581') out$Mstar <- 0.58
        msini <- K2msini.full(par.opt[paste0('K',imax)],Popt[imax],par.opt[paste0('e',imax)],Ms=out$Mstar)
        Mp <- msini$mj/sin(par.opt[paste0('Inc',imax)])
        Mo <- (par.opt[paste0('Mo',imax)]%%(2*pi))*180/pi
        inc <- par.opt[paste0('Inc',imax)]*180/pi
        legend('top',legend=as.expression(c(bquote(P==.(round(Popt[imax]/yr2d))~'yr'),bquote(M==.(round(Mp,1))~M[Jup]),bquote(e==.(round(par.opt['e1'],2))),bquote(I==.(round(inc,1))~'deg'),bquote(omega==.(round(par.opt['omega1']*180/pi,1))~'deg'))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.12,x.intersp=0.2)
        if(astrometry!=8){
            legend('top',legend=as.expression(c(bquote(Omega==.(round(par.opt['Omega1']*180/pi,1))~'deg'),bquote(M[0]==.(round(Mo,1))~'deg'),bquote(logJ==.(round(par.opt['logJ.astro'],1))))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.08,x.intersp=0.2)
        }else{
            legend('top',legend=as.expression(c(bquote(omega==.(round(par.opt['omega1']*180/pi,1))~'deg'),bquote(Omega==.(round(par.opt['Omega1']*180/pi,1))~'deg'),bquote(M[0]==.(round(Mo,1))~'deg'),bquote(logJ[hip]==.(round(par.opt['logJ.hip'],1))),bquote(logJ[gaia]==.(round(par.opt['logJ.gaia'],1))))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.08,x.intersp=0.2)
        }
        text(x=dpmra,y=dpmdec,labels=c('H',ng),col='blue',pos=poss,xpd=NA)
    }
                                        #text(x=dpmra,y=dpmdec,labels=c('Hipparcos','Gaia'),offset=c(2,1),pos=c(4,1),col='black')

    if(!grepl('--',astro.type)){
#####plot position ellipse
        if(astrometry>=7){
            dras <- c()
            ddecs <- c()
            if(grepl('hgca',astro.type)){
                dra <- (data.astrometry[,'ra']-astro0[c(1,3),'ra'])*3.6e6*cos(data.astrometry[,'dec']/180*pi)#mas
                ddec <- (data.astrometry[,'dec']-astro0[c(2,4),'dec'])*3.6e6#mas
                if(Popt[imax]>1000){
                    ylim <- range(dra,dra+data.astrometry[,'ra_error'],dra-data.astrometry[,'ra_error'],planet.dra)
                    plot(data.astrometry[,'epoch_ra'],dra,xlab='BJD',ylab=expression(Delta*alpha*'* [mas]'),ylim=ylim)
                    arrows(data.astrometry[,'epoch_ra'],dra-data.astrometry[,'ra_error'],data.astrometry[,'epoch_ra'],dra+data.astrometry[,'ra_error'],length=0.001,angle=0,code=1,col='grey')
                    text(x=data.astrometry[,'epoch_ra'],y=dra,labels=c('H','G'),col='blue',pos=c(4,2))
                    lines(tsim,planet.dra,col='red')

                    ylim <- range(ddec,ddec+data.astrometry[,'dec_error'],ddec-data.astrometry[,'dec_error'],planet.ddec)
                    plot(data.astrometry[,'epoch_dec'],ddec,xlab='BJD',ylab=expression(Delta*delta*' [mas]'),ylim=ylim)
                    arrows(data.astrometry[,'epoch_dec'],ddec-data.astrometry[,'dec_error'],data.astrometry[,'epoch_dec'],ddec+data.astrometry[,'dec_error'],length=0.001,angle=0,code=1,col='grey')
                    text(x=data.astrometry[,'epoch_dec'],y=ddec,labels=c('H','G'),col='blue',pos=c(4,2))
                    lines(tsim,planet.ddec,col='red')
                }
            }else{
                dra <- (data.astrometry[,'ra']-astro0[,'ra'])*3.6e6*cos(data.astrometry[,'dec']/180*pi)#mas
                ddec <- (data.astrometry[,'dec']-astro0[,'dec'])*3.6e6#mas
                for(j in 1:nastro){
                    ell <- error.ellipse(dra[j],ddec[j],cov.astro[1:2,1:2,j],percent=68)
                    dras <- cbind(dras,ell[,1])
                    ddecs <- cbind(ddecs,ell[,2])
                }
                if(Popt[imax]>1000){
                    plot(dra,ddec,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta~delta~'[mas]'),xlim=rev(range(planet.dra,dras)),ylim=range(planet.ddec,ddecs),pch=20)
                    points(0,0,pch='+')
                    for(j in 1:nastro){
                        lines(dras[,j],ddecs[,j],col='grey')
                    }
                    lines(planet.dra,planet.ddec,col='red')
                    for(k in 1:nastro){
                        i1 <- which.min(abs(tsim-out$astrometry[k,'ref_epoch']))
                        lines(c(planet.dra[i1],dra[k]),c(planet.ddec[i1],ddec[k]),col=ell.col[k])
                    }
                    msini <- K2msini.full(par.opt[paste0('K',imax)],Popt[imax],par.opt[paste0('e',imax)],Ms=out$Mstar)
                    Mp <- msini$mj/sin(par.opt[paste0('Inc',imax)])
                    Mo <- (par.opt[paste0('Mo',imax)]%%(2*pi))*180/pi
                    inc <- par.opt[paste0('Inc',imax)]*180/pi
                    legend('top',legend=as.expression(c(bquote(P==.(round(Popt[imax]/yr2d))~'yr'),bquote(M==.(round(Mp,1))~M[Jup]),bquote(e==.(round(par.opt['e1'],2))),bquote(I==.(round(inc,1))~'deg'),bquote(omega==.(round(par.opt['omega1']*180/pi,1))~'deg'))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.12,x.intersp=0.1)
                    if(astrometry!=8){
                        legend('top',legend=as.expression(c(bquote(Omega==.(round(par.opt['Omega1']*180/pi,1))~'deg'),bquote(M[0]==.(round(Mo,1))~'deg'),bquote(logJ==.(round(par.opt['logJ.astro'],1))))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.08,x.intersp=0.1)
                    }else{
                        legend('top',legend=as.expression(c(bquote(omega==.(round(par.opt['omega1']*180/pi,1))~'deg'),bquote(Omega==.(round(par.opt['Omega1']*180/pi,1))~'deg'),bquote(M[0]==.(round(Mo,1))~'deg'),bquote(logJ[hip]==.(round(par.opt['logJ.hip'],1))),bquote(logJ[gaia]==.(round(par.opt['logJ.gaia'],1))))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.08,x.intersp=0.1)
                    }
                                        #    text(x=dra,y=ddec,labels=c('Hipparcos','Gaia'),offset=c(2,1),pos=c(4,1),col='black')
                    text(x=dra,y=ddec,labels=c('H',ng),col='blue',pos=poss)
                }
            }
        }
    }
}
###add binary motion
nc <- 1
if(any(names(out)=='rel')){
    source('../../pexo/code/timing_function.R')
    source('../../pexo/code/general_function.R')
    source('../../pexo/code/sofa_function.R')
                                        #    ind <- sort(Popt,index.return=TRUE)$ix
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
    planet <- astrometry.rel(par.opt,tt=tsim)
    range.raS <- range(unlist(planet$rel$dra))
    range.decS <- range(unlist(planet$rel$ddec))
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
                tP <- out$rel[[n]][[i]][,1]
                draP <- out$rel[[n]][[i]][,'dra']
                ddecP <- out$rel[[n]][[i]][,'ddec']
                edraP <- out$rel[[n]][[i]][,'era']
                eddecP <- out$rel[[n]][[i]][,'edec']
                points(draP,ddecP,col=cc[nc],pch=20,cex=0.8)
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
                arrows(draS[inds],ddecS[inds],draP,ddecP,length=0.001,angle=0,code=1,col='grey')
                arrows(draP-edraP,ddecP,draP+edraP,ddecP,length=0.01,angle=90,code=3,col=cc[nc])
                arrows(draP,ddecP-eddecP,draP,ddecP+eddecP,length=0.01,angle=90,code=3,col=cc[nc])
                nc <- nc+1
                cs <- c(cs,cc[nc])
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
            msini <- K2msini.full(K,popt,e,Ms=out$Mstar)
            Mp <- msini$mj/sin(par.opt[paste0('Inc',j)])
            legend('top',legend=as.expression(c(bquote(P==.(round(popt/yr2d))~'yr'),bquote(M==.(round(Mp,1))~M[Jup]),bquote(e==.(round(e,2))),bquote(I==.(round(inc,1))~'deg'))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.12,x.intersp=0.2)
            legend('top',legend=as.expression(c(bquote(omega==.(round(omega,1))~'deg'),bquote(Omega==.(round(Omega,1))~'deg'),bquote(M[0]==.(round(Mo,1))~'deg'),bquote(logJ[rel]==.(round(par.opt['logJ.rel'],1))))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.08,x.intersp=0.2)
        }
    }
#dev.off()
}

####add reflex motion
if(length(out$data.ref)>0){
    for(i in out$ins.ref){
        inds <- out$ind.ref[[i]]
        dra <- out$data.ref[inds,'dra']
        era <- out$data.ref[inds,'era']
        ddec <- out$data.ref[inds,'ddec']
        edec <- out$data.ref[inds,'edec']
        ref.sim <- astrometry.ref(par.opt,tt=tsim)$ref
        plot(dra,ddec,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),xlim=rev(range(out$data.ref[inds,'dra'],ref.sim[,'dra'])),ylim=range(out$data.ref[,'ddec'],ref.sim[,'ddec']),main=i,pch=20)
        arrows(dra-era,ddec,dra+era,ddec,length=0.001,angle=0,code=1,col='blue')
        arrows(dra,ddec-edec,dra,ddec+edec,length=0.001,angle=0,code=1,col='blue')
        ref.sim <- astrometry.ref(par.opt,tt=tsim)$ref
        lines(ref.sim[,'dra'],ref.sim[,'ddec'],col='red')
        ii <- sapply(out$data.ref[inds,1],function(t) which.min(abs(tsim-t)))
        segments(ref.sim[ii,'dra'],ref.sim[ii,'ddec'],out$data.ref[inds,'dra'],out$data.ref[inds,'ddec'],col='grey')
    }
}

reflex.sim <- astrometry.epoch(par.opt,tt=tsim)$epoch
####add fit to epoch data
if(length(out$data.epoch)>0){
    ins <- out$ins.epoch
    astro <- astrometry.kepler(pars.kep=par.opt,bases=bases)
    astro.sim <- astrometry.kepler(par.opt,tt=tsim)
    epoch <- astro$epoch
    epoch.sim <- astro.sim$epoch
    for(i in ins){
        pred.res <- loglikelihood(par.opt,prediction=TRUE)
        par.opt1 <- par.opt
        par.opt1[grepl('^K\\d',names(par.opt1))] <- 1e-6
###The following approach may not be good to show the residual relatied to long period companion because the stellar reflex motion is subtraced from the Gaia epoch and then propagate the Gaia epoch to Hipparcos epoch which will result in huge offset if without reflex motion
                                        #        pred.sig <- loglikelihood(par.opt1,prediction=TRUE)
        tepoch <- out$data.epoch[[i]][,'BJD']
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
            dra <- pred.res$res[[paste0('epoch_',i)]]*cpsi+epoch[[i]][,'dra']
            ddec <- pred.res$res[[paste0('epoch_',i)]]*spsi+epoch[[i]][,'ddec']
        }else if(i=='hip2'){
            cpsi <- out$data.epoch[[i]][,'CPSI']
            spsi <- out$data.epoch[[i]][,'SPSI']
            absci <- out$data.epoch[[i]][,'RES']
            eabsci <- out$data.epoch[[i]][,'SRES']
            dra0 <- absci*cpsi
            era <- era0 <- eabsci*cpsi
            ddec0 <- absci*spsi
            edec <- edec0 <- eabsci*spsi
            dra <- pred.res$res[[paste0('epoch_',i)]]*cpsi+epoch[[i]][,'dra']
            ddec <- pred.res$res[[paste0('epoch_',i)]]*spsi+epoch[[i]][,'ddec']
        }else{
            dra <- out$data.epoch[[i]][,'dra']
            ddec <- out$data.epoch[[i]][,'ddec']
            era <- out$data.epoch[[i]][,'era']
            edec <- out$data.epoch[[i]][,'edec']
        }
        dra.sim <- epoch.sim[[i]][,'dra']
        ddec.sim <- epoch.sim[[i]][,'ddec']
        dra.pred <- epoch[[i]][,'dra']
        ddec.pred <- epoch[[i]][,'ddec']

###full fit
        plot(dra,ddec,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),xlim=rev(range(dra,dra.sim[i2])),ylim=range(ddec,ddec.sim[i2]),main=paste(i,'full'),pch=20)
        points(0,0,pch='+')
        ncolor <- 20
        cc <- paletteer_c(palette = "viridis::inferno", n = ncolor)
        jj <- as.factor( as.numeric( cut((tepoch-tepoch[1])%%Popt[imax], ncolor)))
        if(grepl('hip',i)){
            segments(dra-era,ddec-edec,dra+era,ddec+edec,col=tcol('grey',50))
        }else{
            arrows(dra-era,ddec,dra+era,ddec,length=0.001,angle=0,code=1,col='grey')
            arrows(dra,ddec-edec,dra,ddec+edec,length=0.001,angle=0,code=1,col='grey')
        }
        ii <- sapply(tepoch,function(t) which.min(abs(tsim-t)))
        segments(dra.sim[ii],ddec.sim[ii],dra,ddec,col=tcol(cc[jj],50))
        lines(dra.sim,ddec.sim,col='red')
        points(dra.sim[ii],ddec.sim[ii],col=cc[jj],pch=20,cex=0.5)
###add Gaia and Hipparcos data
        if(Popt[imax]>1000 & grepl('hip',i)){
            dra.ref <- (data.astrometry[,'ra']-astro0[,'ra'])*3.6e6*cos(data.astrometry[,'dec']/180*pi)#mas
            ddec.ref <- (data.astrometry[,'dec']-astro0[,'dec'])*3.6e6#mas
            k1 <- which.max(data.astrometry[,1])
#            k1 <- 1:nrow(data.astrometry)
            points(dra.ref[k1],ddec.ref[k1],pch=20,col='black')
            j1 <- sapply(k1,function(i) which.min(abs(tsim-data.astrometry[i,1])))
            for(i1 in k1){
                ell <- error.ellipse(dra.ref[i1],ddec.ref[i1],cov.astro[1:2,1:2,i1],percent=68)
                lines(ell[i1,1],ell[i1,2],col='grey')
            }
#            cat('dra offset:',epoch.sim[[i]][j1,'dra']-dra.ref[k1],'\n')
#            cat('ddec offset:',epoch.sim[[i]][j1,'ddec']-ddec.ref[k1],'\n')
            segments(epoch.sim[[i]][j1,'dra'],epoch.sim[[i]][j1,'ddec'],x1=dra.ref[k1],y1=ddec.ref[k1],col=ell.col[k1])
#            lines(c(planet.dra[i1],dra.ref[1]),c(planet.ddec[i1],ddec.ref[1]),col=ell.col[1])
#            lines(c(planet.dra[i2],dra.ref[2]),c(planet.ddec[i2],ddec.ref[2]),col=ell.col[2])
            if(grepl('hip',i)) text(x=dra.ref[k1],y=ddec.ref[k1],labels=c('H','G')[k1],col='blue',pos=4)
        }
        msini <- K2msini.full(par.opt[paste0('K',imax)],Popt[imax],par.opt[paste0('e',imax)],Ms=out$Mstar)
        Mp <- msini$mj/sin(par.opt[paste0('Inc',imax)])
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
            reflex <- astrometry.epoch(par.opt,tt=tepoch)$epoch
            dra.lin.sim <- dra.sim-reflex.sim[,'dra']
            ddec.lin.sim <- ddec.sim-reflex.sim[,'ddec']
            dra.lin <- dra-reflex[,'dra']
            ddec.lin <- ddec-reflex[,'ddec']
#out$astrometry[imax,'pmdec']*dt.sim+out$astrometry[imax,'parallax']*out$pf.sim[,'dec']+par.opt[paste0('ddec_',i)]
#out$astrometry[imax,'pmra']*dt.sim+out$astrometry[imax,'parallax']*out$pf.sim[,'ra']+par.opt[paste0('dra_',i)]
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
            dra0 <- astro$barycenter[imax,'pmra']*dt+astro$barycenter[imax,'parallax']*out$pf[out$ind.epoch[[i]],'ra']+par.opt[paste0('dra_',i)]
            ddec0 <- out$astrometry[imax,'pmdec']*dt+out$astrometry[imax,'parallax']*out$pf[out$ind.epoch[[i]],'dec']+par.opt[paste0('ddec_',i)]
            dra.nl <- dra-dra0
            ddec.nl <- ddec-ddec0
            dra.nl.sim <- reflex.sim[,'dra']
            ddec.nl.sim <- reflex.sim[,'ddec']
            plot(dra.nl,ddec.nl,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),xlim=rev(range(dra.nl,dra.nl.sim[i2])),ylim=range(ddec.nl,ddec.nl.sim[i2]),main=paste(i,'orbital'),pch=20)
            arrows(dra.nl-era,ddec.nl,dra.nl+era,ddec.nl,length=0.001,angle=0,code=1,col='grey')
            arrows(dra.nl,ddec.nl-edec,dra.nl,ddec.nl+edec,length=0.001,angle=0,code=1,col='grey')
            points(0,0,pch='+')
            ncolor <- 20
            cc <- paletteer_c(palette = "viridis::inferno", n = ncolor)
            jj <- as.factor( as.numeric( cut((tepoch-tepoch[1])%%Popt[imax], ncolor)))
#            segments(dra.nl-era,ddec.nl-edec,dra.nl+era,ddec.nl+edec,col=tcol('grey',50))
            ii <- sapply(out$data.epoch[[i]][,'BJD'],function(t) which.min(abs(tsim-t)))
#            segments(dra.sim[ii],ddec.sim[ii],dra,ddec,col=tcol(cc[jj],50))
            lines(dra.nl.sim,ddec.nl.sim,col='red')
            points(dra.nl.sim[ii],ddec.nl.sim[ii],col=cc[jj],pch=20,cex=0.5)
        }
    }
}

###add GOST fit
if(length(out$gost)>0){
    reflex.gost <- astrometry.epoch(par.opt,tt=out$gost[,'BJD'],bases=bases)$epoch
    kep <- astrometry.kepler(par.opt,bases=bases)
###observed position and proper motion relative to the predicted barycentric position and proper motion
    dra.obs <- ((out$astrometry[out$astro.index,'ra']-kep$barycenter[out$astro.index,'ra']))*cos(out$astrometry[out$astro.index,'dec']/180*pi)*3.6e6
    ddec.obs <- ((out$astrometry[out$astro.index,'dec']-kep$barycenter[out$astro.index,'dec']))*3.6e6
    dpmra.obs <- out$astrometry[out$astro.index,'pmra']-kep$barycenter[out$astro.index,'pmra']
    dpmdec.obs <- out$astrometry[out$astro.index,'pmdec']-kep$barycenter[out$astro.index,'pmdec']
    dplx.obs <- out$astrometry[out$astro.index,'parallax']-kep$barycenter[out$astro.index,'parallax']
###the position and proper motion reconstructed by a linear fit to gost-based observations generated by accounting for reflex motion
    dra.model <- dra.obs+kep$gdr[,1]
    ddec.model <- ddec.obs+kep$gdr[,2]
    dplx.model <- dplx.obs+kep$gdr[,3]
    dpmra.model <- dpmra.obs+kep$gdr[,4]
    dpmdec.model <- dpmdec.obs+kep$gdr[,5]

    dra.obs.range <- range(sapply(1:length(dra.obs), function(k) dra.obs[k]+dpmra.obs[k]*c(min(out$dtg[[k]]),max(out$dtg[[k]]))))
    ddec.obs.range <- range(sapply(1:length(dra.obs), function(k) ddec.obs[k]+dpmdec.obs[k]*c(min(out$dtg[[k]]),max(out$dtg[[k]]))))
    dra.model.range <- range(sapply(1:length(dra.model), function(k) dra.model[k]+dpmra.model[k]*c(min(out$dtg[[k]]),max(out$dtg[[k]]))))
    ddec.model.range <- range(sapply(1:length(dra.model), function(k) ddec.model[k]+dpmdec.model[k]*c(min(out$dtg[[k]]),max(out$dtg[[k]]))))
    dra.range <- range(dra.obs.range,dra.model.range)
    ddec.range <- range(ddec.obs.range,ddec.model.range)
####plot
    plot(reflex.sim[,'dra'],reflex.sim[,'ddec'],xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),type='l',col='grey',xlim=rev(range(reflex.sim[,'dra'],dra.range)),ylim=range(reflex.sim[,'ddec'],ddec.range))
    points(0,0,pch='+')
    jj <- as.factor(((out$gost[,1]-out$gost[1,1])%%max(Popt))/max(Popt))
    cc <- paletteer_c(palette = "viridis::inferno", n = length(jj))
    points(reflex.gost[,'dra'],reflex.gost[,'ddec'],col=cc[jj],pch=20)
    if(length(out$gdrs)==0) out$gdrs <- c('GDR2','GDR3')
    phases <- seq(0,1,by=0.2)
    cl <- paletteer_c(palette = "viridis::inferno", n = length(phases))
    legend('top',legend=phases,pch=20,col=cl[as.factor(phases)],horiz=TRUE,bty='n',inset=c(0,-0.25),xpd=NA,title='Phase')
    ref.epochs <- c()
    cols <- c('darkgreen','blue','red')
#    dra.model.comb <- c()
    for(k in 1:length(out$astro.index)){
        ttsim <- seq(min(out$dtg[[k]]),max(out$dtg[[k]]),by=0.1)#year
        dras.obs <- dra.obs[k]+dpmra.obs[k]*ttsim
        ddecs.obs <- ddec.obs[k]+dpmdec.obs[k]*ttsim
        lines(dras.obs,ddecs.obs,col=cols[k],lty=2)
        points(dra.obs[k],ddec.obs[k],col=cols[k],pch=1)

        dras.model <- dra.model[k]+dpmra.model[k]*ttsim
        ddecs.model <- ddec.model[k]+dpmdec.model[k]*ttsim
        lines(dras.model,ddecs.model,col=cols[k],lty=1)
        points(dra.model[k],ddec.model[k],col=cols[k],pch=20)
    }
#    legend('topright',legend=out$gdrs,col=cols,lty=1,cex=0.8,bty='n',inset=c(-0.23,0),xpd=NA)
    legend('top',legend=out$gdrs,col=cols,lty=1,cex=0.8,bty='n',inset=c(0,-0.1),xpd=NA,horiz=TRUE)
    trefs <- out$astrometry[out$astro.index,'ref_epoch']
    obs <- cbind(dra.obs,ddec.obs,dplx.obs,dpmra.obs,dpmdec.obs)
    eobs <- out$astrometry[out$astro.index,c('ra_error','dec_error','parallax_error','pmra_error','pmdec_error')]
    model <- cbind(dra.model,ddec.model,dplx.model,dpmra.model,dpmdec.model)
    ylabs <- c(expression(Delta*alpha*'* [mas]'),expression(Delta*delta*' [mas]'),'parallax [mas]',expression(mu[alpha]*' [mas/yr]'),expression(mu[delta]*' [mas/yr]'))
    for(j in 1:ncol(obs)){
        plot(trefs,obs[,j],xlab='BJD',ylab=as.expression(ylabs[j]),ylim=range(obs[,j],model[,j]),col=cc[1:nast])
        arrows(trefs,obs[,j]-eobs[,j],trefs,obs[,j]+eobs[,j],length=0.001,angle=0,code=1,col='grey')
        points(trefs,model[,j],col='red')
    }
}
if(plotf){
dev.off()
}
