#library("viridis")
directed.circle <- function(xlim,ylim,clockwise=TRUE,lty=1,col='grey'){
    dx <- abs(diff(xlim))
    dy <- diff(ylim)
    xc <- max(xlim)-0.06*dx
    yc <- min(ylim)+0.06*dy
    a <- 0.05*dx
    b <- 0.05*dy
    theta <- seq(0,2*pi,by=0.01)
    xs <- xc+a*cos(theta)
    ys <- yc+b*sin(theta)
    theta0 <- c(0,pi/2,pi,3*pi/2)
    xa0 <- xc+a*cos(theta0)
    ya0 <- yc+b*sin(theta0)
    dxa <- 0.05*dx*c(0,-1,0,1)
    dya <- 0.05*dy*c(1,0,-1,0)
    if(!clockwise){
        xa1 <- xa0-dxa
        ya1 <- ya0-dya
    }else{
        xa1 <- xa0+dxa
        ya1 <- ya0+dya
    }
    lines(xs,ys,col=col,lty=lty)
    arrows(xa0,ya0,xa1,ya1,length=0.03,code=2,col=col)
}

k5 <- 1
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
#imax <- which.max(Popt)
###add binary motion
nc <- 1
tmp <- astrometry.kepler(par.opt,bases=bases)
astro0 <- tmp$barycenter
reflex.sim <- astrometry.epoch(par.opt,tt=tsim)$epoch
if(target=='UCAC4569-026385'){
    eta <- 0
}else{
    eta <- eta0
}
rel.sim <- astrometry.rel(par.opt,tt=tsim,eta=eta)$rel
ins <- out$ins.epoch
astro <- astrometry.kepler(pars.kep=par.opt,bases=bases)
epoch <- astro$epoch
if(is.null(astro$epoch)){
    epoch.sim <- astrometry.epoch(par.opt,tt=tsim,bases=bases)$epoch
}else{
    epoch.sim <- astrometry.kepler(pars.kep=par.opt,bases=bases,tt=tsim)$epoch
}

#astro.sim <- astrometry.kepler(par.opt,tt=tsim)
if(par.opt[paste0('Inc',j5)]<pi/2){
    clockwise <- FALSE
}else{
    clockwise <- TRUE
}
if(!is.null(ins)){
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
            ts <- out$data.epoch[[i]][,'BJD']
            ts1 <- c(ts[1],ts)
            dt <- diff(ts1)
            index  <- c()
            ii <- 1
            for(j in 1:length(dt)){
                if(dt[j]<0.1){
                    index <- c(index,ii)
                }else{
                    ii <- ii+1
                    index <- c(index,ii)
                }
            }
            w <- 1/eabsci^2
            tbin <- sapply(unique(index),function(i) sum(ts[index==i]*w[index==i])/sum(w[index==i]))
            dra.bin <- sapply(unique(index),function(i) sum(dra[index==i]*w[index==i])/sum(w[index==i]))
            ddec.bin <- sapply(unique(index),function(i) sum(ddec[index==i]*w[index==i])/sum(w[index==i]))
            eabsci.bin <- sapply(unique(index),function(i) sd(absci[index==i]))
            psi <- atan2(spsi,cpsi)
            psi.bin <- sapply(unique(index),function(i) sum(psi[index==i]*w[index==i])/sum(w[index==i]))
            nn <- sapply(unique(index),function(i) length(which(index==i)))
            ii <- which(nn==1)
            if(length(ii)>0){
                tbin <- tbin[-ii]
                dra.bin <- dra.bin[-ii]
                ddec.bin <- ddec.bin[-ii]
                eabsci.bin <- eabsci.bin[-ii]
                psi.bin <- psi.bin[-ii]
            }
            cpsi.bin <- cos(psi.bin)
            spsi.bin <- sin(psi.bin)
            era.bin <- eabsci.bin*cpsi.bin
            edec.bin <- eabsci.bin*spsi.bin
                                        #        era.bin <- sapply(unique(index),function(i) sd(dra[index==i]))
                                        #        edec.bin <- sapply(unique(index),function(i) sd(ddec[index==i]))
            epoch.bin <- astrometry.epoch(par.opt,tt=tbin)$epoch
            res.bin <- sapply(unique(index),function(j) sum(pred.res$res[[paste0('epoch_',i)]][index==j]*w[index==j])/sum(w[index==j]))
                                        #        dra.bin <- res.bin*cpsi.bin+epoch.bin[,'dra']
                                        #        ddec.bin <- res.bin*spsi.bin+epoch.bin[,'ddec']
                                        #        ddec.bin <- sapply(unique(index),function(i) sum(ddec[index==i]*w[index==i])/sum(w[index==i]))
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
        k5 <- k5+1
        screen(ss0[k5])
        par(mar=mm,mgp=mgp)
        Mo <- par.opt[paste0('Mo',j5)]
        Tp <- tmin-(Mo/(2*pi))*Popt
        Omega <- par.opt[paste0('Omega',j5)]
        theta.sim <- atan2(dra.sim,ddec.sim)
        ind <- which.min(abs(theta.sim-Omega))
        Tascend <- tsim[ind]
        par(xpd=NA)
        xlim <- rev(range(dra,dra.sim))
        ylim <- range(ddec,ddec.sim)
                                        #    plot(dra,ddec,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),xlim=rev(range(dra,dra.sim[i2])),ylim=range(ddec,ddec.sim[i2]),pch=ell.pch,col='grey',main='Hipparcos IAD')
        plot(dra,ddec,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),xlim=xlim,ylim=ylim,pch=ell.pch,col='grey',main='',xaxt='n',yaxt='n')
        par(xpd=FALSE)
        magaxis(1:2,cex=size)
        magaxis(3:4,cex=size,labels=FALSE)
####dashed lines
        for(ii in unique(index)){
            jj <- which(index==ii)
            for(j1 in 2:length(jj)){
                                        #            segments(dra[jj[j1-1]],ddec[jj[j1-1]],dra[jj[j1]],ddec[jj[j1]],col='grey',lwd=0.5,lty=2)
            }
        }
        ncolor <- 100

        cc <- paletteer_c(palette = "viridis::inferno", n = length(tbin))
        ##  i3 <- as.factor( as.numeric( cut((tbin-Tascend)%%Popt[j5], length(tsim))))
        i3 <- as.factor( as.numeric( cut((tbin-min(tbin))%%Popt[j5], length(tbin)) ) )
        points(dra.bin,ddec.bin,pch=20,col=cc[i3])
        segments(dra.bin-era.bin,ddec.bin-edec.bin,dra.bin+era.bin,ddec.bin+edec.bin,col=cc[i3],lwd=0.5)#tcol(cc[i3],50)
        phases <- seq(0,1,by=0.2)
        cl <- paletteer_c(palette = "viridis::inferno", n = length(phases))
        legend('top',legend=phases,pch=20,col=cl,horiz=TRUE,bty='n',inset=c(0,-0.18),xpd=NA,title='Phase')
                                        #    legend('top',horiz=TRUE,pch=c(ell.pch,NA,NA),lty=c(NA,1,1),col=c('grey','red','black'),legend=c('Hipparcos IAD','Best fit','Abscissa error'),inset=c(0,-0.15),xpd=NA,bty='n')
                                        #    points(dra.sim,ddec.sim,pch=20,col='darkgrey',cex=0.1)
        lines(dra.sim,ddec.sim,pch=20,col=fit.opt,lwd=2)
        ind <- which.min(abs((tsim-Tascend)%%Popt[j5]))
        segments(0,0,dra.sim[ind],ddec.sim[ind],col='grey',lwd=2)
        points(0,0,pch='+')
        directed.circle(xlim,ylim,clockwise=clockwise,lty=1)
                                        #    points(dra.sim[ind],ddec.sim[ind],pch=15,col='grey')
                                        #    phase <- (tepoch-Tp)%%Popt[j5]
        if(grepl('hip',i)){
                                        #        segments(dra-era,ddec-edec,dra+era,ddec+edec,col='grey',lwd=0.5)
        }else{
            arrows(dra-era,ddec,dra+era,ddec,length=0.001,angle=0,code=1,col='grey')
            arrows(dra,ddec-edec,dra,ddec+edec,length=0.001,angle=0,code=1,col='grey')
        }
        ii <- sapply(tepoch,function(t) which.min(abs(tsim-t)))
###add Gaia and Hipparcos data
        if(Popt[j5]>1000 & grepl('hip',i) & FALSE){
            dra.ref <- (out$astrometry[,'ra']-astro0[,'ra'])*3.6e6*cos(out$astrometry[,'dec']/180*pi)#mas
            ddec.ref <- (out$astrometry[,'dec']-astro0[,'dec'])*3.6e6#mas
            k1 <- which.max(out$astrometry[,1])
                                        #            k1 <- 1:nrow(out$astrometry)
            points(dra.ref[k1],ddec.ref[k1],pch=20,col='black')
            j1 <- sapply(k1,function(i) which.min(abs(tsim-out$astrometry[i,1])))
            for(kk in k1){
                ell <- error.ellipse(dra.ref[kk],ddec.ref[kk],cov.astro[1:2,1:2,kk],percent=68)
                lines(ell[kk,1],ell[kk,2],col='grey')
            }
            segments(epoch.sim[[i]][j1,'dra'],epoch.sim[[i]][j1,'ddec'],x1=dra.ref[k1],y1=ddec.ref[k1],col=ell.col[k1])
            if(grepl('hip',i)) text(x=dra.ref[k1],y=ddec.ref[k1],labels=c('H','G')[k1],col='blue',pos=4)
        }
                                        #    msini <- K2msini.full(par.opt[paste0('K',j5)],Popt[j5],par.opt[paste0('e',j5)],Ms=out$Mstar)
                                        #    Mp <- msini$mj/sin(par.opt[paste0('Inc',j5)])
        Mp <- k2m(par.opt[paste0('K',j5)],Popt[j5],par.opt[paste0('e',j5)],Ms=out$Mstar,Inc=par.opt[paste0('Inc',j5)])$mj
        Mo <- (par.opt[paste0('Mo',j5)]%%(2*pi))*180/pi
        inc <- par.opt[paste0('Inc',j5)]*180/pi
                                        #    legend('top',legend=as.expression(c(bquote(P==.(round(Popt[j5]/yr2d))~'yr'),bquote(M==.(round(Mp,1))~M[Jup]),bquote(e==.(round(par.opt['e1'],2))),bquote(I==.(round(inc,1))~'deg'),bquote(omega==.(round(par.opt['omega1']*180/pi,1))~'deg'))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.12,x.intersp=0.15,cex=0.9)
                                        #    legend('top',legend=as.expression(c(bquote(omega==.(round(par.opt['omega1']*180/pi,1))~'deg'),bquote(Omega==.(round(par.opt['Omega1']*180/pi,1))~'deg'),bquote(M[0]==.(round(Mo,1))~'deg'),bquote(logJ[hip]==.(round(par.opt[grep('J_hip',names(par.opt))],2))),bquote(logJ[gaia]==.(round(par.opt['logJ_gaia'],2))))),bty='n',horiz=TRUE,xpd=TRUE,inset=-0.08,x.intersp=0.15,cex=0.9)
    }
}
###add GOST fit
reflex.gost <- astrometry.epoch(par.opt,tt=out$gost[,'BJD'],bases=bases)$epoch
kep <- astrometry.kepler(par.opt,bases=bases)
###observed position and proper motion relative to the predicted barycentric position and proper motion
dra.obs <- ((out$astrometry[out$astro.index,'ra']-kep$barycenter[out$astro.index,'ra']))*cos(out$astrometry[out$astro.index,'dec']/180*pi)*3.6e6
ddec.obs <- ((out$astrometry[out$astro.index,'dec']-kep$barycenter[out$astro.index,'dec']))*3.6e6
dpmra.obs <- out$astrometry[out$astro.index,'pmra']-kep$barycenter[out$astro.index,'pmra']
dpmdec.obs <- out$astrometry[out$astro.index,'pmdec']-kep$barycenter[out$astro.index,'pmdec']
dplx.obs <- out$astrometry[out$astro.index,'parallax']-kep$barycenter[out$astro.index,'parallax']
era.obs <- out$astrometry[out$astro.index,'ra_error']
edec.obs <- out$astrometry[out$astro.index,'dec_error']
eplx.obs <- out$astrometry[out$astro.index,'parallax_error']
epmra.obs <- out$astrometry[out$astro.index,'pmra_error']
epmdec.obs <- out$astrometry[out$astro.index,'pmdec_error']

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
####plot
##contour of observed position as a function of time given proper motion uncertainty
conPM <- function(dra,ddec,dpmra,dpmdec,era,edec,epmra,epmdec,dt){
##https://byjus.com/jee/equation-of-tangent-to-ellipse/#:~:text=A%20line%20which%20intersects%20the%20ellipse%20at%20a,b%202%5D%20represent%20the%20tangents%20to%20the%20ellipse.
    a <- sqrt(era^2+(epmra*dt)^2)
    b <- sqrt(edec^2+(epmdec*dt)^2)
    m <- dpmdec/dpmra
    c <- sqrt(a^2*m^2+b^2)
    x <- a^2*m/c
    y <- -b^2/c
    x0 <- dra+dpmra*dt
    y0 <- ddec+dpmdec*dt
    x1 <- x0+x
    y1 <- y0+y
    x2 <- x0-x
    y2 <- y0-y
    phi <- seq(0,2*pi,by=0.01)
if(FALSE){
    x3 <- dra+min(a)*cos(phi)
    y3 <- ddec+min(b)*sin(phi)
    x4 <- x0[3]+a[3]*cos(phi)
    y4 <- y0[3]+b[3]*sin(phi)
#    save(list=ls(all=TRUE),file='test0.Robj')
    points(x0[3],y0[3],col='cyan',pch=20)
    lines(x3,y3,col='black')
    lines(x4,y4,col='cyan')
    points(x1,y1,col='red')
    points(x2,y2,col='orange')
}
    return(list(upper=cbind(x1,y1),lower=cbind(x2,y2)))
}
xl <- xu <- yl <- yu <- list()
for(k in 1:length(out$astro.index)){
    ttsim <- seq(min(out$dtg[[k]]),max(out$dtg[[k]]),by=0.1)#year
    con <- conPM(dra.obs[k],ddec.obs[k],dpmra.obs[k],dpmdec.obs[k],era.obs[k],edec.obs[k],epmra.obs[k],epmdec.obs[k],dt=ttsim)
    xl[[k]] <- con$lower[,1]
    xu[[k]] <- con$upper[,1]
    yl[[k]] <- con$lower[,2]
    yu[[k]] <- con$upper[,2]
}
k5 <- k5+1
screen(ss0[k5])
#mm[4] <- 1*margin
par(mar=mm,xpd=NA,mgp=mgp)
xlim <- rev(range(xl,xu,reflex.sim[,'dra'],dra.range))
ylim <- range(yl,yu,reflex.sim[,'ddec'],ddec.range)
plot(reflex.sim[,'dra'],reflex.sim[,'ddec'],xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),type='l',col=fit.opt,xlim=xlim,ylim=ylim,main='',lwd=2,xaxt='n',yaxt='n')#Gaia GOST
par(xpd=FALSE)
magaxis(1:2,cex=size)
magaxis(3:4,cex=size,labels=FALSE)
Mo <- par.opt[paste0('Mo',j5)]
Tp <- tmin-(Mo/(2*pi))*Popt
Omega <- par.opt[paste0('Omega',j5)]
epoch.sim <- astrometry.epoch(par.opt,tt=tsim,bases=bases)$epoch
dra.sim <- epoch.sim[,'dra']
ddec.sim <- epoch.sim[,'ddec']
theta.sim <- atan2(dra.sim,ddec.sim)
ind <- which.min(abs(theta.sim-Omega))
Tascend <- tsim[ind]
ind <- which.min(abs((tsim-Tascend)%%Popt[j5]))
segments(0,0,reflex.sim[ind,'dra'],reflex.sim[ind,'ddec'],col='grey',lwd=2)
points(0,0,pch='+')
directed.circle(xlim,ylim,clockwise=clockwise,lty=1)

#jj <- as.factor(((out$gost[,1]-out$gost[1,1])%%max(Popt))/max(Popt))
#jj <- as.factor(((out$gost[,1]-Tascend)%%Popt[j5])/Popt[j5])
#cc <- paletteer_c(palette = "viridis::inferno", n = length(jj))
cc <- paletteer_c(palette = "viridis::inferno", n = nrow(out$gost))
jj <- as.factor( as.numeric( cut((out$gost[,1]-min(out$gost[,1]))%%Popt[j5], nrow(out$gost))))
#jj <- as.factor( as.numeric( cut((out$gost[,1]-Tascend)%%Popt[j5], length(tsim))))
points(reflex.gost[,'dra'],reflex.gost[,'ddec'],col=cc[jj],pch=ell.pch)
if(length(out$cats)==0) out$cats <- c('GDR2','GDR3')
cols <- c('darkgreen','blue','red')
Ndr <- length(out$cats)
legend('top',legend=out$cats,col=cols[1:Ndr],pch=rep(NA,Ndr),lty=rep(1,Ndr),horiz=TRUE,bty='n',inset=c(0,-0.18),xpd=NA,lwd=8)
legend('top',legend=paste('Fitted',out$cats),col=cols[1:Ndr],pch=rep(20,Ndr),lty=rep(1,Ndr),horiz=TRUE,bty='n',inset=c(0,-0.12),xpd=NA)
ref.epochs <- c()
                                        #    dra.model.comb <- c()
for(k in 1:length(out$astro.index)){
    ttsim <- seq(min(out$dtg[[k]]),max(out$dtg[[k]]),by=0.1)#year
    dras.obs <- dra.obs[k]+dpmra.obs[k]*ttsim
    ddecs.obs <- ddec.obs[k]+dpmdec.obs[k]*ttsim
#    lines(dras.obs,ddecs.obs,col=cols[k],lty=2)
#    points(dra.obs,ddec.obs,col=cols[k],pch=20)
    polygon(c(xl[[k]],rev(xu[[k]])),c(yl[[k]],rev(yu[[k]])),col=tcol(cols[k],50),border=NA)

    dras.model <- dra.model[k]+dpmra.model[k]*ttsim
    ddecs.model <- ddec.model[k]+dpmdec.model[k]*ttsim
    lines(dras.model,ddecs.model,col=cols[k],lty=1,lwd=2)
    points(dra.model[k],ddec.model[k],col=cols[k],pch=20)
}
                                        #    legend('topright',legend=out$cats,col=cols,lty=1,cex=0.8,bty='n',inset=c(-0.23,0),xpd=NA)
#legend('top',legend=out$cats,col=cols,lty=1,cex=0.8,bty='n',inset=c(0,-0.1),xpd=NA,horiz=TRUE)
trefs <- out$astrometry[out$astro.index,'ref_epoch']
obs <- cbind(dra.obs,ddec.obs,dplx.obs,dpmra.obs,dpmdec.obs)
eobs <- out$astrometry[out$astro.index,c('ra_error','dec_error','parallax_error','pmra_error','pmdec_error')]
model <- cbind(dra.model,ddec.model,dplx.model,dpmra.model,dpmdec.model)
ylabs <- c(expression(Delta*alpha*'* [mas]'),expression(Delta*delta*' [mas]'),'parallax [mas]',expression(Delta*mu[alpha]*' [mas/yr]'),expression(Delta*mu[delta]*' [mas/yr]'))

k5 <- k5+1
screen(ss0[k5])
mm[4] <- 1*margin
par(mar=mm,xpd=NA)
jd.ref <- rowSums(time_Cal2JD(cal))
mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
ind <- which.max(mc[,'logpost'])
mc <- rbind(mc[ind,],mc)
jj <- (1:Nsig)[-j5]
if(length(jj)>0) mc[,paste0('K',jj)] <- 0
M0s <- mc[,paste0('Mo',j5)]
ps <- exp(mc[,paste0('per',j5)])
ms <- (M0s+2*pi*((jd.ref-tmin)%%ps)/ps)%%(2*pi)
es <- mc[,paste0('e',j5)]
Es <- kep.mt2(ms,es)
incs <- mc[,paste0('Inc',j5)]
omegas <- mc[,paste0('omega',j5)]
Omegas <- mc[,paste0('Omega',j5)]
Ks <- mc[,paste0('K',j5)]
mstars <- rnorm(nrow(mc),out$Mstar,out$eMstar)
plx <- out$plx
if(star=='UCAC4569-026385'){
    pmfun <- pemfun <- NULL
    if(file.exists('plx_mass.txt')){
        tab <- read.table('plx_mass.txt',header=TRUE)
        pmfun.spec <- approxfun(tab[,'plx'],tab[,'spec_mass'])
        pemfun.spec <- approxfun(tab[,'plx'],0.5*(tab[,'esmass1']+tab[,'esmass2']))
    }
    plx <- out$astrometry[out$iref,'parallax']-mc[,'dplx']
    mstars <- pmfun.spec(plx)
}
if(any(names(par.opt)=='dplx')) plx <- out$plx-par.opt['dplx']
astro <- calc.astro(K=Ks,P=ps,e=es,inc=incs,omega=omegas,Omega=Omegas,E=Es,plx=plx,Mstar=mstars,eta=eta)
mps <- k2m(Ks,ps,es,Ms=mstars,Inc=incs)$ms
#msini <- K2msini.full(Ks,ps,es,Ms=mstars)
#mps <- msini$ms/sin(incs)
xis <- mps/(mstars+mps)
dra.mc <- -astro[,1]/xis
ddec.mc <- -astro[,2]/xis
pa.mc <- atan2(dra.mc,ddec.mc)*180/pi#deg
sep.mc <- sqrt(dra.mc^2+ddec.mc^2)*1e-3#as
cat('Companion position at',paste(unlist(cal),collapse='-'),':dra=',mean(dra.mc)/1e3,'+-',sd(dra.mc)/1e3,'as;ddec=',mean(ddec.mc)/1e3,'+-',sd(ddec.mc)/1e3,'as; or sep=',mean(sep.mc),'+-',sd(sep.mc),'as;pa=',mean(pa.mc),'+-',sd(pa.mc),'deg\n')
getLevel <- function(kk,prob=0.95){
    dx <- diff(kk$x[1:2])
    dy <- diff(kk$y[1:2])
    sz <- sort(kk$z)
    c1 <- cumsum(sz) * dx * dy
    approx(c1, sz, xout = 1 - prob)$y
}

dra.sim <- rel.sim$dra[[paste0('p',j5)]]
ddec.sim <- rel.sim$ddec[[paste0('p',j5)]]
xlim <- rev(range(dra.mc,dra.sim))
ylim <- range(ddec.mc,ddec.sim)
###contour
kk <- kde2d(dra.mc,ddec.mc, n=14)
#kk <- kde2d(dra.mc,ddec.mc, n=20)
#kk <- kde2d(dra.mc,ddec.mc)
#levels <- getLevel(kk,prob=c(0.68,0.95,0.9973))
#levels <- getLevel(kk,prob=c(0.68,0.95))
levels <- getLevel(kk,prob=0.68)
zlim <- c(0.0027,1)#truncate at 3sigma;5sigma: 99.99994 or 1-5sigma=6e-7

par(mar=mm,mgp=mgp)
plot(dra.sim,ddec.sim,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),type='l',col=fit.opt,main=paste0('Companion position on ',paste(unlist(cal),collapse='-')),lwd=2,xlim=xlim,ylim=ylim,xaxt='n',yaxt='n')
par(xpd=FALSE)
magaxis(1:2,cex=size,majorn=3)
magaxis(3:4,cex=size,labels=FALSE,majorn=3)
ind <- which.min(abs((tsim-Tascend)%%Popt[j5]))
segments(0,0,rel.sim$dra$p1[ind],rel.sim$ddec$p1[ind],col='grey',lwd=2)
points(0,0,pch='+')
inds <- sample(1:length(dra.mc),1e5)
#points(dra.mc[inds],ddec.mc[inds],col=tcol('grey',80),pch='.')
nl <- length(levels)
contour(kk, drawlabels=FALSE, nlevels=nl, levels=levels,col='black',add=TRUE,axes=FALSE,lty=nl,lwd=2)
if(FALSE){
    ind <- which(abs(kk$z-levels[1])/kk$z<0.02,arr.ind=TRUE)
    xs <- kk$x[ind[,1]]
    ys <- kk$y[ind[,2]]
    xy <- c(median(xs),median(ys))
    theta <- atan2(ys-median(ys),xs-median(xs))
    ii <- sort(theta,index.return=TRUE)$ix
    polygon(xs[ii],ys[ii],col='red')
}
#filled.contour(kk$x,kk$y,kk$z,color.palette=function(n) hcl.colors(n, "terrain"),levels=sort(levels),add=TRUE,zlim=zlim)
points(dra.mc[1],ddec.mc[1],pch=15,col='black')
directed.circle(xlim,ylim,clockwise=clockwise,lty=1)

###parameter table
if(FALSE){
Ndig <- 1
j <- 1
popt <- exp(par.opt[paste0('per',j)])
K <- par.opt[paste0('K',j)]
e <- par.opt[paste0('e',j)]
inc <- par.opt[paste0('Inc',j)]*180/pi
Mo <- (par.opt[paste0('Mo',j)]%%(2*pi))*180/pi
omega <- par.opt[paste0('omega',j)]*180/pi
Omega <- par.opt[paste0('Omega',j)]*180/pi
msini <- K2msini.full(K,popt,e,Ms=out$Mstar)
Mp <- msini$mj/sin(par.opt[paste0('Inc',j)])
legend('topright',inset=c(-0.7,0),xpd=NA,legend=as.expression(c(bquote(m[c]==.(round(Mp,Ndig))~M[Jup]),bquote(P==.(round(popt/yr2d,Ndig))~'yr'),bquote(e==.(round(e,Ndig))),bquote(I==.(round(inc,Ndig))~'deg'),bquote(omega==.(round(omega,Ndig))~'deg'),bquote(Omega==.(round(Omega,Ndig))~'deg'),bquote(M[0]==.(round(Mo,Ndig))~'deg'),bquote(J[hip]==.(round(exp(par.opt['logJ_hip2']),Ndig))~'mas'),bquote(J[gaia]==.(round(exp(par.opt['logJ_gaia']),Ndig))))),bty='n',title='MAP parameter')
}
