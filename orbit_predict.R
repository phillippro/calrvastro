fname <- 'orbit_predict.pdf'
cat(fname,'\n')
pdf(fname,6,6)
target <- star
dir <- 'orbits/'
#cal <- cbind(2022,7,15)
cal <- cbind(2024,1,1)
jd.ref <- rowSums(time_Cal2JD(cal))
size <- 1.2
#par(mar=c(5,5,1,1),cex=size,cex.lab=size,cex.axis=size)
#par(mar=c(5,5,3,1),cex=size,cex.lab=size,cex.axis=size)
####position at the reference epoch
dras <- c()
ddecs <- c()
#Nsamp <- 1e4
#Nsamp <- 1e3
#Nsamp <- nrow(mc)
Npar <- length(par.opt)
ind.sample <- sample(1:nrow(mc),Nmc)
inds <- c(which.max(mc[,'loglike']),ind.sample)
###MCMC samples of parameters
#pars <- rbind(par.opt,mc[inds,1:Npar])
pars <- mc[ind.sample,1:Npar]
Mstars <- pars[,'Mstar']
mscs <- quantile(Mstars,c(0.16,0.5,0.84))
#Mstars <- Mstars[ind.sample]
#ind.max <- as.integer(which.max(par.opt[grepl('^per',names(par.opt))]))

###settings
nlevel <- 3
cols <- c('red','blue','green','orange','cyan','brown','steelblue','yellow')
#my.cols <- rev(hcl.colors(nlevel))
my.cols <- rev(rainbow(nlevel))
getLevel <- function(kk,prob=0.95){
    dx <- diff(kk$x[1:2])
    dy <- diff(kk$y[1:2])
    sz <- sort(kk$z)
    c1 <- cumsum(sz) * dx * dy
    approx(c1, sz, xout = 1 - prob)$y
}
yrs <- seq(2020,2050,by=5)
#yrs <- c(2023,2025)
#yrs <- 2024
#yrs <- c()
jds <- c(jd.ref,rowSums(time_Yr2jd(yrs)))
#if(length(yrs)>0) jds <- rowSums(time_Yr2jd(yrs))
Pmax <- max(exp(pars[,grepl('per',colnames(pars))]))
#t1 <- min(min(tall),mean(tall)-max(Pmax)/2,jd.ref)
#t2 <- max(max(tall),mean(tall)+max(Pmax)/2,jd.ref)
#out <- e0$out
tall <- out$all$trv.all
t1 <- min(mean(tall)%%24e5-max(Pmax)/2,jd.ref%%24e5)
t2 <- max(mean(tall)%%24e5+max(Pmax)/2,jd.ref%%24e5)
tt <- seq(t1,t2,length.out=1000)
Plow <- 100
ps <- par.opt[grepl('^per',names(par.opt))]
ind.lps <- which(ps>log(Plow))
ind.lps <- ind.lps[sort(ps[ind.lps],index.return=TRUE)$ix]
if(star=='HD42581') ind.lps <- ind.lps[-length(ind.lps)]
####list to save simulation data
ddec.ref.list <- dra.ref.list <- ddec.mc.list <- dra.mc.list <- dra.sim.list <- ddec.sim.list <- dra.list <- ddec.list <- list()

####loop to simulate the orbits
for(k in 1:length(ind.lps)){
    k3 <- as.integer(ind.lps[k])
    P <- exp(pars[,paste0('per',k3)])
    K <- pars[,paste0('K',k3)]
    e <- pars[,paste0('e',k3)]
    omega <- pars[,paste0('omega',k3)]
    Omega <- pars[,paste0('Omega',k3)]
    Inc <- pars[,paste0('Inc',k3)]
    M0 <- pars[,paste0('Mo',k3)]
    ms <- (M0+2*pi*(jd.ref-tmin)/P)%%(2*pi)
    E <- kep.mt2(ms,e)
    T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    alpha0 <- K/sin(Inc)/1e3/4.74047#au/yr
    alpha <- alpha0*out$plx#proper motion in mas/yr
    ##semi-major axis is the astrometric signature in micro-arcsec
    A <- cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(Inc)
    B <- sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(Inc)
    F <- -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(Inc)
    G <- -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(Inc)
    Vx <- -sin(T)
    Vy <- cos(T)+e
    pmraP <- alpha*(B*Vx+G*Vy)
    pmdecP <- alpha*(A*Vx+F*Vy)
###calculate POS
    X <- cos(E)-e
    Y <- sqrt(1-e^2)*sin(E)
    beta0 <- P/365.25*(K*1e-3/4.74047)*sqrt(1-e^2)/(2*pi)/sin(Inc)#au
    beta <- beta0*out$plx#mas
    raP <- beta*(B*X+G*Y)
    decP <- beta*(A*X+F*Y)
    Mp <- k2m(K,P,e,Ms=Mstars,Inc=Inc)$ms
    xi <- Mp/(Mstars+Mp)
    eta <- 1e-3#mas to as
    dra.list[[k]] <- dras <- -raP/xi
    ddec.list[[k]] <- ddecs <- -decP/xi
#####simulated orbit for one companion
    par.opt1 <- par.opt
    if(Nsig>1){
        jj <- (1:Nsig)[-k3]
        par.opt1[paste0('K',jj)] <- 0
    }
#    planet <- astrometry.kepler(par.opt1,tt=tt%%24e5-tmin%%24e5,bases=bases,pa=FALSE)$planet
#    planet <- astrometry.kepler(par.opt1,tt=tt,bases=bases,Pmin=0)$planet
    planet <- astrometry.rel(par.opt1,tt=tt,bases=bases)$rel
    dra.sim.list[[k]] <- dra.sim <- planet$dra[[paste0('p',k3)]]
    ddec.sim.list[[k]] <- ddec.sim <- planet$ddec[[paste0('p',k3)]]

###MCMC
    if(Nmc>0){
        dra.mc <- ddec.mc <- c()
        for(j3 in ind.sample){
            qq <- mc[j3,1:Npar]
            if(Nsig>1){
              jj <- (1:Nsig)[-k3]
              qq[paste0('K',jj)] <- 0
            }
#            planet <- astrometry.kepler(qq,tt=tt%%24e5-tmin%%24e5,bases=bases,pa=FALSE)$planet
            planet <- astrometry.rel(qq,tt=tt,bases=bases,eta=0)$rel
#            msini <- K2msini.full(mc[j3,'K1'],exp(mc[j3,'per1']),mc[j3,'e1'],Ms=rnorm(1,out$Mstar,out$eMstar))
#            mp <- msini$ms/sin(mc[j3,'Inc1'])
            dra.mc <- cbind(dra.mc, planet$dra[[paste0('p',k3)]])
            ddec.mc <- cbind(ddec.mc,planet$ddec[[paste0('p',k3)]])
        }
        dra.mc.list[[k]] <- dra.mc
        ddec.mc.list[[k]] <- ddec.mc
    }
###calculate mass
    Ps <- exp(mc[,paste0('per',k3)])
    Ks <- mc[,paste0('K',k3)]
    es <- mc[,paste0('e',k3)]
    Is1 <- Is <- mc[,paste0('Inc',k3)]
#    Ms <- rnorm(2*nrow(mc),out$Mstar,out$eMstar)
    omegas1  <- omegas <- mc[,paste0('omega',k3)]
    Omegas1 <- Omegas <- mc[,paste0('Omega',k3)]
    Mos1 <- Mos <- mc[,paste0('Mo',k3)]
    omegas1[omegas1>pi] <- omegas1[omegas1>pi]-2*pi
    Mos1[Mos1>pi] <- Mos1[Mos1>pi]-2*pi
    Omegas1[Omegas1>pi] <- Omegas1[Omegas1>pi]-2*pi
    if(sd(Mos1)<sd(Mos)) Mos <- Mos1
    if(sd(omegas1)<sd(omegas)) omegas <- omegas1
    if(sd(Omegas1)<sd(Omegas)) Omegas <- Omegas1
    Tps <- tmin-(Mos/(2*pi))*P

#    msini <- K2msini.full(Ks,Ps,es,Ms=Ms)
#    mps <- msini$mj/sin(Is)
    tmp <- k2m(Ks,Ps,es,Ms=Mstars,Inc=Is,more=TRUE)
    mps <- tmp$mj
    aps <- tmp$a
    if(star=='UCAC4569-026385'){
        Ks <- Ks/1e3#m/s->km/s
        mps <- mps/1048#mj -> ms
    }
    mcs <- quantile(mps,c(0.16,0.5,0.84))
#    aps <- (Ms*(Ps/365.25)^2)^(1/3)#au
    acs <- quantile(aps,c(0.16,0.5,0.84))
    pcs <- quantile(Ps,c(0.16,0.5,0.84))
    pyrcs <- quantile(Ps/365.25,c(0.16,0.5,0.84))
    Kcs <- quantile(Ks,c(0.16,0.5,0.84))
    ecs <- quantile(es,c(0.16,0.5,0.84))
    Ics <- quantile(Is,c(0.16,0.5,0.84))*180/pi
    ocs <- (quantile(omegas,c(0.16,0.5,0.84))%%(2*pi))*180/pi
    Ocs <- (quantile(Omegas,c(0.16,0.5,0.84))%%(2*pi))*180/pi
    Mcs <- (quantile(Mos,c(0.16,0.5,0.84))%%(2*pi))*180/pi
    Tcs <- quantile(Tps,c(0.16,0.5,0.84))-24e5#JD-24e5

    dras <- quantile(mc[,'dra'],c(0.16,0.5,0.84))#mas
    ddecs <- quantile(mc[,'ddec'],c(0.16,0.5,0.84))#mas
    dplxs <- quantile(mc[,'dplx'],c(0.16,0.5,0.84))#mas
    dpmras <- quantile(mc[,'dpmra'],c(0.16,0.5,0.84))#mas
    dpmdecs <- quantile(mc[,'dpmdec'],c(0.16,0.5,0.84))#mas

    p <- letters[-1]
#    plx <- out$astrometry[nrow(out$astrometry),'parallax']
#    eplx <- out$astrometry[nrow(out$astrometry),'parallax_error']
    plxcs <- quantile(out$astrometry[out$iref,'parallax']-mc[,'dplx'],c(0.16,0.5,0.84))

    pyr <- pcs
#    mpacs <- rbind(mpacs,c(star,p[k3],pcs[2],pcs[2]-pcs[1],pcs[3]-pcs[2],Kcs[2],Kcs[2]-Kcs[1],Kcs[3]-Kcs[2],ecs[2],ecs[2]-ecs[1],ecs[3]-ecs[2],Ics[2],Ics[2]-Ics[1],Ics[3]-Ics[2],ocs[2],ocs[2]-ocs[1],ocs[3]-ocs[2],Ocs[2],Ocs[2]-Ocs[1],Ocs[3]-Ocs[2],Mcs[2],Mcs[2]-Mcs[1],Mcs[3]-Mcs[2],dras[2],dras[2]-dras[1],dras[3]-dras[2],ddecs[2],ddecs[2]-ddecs[1],ddecs[3]-ddecs[2],dplxs[2],dplxs[2]-dplxs[1],dplxs[3]-dplxs[2],dpmras[2],dpmras[2]-dpmras[1],dpmras[3]-dpmras[2],dpmdecs[2],dpmdecs[2]-dpmdecs[1],dpmdecs[3]-dpmdecs[2],pyrcs[2],pyrcs[2]-pyrcs[1],pyrcs[3]-pyrcs[2],acs[2],acs[2]-acs[1],acs[3]-acs[2],mcs[2],mcs[2]-mcs[1],mcs[3]-mcs[2],Tcs[2],Tcs[2]-Tcs[1],Tcs[3]-Tcs[2],out$Mstar,out$eMstar,out$eMstar,plx,eplx,eplx))
    mpacs <- rbind(mpacs,c(star,p[k3],pcs[2],pcs[2]-pcs[1],pcs[3]-pcs[2],Kcs[2],Kcs[2]-Kcs[1],Kcs[3]-Kcs[2],ecs[2],ecs[2]-ecs[1],ecs[3]-ecs[2],Ics[2],Ics[2]-Ics[1],Ics[3]-Ics[2],ocs[2],ocs[2]-ocs[1],ocs[3]-ocs[2],Ocs[2],Ocs[2]-Ocs[1],Ocs[3]-Ocs[2],Mcs[2],Mcs[2]-Mcs[1],Mcs[3]-Mcs[2],dras[2],dras[2]-dras[1],dras[3]-dras[2],ddecs[2],ddecs[2]-ddecs[1],ddecs[3]-ddecs[2],dplxs[2],dplxs[2]-dplxs[1],dplxs[3]-dplxs[2],dpmras[2],dpmras[2]-dpmras[1],dpmras[3]-dpmras[2],dpmdecs[2],dpmdecs[2]-dpmdecs[1],dpmdecs[3]-dpmdecs[2],pyrcs[2],pyrcs[2]-pyrcs[1],pyrcs[3]-pyrcs[2],acs[2],acs[2]-acs[1],acs[3]-acs[2],mcs[2],mcs[2]-mcs[1],mcs[3]-mcs[2],Tcs[2],Tcs[2]-Tcs[1],Tcs[3]-Tcs[2],mscs[2],mscs[2]-mscs[1],mscs[3]-mscs[2],plxcs[2],plxcs[2]-plxcs[1],plxcs[3]-plxcs[2]))
    cat('Companion mass: ',paste0(round(mcs[2],3),'_',round(mcs[2]-mcs[1],3),'+',round(mcs[3]-mcs[2],3),' Mjup'),'\n')
    cat('Companion orbital period: ',paste0(round(pcs[2],3),'_',round(pcs[2]-pcs[1],3),'+',round(pcs[3]-pcs[2],3),' day'),'\n')
    cat('Companion semi-major axis: ',paste0(round(acs[2],3),'_',round(acs[2]-acs[1],3),'+',round(acs[3]-acs[2],3),' au'),'\n\n')

###important dates
    yrs.ref <- seq(2020,2030,by=2)
                                        #jds <- c(jd.ref,rowSums(time_Yr2jd(yrs)))
    jds.ref <- rowSums(time_Yr2jd(yrs.ref))
#    planet <- astrometry.kepler(par.opt1,tt=jds.ref%%24e5-tmin%%24e5,bases=bases,pa=FALSE)$planet
    planet <- astrometry.rel(par.opt1,tt=jds,bases=bases)$rel
    dra.ref.list[[k]] <- planet$dra[[paste0('p',k3)]]
    ddec.ref.list[[k]] <- planet$ddec[[paste0('p',k3)]]
}

xlim <- eta*rev(range(unlist(dra.sim.list),unlist(dra.list),unlist(dra.mc.list)))
ylim <- eta*range(unlist(ddec.sim.list),unlist(ddec.list),unlist(ddec.mc.list))
ylim <- range(xlim,ylim)
xlim <- rev(ylim)
#ylim[1] <- ylim[1]-0.2*diff(ylim)
#xlim[1] <- xlim[1]-0.2*diff(xlim)
for(k in 1:length(ind.lps)){
    k3 <- ind.lps[k]
    dra.sim <- dra.sim.list[[k]]
    ddec.sim <- ddec.sim.list[[k]]
    if(k==1){
        plot(dra.sim*eta,ddec.sim*eta,xlab=expression(Delta*alpha*'* [\'\']'),ylab=expression(Delta*delta*' [\'\']'),type='l',xlim=xlim,ylim=ylim,col=cols[k],axes=FALSE)
        magaxis(1:2,frame.plot=TRUE)
        magaxis(3:4,labels=FALSE)
    }else{
        lines(dra.sim*eta,ddec.sim*eta,col=cols[k])
    }
    if(Nmc>0){
#        lines(dra.sim*eta,ddec.sim*eta,col=tcol(cols[k],80))
        for(j3 in 1:Nmc){
            lines(dra.mc.list[[k]][,j3]*eta,ddec.mc.list[[k]][,j3]*eta,col=tcol(cols[k],90))
        }
    }
                                        #points(dra.ref,ddec.ref,pch='+',col='cyan')
    dras <- dra.list[[k]]
    ddecs <- ddec.list[[k]]
    rho <- sqrt(dras^2+ddecs^2)*eta#as
    theta1 <- theta <- atan2(dras,ddecs)*180/pi#deg
    theta1[theta1<0] <- theta1[theta1<0]+360
    if(sd(theta1)<sd(theta)) theta <- theta1
    rho.map <- sqrt(dras[1]^2+ddecs[1]^2)*eta#as
    rhos <- quantile(rho,c(0.16,0.5,0.84))
    thetas <- quantile(theta,c(0.16,0.5,0.84))
    rds <- rbind(rds,c(rhos[1],rhos[2]-rhos[1],rhos[3]-rhos[2],thetas[1],thetas[2]-thetas[1],thetas[3]-thetas[2]))
    if(k==length(ind.lps)){
        if(FALSE){
            text(dra.ref.list[[k]]*eta,ddec.ref.list[[k]]*eta,labels=yrs.ref,col='cyan')
        }
        kk <- kde2d(eta*dras,eta*ddecs, n=20)
                                        #contour(kk,levels=c(0.68,0.95,0.997))
                                        #levels <- getLevel(kk,prob=c(0.68,0.95,0.997))
        levels <- getLevel(kk,prob=c(0.68,0.95))
        contour(kk, drawlabels=FALSE, nlevels=length(levels), levels=levels,col='black',add=TRUE,axes=FALSE,lty=c(2,1),lwd=2)
        points(eta*dras[1],eta*ddecs[1],pch='*',col='gold')

                                        #contour(,levels=c(0.68,0.95,0.997))
        points(0,0,pch='+')
                                        #    if(k3==ind.max) legend('top',inset=-0.2,legend=c(as.expression(bquote(rho[MAP]==.(round(rho.map,2))*' \'\'')),as.expression(bquote('<'*rho*'>='*.(round(mean(rho),2))*' \'\'')),as.expression(bquote(sigma*'('*rho*')='*.(round(sd(rho),2))*' \'\''))),xpd=NA,bty='n',horiz=TRUE)
                                        #    save(list=ls(all=TRUE),file='test.Robj')
        legend('top',inset=-0.2,legend=c(as.expression(bquote(rho==.(round(mean(rho),2))%+-%.(round(sd(rho),2))*second)),as.expression(bquote(theta==.(round(mean(theta),2))%+-%.(round(sd(theta),2))*degree))),xpd=NA,bty='n',horiz=TRUE,x.intersp=0.1)
        legend('topright',bty='n',legend=star)
    }
    dx <- dy <- 0.2*diff(ylim)
    xmin <- xlim[2]
    xmax <- xlim[1]
    ymin <- ylim[1]+0.1*diff(ylim)
    ymax <- ylim[2]
    p0 <- c(xmin,ymin)
    p1 <- c(xmin+dx,ymin)
    p2 <- c(xmin,ymin+dy)
                                        #p0[1] <- p0[1]-0.36
                                        #p1[1] <- p1[1]-0.36
                                        #p2[1] <- p2[1]-0.36
if(FALSE){
    arrows(p0[1],p0[2],p1[1],p1[2],length=0.05,angle=45,code=2,col='grey')
    arrows(p0[1],p0[2],p2[1],p2[2],length=0.05,angle=45,code=2,col='grey')
    text(p1[1],p1[2],labels='E',pos=2,col='grey')
    text(p2[1],p2[2],labels='N',pos=3,col='grey')
}
    loglike.out <- mc[,'loglike']
    p <- data.distr(eta*dras,lp=loglike.out,plotf=FALSE)
    q <- data.distr(eta*ddecs,lp=loglike.out,plotf=FALSE)
    p[1] <- dras[1]
    q[1] <- ddecs[1]
    if(FALSE){
        cat('dra statistics:\n',names(p),'\n')
        cat(round(p,2),'\n\n')
        cat('ddec statistics:\n',names(q),'\n')
        cat(round(q,2),'\n')
    }
}

#if(FALSE) dev.off()
dev.off()
#stop()

