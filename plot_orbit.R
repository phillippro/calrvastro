library("viridis")
source('mcmc_func.R')
m22 <- read.table('../data/combined/mamajek22.txt',header=TRUE)
ms <- as.numeric(m22[,'Msun'])
mg <- as.numeric(m22[,'M_G'])
ind <- which(!is.na(ms) & !is.na(mg))
mrl.m22 <- approxfun(ms[ind],mg[ind])
mlow.m22 <- min(ms[ind])
mup.m22 <- max(ms[ind])

#if(!exists('out')) load('results/HD209100/HD209100_lpspm_photovaryFALSE_relativityFALSE_Niter2200000_Ncores16_ofac2_Nset8_hg123_transit0_P20326_acc0.73_sinI_lnlmax-8245.Robj')
#if(!exists('out')) load('results/GL676A/GL676A_lpspm_photovaryFALSE_relativityFALSE_Niter2200000_Ncores16_ofac2_Nset3_hg123_transit0_P1052d10803d35d4_acc1_sinI_lnlmax-625.Robj')
if(!exists('out')) load('results/HD209100/HD209100_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter110000_dr1FALSE_230617_Nset6_hg123+_Nsig1_P53548_Esd1_astro5TRUE_acc1.2_rvcbarycentric_priorf3_lnpmax-8970_lnlmax-8934.Robj')
star <- target
target <- star
dir <- 'orbits/'
#cal <- cbind(2022,7,15)
cal <- cbind(2019,9,15)
jd.ref <- rowSums(time_Cal2JD(cal))
#yrs.ref <- seq(2018,2026,by=2)
yrs.ref <- 2023
                                          #jds <- c(jd.ref,rowSums(time_Yr2jd(yrs)))
jds.ref <- rowSums(time_Yr2jd(yrs.ref))
f2 <- paste0(dir,target,'_prediction_epoch',paste(cal,collapse='-'),'.pdf')
cat(f2,'\n')
pdf(f2,6,6)
size <- 1.2
#par(mar=c(5,5,1,1),cex=size,cex.lab=size,cex.axis=size)
par(mar=c(5,5,3,1),cex=size,cex.lab=size,cex.axis=size)
####position at the reference epoch
dras <- c()
ddecs <- c()
Npar <- length(par.opt)
mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
Nsamp <- min(nrow(mc),1e3)
ind.sample <- sample(1:nrow(mc),Nsamp)
inds <- c(which.max(mc[,'loglike']),ind.sample)
###MCMC samples of parameters
#pars <- rbind(par.opt,mc[inds,1:Npar])
mstars0 <- mc[,'Mstar']
pars <- mc[ind.sample,1:Npar]
ind.max <- which.max(par.opt[grepl('^per',names(par.opt))])
mstars <- pars[,'Mstar']
Mstar <- pars[ind.max,'Mstar']

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
yrs <- seq(2020,2040,by=5)
#jds <- c(jd.ref,rowSums(time_Yr2jd(yrs)))
jds <- rowSums(time_Yr2jd(yrs))
Pmax <- max(exp(pars[,grepl('per',colnames(pars))]))
tall <- trv.all
t1 <- min(mean(tall)%%24e5-max(Pmax)/2,jd.ref%%24e5)
t2 <- max(mean(tall)%%24e5+max(Pmax)/2,jd.ref%%24e5,jds.ref)
tsim <- tt <- seq(t1,t2,length.out=10000)
Plow <- 1000
Nmc <- 10
ind.lps <- which(par.opt[grepl('^per',names(par.opt))]>log(Plow))
####list to save simulation data
ddec.ref.list <- dra.ref.list <- ddec.mc.list <- dra.mc.list <- dra.sim.list <- ddec.sim.list <- dra.list <- ddec.list <- list()

####loop to simulate the orbits
for(k in 1:length(ind.lps)){
    k3 <- ind.lps[k]
if(FALSE){
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
    Mp <- k2m(K,P,e,Ms=mstars)$ms
    xi <- Mp/(mstars+Mp)
    eta <- 1e-3#mas to as
}
dras <- ddecs <- c()
for(j in 1:nrow(pars)){
    par1 <- pars[j,]
    if(Nsig>1){
        jj <- (1:Nsig)[-k3]
        par1[paste0('K',jj)] <- 0
    }
    planet <- astrometry.rel(par1,tt=jd.ref)
    dras <- c(dras,planet$rel$dra[[paste0('p',k3)]])
    ddecs <- c(ddecs,planet$rel$ddec[[paste0('p',k3)]])
}
    dra.list[[k]] <- dras
    ddec.list[[k]] <- ddecs
#####simulated orbit for one companion
    par.opt1 <- par.opt
    if(Nsig>1){
        jj <- (1:Nsig)[-k3]
        par.opt1[paste0('K',jj)] <- 0
    }
#    planet <- astrometry.kepler(par.opt1,tt=tt%%24e5-tmin%%24e5,bases=bases,pa=FALSE)
    planet <- astrometry.rel(par.opt1,tt=tsim)
#    dra.sim.list[[k]] <- dra.sim <- -planet$rel$dra$p1*(1048*Mstar+mp.opt)/mp.opt
#    ddec.sim.list[[k]] <- ddec.sim <- -planet$rel$ddec$p1*(1048*Mstar+mp.opt)/mp.opt
    dra.sim.list[[k]] <- dra.sim <- planet$rel$dra$p1
    ddec.sim.list[[k]] <- ddec.sim <- planet$rel$ddec$p1
###MCMC
    if(Nmc>0){
        dra.mc <- ddec.mc <- c()
        for(j3 in ind.sample){
            pp <- mc[j3,1:Npar]
            if(Nsig>1){
              jj <- (1:Nsig)[-k3]
              pp[paste0('K',jj)] <- 0
            }
#            planet <- astrometry.kepler(pp,tt=tt%%24e5-tmin%%24e5,bases=bases,pa=FALSE)$planet
            planet <- astrometry.rel(pp,tt=tsim)
            mp <- k2m(mc[j3,'K1'],exp(mc[j3,'per1']),mc[j3,'e1'],Ms=mstars0[j3])$ms
            dra.mc <- cbind(dra.mc, planet$rel$dra$p1)
            ddec.mc <- cbind(ddec.mc,planet$rel$ddec$p1)
        }
        dra.mc.list[[k]] <- dra.mc
        ddec.mc.list[[k]] <- ddec.mc
    }

###important dates
    planet <- astrometry.rel(par.opt1,tt=jds.ref)
#    planet <- astrometry.kepler(par.opt1,tt=jds.ref%%24e5-tmin%%24e5,bases=bases,pa=FALSE)$planet
    dra.ref.list[[k]] <- planet$rel$dra$p1*eta
    ddec.ref.list[[k]] <- planet$rel$ddec$p1*eta
}

if(length(dra.mc.list)==0){
xlim <- eta*rev(range(unlist(dra.sim.list),na.rm=TRUE))
ylim <- eta*range(unlist(ddec.sim.list),na.rm=TRUE)
}else{
xlim <- eta*rev(range(unlist(dra.sim.list),unlist(dra.mc.list),na.rm=TRUE))
ylim <- eta*range(unlist(ddec.sim.list),unlist(ddec.mc.list),na.rm=TRUE)
}
ylim <- range(xlim,ylim)
xlim <- rev(ylim)
#ylim[1] <- ylim[1]-0.2*diff(ylim)
#xlim[1] <- xlim[1]-0.2*diff(xlim)
for(k in 1:length(ind.lps)){
    k3 <- ind.lps[k]
    dra.sim <- dra.sim.list[[k]]
    ddec.sim <- ddec.sim.list[[k]]
    if(k==1){
        plot(dra.sim*eta,ddec.sim*eta,xlab=expression(Delta*alpha*'* [\'\']'),ylab=expression(Delta*delta*' [\'\']'),type='l',xlim=xlim,ylim=ylim,col=cols[k])
    }else{
        lines(dra.sim*eta,ddec.sim*eta,col=cols[k])
    }
    if(Nmc>0){
#        lines(dra.sim*eta,ddec.sim*eta,col=tcol(cols[k],80))
        for(j3 in 1:Nmc){
            lines(dra.mc.list[[k]][,j3]*eta,ddec.mc.list[[k]][,j3]*eta,col=tcol(cols[k],50))
        }
    }
                                        #points(dra.ref,ddec.ref,pch='+',col='cyan')
###best prediction
    points(eta*dras[1],eta*ddecs[1],pch='*',col='blue',cex=3)
###uncertainty of the prediction
    dras <- dra.list[[k]]
    ddecs <- ddec.list[[k]]
    if(k3==ind.max){
    xs <- as.numeric(dra.ref.list[[k]])
    ys <- as.numeric(ddec.ref.list[[k]])
#    points(xs,ys,pch='+',col='darkgreen',cex=0.5)
#    text(xs,ys,labels=yrs.ref,col='darkgreen',pos=1)
    kk <- kde2d(eta*dras,eta*ddecs, n=20)
                                        #contour(kk,levels=c(0.68,0.95,0.997))
                                        #levels <- getLevel(kk,prob=c(0.68,0.95,0.997))
    levels <- getLevel(kk,prob=c(0.68,0.95,0.9974))
#    my.cols <- heat.colors(length(levels))
    cols <- viridis(length(levels))
#    levels <- getLevel(kk,prob=c(0.68,0.95))
    contour(kk, drawlabels=FALSE, nlevels=length(levels), levels=levels,add=TRUE,axes=FALSE,lty=1,col=cols,lwd=2)
                                        #contour(,levels=c(0.68,0.95,0.997))
    points(0,0,pch='+')
    }
    rho <- sqrt(dras^2+ddecs^2)*eta#as
    rho.map <- rho[1]
    pa <- (atan2(dras,ddecs)%%(2*pi))*180/pi
    pa.map <- pa[1]
#    ps <- rd2ps(mean(dras),sd(dras),mean(ddecs),sd(ddecs))
#    if(k3==ind.max) title(main=paste('Separation at epoch: ',paste(as.numeric(cal),collapse='-')),line=2)
    if(k3==ind.max) title(main=paste('epoch: ',paste(as.numeric(cal),collapse='-')),line=2)
    if(k3==ind.max) legend('top',inset=-0.13,legend=c(as.expression(bquote(rho[MAP]==.(round(rho.map,2))*' \'\'')),as.expression(bquote('<'*rho*'>='*.(round(mean(rho),2))*' \'\'')),as.expression(bquote(sigma*'('*rho*')='*.(round(sd(rho),2))*' \'\''))),xpd=NA,bty='n',horiz=TRUE)
    if(k3==ind.max) legend('top',inset=-0.08,legend=c(as.expression(bquote(theta[MAP]==.(round(pa.map,0))*' deg')),as.expression(bquote('<'*theta*'>='*.(round(mean(pa),0))*' deg')),as.expression(bquote(sigma*'('*theta*')='*.(round(sd(pa),0))*' deg'))),xpd=NA,bty='n',horiz=TRUE)
    dx <- dy <- 0.2*diff(ylim)
    xmin <- xlim[2]
    xmax <- xlim[1]
    ymin <- ylim[1]
    ymax <- ylim[2]
    p0 <- c(xmin,ymin)
    p1 <- c(xmin+dx,ymin)
    p2 <- c(xmin,ymin+dy)
    arrows(p0[1],p0[2],p1[1],p1[2],length=0.05,angle=45,code=2,col='grey')
    arrows(p0[1],p0[2],p2[1],p2[2],length=0.05,angle=45,code=2,col='grey')
    text(p1[1],p1[2],labels='E',pos=2,col='grey')
    text(p2[1],p2[2],labels='N',pos=3,col='grey')
    loglike.out <- mc[,'loglike']
    p <- data.distr(eta*dras,lp=loglike.out,plotf=FALSE)
    q <- data.distr(eta*ddecs,lp=loglike.out,plotf=FALSE)
    p[1] <- dras[1]
    q[1] <- ddecs[1]
    cat('dra statistics:\n',names(p),'\n')
    cat(round(p,2),'\n\n')
    cat('ddec statistics:\n',names(q),'\n')
    cat(round(q,2),'\n')
}

dev.off()

if(FALSE){
####mass estimation
f3 <- paste0(dir,target,'_mass_epoch',paste0(cal,collapse='-'),'.pdf')
pdf(f3,8,6)
par(mar=c(5,5,1,5))
mcomp <- Mp*1048
pp <- data.distr(mcomp,plotf=FALSE)
pp[1] <- Mp[1]
hist(mcomp,breaks=20,xlab='Companion mass [Mjup]',ylab='Freq.',main='')
abline(v=pp['med'],col='red')
abline(v=pp[c('xminus.1sig','xplus.1sig')],col='blue')
legend('topright',legend=paste(names(pp),'=',round(pp,1)),bty='n')
dev.off()
}
