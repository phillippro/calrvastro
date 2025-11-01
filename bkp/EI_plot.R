library(RColorBrewer)
if(!exists('mcmc.out')){
    load('/Users/ffeng/Documents/projects/dwarfs/astrometry/results/HD209100/PM_jitter10sd0.2_priormt_mode0_Ndata809_quantifyTRUE_1per1_Nw6_HD209100_snr110v2cut_ind0_1planet_ARMA02_Nsamp1000000_astrometry8parallaxFALSEjitterTRUE_tem1_acc2.9_pretem1P17082.6d_negLmax1612_optpar.Robj')
    load('/Users/ffeng/Documents/projects/dwarfs/astrometry/results/HD209100/PM_jitter10sd0.2_priormt_mode0_Ndata809_quantifyTRUE_1per1_Nw6_HD209100_snr110v2cut_ind0_1planet_ARMA02_Nsamp1000000_astrometry8parallaxFALSEjitterTRUE_tem1_acc2.9_pretem1P17082.6d_negLmax1612.Robj')
    source('mcmc_func.R')
    out <- list()
    out$mcmc.opt <- mcmc.out
    out$par.stat <- par.stat
    out$rel <- list(NA)
    out$Nrel <- 1
    out$Mstar <- 0.754
    out$plx <- 274.8048
    out$Nm <- 0
}
if(tmin<24e5) tmin <- tmin+24e5

####calculate
yrs <- seq(1990,2050,by=0.01)
jds <- time_Yr2jd(yrs)
tsim <- rowSums(jds)
tmp <- astrometry.rel(par.opt,tsim=tsim)
dra.sim <- tmp$rel$raB[[1]]
ddec.sim <- tmp$rel$decB[[1]]

###Monte Carlo
#cal <- cbind(2019,9,1)
#cal <- cbind(2022,9,30)
#cal <- cbind(2022,10,25)
#cal <- cbind(c(2018,10,26),c(2019,07,15))
#cal <- cbind(2018,10,26)
#cal <- cbind(2021,10,1)
cal <- cbind(2022,10,1)
#cal <- cbind(2019,07,15)
jd <- sum(time_Cal2JD(cal))
dras <- c()
ddecs <- c()
#Nsamp <- 1e5
Nsamp <- nrow(mcmc.out)
Npar <- length(par.opt)
inds <- sort(sample(1:nrow(mcmc.out),Nsamp))
#for(k in 1:(Nsamp+1)){
pars <- rbind(par.opt,mcmc.out[inds,1:Npar])
#    if(k==1){
#        par <- par.opt
#    }else{
#        par <- mcmc.out[inds[k-1],1:Npar]
#    }
K <- pars[,'K1']
P <- exp(pars[,'per1'])
e <- pars[,'e1']
omega <- pars[,'omega1']
Omega <- pars[,'Omega1']
Inc <- pars[,'Inc1']
M0 <- pars[,'Mo1']
ms <- (M0+2*pi*(jd-tmin)/P)%%(2*pi)
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
msini <- K2msini.full(K,P,e,Ms=out$Mstar)
Mp <- msini$ms/sin(Inc)
eta <- Mp/(out$Mstar+Mp)
dras <- -raP/eta
ddecs <- -decP/eta
nlevel <- 3
#my.cols <- rev(brewer.pal(nlevel, "RdYlBu"))
#my.cols <- rev(heat.colors(nlevel))
#my.cols <- rev(terrain.colors(nlevel))
#my.cols <- rev(topo.colors(nlevel))
my.cols <- rev(hcl.colors(nlevel))
getLevel <- function(kk,prob=0.95){
    dx <- diff(kk$x[1:2])
    dy <- diff(kk$y[1:2])
    sz <- sort(kk$z)
    c1 <- cumsum(sz) * dx * dy
    approx(c1, sz, xout = 1 - prob)$y
}
pdf('EI_orbit4.pdf',6,6)
size <- 1.2
#par(mar=c(5,5,1,1),cex=size,cex.lab=size,cex.axis=size)
par(mar=c(5,5,3,1),cex=size,cex.lab=size,cex.axis=size)
eta <- 1e-3
xlim <- eta*rev(range(dra.sim,dras))
ylim <- eta*range(dra.sim,ddecs)
ylim <- range(xlim,ylim)
xlim <- rev(ylim)
plot(dra.sim*eta,ddec.sim*eta,xlab=expression(Delta*alpha*'* [arcsecond]'),ylab=expression(Delta*delta*' [arcsecond]'),type='l',xlim=xlim,ylim=ylim)
points(eta*dras,eta*ddecs,pch='.',col=tcol('grey',80))
kk <- kde2d(eta*dras,eta*ddecs, n=20)
#contour(kk,levels=c(0.68,0.95,0.997))
#levels <- getLevel(kk,prob=c(0.68,0.95,0.997))
levels <- getLevel(kk,prob=c(0.68,0.95))
contour(kk, drawlabels=FALSE, nlevels=length(levels), levels=levels,col=my.cols,add=TRUE,axes=FALSE)
#contour(,levels=c(0.68,0.95,0.997))
points(eta*dras[1],eta*ddecs[1],col='red',pch='+')
points(0,0,pch='+')
rho <- sqrt(dras^2+ddecs^2)*eta#as
rho.map <- sqrt(dras[1]^2+ddecs[1]^2)*eta#as
#title(main='epoch: 25-Oct-2022',line=2)
title(main=paste('epoch: ',paste(as.numeric(cal),collapse='-')),line=2)
legend('top',inset=-0.1,legend=c(as.expression(bquote(rho[MAP]==.(round(rho.map,2)))),as.expression(bquote('<'*rho*'>='*.(round(mean(rho),2)))),as.expression(bquote(sigma*'('*rho*')='*.(round(sd(rho),2))))),xpd=NA,bty='n',horiz=TRUE)
dx <- dy <- (ylim[2]-ylim[1])/10
xmax <- xlim[1]
ymin <- ylim[1]
p0 <- c(xmax-3*dx,ymin+dy)
p1 <- c(xmax-1*dx,ymin+dy)
p2 <- c(xmax-3*dx,ymin+3*dy)
arrows(p0[1],p0[2],p1[1],p1[2],length=0.05,angle=45,code=2,col='grey')
arrows(p0[1],p0[2],p2[1],p2[2],length=0.05,angle=45,code=2,col='grey')
text(p1[1],p1[2],labels='E',pos=2,col='grey')
text(p2[1],p2[2],labels='N',pos=3,col='grey')
dev.off()

p <- data.distr(eta*dras,lp=loglike.out,plotf=FALSE)
q <- data.distr(eta*ddecs,lp=loglike.out,plotf=FALSE)
p[1] <- dras[1]
q[1] <- ddecs[1]
cat('dra statistics:\n',names(p),'\n')
cat(round(p,2),'\n\n')

cat('ddec statistics:\n',names(q),'\n')
cat(round(q,2),'\n')
