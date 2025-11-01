target <- star 
dir <- 'orbits/'
f2 <- paste0(dir,target,'_prediction.pdf')
cat(f2,'\n')
pdf(f2,6,6)
size <- 1.2
#par(mar=c(5,5,1,1),cex=size,cex.lab=size,cex.axis=size)
par(mar=c(5,5,3,1),cex=size,cex.lab=size,cex.axis=size)
####position at the reference epoch
cal <- cbind(2021,10,1)
jd.ref <- rowSums(time_Cal2JD(cal))
dras <- c()
ddecs <- c()
Npar <- length(par.opt)
#inds <- c(which.max(mc[,'loglike']),sort(sample(1:nrow(mc),Nsamp)))
inds <- c(which.max(mc[,'loglike']),ind.sample)
#for(k in 1:(Nsamp+1)){
pars <- rbind(par.opt,mc[inds,1:Npar])
#    if(k==1){
#        par <- par.opt
#    }else{
#        par <- mc[inds[k-1],1:Npar]
#    }
ind.max <- which.max(par.opt[grepl('per',names(par.opt))])
K <- pars[,paste0('K',ind.max)]
P <- exp(pars[,paste0('per',ind.max)])
e <- pars[,paste0('e',ind.max)]
omega <- pars[,paste0('omega',ind.max)]
Omega <- pars[,paste0('Omega',ind.max)]
Inc <- pars[,paste0('Inc',ind.max)]
M0 <- pars[,paste0('Mo',ind.max)]
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
msini <- K2msini.full(K,P,e,Ms=out$Mstar)
Mp <- msini$ms/sin(Inc)
eta <- Mp/(out$Mstar+Mp)
dras <- -raP/eta
ddecs <- -decP/eta


####
nlevel <- 3
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
planet <- astrometry.kepler(par.opt,tt=tt%%24e5-tmin%%24e5,bases=bases,pa=FALSE)$planet
dra.sim <- -planet[,'ra']*(1048*Mstar+mp.opt)/mp.opt
ddec.sim <- -planet[,'dec']*(1048*Mstar+mp.opt)/mp.opt
eta <- 1e-3#mas to as
xlim <- eta*rev(range(dra.sim,dras))
ylim <- eta*range(ddec.sim,ddecs)
ylim <- range(xlim,ylim)
xlim <- rev(ylim)
ylim[1] <- ylim[1]-0.2*diff(ylim)
xlim[1] <- xlim[1]-0.2*diff(xlim)
plot(dra.sim*eta,ddec.sim*eta,xlab=expression(Delta*alpha*'* [\'\']'),ylab=expression(Delta*delta*' [\'\']'),type='l',xlim=xlim,ylim=ylim,col='red')
###MCMC
if(Nmc>0){
    mc <- out$mcmc.opt[[paste0('sig',Nkep)]]
    inds <- sample(1:nrow(mc),Nmc)
    Npar <- length(par.opt)
    for(j3 in inds){
        planet <- astrometry.kepler(mc[j3,1:Npar],tt=tt%%24e5-tmin%%24e5,bases=bases,pa=FALSE)$planet
        msini <- K2msini.full(mc[j3,'K1'],exp(mc[j3,'per1']),mc[j3,'e1'],Ms=out$Mstar)
        mp <- msini$ms/sin(mc[j3,'Inc1'])
        dra.sim <- -planet[,'ra']*(Mstar+mp)/mp
        ddec.sim <- -planet[,'dec']*(Mstar+mp)/mp
        lines(dra.sim*eta,ddec.sim*eta,col=tcol('red',80))
    }
}
###important dates
yrs.ref <- seq(2020,2030,by=2)
#jds <- c(jd.ref,rowSums(time_Yr2jd(yrs)))
jds.ref <- rowSums(time_Yr2jd(yrs.ref))
planet <- astrometry.kepler(par.opt,tt=jds.ref%%24e5-tmin%%24e5,bases=bases,pa=FALSE)$planet
dra.ref <- -planet[,'ra']*(1048*Mstar+mp.opt)/mp.opt*eta
ddec.ref <- -planet[,'dec']*(1048*Mstar+mp.opt)/mp.opt*eta
#points(dra.ref,ddec.ref,pch='+',col='cyan')
text(dra.ref,ddec.ref,labels=yrs.ref,col='cyan')

kk <- kde2d(eta*dras,eta*ddecs, n=20)
#contour(kk,levels=c(0.68,0.95,0.997))
#levels <- getLevel(kk,prob=c(0.68,0.95,0.997))
levels <- getLevel(kk,prob=c(0.68,0.95))
contour(kk, drawlabels=FALSE, nlevels=length(levels), levels=levels,col=my.cols,add=TRUE,axes=FALSE)
points(eta*dras[1],eta*ddecs[1],pch='*',col='gold')
#contour(,levels=c(0.68,0.95,0.997))
points(0,0,pch='+')
rho <- sqrt(dras^2+ddecs^2)*eta#as
rho.map <- sqrt(dras[1]^2+ddecs[1]^2)*eta#as
#title(main='epoch: 25-Oct-2022',line=2)
#title(main=paste('Star-companion separation at epoch: ',paste(as.numeric(cal),collapse='-')),line=2)
title(main=paste('Separation at epoch: ',paste(as.numeric(cal),collapse='-')),line=2)
legend('top',inset=-0.1,legend=c(as.expression(bquote(rho[MAP]==.(round(rho.map,2))*' \'\'')),as.expression(bquote('<'*rho*'>='*.(round(mean(rho),2))*' \'\'')),as.expression(bquote(sigma*'('*rho*')='*.(round(sd(rho),2))*' \'\''))),xpd=NA,bty='n',horiz=TRUE)
dx <- dy <- 0.2*diff(ylim)
xmin <- xlim[2]
xmax <- xlim[1]
ymin <- ylim[1]
ymax <- ylim[2]
p0 <- c(xmin,ymin)
p1 <- c(xmin+dx,ymin)
p2 <- c(xmin,ymin+dy)
#p0[1] <- p0[1]-0.36
#p1[1] <- p1[1]-0.36
#p2[1] <- p2[1]-0.36
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
dev.off()

####mass estimation
f3 <- paste0(dir,target,'_mass.pdf')
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
