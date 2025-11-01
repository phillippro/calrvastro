source('mcmc_func.R')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
   f <- args[1]
}else{
   f <- 'HD42581_natural_astro8_Niter10000000_Ncores4_ofac2_Nset5_Esd1_transit0_P39195_acc0.039_lnlmax'
#    f <- 'HD42581_natural_offsetTRUE_Niter3200000_Ncores16_ofac2_Nset5_Esd1_transit0_P78745_acc0.18_lnlmax'
}
star <- gsub('_.+','',f)
secondary <- paste0(star,'B')
fs <- list.files(paste0('results/',star),pattern=paste0(f,'.+Robj'),full.name=TRUE)
if(length(fs)==0) stop(paste0('I cannot find file with pattern:',paste0(f,'.+Robj'),'\n'))
if(!exists('out')) load(fs[1])

mc <- out$mcmc.opt$sig1
Npar <- ncol(mc)-2
Ntry <- 100


####plot
fpdf <- paste0('results/',target,'_prediction.pdf')
cat(fpdf,'\n')
pdf(fpdf,8,8)
par(mfrow=c(2,2),mar=c(4,4,1,1))

for(j in 1:Ntry){
ii <- sample(1:nrow(mc),1)
par.opt <- mc[ii,1:Npar]

yr2d <- 365.25
##astrometry fit
data.astrometry <- out$astrometry
etaH <- etaG <- 1
tsim0 <- seq(data.astrometry[1,1]-2e4,data.astrometry[2,1]+2e4,length.out=1e3)
tsim <- tsim0-tmin
Nt <- length(tsim)
#planet <- astrometry.kepler(par.opt,tt=tsim,bases=bases,more=TRUE,Mstar=Mstar,pa=FALSE)$planet
planet <- astrometry.kepler(par.opt,tt=tsim,bases=bases,more=TRUE,Mstar=Mstar,pa=FALSE)$planet
msini <- K2msini.full(par.opt['K1'],Popt,par.opt['e1'],Ms=data.astrometry[2,'mass'])
Mp <- msini$ms/sin(par.opt['Inc1'])
eta <- Mp/(Mstar+Mp)
planet1 <- -planet/eta

####load data
tab <- read.table('/Users/ffeng/Documents/projects/dwarfs/data/combined/HD42581/HD42581_astrometry.rel',header=TRUE)

####plot
ts <- seq(2000,2030,by=5)
jds <-  rowSums(time_Yr2jd(ts))
inds <- unlist(lapply(jds,function(jd) which.min(abs(tsim0-jd))))
##binary motion
dra <- planet1[,'ra']
ddec <- planet1[,'dec']
dpmra <- planet1[,'pmra']
dpmdec <- planet1[,'pmdec']
xlim <- range(dra,tab[,'dra'])
ylim <- range(ddec,tab[,'ddec'])
plot(dra,ddec,xlab=expression(Delta*alpha*'*[mas]'),ylab=expression(Delta*delta*'[mas]'),type='l',xlim=rev(xlim),ylim=ylim)
ind <- round(nrow(planet)/2)
arrows(dra[ind],ddec[ind],dra[ind+1],ddec[ind+1],length=0.1,angle=30,code=2,col='grey')
points(dra[inds],ddec[inds],pch='+',col='red')
text(dra[inds],ddec[inds],labels=ts,pos=2)
points(0,0,pch='+')
points(tab[1,'dra'],tab[1,'ddec'],col='red')
points(tab[-1,'dra'],tab[-1,'ddec'],col='blue')

xlim <- range(dpmdec)
ylim <- range(dpmdec)
plot(dpmra,dpmdec,xlab=expression(Delta*mu[alpha]*'[mas/yr]'),ylab=expression(Delta*mu[delta]*'[mas/yr]'),type='l')
arrows(dpmra[ind],dpmdec[ind],dpmra[ind+1],dpmdec[ind+1],length=0.1,angle=30,code=2,col='grey')
points(dpmra[inds],dpmdec[inds],pch='+',col='red')
text(dpmra[inds],dpmdec[inds],labels=ts,pos=4)
points(0,0,pch='+')

####reflex motion
dra <- planet[,'ra']
ddec <- planet[,'dec']
dpmra <- planet[,'pmra']
dpmdec <- planet[,'pmdec']
xlim <- range(dra)
ylim <- range(ddec)
plot(dra,ddec,xlab=expression(Delta*alpha*'*[mas]'),ylab=expression(Delta*delta*'[mas]'),type='l',xlim=rev(xlim),ylim=ylim)
ind <- round(nrow(planet)/2)
arrows(dra[ind],ddec[ind],dra[ind+1],ddec[ind+1],length=0.1,angle=30,code=2,col='grey')
points(dra[inds],ddec[inds],pch='+',col='red')
text(dra[inds],ddec[inds],labels=ts,pos=2)
points(0,0,pch='+')

xlim <- range(dpmra)
ylim <- range(dpmdec)
plot(dpmra,dpmdec,xlab=expression(Delta*mu[alpha]*'[mas/yr]'),ylab=expression(Delta*mu[delta]*'[mas/yr]'),type='l')
arrows(dpmra[ind],dpmdec[ind],dpmra[ind+1],dpmdec[ind+1],length=0.1,angle=30,code=2,col='grey')
points(dpmra[inds],dpmdec[inds],pch='+',col='red')
text(dpmra[inds],dpmdec[inds],labels=ts,pos=4)
points(0,0,pch='+')
}
dev.off()
