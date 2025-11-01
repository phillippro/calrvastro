library(paletteer)
library(ggplot2)
source('mcmc_func.R')
source('timing_function.R')
#load('results/HD209100/HD209100_calFALSE_coplanarFALSE_resonanceFALSE_relativityFALSE_Niter2200000_Ncores16_ofac2_Nset5_hg123_transit0_P15663_Esd1_acc1.5_rvcbarycentric_lnlmax-8150.Robj')
#if(!exists('out')) load('results/HD222237/HD222237_lpspm_photovaryFALSE_relativityFALSE_Niter2200000_Ncores16_ofac2_Nset4_hg123_transit0_P22003_acc0.25_sinI_lnlmax-662.Robj')
#load('results/HD22049/HD22049_calFALSE_coplanarFALSE_resonanceFALSE_relativityFALSE_Niter2100000_Ncores16_ofac2_Nset13_hg123_transit0_P2695_Esd1_acc3_rvcbarycentric_lnlmax-4267.Robj')
#load('results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter2200000_sppriorTRUE_dr1FALSE_230617_Nset1_hg123_Nsig1_P886_Esd1_astro5TRUE_acc0.38_rvcbarycentric_lnpmax-223_lnlmax-199.Robj')
load('results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter2200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P879_Esd1_astro5TRUE_acc0.29_rvcbarycentric_priorf3_lnpmax-225_lnlmax-200.Robj')
if(is.null(out$comp.epoch)) out$comp.epoch <- 1
if(is.null(out$cats)) out$cats <- out$gdrs
if(is.null(out$cat.ind)) out$cat.ind <- out$gdr.ind
gdr3.epoch <- c(sum(time_Yr2jd(2016)),rowSums(time_CalHms2JD(rbind(c(2014,07,25,10,30,0),c(2017,05,28,8,44,0)))))
#gost <- read.csv('../data/combined/UCAC4569-026385/UCAC4569-026385_gost.csv')
gost <- read.csv(paste0('../data/combined/',target,'/',target,'_gost.csv'))
inds <- which(gost[,'ObservationTimeAtBarycentre.BarycentricJulianDateInTCB.']>=gdr3.epoch[2] & gost[,'ObservationTimeAtBarycentre.BarycentricJulianDateInTCB.']<=gdr3.epoch[3])
obs <- out$astrometry[out$iref,]
ind.ref <- which.max(obs[,1])
obs['ra'] <- obs[ind.ref,'ra']-par.opt['dra']/3.6e6/cos(obs[ind.ref,'dec']/180*pi)
obs['dec'] <- obs[ind.ref,'dec']-par.opt['ddec']/3.6e6
obs['pmra'] <- obs[ind.ref,'pmra']-par.opt['dpmra']
obs['pmdec'] <- obs[ind.ref,'pmdec']-par.opt['dpmdec']
obs['parallax'] <- obs[ind.ref,'parallax']-par.opt['dplx']
obs['radial_velocity'] <- obs[ind.ref,'radial_velocity']
t <- out$gost[,'BJD']-as.numeric(obs['ref_epoch'])
obss <- obs.lin.prop(obs[c('ra','dec','parallax','pmra','pmdec','radial_velocity')],t,PA=TRUE)
al30.plx <- out$astrometry[out$igdr32,'parallax']*gost[inds,'parallaxFactorAlongScan']
ac30.plx <- out$astrometry[out$igdr32,'parallax']*gost[inds,'parallaxFactorAcrossScan']
al20.plx <- al30.plx[out$gdr.ind[[1]]]
ac20.plx <- ac30.plx[out$gdr.ind[[1]]]
ra30.plx <- al30.plx*sin(out$gost[,'psi'])+ac30.plx*cos(out$gost[,'psi'])
dec30.plx <- -al30.plx*cos(out$gost[,'psi'])+ac30.plx*sin(out$gost[,'psi'])
ra20.plx <- al20.plx*sin(out$gost[out$gdr.ind[[1]],'psi'])+ac20.plx*cos(out$gost[out$gdr.ind[[1]],'psi'])
dec20.plx <- -al20.plx*cos(out$gost[out$gdr.ind[[1]],'psi'])+ac20.plx*sin(out$gost[out$gdr.ind[[1]],'psi'])

al.plx <- obss[,'parallax']*gost[inds,'parallaxFactorAlongScan']
ac.plx <- obss[,'parallax']*gost[inds,'parallaxFactorAcrossScan']
ra3.plx <- al.plx*sin(out$gost[,'psi'])+ac.plx*cos(out$gost[,'psi'])
dec3.plx <- -al.plx*cos(out$gost[,'psi'])+ac.plx*sin(out$gost[,'psi'])
ra2.plx <- ra3.plx[out$gdr.ind[[1]]]
dec2.plx <- dec3.plx[out$gdr.ind[[1]]]
#ra.plx <- obss[,'parallax']*out$gost[,'parf']*sin(out$gost[,'psi'])
#dec.plx <- obss[,'parallax']*out$gost[,'parf']*cos(out$gost[,'psi'])
obss[,'ra'] <- obss[,'ra']+ra3.plx/3.6e6*cos(obss[,'dec']/180*pi)
obss[,'dec'] <- obss[,'dec']+dec3.plx/3.6e6
abscissae <- obss[,'parallax']*out$gost[,'parf']
abscissae0 <- as.numeric(out$astrometry[out$iref,'parallax'])*out$gost[,'parf']
xobs <- (obss[,'ra']-obs['ra'])*cos(obss[,'dec']/180*pi)
yobs <- obss[,'dec']
abs.obs <- (par.opt['dra']+par.opt['dpmra']*t)*sin(out$gost[,'psi'])+(par.opt['ddec']+par.opt['dpmdec']*t)*cos(out$gost[,'psi'])+par.opt['dplx']*out$gost[,'parf']

####model prediction
#tsim <- seq(min(out$[,'ref_epoch']),max(out$[,'ref_epoch']),length.out=1e3)
#obs.sim <- obs.lin.prop(obs[c('ra','dec','parallax','pmra','pmdec','radial_velocity')],tsim-obs['ref_epoch'],PA=TRUE)

if(FALSE){
pdf('test.pdf',16,16)
par(mfrow=c(4,4))
plot(obss[,'ra'],obss[,'dec'],xlab='R.A.',ylab='Decl.',pch=20)
plot(ra30.plx,dec30.plx,xlab='R.A. due to plx',ylab='Decl. due to plx',pch=20)
points(ra20.plx,dec20.plx,col='red',pch=20)
points(ra3.plx,dec3.plx,col='black',pch=20)
points(ra2.plx,dec2.plx,col='red',pch=20)
plot(t,abscissae,xlab='t[day]',ylab='abscissae',pch=20)
plot(t,abscissae,xlab='t[day]',ylab='abscissae',pch=20,ylim=c(-0.8,0.8))
points(t,abscissae0,xlab='t[day]',col='red')
plot(t,obss[,'parallax'],xlab='t[day]',ylab='plx',pch=20)
dev.off()
}

####ggplot
m22 <- read.table('../data/combined/mamajek22.txt',header=TRUE)
m22[m22=='...'|m22=='....'|m22=='.....'|m22=='......'] <- NA
ms <- as.numeric(m22[,'Msun'])
mg <- as.numeric(gsub(':','',m22[,'M_G']))
ind <- which(!is.na(ms) & !is.na(mg))
mrl.m22 <- approxfun(ms[ind],mg[ind])
mlow.m22 <- min(ms[ind])
mup.m22 <- max(ms[ind])
gdr1plx <- 'raw'

reflex.sim <- astrometry.epoch(par.opt,tt=tsim)$epoch

reflex.gost <- astrometry.epoch(par.opt,tt=out$gost[,'BJD'],bases=bases)$epoch
kep <- astrometry.kepler(par.opt,bases=bases)
###observed position and proper motion relative to the predicted barycentric position and proper motion
dra.obs <- ((out$astrometry[out$astro.index,'ra']-kep$barycenter[out$astro.index,'ra']))*cos(out$astrometry[out$astro.index,'dec']/180*pi)*3.6e6
ddec.obs <- ((out$astrometry[out$astro.index,'dec']-kep$barycenter[out$astro.index,'dec']))*3.6e6
dpmra.obs <- out$astrometry[out$astro.index,'pmra']-kep$barycenter[out$astro.index,'pmra']
dpmdec.obs <- out$astrometry[out$astro.index,'pmdec']-kep$barycenter[out$astro.index,'pmdec']
dplx.obs <- out$astrometry[out$astro.index,'parallax']-kep$barycenter[out$astro.index,'parallax']

###simulated position and proper motion relative to the barycenter motion
dra.model <- ddec.model <- dplx.model <- dpmra.model <- dpmdec.model <- c()
mc <- out$mcmc.opt$sig1
nmc <- nrow(mc)
#Nsamp <- 1e4
Nsamp <- 1e3
inds <- sort(sample(1:nmc,Nsamp))
for(j in inds){
#    if(j%%100==0) cat(j,'\n')
    sim <- -astrometry.kepler(mc[j,],bases=bases)$cats
    dra.model <- rbind(dra.model,sim[,1]+dra.obs)
    ddec.model <- rbind(ddec.model,sim[,2]+ddec.obs)
    dplx.model <- rbind(dplx.model,sim[,3]+dplx.obs)
    dpmra.model <- rbind(dpmra.model,sim[,4]+dpmra.obs)
    dpmdec.model <- rbind(dpmdec.model,sim[,5]+dpmdec.obs)
}

###the position and proper motion reconstructed by a linear fit to gost-based observations generated by accounting for reflex motion
obs <- cbind(dra.obs,ddec.obs,dplx.obs,dpmra.obs,dpmdec.obs)
eobs <- out$astrometry[out$astro.index,c('ra_error','dec_error','parallax_error','pmra_error','pmdec_error')]
model <- cbind(dra.model,ddec.model,dplx.model,dpmra.model,dpmdec.model)
ylabs <- c(expression(Delta*alpha*'* [mas]'),expression(Delta*delta*' [mas]'),'parallax [mas]',expression(Delta*mu[alpha]*' [mas/yr]'),expression(Delta*mu[delta]*' [mas/yr]'))
nast <- length(out$astro.index)
trefs <- out$astrometry[out$astro.index,'ref_epoch']
ts <- trefs
dat <- c()
pred <- c()
#cn <- c(expression(Delta*alpha*'*'),expression(Delta*delta),expression(tilde(omega)),expression(mu[alpha]),expression(mu[delta]))
#cn <- c("Delta*alpha*'*'^GDR2","Delta*delta^GDR2","tilde(omega)^GDR2","mu[alpha]^GDR2","mu[delta]^GDR2","Delta*alpha*'*'^GDR3","Delta*delta^GDR3","tilde(omega)^GDR3","mu[alpha]^GDR3","mu[delta]^GDR3")
#cn <- cbind(c("Delta*alpha['*'][2]","Delta*delta[2]","tilde(omega)[2]","mu[alpha][2]","mu[delta][2]"),c("Delta*alpha['*'][3]","Delta*delta[3]","tilde(omega)[3]","mu[alpha][3]","mu[delta][3]"))
cn <- as.character(t(cbind(c("Delta*alpha['*'][2]","Delta*delta[2]","tilde(omega)[2]","mu[alpha][2]","mu[delta][2]"),c("Delta*alpha['*'][3]","Delta*delta[3]","tilde(omega)[3]","mu[alpha][3]","mu[delta][3]"))))
nn0 <- as.character(cn)
nn <- as.expression(cn)
gdrs <- c()
for(j in 1:ncol(obs)){
    inds <- 1:nrow(obs)
    if(j>2 & out$Nastro==3 & gdr1plx=='synthetic'){
        inds <- (1:nrow(obs))[-1]
    }
    dat <- rbind(dat,cbind(trefs[inds],obs[inds,j],eobs[inds,j]))
    gdrs <- c(gdrs,c('GDR2','GDR3'))
    pred <- rbind(pred,cbind(trefs[inds],model[inds,j]))
}
df <- data.frame(t=dat[,1],x=dat[,2],ex=dat[,3],xmin=dat[,2]-dat[,3],xmax=dat[,2]+dat[,3],Data=gdrs,xp=pred[,2],Parameter=nn0)
###add model prediction
dras <- as.numeric(dra.model)
ddecs <- as.numeric(ddec.model)
dplxs <- as.numeric(dplx.model)
dpmras <- as.numeric(dpmra.model)
dpmdecs <- as.numeric(dpmdec.model)
gdr.model <- c(rep('GDR2',Nsamp),rep('GDR3',Nsamp))
xs <- c(dras,ddecs,dplxs,dpmras,dpmdecs)
ns <- as.character(t(replicate(Nsamp,cn)))
df.model <- data.frame(Parameter=ns,x=xs,Data=as.character(replicate(5,gdr.model)))

###
#f1 <- ggplot(df.model,aes(x=Parameter,y=x))#col=GDR
f1 <- ggplot(df.model,aes(x=Parameter,y=x))
#f1 <- f1+geom_boxplot(aes(col=Data),outlier.shape = NA)
#f1 <- f1+geom_violin(aes(col=Data))
f1 <- f1+geom_boxplot(outlier.shape = NA)
#f1 <- f1+geom_violin(col='black')
f1 <- f1 + geom_pointrange(df,mapping=aes(x = Parameter,y = x,ymin=x-ex,ymax=x+ex,col=Data),size=1.2,alpha = 0.8)
if(target=='HD22049'){
    ylim <- c(min(dat[,2])-0.2,max(dat[,2])+0.4)
    f1 <- f1+annotate("text", cn, rep(ylim[1]-0.3,10), label = cn, parse = TRUE ,size = 6)
    f1 <- f1+coord_cartesian(ylim = ylim, clip = "off")
}
if(target=='HD209100'){
    ylim <- c(min(dat[,2])-0.4,max(dat[,2])+0.3)
    f1 <- f1+annotate("text", cn, rep(ylim[1]-1,10), label = cn, parse = TRUE ,size = 6)
    f1 <- f1+coord_cartesian(ylim = ylim, clip = "off")
}
if(target=='UCAC4569-026385'){
    ylim <- c(min(dat[,2])-0.7,max(dat[,2])+0.4)
    f1 <- f1+annotate("text", cn, rep(ylim[1]-0.16,10), label = cn, parse = TRUE ,size = 6)
    f1 <- f1+coord_cartesian(ylim = ylim, clip = "off")
}
if(target=='HD222237'){
    ylim <- c(min(dat[,2])-0.2,max(dat[,2])+0.2)
    f1 <- f1+annotate("text", cn, rep(ylim[1]-0.5,10), label = cn, parse = TRUE ,size = 6)
    f1 <- f1+coord_cartesian(ylim = ylim, clip = "off")
}
f1 <- f1+theme(axis.text.x = element_blank(),text=element_text(size = 20),axis.text.y = element_text (size=15))#axis.title.x = element_text (margin=margin(t=40))
f1 <- f1+labs(x='',y='Astrometric offset [mas or mas/yr]',size=2)
if(target=='HD22049') f1 <- f1+ggtitle(expression(epsilon*' Eridani b'))
if(target=='HD209100') f1 <- f1+ggtitle(expression(epsilon*' Ind A b'))
if(target=='HD222237') f1 <- f1+ggtitle('HD 222237')
if(target=='HD22049' | target=='HD209100') dd <- '/Users/ffeng/Documents/projects/dwarfs/paper/nearest_jupiters/arxiv/'
if(target=='HD222237') dd <- '/Users/ffeng/Documents/projects/dwarfs/paper/HD222237/'
if(target=='UCAC4569-026385') dd <- '/Users/ffeng/Documents/projects/dwarfs/paper/BH/'
fout <- paste0(dd,'model_fit_',target,'.pdf')
cat(fout,'\n')
ggsave(fout,width=8,height=8)
