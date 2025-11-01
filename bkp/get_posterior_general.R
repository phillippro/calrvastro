source('mcmc_func.R')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    f3 <- args[1]
    showtype <- args[2]
    masstype <- args[3]
}else{
#f3 <- 'results/HD4113/HD4113_lpspm_photovaryFALSE_relativityFALSE_Niter2100000_Ncores16_ofac2_Nset5_hg123_transit0_P39813d527_acc0.15_sinI_lnlmax-924.Robj'
#f3 <- 'results/HR8799/HR8799_calFALSE_coplanarFALSE_resonanceFALSE_staTRUE_relativityFALSE_Niter2100000_sppriorFALSE_dr1FALSE_230617_Nset0_hg123+_Nsig5_P290919d132616d59697d21724d3295_Esd1_astro5TRUE_acc0.37_rvcbarycentric_lnpmax-1692_lnlmax-1595.Robj'
#f3 <- 'HR8799/HR8799_calFALSE_coplanarFALSE_resonanceFALSE_staTRUE_relativityFALSE_Niter2200000_sppriorFALSE_dr1FALSE_230617_Nset0_hg123+_Nsig4_P258166d125784d43294d19635_Esd1_astro5TRUE_acc0.25_rvcbarycentric_lnpmax-1654_lnlmax-1576.Robj'
#f3 <- 'results/HD222237/HD222237_lpspm_photovaryFALSE_relativityFALSE_Niter2200000_Ncores16_ofac2_Nset4_hg123_transit0_P22003_acc0.25_sinI_lnlmax-662.Robj'
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter2200000_sppriorFALSE_dr1FALSE_230617_Nset1_hg123_Nsig1_P884_Esd1_astro5TRUE_acc0.28_rvcbarycentric_lnpmax-224_lnlmax-199.Robj'
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_sppriorFALSE_dr1FALSE_230617_Nset1_hg123_Nsig1_P887_Esd1_astro5TRUE_acc0.3_rvcbarycentric_lnpmax-224_lnlmax-199.Robj'
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter2100000_sppriorFALSE_dr1FALSE_230617_Nset1_hg123_Nsig1_P885_Esd1_astro5TRUE_acc0.25_rvcbarycentric_lnpmax-224_lnlmax-199.Robj'
########new results
####calibrate=FALSE
#f3 <- 'results/GaiaBH2/GaiaBH2_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P1310_Esd1_astro5TRUE_acc4.6_rvcbarycentric_priorf1_lnpmax-243_lnlmax-218.Robj'
#f3 <- 'results/GaiaBH2/GaiaBH2_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P1321_Esd1_astro5TRUE_acc15_rvcbarycentric_priorf2_lnpmax-244_lnlmax-219.Robj'
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P885_Esd1_astro5TRUE_acc0.28_rvcbarycentric_priorf1_lnpmax-222_lnlmax-198.Robj'
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P885_Esd1_astro5TRUE_acc0.34_rvcbarycentric_priorf2_lnpmax-222_lnlmax-198.Robj'
###calibrate=TRUE
#f3 <- 'results/GaiaBH2/GaiaBH2_calTRUE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P1333_Esd1_astro5TRUE_acc8.8_rvcbarycentric_priorf2_lnpmax-243_lnlmax-219.Robj'
#f3 <- 'results/GaiaBH2/GaiaBH2_calTRUE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P1324_Esd1_astro5TRUE_acc4.1_rvcbarycentric_priorf1_lnpmax-243_lnlmax-218.Robj'
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calTRUE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P885_Esd1_astro5TRUE_acc0.31_rvcbarycentric_priorf2_lnpmax-223_lnlmax-198.Robj'
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calTRUE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P886_Esd1_astro5TRUE_acc0.33_rvcbarycentric_priorf1_lnpmax-222_lnlmax-198.Robj'

###prior 3
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P886_Esd1_astro5TRUE_acc0.39_rvcbarycentric_priorf3_lnpmax-223_lnlmax-198.Robj'
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calTRUE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P886_Esd1_astro5TRUE_acc0.27_rvcbarycentric_priorf3_lnpmax-223_lnlmax-198.Robj'

#f3 <- 'results/GaiaBH2/GaiaBH2_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P1322_Esd1_astro5TRUE_acc3.6_rvcbarycentric_priorf3_lnpmax-243_lnlmax-219.Robj'
#f3 <- 'results/GaiaBH2/GaiaBH2_calTRUE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P1309_Esd1_astro5TRUE_acc3.6_rvcbarycentric_priorf3_lnpmax-243_lnlmax-219.Robj'

###prior 0
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter2200000_sppriorTRUE_dr1FALSE_230617_Nset1_hg123_Nsig1_P886_Esd1_astro5TRUE_acc0.38_rvcbarycentric_lnpmax-223_lnlmax-199.Robj'
#f3 <- 'results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1100000_sppriorTRUE_dr1FALSE_230617_Nset1_hg123_Nsig1_P885_Esd1_astro5TRUE_acc0.35_rvcbarycentric_lnpmax-224_lnlmax-201.Robj'
#f3 <- 'results/HD209100/HD209100_calFALSE_coplanarFALSE_resonanceFALSE_relativityFALSE_Niter2200000_Ncores16_ofac2_Nset5_hg123_transit0_P15663_Esd1_acc1.5_rvcbarycentric_lnlmax-8150.Robj'#rv+hg23
#f3 <- 'results/HD209100/HD209100_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter2200000_dr1FALSE_230617_Nset5_hg123+_Nsig1_P13993_Esd1_astro5TRUE_acc1.7_rvcbarycentric_priorf2_lnpmax-8195_lnlmax-8167.Robj'#rv+hg23+psf1
f3 <- 'results/HD209100/HD209100_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset5_hg123+_Nsig1_P20781_Esd1_astro5TRUE_acc1.2_rvcbarycentric_priorf2_lnpmax-8494_lnlmax-8464.Robj'#rv+hg23+psf2(cc1)
#showtype <- 'simple'
    showtype <- 'complex'
#    masstype <- 'gra'
#    masstype <- 'spec'
    masstype <- 'evo'
}
load(f3,envir=ee <- new.env())
out <- ee$out
Nsig <- ee$Nsig
tmin <- ee$tmin
mc <- out$mcmc.opt[[paste0('sig',Nsig)]]

if(grepl('GaiaBH2',f3)){
    tab <- read.table('plx_mass_GBH2.txt',header=TRUE)
    gbh2.mgra <- approxfun(tab[,'plx'],tab[,'gramass'],rule=2)
    gbh2.emgra <- approxfun(tab[,'plx'],0.5*(abs(tab[,'err_low'])+abs(tab[,'err_up'])),rule=2)
}

if(grepl('UCAC4569-026385',f3)){
    tab <- read.table('plx_mass_UCAC4569-026385.txt',header=TRUE)
    g3425.mspec <- approxfun(tab[,'plx'],tab[,'spec_mass'],rule=2)
    g3425.mevo <- approxfun(tab[,'plx'],tab[,'evo_mass'],rule=2)
    g3425.emspec <- approxfun(tab[,'plx'],0.5*(abs(tab[,'esmass1'])+abs(tab[,'esmass2'])),rule=2)
    g3425.emevo <- approxfun(tab[,'plx'],0.5*(abs(tab[,'eemass1'])+abs(tab[,'eemass2'])),rule=2)
}

#Npar <- ncol(mc)-2
###modify period
indp <- grep('^per\\d',colnames(mc))
mc[,indp] <- exp(mc[,indp])
c0 <- colnames(mc)
c0[indp] <- paste0('Pd',1:Nsig)
colnames(mc) <- c0
plxs <- out$astrometry[nrow(out$astrometry),'parallax']-mc[,'dplx']
ds <- 1/plxs#kpc

#if(target=='UCAC4569-026385'){
#    out$Mstar <- 2.8
#    out$eMstar <- 0.2
#    mstars <- rnorm(nrow(mc),out$Mstar,out$eMstar)
#}

###calculate mass
mcs <- c()
Tps <- c()
if(any(colnames(mc)=='Mstar')){
    mstars <- mc[,'Mstar']
}else{
    mstars <- rnorm(nrow(mc),out$Mstar,out$eMstar)
}
if(grepl('GaiaBH2',f3) & masstype=='gra'){
    mstars <- gbh2.mgra(plxs)+rnorm(length(plxs),0,gbh2.emgra(plxs))
}
if(grepl('UCAC4569-026385',f3) & masstype=='spec'){
    mstars <- g3425.mspec(plxs)+rnorm(length(plxs),0,g3425.emspec(plxs))
}
if(grepl('UCAC4569-026385',f3) & masstype=='evo'){
    mstars <- g3425.mevo(plxs)+rnorm(length(plxs),0,g3425.emevo(plxs))
}
inds <- which(mstars<0)
if(length(inds)>0)   mstars[inds] <- sample(mstars[-inds],length(inds))

for(j in 1:Nsig){
    Mc <- k2m(mc[,paste0('K',j)],mc[,paste0('Pd',j)],mc[,paste0('e',j)],Ms=mstars,Inc=mc[,paste0('Inc',j)])$ms
    Tp <- M02Tp(M0=mc[,paste0('Mo',j)],T0=tmin,P=mc[,paste0('Pd',j)])
    mcs <- cbind(mcs,Mc)
    Tps <- cbind(Tps,Tp)
}
colnames(Tps) <- paste0('Tp',1:Nsig)
mstars <- t(t(mstars))
colnames(mcs) <- paste0('Mc',1:Nsig)
colnames(mstars) <- 'Mstar'
#cc <- c(paste0('Pd',1:Nsig),paste0('K',1:Nsig),paste0('e',1:Nsig),paste0('Inc',1:Nsig),paste0('omega',1:Nsig),paste0('Omega',1:Nsig),paste0('Mo',1:Nsig))
if(showtype=='simple'){
    dat <- cbind('parallax'=plxs,'distance'=ds,Tps,mc)
}else{
    dat <- cbind(mcs,mstars,'parallax'=plxs,'distance'=ds,Tps,mc)
}
if(grepl('priorf',f3)){
    priortype <- gsub('_.+','',gsub('.+_priorf','',f3))
}else{
    priortype <- 0
}
calibrate <- gsub('_.+','',gsub('.+_cal','',f3))
#fout <- paste0('coplanar',coplanar,'_lowe',lowe,'_resonance',resonance,'_nsig',nsig,'.txt')
target <- gsub('\\/.+','',gsub('results\\/','',f3))
fout <- paste0('results/',target,'/',target,'_Nsig',Nsig,'_prior',priortype,'_mass',masstype,'_calib',calibrate,'.txt')
cat(fout,'\n')
jj <- sample(1:nrow(dat),1e4)
write.table(dat[jj,],file=fout,quote=FALSE,row.names=FALSE)

#fpdf <- gsub('.Robj','_simple.pdf',f3)
fpdf <- gsub('txt','pdf',fout)
cat(fpdf,'\n')
pdf(fpdf,8,8)
par(mfrow=c(2,2))
p1 <- data.distr(ds,lp=mc[,'logpost'],xlab='Distance [kpc]')
p2 <- data.distr(mc[,'Inc1']*180/pi,lp=mc[,'logpost'],xlab='I [deg]')
p3 <- data.distr(mstars,xlab='Mstar [Msun]')
p4 <- data.distr(mcs,lp=mc[,'logpost'],xlab='Mc [Msun]')
dev.off()
