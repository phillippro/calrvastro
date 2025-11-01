f5 <- 'results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter2200000_sppriorTRUE_dr1FALSE_230617_Nset1_hg123_Nsig1_P885_Esd1_astro5TRUE_acc0.33_rvcbarycentric_priorfTRUE_lnpmax-224_lnlmax-198.Robj'
load(f5)
source('mcmc_func.R')
mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
incs <- mc[,'Inc1']
plxs <- out$astrometry[2,'parallax']-mc[,'dplx']
mstars <- pmfun.spec(plxs)
mcs <- k2m(mc[,'K1'],exp(mc[,'per1']),mc[,'e1'],Ms=mstars,Inc=incs)$ms


fout <- gsub('\\.Robj','_mass.pdf',f5)
cat(fout,'\n')
pdf(fout,8,8)
par(mfrow=c(2,2))
#hist(plxs,xlab='Parallax [mas]')
#hist(incs*180/pi,xlab='Inclination [deg]')
#hist(mstars,xlab='Primary Mass [Msun]')
#hist(mcs,xlab='Companion Mass [Msun]')
data.distr(plxs,mc[,'logpost'],xlab='Parallax [mas]')
data.distr(incs*180/pi,mc[,'logpost'],xlab='Inclination [deg]')
data.distr(mstars,mc[,'logpost'],xlab='Primary Mass [Msun]')
data.distr(mcs,mc[,'logpost'],xlab='Companion Mass [Msun]')
dev.off()

