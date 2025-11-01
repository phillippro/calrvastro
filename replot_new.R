library(viridis)
library(paletteer)
library(magicaxis)
library(orthopolynom)
source('mcmc_func.R')
source('periodograms.R')
source('periodoframe.R')
source('general_function.R')
source('sofa_function.R')
source('timing_function.R')
source('constants.R')
source('astrometry_function.R')
args <- commandArgs(trailing=TRUE)
if(length(args)>0){
   f <- args[1]
}else{
#   f <- 'results/WASP22/WASP22_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityTRUE_Niter1200000_dr1FALSE_230617_Nset3_hg123_Nsig2_P4d2902_Esd1_astro5TRUE_acc0.54_rvcbarycentric_priorf3_lnpmax-365_lnlmax-326.Robj'
#    f <- 'results/HAT-P-7/HAT-P-7_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityTRUE_Niter1100000_dr1FALSE_230617_Nset4_hg123_Nsig2_P2d16296_Esd1_astro5TRUE_acc0.2_rvcbarycentric_priorf1_lnpmax-1131_lnlmax-1086.Robj'
#     f <- 'results/WD0133-11/WD0133-11_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityTRUE_Niter1400000_dr1FALSE_230617_Nset0_hg123+_Nsig1_P8657_Esd1_astro5TRUE_acc0.0064_rvcbarycentric_priorf1_lnpmax-1915_lnlmax-1893.Robj'
#    f <- 'results/GJ65A/GJ65A_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1000_dr1FALSE_230617_Nset1_hg123+_Nsig2_P9568d156_Esd1_astro5TRUE_acc3.4_rvcbarycentric_priorf2_lnpmax-30_lnlmax12.Robj'
#    f <- 'results/GaiaDR31903530403238236032/GaiaDR31903530403238236032_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityTRUE_Niter130000_jitteruniform_240715_Nset0_hg123_Nsig1_P8740_Esd1_astro5TRUE_acc0.0015_rvcbarycentric_priorf1_lnpmax-230_lnlmax829.Robj'
    f <- 'results/GaiaDR31903530403238236032/GaiaDR31903530403238236032_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityTRUE_Niter11000_jitteruniform_240715_Nset0_hg123_Nsig1_P7119_Esd1_astro5TRUE_acc1.2_rvcbarycentric_priorf1_lnpmax-1944_lnlmax-1923.Robj'
}
load(f)
if(!exists('Prefs')) Prefs <- Pref
mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
if(any(colnames(mc)=='dplx')){
    plx.opt <- out$astrometry[out$iref,'parallax']-par.opt['dplx']
    plxs <- out$astrometry[out$iref,'parallax']-mc[,'dplx']
}else{
    plx.opt <- out$astrometry[out$iref,'parallax']
    plxs <- rnorm(nrow(mc),plx.opt,out$astrometry[out$iref,'parallax_error'])
}
if(target=='UCAC4569-026385'){
    mstar <- pmfun.spec(plx.opt)
    emstar1 <- pemfun.spec(plx.opt)#intrinsic uncertainty of the model
    emstar2 <- sd(pemfun.spec(plxs))#uncertainty due to MCMC
    emstar <- sqrt(emstar1^2+emstar2^2)
    mstars <- rnorm(2*nrow(mc),mstar,emstar)
    mstars <- mstars[mstars>0][1:nrow(mc)]
}else{
    if(any(names(par.opt)=='Mstar')){
        mstar <- par.opt['Mstar']
        emstar <- sd(mc[,'Mstar'])
        mstars <- mc[,'Mstar']
    }else{
        mstar <- out$Mstar
        emstar <- out$eMstar
        mstars <- rnorm(2*nrow(mc),mstar,emstar)
        mstars <- mstars[mstars>0][1:nrow(mc)]
    }
}
ofac <- 0.1
source('generate_figures.R')
