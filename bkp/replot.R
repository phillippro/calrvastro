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
library(foreach)
source('astrometry_function.R')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    fs <- args[1]
}else{
#    fs <- c()
#    fs <- 'HR8799_lpspm_photovaryFALSE_relativityFALSE_Niter2100000_Ncores16_ofac2_Nset0_hg123+_transit0_P163729d123441d31574d19191_acc0.056_sinI_lnlmax-2770'
#    fs <- 'TOI1830_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1200000_dr1FALSE_230617_Nset3_hg123_Nsig2_P10d659_Esd1_astro5TRUE_acc0.54_rvcbarycentric_priorf2_lnpmax'
#    fs <- 'TOI1830_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1100000_dr1FALSE_230617_Nset3_hg123_Nsig2_P10d658_Esd1_astro5TRUE_acc0.16_rvcbarycentric_priorf2_lnpmax'
     fs <- 'HD209100_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1100000_dr1FALSE_230617_Nset5_hg123+_Nsig1_P16761_Esd1_astro5TRUE_acc2.3_rvcbarycentric_priorf2_lnpmax-8200_lnlmax-8172'
}
m22 <- read.table('../data/combined/mamajek22.txt',header=TRUE)
ms <- as.numeric(m22[,'Msun'])
mg <- as.numeric(m22[,'M_G'])
ind <- which(!is.na(ms) & !is.na(mg))
mrl.m22 <- approxfun(ms[ind],mg[ind])
mlow.m22 <- min(ms[ind])
mup.m22 <- max(ms[ind])
stars <- list.dirs('results',recursive=F,full.names=FALSE)

#MPonly <- TRUE
MPonly <- FALSE
moon <- ''
#for(star in stars){
#foreach(star = stars, .combine='rbind') %dopar% {
#fs <- list.files(paste0('results/',star),pattern='',full.name=TRUE)
#ind <- which.max(as.numeric(gsub('.+_lnlmax|\\.Robj','',fs)))
#fs <- fs[ind]
for(f in fs){
file <- gsub('\\+','\\\\+',f)
cat(file,'\n')
    if(!file.exists(file)){
        target <- gsub('.+\\/|_.+','',file)
        f <- list.files(paste0('results/',target),pattern=file,full.name=TRUE)
        file <- f[grepl('Robj',f)][1]
#        file <- gsub('pdf','Robj',f)
    }
    cat('file:',file,'\n')
    load(file)
if(!exists('mstars')){
 if(any(names(par.opt)=='Mstar')){
         mstar <- par.opt['Mstar']
         emstar <- sd(mc[,'Mstar'])
         mstars <- mc[,'Mstar']
     }else{
         mstar <- out$Mstar
         emstar <- out$eMstar
         mstars <- rnorm(nrow(mc),mstar,emstar)
     }
}
#load(file)
    if(!exists('bases')) bases <- rep('natural',length(Popt))
    if(!exists('binary')){
        if(grepl('HD128620',target)){
	    binary <- TRUE
	}else{
	    binary <- FALSE
	}
    }
    fname <- paste0(fname,'_new')
    source('generate_figures.R')
}
#source('calcres.R')
