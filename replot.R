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
#    fs <- 'results/HD128620/HD128620_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityTRUE_Niter110000_jitteruniform_240715_Nset1_Nrv79_hg123+_Nsig1_P27796_Esd1_dtyperaw_acc9.8_rvcbarycentric_priorf1_lnpmax-743_lnlmax-724.Robj'
#    fs <- 'results/HD128620/HD128620_calTRUE_coplanarFALSE_resonanceFALSE_staFALSE_relativityTRUE_Niter110000_jitteruniform_240715_Nset1_Nrv76_hg123+_Nsig1_P27666_Esd1_dtyperaw_acc14_rvcbarycentric_priorf1_lnpmax-1043_lnlmax-1021.Robj'
    fs <- 'HD128620_calTRUE_coplanarFALSE_resonanceFALSE_staFALSE_relativityTRUE_Niter110000_jitteruniform_240715_Nset1_Nrv76_hg123+_Nsig1_P27666_Esd1_dtyperaw_acc14_rvcbarycentric_priorf1_lnpmax-1043_lnlmax-1021'
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
if(!grepl('Robj',fs)){
file <- gsub('\\+','\\\\+',f)
cat(file,'\n')
    if(!file.exists(file)){
        target <- gsub('.+\\/|_.+','',file)
        f <- list.files(paste0('results/',target),pattern=file,full.name=TRUE)
        file <- f[grepl('Robj',f)][1]
#        file <- gsub('pdf','Robj',f)
    }
}else{
file <- fs
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
source('mcmc_func.R')
source('generate_figures.R')
}
#source('calcres.R')
