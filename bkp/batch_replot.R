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
Ncores <- 16
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    fs <- args[1]
}else{
    fs <- c()
    fs <- 'HD125607_lpspm_photovaryFALSE_relativityFALSE_Niter1200000_Ncores16_ofac2_Nset1_htg_transit0_P579_acc16_sinI_lnlmax'
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
foreach(star = stars, .combine='rbind') %dopar% {
fs <- list.files(paste0('results/',star),pattern='',full.name=TRUE)
ind <- which.max(as.numeric(gsub('.+_lnlmax|\\.Robj','',fs)))
file <- gsub('\\+','\\\\+',fs[ind])
cat(file,'\n')
    if(!file.exists(file)){
        target <- gsub('.+\\/|_.+','',file)
        f <- list.files(paste0('results/',target),pattern=file,full.name=TRUE)
        file <- f[grepl('Robj',f)][1]
#        file <- gsub('pdf','Robj',f)
    }
    cat('file:',file,'\n')
    load(file)
#load(file)
    if(!exists('bases')) bases <- rep('natural',length(Popt))
    if(!exists('binary')){
        if(grepl('HD128620',target)){
	    binary <- TRUE
	}else{
	    binary <- FALSE
	}
    }
#    fname <- paste0(fname,'_new')
    source('generate_figures.R')
}
#source('calcres.R')
