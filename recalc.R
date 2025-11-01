library(magicaxis)
source('mcmc_func.R')
source('periodograms.R')
source('periodoframe.R')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    f <- args[1]
}else{
#    f <- 'HIP22762_natural_Nmax6_Niter2000000_Ncores16_ofac2_Nset5_Esd0.1_transit0_P17d106_acc0.0052_lnlmax'
    f <- 'HD190007_natural_Nmax6_Niter2000000_Ncores8_ofac2_Nset2_Esd0.2_transit0_P12d29_acc0.33_lnlmax'
}
files <- f
##files <- c('HIP12961/HIP12961_Nmax6_Niter1000000_Ncores10_ofac2_Esd0.1_P57_acc0.13_lnlmax-268.Robj')
for(file in files){
    ##if(!file.exists(file)) file <- paste0('results/',file)
    if(!file.exists(file)){
        target <- gsub('.+/|_.+','',file)
        file <- list.files(paste0('results/',target),pattern=file,full.name=TRUE)
        file <- file[grepl('Robj',file)][1]
    }
    cat('file:',file,'\n')
    load(file)
    Ncores <- 8
    if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
    ofac <- 1
    save.memory <- TRUE
    ##BFPs for data and proxies
    fmax <- 1/1.01
    fmin <- 1/1e5
    cat('\nStep 6: Calculate BFPs\n')
    source('calBFPs_parallel.R')
    ##calculate MP and WP
    cat('\nStep 7: Calculate MP and WP\n')
    if(Nsig>0){
        source('MP_WP.R')
    }
    ##calculate MP and WP
    cat('\nStep 8: Save data\n')
    source('save_output.R')
    ##plot
    cat('\nStep 9: plot results\n')
    source('generate_figures.R')
}
