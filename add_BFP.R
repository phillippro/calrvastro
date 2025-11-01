library(magicaxis)
library(RColorBrewer)
source('mcmc_func.R')
source('periodograms.R')
source('periodoframe.R')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    fs <- args[1]
}else{
    fs <- c()
    fs <- c(fs,'HD216520_natural_offsetTRUE_Niter2000000_Ncores8_ofac2_Nset2_Esd0.2_transit0_P35d7000d155_acc0.11_lnlmax')
}
#MPonly <- TRUE
MPonly <- FALSE
alpha.dP <- 1e-6
files <- fs
for(file in files){
    file0 <- file
    if(!file.exists(file)){
        target <- gsub('.+\\/|_.+','',file)
        file <- list.files(paste0('results/',target),pattern=file0,full.name=TRUE)
        file <- file[grepl('Robj',file)][1]
    }
    if(!file.exists(file)){
        file <- list.files(paste0('/malibu/ffeng/',target),pattern=file0,full.name=TRUE)
        file <- file[grepl('Robj',file)][1]
    }
    cat('file:',file,'\n')
    if(!exists('out') | target!=gsub('.+\\/|_.+','',file)) load(file)
                                        #    Nsig <- 1
    if(target=='HD239960') Nsig <- 3

                                        #load(file)
    if(!exists('bases')) bases <- rep('natural',length(Popt))
    source('global_setting.R')

    if(target=='HD216520'){
        dd <- read.table(paste0('../data/combined/',target,'/',target,'_APF.vels'),header=TRUE)
        ns <- cbind('APF',colnames(dd)[4:ncol(dd)],'MA')
    }
    for(j in 1:nrow(ns)){
        y <- dd[,3+j]
        inds <- which(abs(y-mean(y))<3*sd(y))
        y <- scale(y[inds])
        t <- dd[inds,1]-min(dd[inds,1])
        ey <- rnorm(length(inds),0.1,0.001)
        per <- BFP(t,y,ey,Nma=0,Nar=0,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,renew=TRUE)
        tmp <- list()
        tmp[[ns[j,3]]] <- per
        out$BFP[[ns[j,1]]][[ns[j,2]]] <- tmp
    }
    source('save_output.R')
}
