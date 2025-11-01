##I have added the if(noise.model=='GP') in section 2.3
library(MASS)
#library(coda)
library(foreach)
library(doMC)
library(doParallel)
library(parallel)
#library(SnowballC)
#library(Rmpi)
#library(snow)
####
source('mcmc_func.R',local=TRUE)
source('periodograms.R',local=TRUE)
source('periodoframe.R',local=TRUE)
#########examples of execution
#####****Rscript mcmc_red.R GJ14_SOPHIE_TERRA_1AP1w0ro2ab6ap 1000 1 1 0.2 1e-3 200 ARMA01 logP 4 1 data D 1 mt 0
##########################################################################
####part I: preparation
##########################################################################
cat('\n set up!\n\n')
period.optimize <- FALSE
#period.optimize <- TRUE
source('set_up.R',local=TRUE)#set up parameters and read data 
##simulating data sets
IDsim <- 'nosim'
cat('\n prepare parameters!\n\n')
source('prepare_par.R',local=TRUE)##set parameters prior, initial values, initial covariance
cat('Pmin=',Pmin,'\n')
cat('Pmax=',Pmax,'\n')
Pmax0 <- Pmax
Pmin0 <- Pmin
##########################################################################
####part II: run MCMC and save the results
##########################################################################
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
if(Nbin.per>1){
   Nr <- 1
}else{
   Nr <- 1
}
for(nr in 1:Nr){
    t1 <- proc.time()
#    if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
    cat('Period range:',Pmin,'d to',Pmax,'d\n\n')
###################
######Find signals
###################
#    cat('node:',Sys.info()[['nodename']],'\n\n')
#    source('set_up.R',local=TRUE)#set up parameters and read data 
    cat('\nFind',Np,'signals!\n\n')
    if(Np>1){
        fixP <- TRUE
    }else{
        fixP <- FALSE
    }
    quantify <- FALSE
    cat('\n Run MCMC!\n\n')
    if(Niter>1e4 & Np>0){
#    if(Niter>1e7){
       chain.type <- 'adapt'
#       chain.type <- 'parallel'
    }else{
       chain.type <- 'parallel'
    }
    cat('chain.type=',chain.type,'\n')
    source('prepare_par.R',local=TRUE)
#    cat('Pini=',Pini,'d\n')
#    source('mcmc_func.R',local=TRUE)
    source('gen_mcmc.R',local=TRUE)
    cat('\n convergency test!\n\n')
    source('convergence_test.R',local=TRUE)
    cat('\n generating figures!\n\n')
    source('figure_gen.R',local=TRUE)
###################
##quantify signal
###################
#    if(conv & Np>0){
#    if(Np>0){
    if(Np>=0){
        for(kk in 1:1){
        quantify <- TRUE
        cat('\n quantify',Np,'signals',kk,'times!\n\n')
        startvalue <- par.opt
        cat('par.opt=',par.opt,'\n')
        #cov.start <- cov.start*1e-3
        tem0 <- tem
        tem <- 1
###output
        cat('\n Run MCMC!\n\n')
#        Nbin.per <- nbin.per <- 1
        reprepare <- FALSE
        fixP <- FALSE
        source('prepare_par.R',local=TRUE)
        source('mcmc_func.R',local=TRUE)
        chain.type <- 'talk'
        source('gen_mcmc.R',local=TRUE)
        cat('\n convergency test!\n\n')
        source('convergence_test.R',local=TRUE)
        source('figure_gen.R',local=TRUE)
        if(conv & Nc2p>0) break()
        }
        tem <- tem0
        if(!conv) break()
#####renew parameters for next MCMC run
        if(Np>0){
            source('isolate_signals.R',local=TRUE)
        }
        Np <- Np+1
    }else{
        break()
    }
}

