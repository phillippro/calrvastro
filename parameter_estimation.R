library(MASS)
#library(coda)
library(foreach)
library(doMC)
library(doParallel)
library(parallel)
#library(Matrix)
rm(list=ls())
source('mcmc_func.R')
source('periodograms.R')
source('periodoframe.R')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    estimation <- as.logical(args[1])
    Nbin.per1 <- as.integer(args[2])
    nbin.per1 <- as.integer(args[3])
    Niter1 <- as.integer(args[4])
    add.par0 <- as.character(args[5])
    trans.id <- as.character(args[6])
    Ncores0 <- as.integer(args[7])
    fs <- as.character(args[8:length(args)])
}else{
     estimation <- TRUE
     Nbin.per1 <- 1
     nbin.per1 <- 1
     Niter1 <- 1e3
     add.par0 <- 'NA'#E,...
     trans.id <- 'self'###'self','1AP1','C1AP1','CCF','1AP1+KECK','CCF+KECK','1AP1+KECK','C1AP1+KECK'
     Ncores0 <- 4
     fs <- c('keppure_priormt_poly20_Ndata831_quantifyTRUE_1per1_Nw1_HD209100_HARPS_modified_ind0_2planet_ARMA02_Nsamp4000000_tem1_acc0.13_pretem1P17.9d6128.3d_negLmax1295')
}
####example:
##Rscript parameter_estimation.R TRUE 1 1 1000 QP mixFALSE_mikko_data_quantifyTRUE_1per1_Pmin1Pmax3604_Nw1ABTRUE_DtypeDforce6_modedata_HD10700_HARPS_CCF_ro_nosim_1planet_ARMAp0q4_beta8.19_Nsamp5000000_tem1_acc0.085
#Nbin.per0 <- 1
#nbin.per0 <- 1
f1 <- gsub('_ind0.+','',fs)
f1 <- gsub('.+Nw\\d_','',f1)
target <- f1[1]
ins0 <- c('HARPS','HARPN','SOPHIE','KECK','resGP2_sim')
cat('f1[1]=',f1[1],'\n')
star <- gsub('_.+','',f1[1])
if(grepl('MP_Mus',f1[1])){
   star <- target <- 'MP_Mus'
}
if(grepl('PlSystem',f1[1])){
   star <- target <- 'PlSystem2'
}
cat('star=',star,'\n')
cat('target=',target,'\n')
###
folder0 <- paste0('/car-data/ffeng/dwarfs/output/',target,'/')
if(!file.exists(folder0)){
    folder0 <- paste0('/car-data/ffeng/dwarfs/output/',gsub('_.+','',target),'/')
}
if(!file.exists(folder0)){
    folder0 <- paste0('../output/',star,'/')
}
cat('folder0=',folder0,'\n')
cat('folder0=',folder0,'\n')
newpar1 <- c()
Np1 <- 0
#Pbins <- c()
for(i in 1:length(fs)){
    file <- fs[i]
    if(!grepl('Robj',file)){
          file <- paste0(file,'.Robj')
    }
    fname <- paste0(folder0,file)
    cat('fname=',fname,'\n')
    if(!file.exists(fname)){
       folder0 <- gsub(target,star,folder0)
       fname <- paste0(folder0,file)
    }
    if(!file.exists(fname) & !grepl('output_bkp',folder0)){
       folder0 <- gsub('output','output_bkp',folder0)
       fname <- paste0(folder0,file)
    }
    cat('file name:',fname,'\n')
    load(fname)
###extract kep.type
    kep.type0 <- kep.type
####extract prior.type
    prior.type0 <- prior.type
    Nw0 <- Nw
    Npoly0 <- Npoly
    Npoly.sub0 <- Npoly.sub
#    if(exists('P.bin')){
#         Pbins <- cbind(Pbins,P.bin)
#    }
    Np1 <- Np1+Np
    ind <- which.max(post.out)
#    ind <- which.max(loglike.out)
    par.opt <- mcmc.out[ind,]
    if(Np>0){
        newpar1<- c(newpar1,par.opt[1:(Np*Nkeppar)])
    }
}
if(prior.type0!='e0'){
   if(grepl('E',add.par0)){
       add.par0 <- gsub('E','',add.par0)
   }
}
if(grepl('E',add.par0)){
        Nt <- Np1*5
	tmp <- rep(NA,Nt)
	for(j in 1:Nt){
	    if(j%%5==1 | j%%5==2 | j%%5==4){
	        np <- ceiling(j/5)-1
		if(j%%5==4){
            	    ind <- 3
		}else{
                    ind <- j%%5
		}
	        tmp[j] <- newpar1[np*3+ind]
	    }else if(j%%5==3){
	        tmp[j] <- 0.1
	    }else{
                tmp[j] <- 0
	    }
	}
        newpar1 <- tmp
	prior.type0 <- 'mt'
}
kep.type0 <- 'pure'
if(grepl('AK',add.par0)){
    Nkeppar <- length(newpar1)/Np
    Nkeppar0 <- Nkeppar+2
    ind.old <- 1:Nkeppar
    ind.new <- (Nkeppar+1):Nkeppar0
    Nt <- Nkeppar0*Np1
    tmp <- rep(NA,Nt)
    for(j in 1:Nt){
        if(any(j%%Nkeppar0==ind.old)){
            np <- ceiling(j/Nkeppar0)-1
            ind <- j%%Nkeppar0
            tmp[j] <- newpar1[np*Nkeppar+ind]
        }else if(j%%Nkeppar0==Nkeppar+1){
            tmp[j] <- log(5000)#logtau
        }else{
            tmp[j] <- 1000#ta
        }
    }
    newpar1 <- tmp
    kep.type0 <- 'AK'
}
if(grepl('QP',add.par0)){
#    Nkeppar <- length(newpar1)/Np
    Nkeppar <- length(newpar1)/Np1
    Nkeppar0 <- Nkeppar+2
    ind.old <- 1:Nkeppar
    ind.new <- (Nkeppar+1):Nkeppar0
    Nt <- Nkeppar0*Np1
    tmp <- rep(NA,Nt)
    for(j in 1:Nt){
        if(any(j%%Nkeppar0==ind.old)){
            np <- ceiling(j/Nkeppar0)-1
            ind <- j%%Nkeppar0
            tmp[j] <- newpar1[np*Nkeppar+ind]
        }else if(j%%Nkeppar0==Nkeppar+1){
            tmp[j] <- log(5000)#logtq
        }else{
            tmp[j] <- 1#aq
        }
    }
    newpar1 <- tmp
    kep.type0 <- 'QP'
}
if(grepl('AKQP',add.par0)){
    kep.type0 <- 'AKQP'
}
if(Np1>0){
   par.opt0 <- c(newpar1,par.opt[(Np*Nkeppar+1):length(par.opt)])
}else{
   par.opt0 <- c(newpar1,par.opt[1:length(par.opt)])
}
par.opt2 <- par.opt <- par.opt0
epoch.type <- 'duration'
Nepoch <- 1
period.optimize <- TRUE
if(grepl('mix',file)){
    f1 <- gsub('mix','',file)
    mix<- as.logical(gsub('_.+','',f1))
}else if(!grepl('Nw1',file)){
    mix<- TRUE
}else{
    mix <- FALSE
}
if(!grepl(paste0('^',target),file)){
    f1 <- gsub('.+modedata_','',file)
}else{
    f1 <- file
}
#id <- gsub('HD020794','HD20794',id)
if(!exists('ids')){
    ids <- id
}
if(!exists('noise.models')){
   if(noise.model=='ARMA'){
      noise.models <- paste0(noise.model,p,q)
   }else{
      noise.models <- noise.model
   }
}
ids0 <- ids
noise.models0 <- noise.models
###
par.opt <- par.opt2
rm(par.opt2)
f1 <- gsub('.+force','',file)
Ndata0 <- length(trv[!is.na(trv)])
Inds <- as.integer(gsub('_.+','',gsub('.+ind','',file)))
if(grepl('HD020794',id)){
    ids0 <- gsub('HD020794','HD20794',ids0)
}
######find the optimal initial period
var.keep <- c('Npoly0','Ncores0','noise.models0','kep.type0','Nw0','fs','prior.type0','Niter1','Nbin.per1','nbin.per1','tem','inicov','ids0','par.opt','estimation','data.version','mode','Nepoch','period.optimize','epoch.type','Np1','P.bin','pmax.res','res.glst.all','trv.all','RV.all','eRV.all','res.all','res','res.glst','eRV','RV','trv','Npoly.sub0')
rm(list= ls()[!(ls() %in% var.keep)])
source('mcmc_func.R',local=TRUE)
source('periodograms.R',local=TRUE)
commandArgs <- function(trailingOnly=TRUE) c()
source('set_up.R')#set up parameters and read data
IDsim <- 'nosim'
##########################################################################
####parameter estimation
##########################################################################
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
cat('Ncores=',Ncores,'\n')
if(!estimation){
    Np1 <- Np1+1
    Np <- Np1
    fixP <- TRUE
    quantify <- FALSE
    Nbin.per <- Nbin.per1
    nbin.per <- nbin.per1
    cat('Np=',Np,'\n\n')
    source('prepare_par.R',local=TRUE)##general plows
    source('isolate_signals.R',local=TRUE)
    source('prepare_par.R',local=TRUE)
    if(Niter<1e5){
        chain.type <- 'parallel'
    }else{
#	chain.type <- 'adapt'
	chain.type <- 'parallel'
    }
    tem <- 1
    cat('Pini=',Pini,'d\n')
    source('gen_mcmc.R',local=TRUE)
    cat('\n convergency test!\n\n')
    source('convergence_test.R',local=TRUE)
    source('figure_gen.R',local=TRUE)
    rm(plows)
}
Np <- Np1
fixP <- FALSE
quantify <- TRUE
tem0 <- tem
tem <- 1
Nbin.per <- nbin.per <- 1
source('prepare_par.R',local=TRUE)##set parameters prior, initial values, initial covariance
chain.type <- 'talk'
if(Niter<1e5){
   chain.type <- 'parallel'
}
source('gen_mcmc.R',local=TRUE)
cat('\n convergency test!\n\n')
source('convergence_test.R',local=TRUE)
source('figure_gen.R',local=TRUE)
