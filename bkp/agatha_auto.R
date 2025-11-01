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
time.start <- proc.time()
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    target <- as.character(args[1])
    Niter0 <- as.integer(args[2])
    Esd <- as.numeric(args[3])
    Ncores <- as.integer(args[4])
    coplanar <- as.character(args[5])
    astro.type <- as.character(args[6])
    resonance <- as.logical(args[7])
    relativity <- as.logical(args[8])
    calibrate <- as.logical(args[9])
    priorf <- as.integer(args[10])
}else{
#    target <-'UCAC4569-026385'
#    target <-'WD1202-232'
#    target <-'GaiaDR33221450528288877440'
#    target <-'GaiaDR31992469104239732096'
#    target <-'HD259440'
#    target <-'GaiaDR31903530403238236032'
#    target <-'HIP14754'
#    target <-'HD211847'
#    target <- 'HD154791'
#    target <- 'HIP56662'
#    target <- 'GJ440'
#    target <- 'WD0643-16'
#    target <- 'KIC9472174'
#    target <- 'GaiaDR33971777500966832384'
#    target <- 'GaiaDR31191504471436192512'
#    target <-'KIC12400729'
#    target <- 'GaiaDR35084805635638179584'
#    target <- 'GaiaDR36911950900211768704'
#    target <- 'HD118203'
#    target <- 'HD83443'
#    target <- 'WASP-53'
#    target <- 'WASP-41'
#    target <- 'GaiaDR36911950900211768704'
#    target <- 'GaiaDR34788741548375134336'
#    target <- 'K2-99'
#    target <- 'HD13908'
#    target <- 'HD37605'
#    target <- 'Kepler-432'
#    target <- 'Kepler-419'
#    target <- 'HATS-59'
#    target <- 'Kepler-432'
#    target <- 'TOI-2537'
#    target <- 'WD2149+021'
    target <- 'GD66'
#    target <- 'Pr0211'
#    target <- "WASP-81"
#    target <- 'HD87646'
#    target <- 'K2-99'
#    target <- 'CoRoT-20'
#    target <-'HD118203'
#    target <-'HD187123'
#    target <-'HD83443'
#    target <- 'HD87646'
#    target <-'HIP12469'
#    target <-'CPD-632495'
#    target <- 'TYC3588-11669-1'
#    target <-'HD165341A'
#    target <-'WASP4'
#    target <-'WASP22'
#    target <-'HAT-P-7'
#    target <-'TOI-1227'
#    target <-'HD201091'
#    target <-'KOI-12'
#    target <-'TOI-1227'
#    target <-'Kepler-396'
#    target <-'HD108236'
#    target <-'HAT-P-13'
#    target <-'PDS70'
#    target <- '2MASSJ02000451-0229333'
#    target <- '2MASSJ04412976-7439249'
#    target <- 'HIP32791'
#    target <- 'HIP89039'
#    target <-'WD0133-11'
#    target <-'HD41004'
#    target <-'HD28185'
#    target <-'GJ65A'
#    target <-'GJ65B'
#    target <-'GJ680'
#   target <-'HIP32791'
#   target <-'HIP89039'
#   target <-'GD133'
#    target <- 'HIP89039'
#    target <-'2MASSJ12451043+1217401'
#    target <-'WD1202-232'
#    target <-'HAT-P-7'
#    target <-'HD201091'
#    target <-'HD209100'
#    target <-'HD209100bin'
#    target <-'HD30219'
#    target <-'HD209100bin'
#    target <-'XO3'
#    target <-'WASP22'
#    target <-'TOI1830'
#    target <- 'HD165341A'
#    target <- 'HD253754'
#   target <-'HR8799'
#   target <-'GJ234'
#   target <-'PDS70'
#   target <-'HD209100'
#   target <-'HD22049'
#    Niter0 <- 1e4
    Niter0 <- 1e5
    Esd <- 1
    Ncores <- 8
#    Ncores <- 1
    coplanar <- FALSE
#    coplanar <- TRUE
    astro.type <- 'hg123'
#    astro.type <- 'hg123'
    resonance <- FALSE
#   resonance <- TRUE
#    relativity <- TRUE
    relativity <- FALSE
#    calibrate <- TRUE
    calibrate <- FALSE
    priorf <- 0
#    priorf <- 1
}
ruweDR <- 3
###if ruweDR==2 or 3, we will use ruwe DR2 or DR3 to constrain the orbit
#ruweDR <- 3
#jitter <- 'gauss'
jitter <- 'uniform'
#jitter <- 'fixed'
gmag <- 9.63#this value is to determine which calibration parameters to be used for calibrating astrometric catalogs
###for WD targets, priorf=0: fixed par; priorf=1: gaussian prior; priorf=2: uniform prior
Nfix <- 0
#priorf <- FALSE
stability <- FALSE
#stability <- TRUE
if(target=='HR8799' | target=='PDS70') stability <- TRUE
####targets without red noise modeling
target.except  <- c('GJ245.1','HD24040','GJ559A','HD129814','HD218935','HD26736','HD26990','HD39094','HD4203','HD4208','HD48122','HD55696','HD66428','GJ1089')
moon <- ''
#    moon <- 'b1c2'
#   basis <- 'linear2'##transit basis for Keplerian parameters
#    photovary <- TRUE
gammaf <- FALSE
#gammaf <- TRUE
astro5 <- TRUE
par.fix <- NULL
#rel.type <- 'sp'#PA and separation
rel.type <- 'rd'#dRA and dDEC
#if(any(target==c('CPD-632495','TYC3588-11669-1'))) par.fix <- c('per1','e1','omega1','Mo1','Inc1','Omega1')
if(any(target==c('CPD-632495','TYC3588-11669-1'))) par.fix <- c('per1','e1','omega1','Mo1')
if(any(target==c('WD1202-232','WD2105-82','PDS70','WD0643-16')) & priorf==0) par.fix <- c('Mstar')
#if(any(target==c('HD259440','HD30219'))) par.fix <- c('per1')
if(any(target==c('HD30219','HATS-59','CoRoT-20'))) par.fix <- c('per1')
#astro5 <- FALSE
rvc.type <- 'barycentric'#companion RV type
dP <- 0
if(FALSE) source('check_transit.R')
#use.gdr1 <- TRUE
use.gdr1 <- FALSE
laplace <- c(1,1/2,1/4,1/8,1/32)
if(target=='PDS70') laplace <- c(1,2,1/2)
photovary <- FALSE
Nmax <- 2
etaH <- etaG <- rep(NA,10)
bc <- FALSE
if(any(target==c('TYC3588-11669-1','HD253754','GaiaBH2','UCAC4569-026385','HD150554','HD13507','HD133621','HD239960','HD142','HD10790','HIP35305','GJ9476','CPD-632495','HD201091','WD2105-82','WD1202−232','PDS70','HD259440','GJ65A','GJ65B','HD41004','GJ680','HIP32791','GaiaDR31903530403238236032','WD0643-16','KIC5307780','KIC7691553','KIC9161428','KIC7273033','KIC12400729','KIC9947924','KIC8669092','WD2149+021'))) etaH[1] <- etaG[1] <- 0
if(any(target=='HD41004')){
    dg <- 3.7
    etaG[2] <- (1+10^(0.074*dg))/(10^(0.4*dg)-10^(0.074*dg))
    etaH[2] <- 0
}
basis <- 'natural'
sim <- FALSE
#Tmin <- 1000#minimum orbital period for Gaia and Hipparcos catalog data to be considered instataneous
Tmin <- 1e9
barycor <- TRUE
#barycor <- FALSE
Nsim <- 1e4
tycf <- FALSE
#binary <- 'HD128621'
#binary <- 'HD41004B'
binary <- ''
if(target=='GJ65A') binary <- 'GJ65B'
if(target=='GJ65B') binary <- 'GJ65A'
#if(target=='HD41004') binary <- 'HD41004B'
if(is.na(as.numeric(moon))) moon <- ''
###out is a list store most important variables for fitting and output
out <- list()
ofac <- 0.5
#source('adjust_astrometry_file.R',local=TRUE)
astrometry <- 0
Nastro <- 0
if(Nfix==0){
    Nmin <- 0
    Nmax <- 8
}else{
    Nmax <- Nmin <- Nfix
}
fpar <- paste0('pars/',target,'.par')
if(target=='HR8799'){
   fpar <- paste0('pars/',target,'_c',coplanar,'_l',as.logical(Esd<1),'_r',resonance,'_n',Nfix,'.par')
   if(!file.exists(fpar)){
      fpar <- paste0('pars/',target,'_n',Nfix,'.par')
   }
   cat(fpar,'\n')
}
if(file.exists(fpar)){
    np <- length(grep('^per|^dP',read.table(fpar)[,1]))
    Nmin <- Nmax <- np
}
if(target=='HD190360') Nmin <- Nmax <- 1
cat('Nmin=',Nmin,';Nmax=',Nmax,'\n')

###If relativity astrometry for multiple planets are analyzed, there is an option to specify the order of the relativity astrometry data files.
out$rel.order <- 1:Nmax# default order
#out$rel.order <- c(1)
basis2 <- 'natural'
#basis2 <- 'linear1'#non-transit basis for Keplerian parameters
cov.astro <- NULL
target <- gsub(' ','',target)
###low cadence for test
if(Niter0<1e6){
    ofac <- 0.1
}
if(Niter0<=1e4){
    ofac <- 0.01
}
#star.except <- c('GJ551','GJ699')
star.except <- c('GJ551')
offset <- TRUE
if(any(target==star.except)) offset <- FALSE
#save.memory <- FALSE
save.memory <- TRUE
ind.transit <- 0#'auto'
#ind.transit <- 0
cat('target:',target,';Nmax=',Nmax,';Niter=',Niter0,';Ncores=',Ncores,';ofac=',ofac,';Esd=',Esd,';basis=',basis,'priorf=',priorf,'\n')

t00 <- proc.time()
cat('start time:',t00,'\n')

Niter <- Niter0
##output data structure as a list
out$relativity <- relativity
##register cores
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
#trace(parallel:::sendMaster, at = 3L, tracer = quote({ str(list(what = what)) }))

###derive moon parameters
ss <- unlist(strsplit(moon,''))
if(length(ss)%%2!=0) stop('Error: wrong moon argument format!')
out$moon <- lapply(1:10,function(i) c(0,0,0))#default moon list
out$Nm <- 0
if(moon!=''){
    Ns <- as.integer(ss[(1:length(ss))%%2==0])
    ps <- ss[(1:length(ss))%%2!=0]
###sorting moon according to their associated planets
    ix <- sort(ps,index.return=TRUE)$ix
    ps <- ps[ix]
    Ns <- Ns[ix]
    out$Nm <- sum(Ns)
    Ncum <- cumsum(Ns)
}
if(out$Nm>0){
    for(k in 1:length(ps)){
        p <- ps[k]
        ind <- match(p,c('b','c','d','e','f'))
        if(k==1){
            ind1 <- 1
            ind2 <- Ncum[1]
        }else{
            ind1 <- Ncum[k-1]+1
            ind2 <- Ncum[k]
        }
        out$moon[[ind]] <- c(Ns[k],ind1,ind2)
    }
}

##load data
##set up
cat('Step 1: Global settings\n')
source('global_setting.R')
cat('\nStep 2: Load data & calibration\n')
source('load_data.R')
###calibrate astrometric data
source('astro_calibration.R')
##barycentric correction & best noise model
cat('\nStep 3: PEXO-based barycentric correction\n')
if(barycor | length(out$tiall)>0){
#if(barycor){
    cat('Barycentric correction using corrected astrometry!')
    source('barycor.R')
}

cat('\nStep 4: Noise model comparison\n')
if(length(out$ins.rv)>0){
    source('model_selection.R')
}else{
    nqp <- rep(0,3)
}
##0 to Nmax-planet model fit using DRAM
cat('\nStep 5: Signal selection through DRAM\n')
source('dram.R')

##analysis of the MCMC results
cat('\nStep 6: Analysis of MCMC results\n')
source('analysis_mcmc.R')

###simulation
sim <- TRUE
ts <- c(tmin,tmax,out$tiall[,1])
if(length(out$astrometry)>0){
    ts <- c(ts,out$astrometry[,'ref_epoch'])
}
if(length(out$data.ref)>0){
    ts <- c(ts,out$data.ref[,1])
}
if(length(out$data.epoch)>0){
    for(i in out$ins.epoch){
        ts <- c(ts,out$data.epoch[[i]][,"BJD"])
    }
}
if(any(names(out)=='rel')){
    tt <- c()
    for(n1 in names(out$rel)){
        for(n2 in names(out$rel[[n1]])){
            tt <- c(tt,out$rel[[n1]][[n2]][['BJD']])
        }
    }
    ts <- c(ts,tt)
}
if((max(ts)-min(ts))<max(Popt)){
    dT <- max(Popt)-(max(ts)-min(ts))
    ts <- c(ts,min(ts)-dT)
}
if(any(names(out)=='gost')) ts <- c(ts,out$gost[,'BJD'])
tsim <- seq(min(ts),max(ts),length.out=Nsim)
#tsim <- trv.all
if(barycor){
    if(!is.null(out$tiall)){
        source('barycor.R')
    }else{
        source('barycor_hipgaia.R')
    }
}

mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
##save in case of later errors
source('save_output.R')

##BFPs for data and proxies
cat('\nStep 6: Calculate BFPs\n')
#t1 <- proc.time()
#if(out$Nrv>0){
if(!is.null(out$ins.rv)){
    source('calBFPs_parallel.R')
}

##calculate MP and WP
cat('\nStep 7: Calculate MP and WP\n')
#if(Nsig>0 & out$Nrv>0){
#if(Nsig>0 & !is.null(out$ins.rv)){
if(FALSE){
    source('MP_WP.R')
}

##calculate MP and WP
cat('\nStep 8: Save data\n')
source('save_output.R')

###add plxs
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

##plot
cat('\nStep 9: plot results\n')
source('generate_figures.R')

t11 <- proc.time()
cat('end time:',t11,'\n')
dth <- (t11-t00)[3]/3600
cat('Total time consumed:',dth,'hour\n')
