library("plotrix")
source('OrbitFunction.R')
source('../../pexo/code/timing_function.R')
source('../../pexo/code/general_function.R')
source('../../pexo/code/sofa_function.R')
#if(!exists('par0s') | TRUE){
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    version <- as.integer(args[1])
}else{
    version <- 0
}
set.seed(100)
fs <- c('HD131664_abs1_astro8_Niter5000000_Ncores8_ofac2_Nset2_Esd1_transit0_P1984_acc0.51_lnlmax',
'HD16160_natural_offsetTRUE_Niter15000000_Ncores4_ofac2_Nset6_Esd1_transit0_P30213_acc0.17_lnlmax',
'HD161797_natural_astro8_Niter11000000_Ncores4_ofac2_Nset3_Esd1_transit0_P28263_acc0.07_lnlmax',
'HD190360_natural_astro8_Niter10000000_Ncores4_ofac2_Nset5_Esd1_transit0_P2851_acc1.8_lnlmax',
'HD190406_abs1_astro8_Niter4000000_Ncores4_ofac2_Nset3_Esd1_transit0_P22370_acc2.9_lnlmax',
'HD39587_natural_astro8_Niter15000000_Ncores4_ofac2_Nset3_Esd1_transit0_P5155_acc1.6_lnlmax',
'HD4747_natural_astro8_Niter10000000_Ncores4_ofac2_Nset3_Esd1_transit0_P12649_acc2.1_lnlmax',
'HIP2552_natural_astro8_Niter10000000_Ncores4_ofac2_Nset1_Esd1_transit0_P5655_acc0.03_lnlmax'
)
etaH <- etaG <- 1
fs <- sort(fs)
###modifiy input file names
stars <- gsub('_.+','',fs)
ff <- c()
for(k in 1:length(stars)){
    star  <- stars[k]
    ff <- c(ff,list.files(path=paste0('results/',star),pattern=paste0(fs[k],'.+simple.Robj'),full.name=TRUE)[1])
}
if(version==1){
      fs <- ff[1:4]
      stars <- stars[1:4]
}else if(version==2){
      fs <- ff[5:8]
      stars <- stars[5:8]
}else{
      fs <- ff
}
planets <- read.csv('../data/code/SigClassify/ranking_complex2.csv')
n0 <- gsub(' ','',as.character(planets[,'Name']))
n1 <- gsub(' ','',as.character(planets[,'ID']))
n2 <- gsub(' ','',as.character(planets[,'StarKnown']))
eMstars <- Mstars <- c()
for(s in stars){
    s1 <- gsub('GJ','GL',s)
    ind <- which(n0==s | n1==s | n2==s | n0==s1 | n1==s1 | n2==s1)
    if(s=='HD131664'){
        Mstar <- 1.060
        eMstar <- 0.129
    }
    if(s=='HD190406'){
        Mstar <- 1.080
        eMstar <- 0.137
    }
    if(s=='HD16160'){
        Mstar <- 0.780
        eMstar <- 0.091
    }
    if(s=='HD161797'){
###First Results from the Hertzsprung SONG Telescope: Asteroseismology of the G5 Subgiant Star {\ensuremath{\mu}} Herculis by grundahl17
        Mstar <- 1.11
        eMstar <- 0.01
    }
    if(s=='HD182488'){
        Mstar <- 0.930
        eMstar <- 0.113
    }
    if(s=='HD190360'){
        Mstar <- 0.980
        eMstar <- 0.120
    }
    if(s=='HD39587'){
        Mstar <- 1.100
        eMstar <- 0.134
    }
    if(s=='HD42581'){
        Mstar <-  0.544
        eMstar <-  0.041
    }
    if(s=='HD4747'){
        Mstar <- 0.910
        eMstar <- 0.114
    }
    if(s=='HIP2552'){
        Mstar <- 0.530
        eMstar <- 0.083
    }
    if(s=='HD131664'){
        Mstar <- 1.060
        eMstar <- 0.129
    }
    if(s=='HD190406'){
        Mstar <- 1.080
        eMstar <- 0.137
    }
    cat(s,'\n')
    cat('stellar mass=',Mstar,'Msun\n')
    cat('stellar mass error=',eMstar,'Msun\n\n')
    Mstars <- c(Mstars,Mstar)
    eMstars <- c(eMstars,eMstar)
}

###collect all data needed for the plot
Nmc <- 0#Number of Monte Carlo orbits
ins.tot <- list()
res.sig <- res.tot <- list()
eRV.tot <- list()
RV.tot <- list()
nqp.tot <- JD.tot <- list()
Nrv.tot <- cov.tot <- astrometry.tot <- list()
id.tot <- list()
par0.tot <- par1.tot<- list()
ffs <- fs
starss <- stars
Mstarss <- Mstars
eMstarss <- eMstars
dras <- ddecs <- dpmras <- dpmdecs <- c()
for(j3 in 1:length(ffs)){
    star <- starss[j3]
    cat('load ',ffs[j3],'\n')
#    load(ffs[j3],env=e0<-new.env())
    load(ffs[j3])
    star <- starss[j3]
    source('mcmc_func1.R')
    par1 <- par0 <- par.opt
    par0['K1'] <- 0
    par1['logJ.gaia'] <- log(0.743)
    out$Mstar <- Mstarss[j3]
    out$eMstar <- eMstarss[j3]
    out$plx <- out$astrometry$parallax[2]
    out$Nm  <- 0
    out$Nrv <- nrow(out$all)
    out$Nastro <- 2
    out$relativity <- FALSE
    out$Nrvc <- 0
    etaG <- etaH <- 1
    ll1 <- loglikelihood(par.opt)$llastro
   plx0 <- out$plx
#    out$plx <- 1000/(1000/plx0-out$astrometry[1,'radial_velocity']/4.74047/206265)
#    out$plx <- plx0+0.029#correct plx
#    ll0 <- loglikelihood(par0)$llastro
#    ll2 <- loglikelihood(par.opt)$llastro
#    out$astrometry[,'ref_epoch'] <- out$astrometry[,'ref_epoch']+c(0.13,0.14)
    msini <- K2msini.full(par.opt['K1'],exp(par.opt['per1']),par.opt['e1'],Ms=out$Mstar)
    Mp <- msini$ms/sin(par.opt['Inc1'])
if(star=='HD131664'){
    teff1 <- 6000
    teff2 <- 1200
}
if(star=='HD16160'){
    teff1 <- 5500
    teff2 <- 2850
}
if(star=='HD161797'){
    teff1 <- 5560
#https://arxiv.org/abs/1701.03365
    teff2 <- 3100
}
if(star=='HD190360'){
    teff1 <- 5781
#https://arxiv.org/abs/1511.03197
    teff2 <- 149.5
}
if(star=='HD39587'){
    teff1 <- 6000.0
#
    teff2 <- 3030
}
if(star=='HD190406'){
    teff1 <- 6000
    teff2 <- 1680
}
if(star=='HD4747'){
    teff1 <- 5500.0
    teff2 <- 1700
}
if(star=='HIP2552'){
    teff1 <- 4000
    teff2 <- 3100
}
###non-cooling-model-mass from: http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
cat('star=',star,'\n')
cat('teff1=',teff1,'\n')
cat('teff2=',teff2,'\n')
    zeta <- relFlux(teff1,teff2)
    zetaG <- zeta$fG
    zetaH <- zeta$fH
    cat('zetaG=',zetaG,'\n')
    cat('zetaH=',zetaH,'\n')
    alphaG <- (Mp/(out$Mstar+Mp))^2.5*zetaG
    alphaH <- (Mp/(out$Mstar+Mp))^2.5*zetaH

    cat('alphaG=',alphaG,'\n')
    cat('alphaH=',alphaH,'\n')
#    cat('alpha*out$astrometry[,parallax]=',alpha*out$astrometry[,'parallax'],'\n')
#    cat('out$astrometry[,parallax_error]=',out$astrometry[,'parallax_error'],'\n')
    etaG <- 1- alphaG
    etaH <- 1- alphaH
#    etaG <- etaH <- 1
    ll2 <- loglikelihood(par.opt)$llastro
    out$plx <- plx0
#    out$plx <- plx0+0.029
#    dll <- ll1-ll0##likelihood increasement
    dll <- ll2-ll1##likelihood increasement
    bic <- 2*(dll-log(8))#BIC
    lnbf <- bic/2#lnBF
    cat('dLogLike=',dll,'\n\n')
#    cat('BIC=',bic,'\n')
#    cat('lnBF=',round(lnbf),'\n')

####simulation at Gaia starting and end epochs
    tsim0 <- rowSums(time_Cal2JD(rbind(c(2014,7,25),c(2016,5,23))))
    tsim <- tsim0-tmin
    Nt <- length(tsim)
    planet <- astrometry.kepler(par.opt,tt=tsim,bases=bases)$planet
    dras <- rbind(dras,planet[,'ra'])
    ddecs <- rbind(ddecs,planet[,'dec'])
    dpmras <- rbind(dpmras,planet[,'pmra'])
    dpmdecs <- rbind(dpmdecs,planet[,'pmdec'])
}
pdf('test.pdf',8,8)
par(mfrow=c(2,2))
dpos <- sqrt((dras[,1]-dras[,2])^2+(ddecs[,1]-ddecs[,2])^2)
dpm <- sqrt((dpmras[,1]-dpmras[,2])^2+(dpmdecs[,1]-dpmdecs[,2])^2)
ruwe <- c(1.34629,1.1768,1.3787,0.83358,0.94779,1.03412,0.86374,4.23935)
plot(dpos,ruwe,main=paste0('r=',cor(ruwe,dpos)))
plot(dpm,ruwe,main=paste0('r=',cor(ruwe,dpm)))
plot(dpos,ruwe,log='y',main=paste0('r=',cor(log(ruwe),dpos)))
plot(dpm,ruwe,log='y',main=paste0('r=',cor(log(ruwe),dpm)))
dev.off()

