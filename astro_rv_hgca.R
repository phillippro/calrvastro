library("plotrix")
source('OrbitFunction.R')
#if(!exists('par0s') | TRUE){
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    version <- as.integer(args[1])
}else{
    version <- 1
}
set.seed(100)
fs <- c('HD131664_fix1_relativityFALSE_Niter1100000_Ncores8_ofac2_Nset3_hg2--_transit0_P1984_acc11_sinI_lnlmax',
'HD16160_fix1_relativityFALSE_Niter1100000_Ncores8_ofac2_Nset8_hg2--_transit0_P27058_acc3_sinI_lnlmax',
'HD161797_fix1_relativityFALSE_Niter1100000_Ncores8_ofac2_Nset5_hg2--_transit0_P30059_acc0.19_sinI_lnlmax',
'HD190360_fix1_relativityFALSE_Niter1100000_Ncores8_ofac2_Nset6_hg2--_transit0_P2853_acc12_sinI_lnlmax',
'HD190406_fix1_relativityFALSE_Niter1000000_Ncores8_ofac2_Nset5_hg2--_transit0_P22879_acc7.6_sinI_lnlmax',
'HD39587_fix1_relativityFALSE_Niter1000000_Ncores8_ofac2_Nset5_hg2--_transit0_P5153_acc15_sinI_lnlmax',
'HD4747_fix1_relativityFALSE_Niter1000000_Ncores8_ofac2_Nset3_hg2--_transit0_P12028_acc17_sinI_lnlmax',
'HIP2552_fix1_relativityFALSE_Niter1100000_Ncores8_ofac2_Nset1_hg2--_transit0_P5500_acc0.13_sinI_lnlmax'
)
plotf <- FALSE
etaH <- etaG <- 1
fs <- sort(fs)
###modifiy input file names
stars <- gsub('_.+','',fs)
ff <- c()
for(k in 1:length(stars)){
    star  <- stars[k]
    ff <- c(ff,list.files(path=paste0('results/',star),pattern=paste0(fs[k],'.+Robj'),full.name=TRUE)[1])
}
fs <- ff
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
for(j in 1:length(fs)){
    star <- stars[j]
    cat('load ',fs[j],'\n')
    load(fs[j],env=e0<-new.env())
    if(Nmc > 0){
        par1.tot[[j]] <- t(rbind(t(e0$out$par.stat$sig1[1,]),e0$out$mcmc.opt$sig1[sample(1:e0$Niter,Nmc),]))
    }else{
        par1.tot[[j]] <- e0$out$par.stat$sig1[1,]
    }
#    cat('Mstar=',e0$Mstar,'\n')
#    Mstar <- Mstars[j]
#    Mstars[j] <- e0$Mstar
#    eMstar <- eMstars[j]
    cov.tot[[j]] <- e0$cov.astro
    ins <- e0$ins
    ins.tot[[j]] <- ins
    JD.tot[[j]] <- RV.tot[[j]] <- eRV.tot[[j]] <- res.tot[[j]] <- res.sig[[j]] <- nqp.tot[[j]] <- list()
    Nrv <- length(ins)
    Nrv.tot[[j]] <- Nrv
    cat('Nrv=',Nrv,'\n')
    cat('ins=',ins,'\n')
    for(k in 1:Nrv){
        t <- e0$out[[ins[k]]]$RV[,1]
        if(max(t)<24e5) t <- t+24e5
        JD.tot[[j]][[ins[k]]] <- t
        nqp.tot[[j]][[ins[k]]] <- e0$out[[ins[k]]]$noise$nqp
        RV.tot[[j]][[ins[k]]] <- e0$out[[ins[k]]]$RV[,2]
        eRV.tot[[j]][[ins[k]]] <- e0$out[[ins[k]]]$RV[,3]
        res.tot[[j]][[ins[k]]] <- e0$out$res.all$sig1[[ins[k]]]
        res.sig[[j]][[ins[k]]] <- e0$out$res.sig$sig1[[ins[k]]]
    }
    astrometry.tot[[j]] <- e0$out$astrometry
}
source('astro_rv_plot.R')


