####This code is to calculate the observer's barycentric velocity at each epoch of the RV data.
#1. get the parameters and astrometry information needed
ast <- out$astrometry[nrow(out$astrometry),c('ref_epoch','ra','dec','parallax','pmra','pmdec','radial_velocity')]
if(!is.null(out$tiall)>0){
    ts <- out$tiall[,'bjd']
    instrs <- out$tiall[,'instr']
    types <- out$tiall[,'type']
    its <- toupper(paste0(instrs,'_',types))
    it <- unique(its)
    peh <- array(NA,dim=c(nrow(out$tiall),3))
    for(kk in 1:length(it)){
        i <- gsub('_.+','',it[kk])
        if(grepl('AAT',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(-31.27704,149.0661,1.164)))
        }
        if(grepl('MINERVA',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(-27.797656941050338,151.8556103281065,0.682)))
        }
        if(grepl('PFS|MIKE',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(-29.013983,-70.692633,2.40792)))
        }
        if(grepl('ALMA',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(-23.026,-67.754,5.050)))
        }
        if(grepl('XINGLONG|LAMOST|lamost',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(40.389831774, 117.57166438,0.9)))
        }
        if(grepl('SMITH|HJS|MCD|McDonald|HPF',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(30.6718, -104.022,2.070)))
        }
        if(grepl('AJT',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(50.98014,11.71124,0.3507)))
        }
        if(grepl('HDS|SUBARU',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(19.701609,-155.08957,4.139)))
        }
        if(grepl('TRES',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(31.68895,-110.88460,2.630)))
        }
        if(grepl('WIYN|NEID|CFA',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(31.9575, -111.601,2.096)))
        }
        if(grepl('HIDES',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(34.57389,133.5965,0.359)))
        }
        if(grepl('HRS|TULL',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(30.681,-104.015,2.026)))
        }
        if(grepl('KECK|HIRES',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(19.82636,-155.47501,4.145)))
        }
        if(grepl('UVES|ESPRESSO|VLT|FORS',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(-24.625407,-70.403015,2.648)))
        }
        if(grepl('LBA',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(-33.000000,148.261667,0.392)))
        }
        if(grepl('HARPN|TV20|HARPSN|HARPS-N|SARG',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(28.754,-17.88814,2.37)))
        }
        if(grepl('BOES|SENS',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(36.16483,128.97670,1.146)))
        }
        if(grepl('SOPHIE|ELODIE',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(43.92944,5.7125,0.65)))
        }
        if(grepl('CARM|CARMENES|CAR',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(37.220791,-2.546847,2.168)))
        }
        if(grepl('EXPRES',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(35.202778,-111.664444,2.210)))
        }
        if(grepl('APO|APOGEE',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(32.780278,-105.820278,2.788)))
        }
                                        #        if(grepl('APF|LICK|C07|C14|C98',i)){
        if(grepl('APF|LICK',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(37.3425,-121.63825,1.274)))
        }
        if(grepl('HARPS03',i)){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(28.76206,-17.87764,2.336)))
        }
        if(grepl('C98|C07|C14|COR|HARPS$|HARPS03|HARPSP|LC|VLC|^CES$|CE$|FEROS',i) | (i=='COMB' & target=='GaiaBH2')){
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(-29.2584,-70.7345,2.4)))
        }
        if(grepl('^ES$|CHIRON|CTIO',i)){
            ##https://arxiv.org/pdf/1210.1616.pdf
            inds <- which(its==it[kk])
            peh[inds,] <- t(replicate(length(inds),c(-30.16928,-70.80679, 2.2419)))
        }
    }
    colnames(peh) <- c('phi','elong','height')
    deg2rad <- pi/180
    peh[,1:2] <- peh[,1:2]*deg2rad
    ra0 <- out$astrometry[out$iref,'ra']/180*pi#rad
    dec0 <- out$astrometry[out$iref,'dec']/180*pi#rad
    pqu0 <- cbind(c(-sin(ra0),cos(ra0),0),c(-sin(dec0)*cos(ra0),-sin(dec0)*sin(ra0),cos(dec0)),c(cos(dec0)*cos(ra0),cos(dec0)*sin(ra0),sin(dec0)))
    tpos <- out$astrometry[out$iref,'ref_epoch']#jd
    ind.nna <- which(!is.na(peh[,1]))
    ind.na <- which(is.na(peh[,1]))
    xyz <- peh
#    xyz[ind.nna,] <- sofa_Gd2gc(1,elong=peh[ind.nna,'elong'],phi=peh[ind.nna,'phi'],height=peh[ind.nna,'height']*1e3)/1e3#km
    xyz[ind.nna,] <- t(sapply(ind.nna,function(j) sofa_Gd2gc(1,elong=peh[j,'elong'],phi=peh[j,'phi'],height=peh[j,'height']*1e3)/1e3))
####for epochs without ephemerides we assume they are observed at the PFS/Megallen location because the parallax is only modified by about 1% for satellites located at L2 and it does not infuence the target system roemer delay
    if(length(ind.na)>0){
        peh[ind.na,] <- t(replicate(length(ind.na),c(-29.013983,-70.692633,2.40792)))
        xyz[ind.na,] <- t(sapply(ind.na,function(j) sofa_Gd2gc(1,elong=peh[j,'elong'],phi=peh[j,'phi'],height=peh[j,'height']*1e3)/1e3))
    }
# xyz[ind.na,] <- 0
    colnames(xyz) <- c('xtel','ytel','ztel')
    Par <- list(stars=target,secondary="NA",TaiType = 'instant',TtType='BIPM',TtTdbMethod='eph',DE=430,xgeoOff=0,ygeoOff=0,zgeoOff=0,vxgeoOff=0,vygeoOff=0,vzgeoOff=0,binary=FALSE,mode='emulate',SBscaling=FALSE,ra=out$astrometry[out$iref,'ra']*deg2rad,dec=out$astrometry[out$iref,'dec']*deg2rad,plx=out$astrometry[out$iref,'parallax'],pmra=out$astrometry[out$iref,'pmra'],pmdec=out$astrometry[out$iref,'pmdec'],rv=out$astrometry[out$iref,'radial_velocity'],tpos=tpos,pqu=pqu0,CompareT2=FALSE,near=FALSE,PlanetShapiro=FALSE,unit='TDB',Np=0)
    delat <- read.csv('tai_utc.csv',header=TRUE)
    Par$ObsInfo <- Data <-  data.frame(out$tiall,xyz,type='rv',peh,star=target,ObsType='ground',pmb=1013.25)
    dir <- 'data/'
    ##source('update_bipm.R')
    source('pexo_load_data.R')
####find UTC, TDB and the correct Earth's velocity interactively
    if(!sim){
        utc <- utc0 <- bjd0 <- cbind(floor(ts),ts%%1)
    }else{
        utc <- utc0 <- bjd0 <- cbind(floor(tsim),tsim%%1)
    }
#    if(dtype!='raw'){
    if(dtype=='raw'){
        Ntry <- 10
        for(j in 1:Ntry){
            OutObs <- time_Utc2tb(utc,Par)

            dt <- bjd0-bjd
            dt <- rowSums(utc-utc0)
            if(max(abs(dt))<1e-6) break()
            utc0 <- utc
        }
    }else{
        OutObs <- time_Utc2tb(utc,Par)
    }
###don't need so precise direction of OT
    OutTime <- time_Ta2te(OutObs,Par,fit=FALSE)
    out$vST <- OutTime$vST
    ZgO.geo <- 0.00887/2/1000/sqrt(rowSums(OutObs$GO[,1:3]^2))
    ZgO.sun <- 2.95/2/sqrt(rowSums(OutObs$SO[,1:3]^2))
#    out$SO <- OutObs$SO
    out$rSO <- OutObs$SO[,1:3]
    out$vSO <- OutObs$SO[,4:6]
    out$ZgO <- ZgO.geo+ZgO.sun
    out$ZsO <- rowSums(OutObs$SO[,4:6]^2)/CKMPS^2/2
    out$vSB <- OutTime$vSB
    out$rSB <- OutTime$rSB
    out$tauT <- out$BJDtdb <- rowSums(OutTime$BJDtdb)
    out$BJDtcb <- rowSums(OutTime$BJDtcb)
    out$RoemerS <- OutTime$RoemerSolar/3600/24#day
    out$ZgsO <- 1-1/OutObs$dTCB.dTT#positive values
#    lensing <- rv_LenSolar(OutTime$OL,OutTime$rOT,OutTime$vST,OutObs$SO[,4:6]/auyr2kms,g=1,LenRVmethod='T2')
#    out$ZlO <- lensing$Zall
    out$ZlO <- 0
    cat('Roemer delay:',head(out$RoemerS),'day\n')
    out$uOT <- OutTime$uOT
    out$uOB <- OutTime$uOB
    ##    du0 <- out$du <- cbind(ra=c(-cos(dec0)*sin(ra0),cos(dec0)*cos(ra0),0),dec=c(-sin(dec0)*cos(ra0),-sin(dec0)*sin(ra0),cos(dec0)))#error; before 20220419
    pqu0 <- out$pqu <- cbind(ra=c(-sin(ra0),cos(ra0),0),dec=c(-sin(dec0)*cos(ra0),-sin(dec0)*sin(ra0),cos(dec0)),u=c(cos(dec0)*cos(ra0),cos(dec0)*sin(ra0),sin(dec0)))#corrected; after 20220419
    if(!exists('sim')) sim <- FALSE
    ##    plx.vec <- -gen_CalUnit(t(sapply(1:nrow(OutObs$SO), function(i) OutObs$SO[i,1:3]-as.numeric(OutObs$SO[i,1:3]%*%pqu0[,3])*pqu0[,3])))
    uOS <- -gen_CalUnit(OutObs$SO[,1:3]/au2km)
    plx.vec <- t(sapply(1:nrow(uOS), function(i) uOS[i,]-as.numeric(uOS[i,]%*%pqu0[,3])*pqu0[,3]))
    plx.ra <- sapply(1:nrow(OutObs$SO),function(i) plx.vec[i,]%*%pqu0[,1])
    plx.dec <- sapply(1:nrow(OutObs$SO),function(i) plx.vec[i,]%*%pqu0[,2])
#    if(!sim){
        vSO0 <- out$vSO <- OutObs$SO[,4:6]
        rSB0 <- out$rSB <- OutTime$rSB
        out$pf <- data.frame(ra=plx.ra,dec=plx.dec)
#    }else{
#        out$vSO.sim <- OutObs$SO[,4:6]
#        out$rSB.sim <- OutTime$rSB
#        out$pf.sim <- data.frame(ra=plx.ra,dec=plx.dec)
#    }
    out$rvSO <- rowSums(out$vSO*out$uOT*1e3)#m/s; RV due to Earth's motion
}else{
#    out$vSO <- vSO0
#    out$rSB <- rSB0
#    out$pqu <- pqu0
    out$ind.all <- out$vSO <- out$rSB <- out$pqu <-  out$pf <- NULL
}
