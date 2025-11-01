####This code is to calculate the observer's barycentric velocity at each epoch of the RV data.
#1. get the parameters and astrometry information needed
ast <- out$astrometry[nrow(out$astrometry),c('ref_epoch','ra','dec','parallax','pmra','pmdec','radial_velocity')]
deg2rad <- pi/180
peh <- xyz <- t(replicate(length(tsim),c(0,0,0)))
colnames(xyz) <- c('xtel','ytel','ztel')
colnames(peh) <- c('phi','elong','height')
ra0 <- out$astrometry[out$iref,'ra']/180*pi#rad
dec0 <- out$astrometry[out$iref,'dec']/180*pi#rad
pqu0 <- cbind(c(-sin(ra0),cos(ra0),0),c(-sin(dec0)*cos(ra0),-sin(dec0)*sin(ra0),cos(dec0)),c(cos(dec0)*cos(ra0),cos(dec0)*sin(ra0),sin(dec0)))
tpos <- out$astrometry[out$iref,'ref_epoch']#jd
Par <- list(stars=target,secondary="NA",TaiType = 'instant',TtType='BIPM',TtTdbMethod='eph',DE=430,xgeoOff=0,ygeoOff=0,zgeoOff=0,vxgeoOff=0,vygeoOff=0,vzgeoOff=0,binary=FALSE,mode='emulate',SBscaling=FALSE,ra=out$astrometry[out$iref,'ra']*deg2rad,dec=out$astrometry[out$iref,'dec']*deg2rad,plx=out$astrometry[out$iref,'parallax'],pmra=out$astrometry[out$iref,'pmra'],pmdec=out$astrometry[out$iref,'pmdec'],rv=out$astrometry[out$iref,'radial_velocity'],tpos=tpos,pqu=pqu0,CompareT2=FALSE,near=FALSE,PlanetShapiro=FALSE,unit='TDB',Np=0)
delat <- read.csv('tai_utc.csv',header=TRUE)
Par$ObsInfo <- Data <-  data.frame(xyz,type='rv',peh,star=target,ObsType='space',pmb=1013.25)
dir <- 'data/'
source('pexo_load_data.R')
####find UTC, TDB and the correct Earth's velocity interatively
utc <- utc0 <- bjd <- cbind(floor(tsim),tsim%%1)
###don't need so precise direction of OT
OutObs <- time_Utc2tb(utc,Par)
OutTime <- time_Ta2te(OutObs,Par,fit=FALSE)
Roemer <- OutTime$RoemerOrder$Roemer1/3600/24#day
cat('Roemer delay:',head(Roemer),'day\n')
u <- OutTime$uOT
pqu0 <- out$pqu <- cbind(ra=c(-sin(ra0),cos(ra0),0),dec=c(-sin(dec0)*cos(ra0),-sin(dec0)*sin(ra0),cos(dec0)),u=c(cos(dec0)*cos(ra0),cos(dec0)*sin(ra0),sin(dec0)))#corrected; after 20220419
uOS <- -gen_CalUnit(OutObs$SO[,1:3])
plx.vec <- t(sapply(1:nrow(uOS), function(i) uOS[i,]-as.numeric(uOS[i,]%*%pqu0[,3])*pqu0[,3]))
plx.ra <- sapply(1:nrow(OutObs$SO),function(i) plx.vec[i,]%*%pqu0[,1])
plx.dec <- sapply(1:nrow(OutObs$SO),function(i) plx.vec[i,]%*%pqu0[,2])
out$vSO.sim <- OutObs$SO[,4:6]
out$rSB.sim <- OutTime$rSB
out$pf.sim <- data.frame(ra=plx.ra,dec=plx.dec)
