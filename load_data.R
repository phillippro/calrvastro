ff <-paste0('../data/combined/',target)

pmfun <- pemfun <- NULL
if(file.exists('plx_mass.txt')){
    tab <- read.table('plx_mass.txt',header=TRUE)
    pmfun.spec <- approxfun(tab[,'plx'],tab[,'spec_mass'])
    pmfun.evo <- approxfun(tab[,'plx'],tab[,'evo_mass'])
    pemfun.spec <- approxfun(tab[,'plx'],0.5*(tab[,'esmass1']+tab[,'esmass2']))
    pemfun.evo <- approxfun(tab[,'plx'],0.5*(tab[,'eemass1']+tab[,'eemass2']))
#    plxs <- seq(min(tab[,'plx']),max(tab[,'plx']),by=1e-3)
#    mspec <- pmfun.spec(plxs)
#    mevo <- pmfun.evo(plxs)
#    ind <- which.min(abs(mspec-mevo))
#    inds <- which(abs(mspec-mevo)/mspec<0.1)
}
####https://link.springer.com/article/10.1007/s10509-022-04066-1#data-availability
###http://adsabs.harvard.edu/abs/2013ApJS..208....9P
###http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
m22 <- read.table('../data/combined/mamajek22.txt',header=TRUE,check.names=TRUE)
m22[m22=='...'|m22=='....'|m22=='.....'|m22=='......'] <- NA
ms <- as.numeric(m22[,'Msun'])
mg <- as.numeric(gsub(':','',m22[,'M_G']))
ind <- which(!is.na(ms) & !is.na(mg))
mrl.m22 <- approxfun(ms[ind],mg[ind])
mlow.m22 <- min(ms[ind])
mup.m22 <- max(ms[ind])

fs <- list.files(path=ff,pattern=paste0(target,'|',binary),full.name=TRUE)
fs <- fs[!file.info(fs)$isdir]
fs <- fs[!grepl('~$|bin15min|bin30min|bkp|pdf|R$',fs)]
if(grepl('\\+',astro.type)){
    frel <- fs[grepl('rel',fs)]
}else{
    frel <- c()
}
fs0 <- fs
fs <- fs[!grepl('hg123$|g123$|htg$|hg3$|astro$|g23$|rel|hgca',fs)]


if(grepl('rel',astro.type)){
      fs <- c(fs,fs0[grep('rel',fs0)])
}else if(grepl('g23',astro.type)){
      fs <- c(fs,fs0[grep('g23$',fs0)])
}else if(grepl('htg',astro.type)){
      fs <- c(fs,fs0[grep('htg$',fs0)])
}else if(grepl('astro',astro.type)){
      fs <- c(fs,fs0[grep('astro$',fs0)])
}else if(grepl('hg3',astro.type)){
      fs <- c(fs,fs0[grep('hg3$',fs0)])
}else if(grepl('hgca2',astro.type)){
      fs <- c(fs,fs0[grep('hgca2$',fs0)])
}else if(grepl('hgca3',astro.type)){
      fs <- c(fs,fs0[grep('hgca3$',fs0)])
}else if(grepl('htg',astro.type)){
      fs <- c(fs,fs0[grep('htg$',fs0)])
}else if(grepl('hg123',astro.type)){
      fs <- c(fs,fs0[grep('hg123$',fs0)])
}else if(grepl('^g123',astro.type)){
      fs <- c(fs,fs0[grep('_g123$',fs0)])
}
fs <- unique(c(fs,frel))
#fs <- fs[!grepl('~$|bin15min|bin30min|bkp|csv|pdf|R$|astro$|rel$',fs)]
if(target=='HD42581'){
    fs <- fs[!grepl('LICK',fs)]
}
hip.epoch <- sum(time_Yr2jd(1991.25))
gdr1.epoch <- c(sum(time_Yr2jd(2015)),rowSums(time_CalHms2JD(rbind(c(2014,07,25,10,30,0),c(2015,09,16,0,0,0)))))
gdr2.epoch <- c(sum(time_Yr2jd(2015.5)),rowSums(time_CalHms2JD(rbind(c(2014,07,25,10,30,0),c(2016,5,23,11,35,0)))))
gdr3.epoch <- c(sum(time_Yr2jd(2016)),rowSums(time_CalHms2JD(rbind(c(2014,07,25,10,30,0),c(2017,05,28,8,44,0)))))
if(target=='HD197481') fs <- fs[!grepl('CSHELL',fs)]
out$Nastro <- 0
out$cov.epoch <- list()
#ins.ref <-
comp.epoch <- ins.epoch <- ins.rv <- c()
timetypes <- c()
data.binary <- c()
ind.binary <- list()
if(binary!='') out$binary <- list()
#if(!any(grepl('hg3',fs)) & any(grepl('astro',fs))) astro.type <- 'astro'
stars <- c()
catalogs <- c()
for(f in fs){
    star <- gsub('_.+','',gsub('.+\\/','',f))
    stars <- c(stars,star)
#    fin <- paste0(ff,f)
    fin <- f
    f <- gsub('.+\\/','',f)
    if(grepl('csv$',fin)){
        tab <- read.csv(fin,header=TRUE)
        if(grepl('gost',fin)){
            tab <- tab[,c('ObservationTimeAtBarycentre.BarycentricJulianDateInTCB.','scanAngle.rad.','parallaxFactorAlongScan','parallaxFactorAcrossScan')]
            colnames(tab) <- c('BJD','psi','parf','parx')
            inds <- which(tab[,'BJD']>=gdr3.epoch[2] & tab[,'BJD']<=gdr3.epoch[3])
            tab <- tab[inds,]
        }
    }else{
#        tab <- read.table(fin,sep=' ',numerals='no.loss')
#        if(class(tab[,1])=='factor' | class(tab[,1])=='character') tab <- read.table(fin,header=TRUE,check.name=FALSE,sep=' ',numerals='no.loss')
        tab <- try(read.table(fin,header=TRUE),TRUE)
        if(class(tab)=='try-error'){
            tab <- try(read.table(fin,header=TRUE,sep=' '),TRUE)
        }
        if(grepl('hipgaia',fin) & !use.gdr1){
            ii <- which(tab[,'ref_epoch']==gdr1.epoch[1])
            if(length(ii)>0) tab <- tab[-ii,]
        }
    }

#    if(!grepl('tim\\d$',fin)){
    if(max(tab[,1])<24e5 & max(tab[,1])>5e4){
        tab[,1] <- tab[,1]+24e5
    }else if(max(tab[,1])<5e4){
        tab[,1] <- tab[,1]+245e4
    }

    ts <- c()
    for(j in 1:ncol(tab)) ts <- c(ts,class(tab[,j]))
    ind <- which(ts!='factor')
    tab <- tab[,ind]
    if(nrow(tab)>0){
        cat('Read',fin,'\n')
        ##assume the outliers in data are removed
        if(grepl('tim\\d$',fin)){
            inds <- sort(tab[,1],index.return=TRUE)$ix
            tab <- tab[inds,]
            ncomp <- as.integer(gsub('.+tim','',fin))#companion number
            if(!any(names(out)=='timing')){
                    out$timing <- list()
            }
            np <- gsub('.+\\.tim','',fin)
            timetype <- gsub('.+_|\\.tim.+','',fin)#occultation or midtransit
            if(!any(names(out$timing)==paste0('p',np))) out[['timing']][[paste0('p',np)]] <- list()
            if(!grepl('pulsation|eb',timetype)){
                    out[['timing']][[paste0('p',np)]][[timetype]] <- tab[,c('BJD','eBJD')]
            }else{
                tim <- tab[,c('BJD','dt','et')]
                if(median(diff(sort(tim[,1])))<1 & nrow(tim)>1e3){
                    tim <- wtb(tim[,1],tim[,2],tim[,3],dt=1)
                    colnames(tim) <- c('BJD','dt','et')
                }
                out[['timing']][[paste0('p',np)]][[timetype]] <- tim
            }
            timetypes <- c(timetypes,timetype)
        }else if(!grepl('error|photo|par|astro|hip|rel|hg3|hg123|g123|gost|g23|hgca|htg|\\d.+AP\\.dat|rv\\d|csv|epoch|abs',f)){
            ##sort data; default first column: time
            inds <- sort(tab[,1],index.return=TRUE)$ix
            tab <- tab[inds,]
            if(max(tab[,1])< 24e5) tab[,1] <- tab[,1]+2400000
#            tab[,2] <- tab[,2]-mean(tab[,2])
#            if(offset) tab[,2] <- tab[,2]-tab[1,2]
            if(star==target){
                i <- gsub(paste0(target,'_|\\..+'),'',f)
                ins.rv <- c(ins.rv,i)
                out[[i]][['RV']] <- tab[,1:3]
                out[[i]][['data']] <- tab
            }else if(star==binary){
                instr <- gsub(paste0(binary,'_|\\..+'),'',f)
                out[[star]][[instr]] <- tab[,1:3]
                if(length(data.binary)==0){
                    nb <- 0
                }else{
                    nb <- nrow(data.binary)
                }
                ind.binary[[instr]] <- nb+1:nrow(tab)
                tmp <- tab[,1:3]
                colnames(tmp) <- c('BJD','RV','eRV')
                data.binary <- rbind(data.binary,data.frame(tmp,ins=rep(instr,nrow(tab))))
            }
            ##prepare indices
            Nact <- 0
            if(ncol(tab)>3){
                tmp <- tab[,4:ncol(tab),drop=FALSE]
                ferr <- gsub(i,paste0(i,'_error'),f)
                ind.err <- grep(ferr,fs)
                err <- NULL
                if(length(ind.err)>0){
#                    err <- read.table(paste0(ff,fs[ind.err]),header=TRUE,check.name=FALSE,sep=' ',numerals='no.loss')
                    err <- read.table(fs[ind.err],header=TRUE,check.name=FALSE,sep=' ',numerals='no.loss')
                }
                ns <- colnames(tmp)
                ns.act <-ns[!grepl('error',ns)]
                for(n in ns.act){
                    ind <- grep(paste0(n,'_error'),ns)
                    act <- tab[,n]
                    if(length(act)>1){
                        if(sd(act,na.rm=TRUE)>0){
                            if(grepl('AP',n)){
                                ind.keep <- which((act-mean(act))<xi*sd(act) &  !is.na(act))
                            }else{
                                ind.keep <- which((act-mean(act))<xi*sd(act) & act!=-1 & !is.na(act))
                            }
                            if(length(ind.keep)>0){
                                if(length(ind)==0 & is.null(err)){
                                    act <- scale(act[ind.keep])
                                    dact <- rnorm(length(ind.keep),0.1,0.001)
                                }else if(is.null(err)){
                                    dact <- tmp[ind.keep,ind]
                                    if(sd(dact,na.rm=TRUE)==0 | all(is.na(dact))){
                                        dact <- rnorm(length(ind.keep),0.1,0.001)
                                    }else{
                                        dact <- dact/sd(act[ind.keep])
                                    }
                                    act <- scale(act[ind.keep])
                                }else{
                                    dact <- rnorm(length(ind.keep),0.1,0.001)
                                    if(any(grepl(n,colnames(err)))){
                                        if(!all(is.na(err[ind.keep,n]))){
                                            if(sd(err[ind.keep,n],na.rm=TRUE)>0){
                                                ind.na <- which(is.na(err[ind.keep,n]))
                                                if(length(ind.na)==0){
                                                    dact <- err[ind.keep,n]/sd(act[ind.keep])
                                                }else{
                                                    dact <- rep(0,length(ind.keep))
                                                    dact[-ind.na] <- err[-ind.na,n]/sd(act[-ind.na])
                                                    dact[ind.na] <- rep(mean(dact[-ind.na]),length(ind.na))
                                                }
                                            }
                                        }
                                    }
                                    act <- scale(act[ind.keep])
                                }
                                data.act <- cbind(tab[ind.keep,1],act,dact)
                                out[[i]][['activity']][[n]] <-data.act
                            }
                        }
                    }
                }
            }
            out[[i]][['activity']][['window']] <- cbind(tab[,1],rep(1,nrow(tab)),rep(0.1,nrow(tab)))
        }else{
            if(grepl('\\d.+AP',f)){
                n <- gsub(paste0(target,'_|\\..+'),'',f)
                out[[n]] <- tab
            }
            if(grepl('photo',f)){
                nam <- gsub('.+_|\\..+','',f)
	        if(ncol(tab)<3) tab <- cbind(tab,sd(tab[,2])*0.01)
                out[['photometry']][[nam]] <- tab
                                        #if(!any(grepl('GP',noise.types))) noise.types <- c(noise.types,'GP011')
            }
            astrof <- FALSE
            if(grepl('astro',astro.type)) astrof <- grepl('astro$',f)
            if(grepl('hg3',astro.type)) astrof <- grepl('hg3$',f)
            if(grepl('htg',astro.type)) astrof <- grepl('htg$',f)
            if(grepl('g23',astro.type)) astrof <- grepl('g23$',f)
            if(grepl('hgca2',astro.type)) astrof <- grepl('hgca2$',f)
            if(grepl('hgca3',astro.type)) astrof <- grepl('hgca3$',f)
            if(grepl('hg123',astro.type)) astrof <- grepl('hg123$',f)
            if(grepl('^g123',astro.type)) astrof <- grepl('_g123$',f)
            if(grepl('.+abs',f)){
#                ii <- gsub(paste0(target,'_|\\..+'),'',f)
                ii <- gsub('.+_|\\..+','',f)
                out$data.epoch[[ii]] <- tab
                if(grepl('hip1',f)){
                    orbits <- unique(tab[,'orbit'])
                    cov.epoch <- array(0,dim=c(nrow(tab),nrow(tab)))
                    for(j in 1:length(orbits)){
                        indF <- which(tab[,'orbit']==orbits[j] & (tab[,'solution']=='F' | tab[,'solution']=='f'))
                        indN <- which(tab[,'orbit']==orbits[j] & (tab[,'solution']=='N' | tab[,'solution']=='n'))
                        if(length(indF)==1 & length(indN)==1){
#                            Tepoch[indF,indN] <- cbind(cov.epoch,rbind(c(1/tab[indF,'ev'],0),c(-tab[indF,'qNF']/tab[indF,'ev']/sqrt(1-tab[indF,'qNF']^2),1/tab[indN,'ev']/sqrt(1-tab[indF,'qNF']^2))))
                            cov.epoch[c(indF,indN),c(indF,indN)] <- rbind(c(tab[indF,'ev']^2,tab[indF,'qNF']*tab[indF,'ev']*tab[indN,'ev']),c(tab[indF,'qNF']*tab[indF,'ev']*tab[indN,'ev'],tab[indN,'ev']^2))
                        }else if(length(indF)==1){
                            cov.epoch[indF,indF] <- tab[indF,'ev']^2
                        }else if(length(indN)==1){
                            cov.epoch[indN,indN] <- tab[indN,'ev']^2
                        }
                    }
                    out$cov.epoch[[ii]] <- cov.epoch
#                    out$astro.epoch <- read.table(paste0('../data/combined/hip1/',target,'_hip1_cat.txt'),header=TRUE)
                }
                ins.epoch <- c(ins.epoch,ii)
                if(grepl('abs2',f)){
                    comp.epoch <- c(comp.epoch,2)
                }else{
                    comp.epoch <- c(comp.epoch,1)
                }
            }else if(astrof){
                if(any(diff(tab[,1])==0)){
                    tab <- tab[-which(diff(tab[,1])==0),]
                }
                tab <- tab[!is.na(tab[,'parallax']),]
###remove Gaia DR1 targets without parallax
                if(target=='HIP111958'){
                    tab[,'radial_velocity'] <- -4.84
                    tab[,'radial_velocity_error'] <- 0.24
                }
                if(target=='WD0643-16'){
###The gravitational redshift of Sirius B
##https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.2361J/abstract
                    tab[,'radial_velocity'] <- -7.69
                    tab[,'radial_velocity_error'] <- 0.01
                    tab[,'mass'] <- 1.018
                    tab[,'mass_error'] <- 0.011
                }
                rv <- tab[,'radial_velocity']
                erv <- tab[,'radial_velocity_error']
		ind.na <- which(is.na(rv))
		if(length(ind.na)>=2){
                    print('Warning: No valide systematic RV in astrometry data and assign RV=0!')
                    rv <- erv <- tab[,'radial_velocity'] <- tab[,'radial_velocity_error'] <- 0
                }
		if(length(ind.na)==1){
			rv[ind.na] <- rv[-ind.na]
			erv[ind.na] <- erv[-ind.na]
		}
		tab[,'radial_velocity'] <- rv
		tab[,'radial_velocity_error'] <- erv
                Mstar <- tab[2,'mass']
                if(is.na(Mstar)){
                    Mstar <- 1
                    eMstar <- 1
                }else{
                if(tab[nrow(tab),'mass.upper']<tab[nrow(tab),'mass']){
                    eMstar <- 0.5*(tab[nrow(tab),'mass.upper']+tab[nrow(tab),'mass.lower'])
                }else{
                    eMstar <- 0.5*(tab[nrow(tab),'mass.upper']-tab[nrow(tab),'mass.lower'])
                }#Msol
                }
                if(star==target){
                    out$astrometry <- tab
                    out$Mstar <- Mstar
                    out$eMstar <- eMstar
                }else if(star==binary){
                    out$binary$astrometry <- tab
                    out$binary$Mstar <- Mstar
                    out$binary$eMstar <- eMstar
                }
		if(is.na(Mstar)) cat('No stellar mass available in astrometry data!\n')
                Nastro <- nrow(tab)
                out$Nastro <- Nastro
                ##derive covariance matrix of parallax and proper motion
                cov.astro <- array(0,dim=c(6,6,Nastro))
                for(i in 1:Nastro){
                    cov1 <- array(0,dim=c(5,5))
                    if(is.na(tab[i,'parallax'])){
                        tab[i,c('parallax','pmra','pmdec','parallax_error','pmra_error','pmdec_error','ra_dec_cov','ra_parallax_cov','ra_pmra_cov','ra_pmdec_cov','dec_parallax_cov','dec_pmra_cov','dec_pmdec_cov','parallax_pmra_cov','parallax_pmdec_cov','pmra_pmdec_cov')] <- c(tab[Nastro,c('parallax','pmra','pmdec')],1e3,1e3,1e3,rep(0,10))
                    }
                    diag(cov1) <- as.numeric(tab[i,c('ra_error','dec_error','parallax_error','pmra_error','pmdec_error')]^2)
                    cov1[upper.tri(cov1)] <- as.numeric(tab[i,c('ra_dec_cov','ra_parallax_cov','dec_parallax_cov','ra_pmra_cov','dec_pmra_cov','parallax_pmra_cov','ra_pmdec_cov','dec_pmdec_cov','parallax_pmdec_cov','pmra_pmdec_cov')])
                    cov1[lower.tri(cov1)] <- t(cov1)[lower.tri(cov1)]
                    cov.astro[1:5,1:5,i] <- cov1
                    erv <- as.numeric(tab[i,'radial_velocity_error'])
                    if(is.na(erv)) erv <- 0.1*as.numeric(tab[i,'radial_velocity'])
                    cov.astro[6,6,i] <- erv^2
                }
                if(star==target){
                    out$cov.astro <- cov.astro
                }else if(star==binary){
                    out[[star]]$cov.astro <- cov.astro
                }
                Nkeppar <- Nkeppar+2
            }else if(grepl('hip2',f)){
#                astrometry <- read.table(paste0(ff,f),header=TRUE,sep=' ',numerals='no.loss')
                out$astrometry <- read.table(fin,header=TRUE,numerals='no.loss')
                catalogs <- c(catalogs,'hip2')
            }else if(grepl('gost',f)){
                if(star==target){
                    out$gost <-  tab
                }else if(star==binary){
                    out[[star]]$gost <- tab
                }
            }
            if(grepl('\\.rel',f)){
#                tab <- read.table(paste0(ff,f),header=TRUE,sep=' ',numerals='no.loss')
#                tab <- read.table(f,header=TRUE,sep=' ',numerals='no.loss')
                tab <- read.table(fin,header=TRUE,numerals='no.loss')
                cn <- colnames(tab)
                if(any(grepl('PA|pa|rho|sep',cn))){
                    rel.type <- 'sp'
                    cn <- gsub('PA','pa',cn)
                    cn <- gsub('rho','sep',cn)
                }else{
                    rel.type <- 'rd'
                    cn <- gsub('edra','era',cn)
                    cn <- gsub('edpmra','epmra',cn)
                    cn <- gsub('eddec','edec',cn)
                    cn <- gsub('edpmdec','epmdec',cn)
                    cn <- gsub('edparallax','eparallax',cn)
                }
                colnames(tab) <- cn
                sigfit <- as.integer(gsub('.+rel','',f))
                ii <- gsub(paste0(target,'_|.rel.+'),'',f)
                index <- match(sigfit,out$rel.order)
                if(!is.na(index) & (sigfit<=4 | astro5)){
                    if(!any(names(out)=='rel')){
                        out[['rel']] <- list()
#                        out[['cov']] <- list()
                    }
                    n <- paste0('p',index)
                    if(!any(names(out$rel)==n)){
                        out[['rel']][[n]] <- list()
#                        out[['cov']][[n]] <- list()
                    }
                    if(ii=='edr3' | ii=='dr2'){
                        qq <- tab
                    }else{
                        if(rel.type=='sp'){
                            if(ncol(tab)>5){
                                tab[,6] <- 0
                            }else{
                                tab <- cbind(tab,0)
                            }
                            colnames(tab)[6] <- 'cov.ra.dec'
                        }
                        if(ncol(tab)==5){
                            qq <- cbind(tab,0)
                        }else if(ncol(tab)==6){
                            if(is.numeric(tab[,6])){
                                if(colnames(tab)[6]=='cor.ra.dec' | colnames(tab)[6]=='cor'){
                                    cov <- tab[,6]*tab[,3]*tab[,5]
                                }else{
                                    cov <- tab[,6]
                                }
                                qq <- cbind(tab[,1:5],cov)
                            }else{
                                qq <- cbind(tab[,1:5],0)
                            }
                        }
                        if(rel.type=='rd'){
                            colnames(qq) <- c('BJD','dra','era','ddec','edec','cov.ra.dec')
                        }else{
                            colnames(qq) <- c('BJD','sep','esep','pa','epa','cov')
                        }
                    }
                    out$rel[[n]][[ii]] <- qq

                    if(grepl('edr3',f)){
                        cov.rel <- list()
                        for(k3 in 1:nrow(tab)){
                            cov <- vec2cov(qq[k3,])
                            cov <- cov[,-3]
                            cov <- cov[-3,]
                            cov.rel[[k3]] <- cov
                        }
                        out$cov[[n]][[ii]] <- cov.rel
                    }
                    if(!any(names(out)=='astrometry')){
                        fin <- '../data/code/SigClassify/common_withSig1.csv'
                        if(file.exists(fin)){
                            dat <- read.csv(fin,header=TRUE)
                            ind <- which(target==dat[,'target'])
                            if(length(ind)>0){
#                                out$plx <- dat[ind[1],'parallax']
                                out$Mstar <- dat[ind[1],'Ms']
                                out$eMstar <- 0.1*dat[ind[1],'Ms']
                            }
                        }
                    }
                }
            }
            if(grepl('rv\\d',fin)){
                tab <- read.table(fin,header=TRUE,sep=' ',numerals='no.loss')
                sigfit <- as.integer(gsub('.+rv','',f))
                index <- match(sigfit,out$rel.order)
                if(!is.na(index)){
                    if(!any(names(out)=='rvc')){
                        out[['rvc']] <- list()
                    }
                    instr <- gsub(paste0('.+_|\\..+'),'',f)
                    n <- paste0('p',index)
                    if(!any(names(out$rvc)==n)) out$rvc[[n]] <- list()
                    out$rvc[[n]] <- data.frame(tab,ins=rep(instr,nrow(tab)))
                }
            }
            if(grepl('par',f)){
                out[['parameter']] <- tab
            }
        }
    }
}


out$timetypes <- timetypes
out$plx <- NA
for(star in unique(stars)){
    if(star==target){
        astrometry <- out$astrometry
        gost <- out$gost
    }else{
        astrometry <- out[[star]]$astrometry
        gost <- out[[star]]$gost
    }
    if(length(astrometry)>0){
        ihip <- which(abs(astrometry[,'ref_epoch']-hip.epoch)<1)#reference epoch
        ityc <- c()
        if(any(colnames(astrometry)=='catalog')){
            if(tycf) ityc <- which(astrometry[,'catalog']=='TYC'|astrometry[,'catalog']=='tycho2')
            if(length(ityc)>0){
                ihip <- c()
####don't treat TYCHO catalog data as epoch data
                                        #            data.epoch[['TYC']] <- astrometry[ityc,]
                                        #            ins.epoch <- c(ins.epoch,'TYC')
                                        #            era <- astrometry[ityc,'ra_error']
                                        #            edec <- astrometry[ityc,'dec_error']
                                        #            cov.radec <- astrometry[ityc,'ra_dec_cov']
                                        #            cov.epoch[['TYC']] <- rbind(c(era^2,cov.radec),c(cov.radec,edec^2))
            }
        }
        igdr1 <- which(abs(astrometry[,'ref_epoch']-gdr1.epoch[1])<1)##don't use GDR1 because it is correlated with Hipparcos***
        igdr2 <- which(abs(astrometry[,'ref_epoch']-gdr2.epoch[1])<1)
        igdr3 <- which(abs(astrometry[,'ref_epoch']-gdr3.epoch[1])<1)
        iref <- which.max(astrometry[,'ref_epoch'])#reference epoch
        if(star==target) out$plx <- astrometry[iref,'parallax']
####remove GDR1
        if(length(igdr1)>0){
            if(is.na(astrometry[igdr1,'parallax'])){
                astrometry[igdr1,c('parallax','pmra','pmdec','parallax_error','pmra_error','pmdec_error','ra_dec_cov','ra_parallax_cov','ra_pmra_cov','ra_pmdec_cov','dec_parallax_cov','dec_pmra_cov','dec_pmdec_cov','parallax_pmra_cov','parallax_pmdec_cov','pmra_pmdec_cov')] <- c(astrometry[igdr3,c('parallax','pmra','pmdec')],1e3,1e3,1e3,rep(0,10))
            }
        }
        dtgs <- a1s <- a2s <- a3s <- a4s <- a5s <- list()
        astro.index <- index <- 1
        astro.gost <- c()
        if(nrow(astrometry)>1){
            if(length(ihip)>0){
                index <- (1:nrow(astrometry))[-ihip]
            }else{
                index <- 1:nrow(astrometry)
            }
            if(all(is.na(astrometry[,'radial_velocity']))) astrometry[,'radial_velocity'] <- 0

            for(j in index){
                dra <- (astrometry[j,'ra']-astrometry[iref,'ra'])*cos(astrometry[j,'dec']/180*pi)*3.6e6
                ddec <- (astrometry[j,'dec']-astrometry[iref,'dec'])*3.6e6
                astro.gost <- rbind(astro.gost,c(dra,ddec,unlist(astrometry[j,c('parallax','pmra','pmdec')])))
            }
            colnames(astro.gost) <- c('dra','ddec','parallax','pmra','pmdec')
            astro.index <- index
            if(length(gost)>0){
###if both Gaia and Hipparcos have (simulated) epoch data, then set Nastro=0, or no need to approximate the reference astrometry as instatanous astrometry
                if(length(out$data.epoch)>0) out$Nastro <- 0
                ind.hip <- NA
                ind.gdr1 <- which(gost[,'BJD']>=gdr1.epoch[2] & gost[,'BJD']<=gdr1.epoch[3])
                ind.gdr2 <- which(gost[,'BJD']>=gdr2.epoch[2] & gost[,'BJD']<=gdr2.epoch[3])
                ind.gdr3 <- which(gost[,'BJD']>=gdr3.epoch[2] & gost[,'BJD']<=gdr3.epoch[3])
                index <- list()
                for(k in 1:length(astro.index)){
                    j <- astro.index[k]
                    jj <- which.min(abs(astrometry[j,'ref_epoch']-c(hip.epoch[1],gdr1.epoch[1],gdr2.epoch[1],gdr3.epoch[1])))
                    if(jj==1){
                        index[[k]] <- ind.hip
                        if(length(ityc)>0)    catalogs <- c(catalogs,'TYC')
                        if(length(ihip)>0)    catalogs <- c(catalogs,'HIP')
                    }
                    if(jj==2){
                        index[[k]] <- ind.gdr1
                        catalogs <- c(catalogs,'GDR1')
                    }
                    if(jj==3){
                        index[[k]] <- ind.gdr2
                        catalogs <- c(catalogs,'GDR2')
                    }
                    if(jj==4){
                        index[[k]] <- ind.gdr3
                        catalogs <- c(catalogs,'GDR3')
                    }
                    if(jj>1){
                        dtgs[[k]] <- dtg <- (gost[index[[k]],'BJD']-astrometry[j,'ref_epoch'])/365.25
                        a1s[[k]] <- sin(gost[index[[k]],'psi'])
                        a2s[[k]] <- cos(gost[index[[k]],'psi'])
                        a3s[[k]] <- gost[index[[k]],'parf']
                        a4s[[k]] <- sin(gost[index[[k]],'psi'])*dtg
                        a5s[[k]] <- cos(gost[index[[k]],'psi'])*dtg
                    }else{
                        dtgs[[k]] <- (astrometry[jj,'ref_epoch']-astrometry[j,'ref_epoch'])/365.25
                        a1s[[k]] <- a2s[[k]] <- a3s[[k]] <- a4s[[k]] <- a5s[[k]] <- NA
                    }
                }
            }
        }
        if(star==target){
            out$iref <- iref
            out$ityc <- ityc
            out$ihip <- ihip
            out$igdr1 <- igdr1
            out$igdr2 <- igdr2
            out$igdr3 <- igdr3
            out$catalogs <- catalogs
            out$cats <- out$catalogs[!grepl('hip',out$catalogs)]
            out$cat.ind <- index
            out$dtg <- dtgs
            out$a1 <- a1s
            out$a2 <- a2s
            out$a3 <- a3s
            out$a4 <- a4s
            out$a5 <- a5s
            out$astro.gost <- astro.gost
            out$astro.index <- astro.index
        }else if(star==binary){
            out[[star]]$iref <- iref
            out[[star]]$ityc <- ityc
            out[[star]]$ihip <- ihip
            out[[star]]$igdr1 <- igdr1
            out[[star]]$igdr2 <- igdr2
            out[[star]]$igdr3 <- igdr3
            out[[star]]$catalogs <- catalogs
            out[[star]]$cats <- out$catalogs[!grepl('hip',out$catalogs)]
            out[[star]]$cat.ind <- index
            out[[star]]$dtg <- dtgs
            out[[star]]$a1 <- a1s
            out[[star]]$a2 <- a2s
            out[[star]]$a3 <- a3s
            out[[star]]$a4 <- a4s
            out[[star]]$a5 <- a5s
            out[[star]]$astro.gost <- astro.gost
            out[[star]]$astro.index <- astro.index
        }
    }
}

if(!any(names(out)=='plx')){
    cmd <- paste('python get_cor_from_simbad.py',target)
    cat('cmd')
    system(cmd)
    if(!file.exists('test.txt')){
        cmd <- paste('python2.7 get_cor_from_simbad.py',target)
        cat('cmd')
        system(cmd)
    }
    if(!file.exists('test.txt')) stop('I cannot find the paralax information for relative astrometry!')
                                        #                            out$plx <- read.table('test.txt',header=TRUE)[1,3]
    system('rm test.txt')
}
if(target=='HD128620'){
    out$plx <- 750.81
    out$eplx <- 0.38
    out$Mstar <- 1.0788
    out$eMstar <- 0.0029
}

out$Nrvc <- length(out$rvc)
out$data.binary <- data.binary
out$ind.binary <- ind.binary
out$ins.binary <- unique(data.binary[,'ins'])
out$ins.rv <- ins.rv
out$ins.epoch <- ins.epoch
out$comp.epoch <- comp.epoch

###number of astrometry data sets
astro.ins <- c()
astro.par <- 0
if(length(out$data.epoch)>0){
    astro.ins <- c(astro.ins,out$ins.epoch)
    astro.par <- max(astro.par,5)
}
if(length(out$gost)>0){
    astro.ins <- c(astro.ins,'gaia')
    astro.par <- max(astro.par,5)
}
if(length(out$astrometry)>0){
    if(length(out$ihip)>0 & !any(grepl('hip',astro.ins))){
        astro.ins <- c(astro.ins,'hip')
        if(out$Nastro>0) astro.par <- max(astro.par,4)
    }
    if(length(out$ityc)>0 & !any(grepl('TYC',astro.ins))){
        astro.ins <- c(astro.ins,'TYC')
        if(out$Nastro>0) astro.par <- max(astro.par,4)
    }
    if(length(out$igdr1)+length(out$igdr2)+length(out$igdr3)>0){
        astro.ins <- c(astro.ins,'gaia')
        if(out$Nastro>0) astro.par <- max(astro.par,4)
    }
}
if(astro.par>0 & relativity>0){
    astro.par <- max(astro.par,5)
}
out$astro.par <- astro.par
out$astro.ins <- unique(astro.ins)
if(any(grepl('hip|tyc',names(out$data.epoch)))) out$ind.astro <- 0
out$Irvc <- as.integer(gsub('p','',names(out[['rvc']])))
cat('out$plx=',out$plx,'\n')
if(!any(names(out)=='Mstar')){
   out$Mstar <- 1
   out$eMstar <- 0.1
   if(target=='HR8799') out$Mstar <- 1.48516057021498
   if(target=='GJ559A') out$Mstar <- 1.1
}
mm <- read.table('mstars.txt',header=TRUE)
if(any(mm[,1]==target)){
    ii <- which(target==mm[,1])
    out$Mstar <- mm[ii,2]
    out$eMstar <- mm[ii,3]
}
cat('out$Mstar=',out$Mstar,'\n')
cat('out$eMstar=',out$eMstar,'\n')
rel.name <- c()
rel.ins <- c()
if(any(names(out)=='rel')){
    out[['indrel']] <- list()
    out[['typerel']] <- list()
    out[['Nrel']] <- list()#cummulative number of relative astrometry data
    Nrel <- trel <- c()
    for(j in 1:Nmax){
        n <- paste0('p',j)
        if(!any(names(out$rel)==n)){
            out$rel[n] <- list(NULL)
            if(grepl('edr3.+rel',f)){
                out$cov[n] <- list(NULL)
                eta0 <- rep(NA,10)
                ss <- unlist(strsplit(f,''))
                ind <- as.integer(ss[length(ss)])
                eta0[ind] <- 0
            }
        }
        if(!is.null(out$rel[[n]])){
            tt <- c()
            inss <- names(out$rel[[n]])
            out$indrel[[n]] <- list()
            out$typerel[[n]] <- list()
            for(i in inss){
                rel.ins <- c(rel.ins,i)
                ts <- out$rel[[n]][[i]][,1]
                rel.name <- c(rel.name,paste0(n,'.',i))
                N0 <- length(trel)+1
                trel <- c(trel,ts)
                out$indrel[[n]][[i]] <- N0:length(trel)
                reltype <- 'rd'
                if(any(grepl('sep',colnames(out$rel[[n]][[i]])))) reltype <- 'sp'
                out$typerel[[n]][[i]] <- reltype
            }
        }else{
            out$indrel[j] <- list(NULL)
            out$typerel[j] <- list(NULL)
        }
        Nrel <- c(Nrel,length(trel))
    }
    out$trel <- trel
    out$Nrel <- Nrel
}
out$rel.ins <- unique(rel.ins)
###data property
Nset <- length(out$ins.rv)
trv.all <- RV.all <- eRV.all <- c()
index <- 0
ins.all <- trv.all <- RV.all <- eRV.all <- NULL
out$Nrv <- 0
if(Nset>0){
    for(i3 in 1:Nset){
        trv.all <- c(trv.all,out[[out$ins.rv[i3]]]$RV[,1])
        index <- c(index,length(trv.all))
        RV.all <- c(RV.all,out[[out$ins.rv[i3]]]$RV[,2])
        eRV.all <- c(eRV.all,out[[out$ins.rv[i3]]]$RV[,3])
        ins.all <- c(ins.all,rep(out$ins.rv[i3],length(out[[out$ins.rv[i3]]]$RV[,3])))
    }
#    ind <- sort(trv.all,index.return=TRUE)$ix; sorting
    ind <- 1:length(trv.all)#no sorting
    for(i3 in 1:Nset){
        out[[out$ins.rv[i3]]]$index <- match((index[i3]+1):index[i3+1],ind)
    }
    trv.all <- trv.all[ind]
    RV.all <- RV.all[ind]
    eRV.all <- eRV.all[ind]
    ins.all <- ins.all[ind]
    out[['all']] <- data.frame(trv.all,RV.all,eRV.all,ins.all)
    Ndata <- length(trv.all)
    rmin <- min(RV.all)
    rmax <- max(RV.all)
    tmin <- min(trv.all)
    tmax <- max(trv.all)
    out$Nrv <- length(trv.all)
}else{
    if(out$Nrvc>0){
        tt <- unlist(lapply(1:length(out$rvc),function(i) out$rvc[[i]][,1]))
        out$tmin <- tmin <- min(tt)
        out$tmax <- tmax <- max(tt)
    }else if(length(out$astrometry)>0){
        out$tmin <- tmin <- out$astrometry[1,1]
        out$tmax <- tmax <- out$astrometry[2,1]
    }else if(any(names(out)=='rel')){
        tt <- unlist(lapply(1:length(out$rel),function(i) out$rel[[i]][,1]))
        out$tmin <- tmin <- min(tt)
        out$tmax <- tmax <- max(tt)
    }
    rmin <- -1e3
    rmax <- 1e3
}
if(length(out$astrometry)>0){
    out$tref <- out$astrometry[iref,'ref_epoch']
}else{
    out$tref <- tref <- gdr3.epoch[1]
}
tspan <- tmax-tmin
if(tmax<2400000 & Ntr>0) Tc <- Tc%%2400000
Pmin <- 0.1
Pmax <- max(1e7,10*(tmax-tmin))
####limit the period to specific ranges for the following stars
yr2d <- 365.25
cat('Pmin=',Pmin,'\n')
cat('Pmax=',Pmax,'\n')

fmin0 <- fmin <- 1/Pmax
fmax0 <- fmax <- 1/Pmin

###whether exist pars/target.par
#fs <- list.files(path='pars',pattern=paste0('^',target,'.+par$'),full.name=TRUE)
out$replacePar <- FALSE
fin  <- paste0('pars/',target,'.par')
if(file.exists(fin)){
    out$replacePar <- TRUE
    tmp <- read.table(fin,numerals='no.loss')
    pp <- as.numeric(tmp[,2])
    names(pp) <- as.character(tmp[,1])
    out$IniPar <- pp
}

NormMass <- FALSE
out$mp <- NA
#if(out$Nrv==0 & out$Nm==0 & out$Nastro==0){
if(out$Nrv==0 & out$Nm==0 & out$Nastro==0 & is.null(out$data.epoch)){
    NormMass <- TRUE
    if(target=='HR8799') out$mp <- rbind(c(5,7,7,10,10,10),c(2,3,3,7,7,7))#Mjup
    if(target=='HD163296') out$mp <- rbind(c(2.18),c(0.5))#Mjup
}
out$NormMass <- NormMass
nepoch <- out$Nastro
if(grepl('hgca',astro.type)){
    nepoch <- 4
}

###change epoch data
out$all.epoch <- c()
ttmp <- c()
if(length(out$data.epoch)>0){
    ind <- which(!grepl('hip',out$ins.epoch))
#    ind <- 1:length(out$ins.epoch)
    off.epoch <- list()
    if(length(ind)>0){
        for(j in ind){
            i <- out$ins.epoch[j]
            c <- out$comp.epoch[j]
#            pos.dif <- sapply(1:nrow(out$data.epoch[[i]]),function(j) AstroDiff(out$astrometry[out$iref,],out$data.epoch[[i]][j,]))
#            cosdec <- cos(out$data.epoch[[i]][,'dec']/180*pi)
#            tmp <- data.frame(out$data.epoch[[i]],dra=(out$data.epoch[[i]][,'ra']-out$data.epoch[[i]][1,'ra'])*3.6e6*cosdec,ddec=(out$data.epoch[[i]][,'dec']-out$data.epoch[[i]][1,'dec'])*3.6e6,instr=i,comp=c)
            cosdec <- cos(out$astrometry[out$iref,'dec']/180*pi)
            tmp <- data.frame(out$data.epoch[[i]],dra=(out$data.epoch[[i]][,'ra']-out$astrometry[out$iref,'ra'])*3.6e6*cosdec,ddec=(out$data.epoch[[i]][,'dec']-out$astrometry[out$iref,'dec'])*3.6e6,instr=i,comp=c)
            out$data.epoch[[i]] <- tmp
            ttmp <- rbind(ttmp,tmp)
        }
    }
    out$all.epoch <- ttmp
}
out$tiall <- c()
tmp <- c()
if(length(out$ins.rv)>0){
    tmp0 <- data.frame(bjd=out$all[,1],instr=out$all$ins.all,type='rv')
#    tmp0 <- tmp0[sort(tmp0[,1],index.return=TRUE)$ix,]
    tmp <- rbind(tmp,tmp0)
}
if(length(out$astrometry)>0){
    epochs <- c(hip.epoch,gdr1.epoch[1],gdr2.epoch[1],gdr3.epoch[1])
    cat.name <- c('hip','GDR1','GDR2','GDR3')
    ii <- match(floor(out$astrometry[,1]),floor(epochs))
    if(length(which(!is.na(ii)))==0){
        cat('there is no known catalogs matching any astrometry catalog data!')
        stop()
    }
    tmp0 <- data.frame(bjd=out$astrometry[,1],instr=cat.name[ii],type='catalog')
    tmp <- rbind(tmp,tmp0)
}
if(length(out$timing)>0){
    tts <- c()
    for(n1 in names(out$timing)){
      for(n2 in names(out$timing[[n1]])){
        tmp0 <- data.frame(bjd=out$timing[[n1]][[n2]][,1],instr=paste0(n1,'_',n2),type='timing')
        tmp <- rbind(tmp,tmp0)
        tts <- c(tts,out$timing[[n1]][[n2]][,1])
      }
    }
    out$ttv0 <- min(tts)
}
if(length(out$data.epoch)>0){
    if(!is.null(out$all.epoch)){
        tmp0 <- data.frame(bjd=out$all.epoch[,1],instr=out$all.epoch[,'instr'],type=paste0('epoch',out$all.epoch[,'comp']))
        tmp <- rbind(tmp,tmp0)
    }
    if(any(grepl('hip',names(out$data.epoch)))){
        ii <- grep('hip',names(out$data.epoch))
        for(j in ii){
            i <- names(out$data.epoch)[j]
            tmp0 <- data.frame(bjd=out$data.epoch[[i]][,1],instr=i,type=paste0('epoch',out$comp.epoch[j]))
            tmp <- rbind(tmp,tmp0)
        }
    }
}
out$ind.rel <- c()
if(!is.null(out$rel)){
    nps <- gsub('p','',names(out$rel))
    for(p in nps){
        rels <- out$rel[[paste0('p',p)]]
        inss <- names(rels)
        for(ins in inss){
            tmp0 <- data.frame(bjd=rels[[ins]][,1],instr=ins,type=paste0('rel',p))
            out$ind.rel <- c(out$ind.rel,nrow(tmp)+(1:nrow(tmp0)))
            tmp <- rbind(tmp,tmp0)
        }
    }
}
if(length(out$data.binary)>0){
    tmp0 <- data.frame(bjd=out$data.binary[,1],instr=out$data.binary[,'ins'],type='binary')
    tmp <- rbind(tmp,tmp0)
    out$ind.binary <- list()
    for(i in unique(out$data.binary[,'ins'])){
        out$ind.binary[[i]] <- which(out$data.binary[,'ins']==i)
    }
}

if(length(out$gost)>0){
    tmp0 <- data.frame(bjd=out$gost[,1],instr='Gaia',type='gost')
    tmp <- rbind(tmp,tmp0)
}

if(length(tmp)>0){
    colnames(tmp) <- c('bjd','instr','type')
    out$tiall <- tmp
    out$indP <- which(tmp[,'type']!='binary')
    if(binary!=''){
        out$indC <- which(tmp[,'type']=='binary')
    }else{
        out$indC <- NULL
    }
}

out$ind.all <- list()
for(t in unique(out$tiall[,'type'])){
    out$ind.all[[t]] <- list()
    for(i in unique(out$tiall[out$tiall[,'type']==t,'instr'])){
        cat('t=',t,'; i=',i,'; range:',range(which(out$tiall[,'instr']==i & out$tiall[,'type']==t)),'\n')
        out$ind.all[[t]][[i]]  <- which(out$tiall[,'instr']==i & out$tiall[,'type']==t)
    }
}
out$rvs <- 0
if(length(out$astrometry)>0) out$rvs <- out$astrometry[out$iref,'radial_velocity']*1e3
###RUWE-related quantities
if(ruweDR>1 & nrow(out$astrometry)>1){
   source('ruwe_estimate.R')
}else{
    ruweDR <- 0
}
