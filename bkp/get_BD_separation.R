source('mcmc_func.R')
Ncores <- 48
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
info <- read.csv('../data/code/SigClassify/common_withSig1.csv',sep=',')
tab <- read.table('companion_files.txt')[,1]
tab <- tab[grep('^HD74014',tab)]
#tab <- tab[length(tab)-(10:0)]
#tab <- tab[1:5]
mlow <- 5
mup <- 120
##error.type <- 'meansd'
error.type <- 'med1sig'
bd <- c()
types <- c()
cal <- cbind(2022,5,1)
jd <- sum(time_Cal2JD(cal))
targets <- gsub('\\/.+','',tab)
#out.bd <- foreach(i3 = 1:length(tab), .combine='rbind') %dopar% {
for(i3 in 1:length(tab)){
    star <- targets[i3]
    f <- tab[i3]
    cat('\n',star,'\n')
    fobj <- gsub('pdf','Robj',f)
    fobj1 <- gsub('.+\\/','',fobj)
    fs <- list.files(path=paste0('results/',star))
    if(!any(fs==fobj1)){
	cat('scp data file from HPC!\n')
	cmd <- paste0('rsync -azP tdlffb@login.hpc.sjtu.edu.cn:Documents/projects/dwarfs/rvastro/results/',star,'/',fobj1,' results/',star,'/')
	system(cmd)
        fs <- list.files(path=paste0('results/',star))
    }
#stop()
    if(any(fs==fobj1)){
        cat('AstroCluster!\n')
        load(paste0('results/',star,'/',fobj1), envir=e <- new.env())
        out <- e$out
        cat('0out whether exists:',exists('out'),'\n')
        tmin <- e$tmin
        par.opt <- e$par.opt
        etaH <- e$etaH
        etaG <- e$etaG
        source('mcmc_func.R')
#        load(paste0('results/',star,'/',fobj1))
        if(out$Nastro==0){
            cat('No astrometry data found!\n')
        }else{
            shortP <- FALSE
            if(grepl('shortPTRUE',fobj1)){
		shortP <- TRUE
            }
            if(any(names(out)=='Mstar')){
                Mstar <- out$Mstar
            }else{
                Mstar <- out$astrometry[1,'Mstar']
            }
            eMstar <- 0.5*(out$astrometry[1,'mass.lower']+out$astrometry[1,'mass.upper'])
            if(any(names(out)=='eMstar')){
                if(out$eMstar>0) eMstar <- out$eMstar
            }

        cat('1out whether exists:',exists('out'),'\n')


            if(length(eMstar)==0 | is.null(eMstar)) eMstar <- 0.1*Mstar
            Np <- length(grep('per',names(par.opt)))
            for(k in 1:Np){
                msini <- K2msini.full(par.opt[paste0('K',k)],exp(par.opt[paste0('per',k)]),par.opt[paste0('e',k)],Ms=Mstar)
                Mp <- msini$ms/sin(par.opt[paste0('Inc',k)])
                Mcomp <- msini$mj/sin(par.opt[paste0('Inc',k)])#companion mass in jupiter mass
                if(Mcomp>mlow & Mcomp<mup){
                    a <- msini$a#au
####MC -> parameter uncertainty
                    mc <- out$mcmc.opt[[paste0('sig',Np)]]
                    mc <- mc[sample(1:nrow(mc),1e4),]
                    mstars <- rnorm(nrow(mc),Mstar,eMstar)
                    mstars[mstars<0] <- median(mstars)
                    ps <- exp(mc[,paste0('per',k)])
                    msinis <- K2msini.full(mc[,paste0('K',k)],exp(mc[,paste0('per',k)]),mc[,paste0('e',k)],Ms=mstars)
                    mps <- msinis$mj/sin(mc[,paste0('Inc',k)])
                    as <- msinis$a
                    if(error.type!='meansd'){
                        aa <- data.distr(as,lp=mc[,'loglike'],plotf=FALSE)[c('med','xminus.1sig','xplus.1sig')]
                        pp <- data.distr(ps,lp=mc[,'loglike'],plotf=FALSE)[c('med','xminus.1sig','xplus.1sig')]
                        mm <- data.distr(mps,lp=mc[,'loglike'],plotf=FALSE)[c('med','xminus.1sig','xplus.1sig')]
                    }else{
                        aa <- data.distr(as,lp=mc[,'loglike'],plotf=FALSE)[c('mean','xsd')]
                        pp <- data.distr(ps,lp=mc[,'loglike'],plotf=FALSE)[c('mean','xsd')]
                        mm <- data.distr(mps,lp=mc[,'loglike'],plotf=FALSE)[c('mean','xsd')]
                    }
                    orb <- par.opt[paste0(c("per", "K", "e",  "omega", "Mo", "Inc",   "Omega"),k)]
                    Tp <- M02Tp(orb[5],tmin,exp(orb[1]))
                    orb[1] <- exp(orb[1])/365.25#yr
                    orb[5] <- Tp
                    eta <- Mp/(Mstar+Mp)#
                    par.opt1 <- par.opt
                    if(Np>1){
                        jj <- (1:Np)[-k]
                        par.opt1[paste0('K',jj)] <- 0
                    }
                    cat('3out whether exists:',exists('out'),'\n')
                    if(!any(names(out)=='indrel')){
                        indrel <- list()
                    }else{
                        indrel <- out$indrel
                    }
                    astro <- astrometry.rel(par.opt1,tt=jd%%24e5-tmin%%24e5)
                    rho <- sqrt(astro$ra^2+astro$dec^2)/eta#[mas]
                    x <- astro$ra
                    y <- astro$dec
                    pa <- atan2(x,y)*180/pi#deg
                    xs <- ys <- rhos <- pas <- c()
                    for(k3 in 1:nrow(mc)){
                        par <- mc[k3,1:length(par.opt)]
                        if(Np>1){
                            jj <- (1:Np)[-k]
                            par[paste0('K',jj)] <- 0
                        }
                        astro <- astrometry.rel(par,tt=jd%%24e5-tmin%%24e5)
                        rhos <- c(rhos,sqrt(astro$ra^2+astro$dec^2)/eta)
                        xs <- c(xs,astro$ra)
                        ys <- c(ys,astro$dec)
                    }
                    pas <- atan2(xs,ys)
                    if(error.type!='meansd'){
                        xx <- data.distr(xs,lp=mc[,'loglike'],plotf=FALSE)[c('med','xminus.1sig','xplus.1sig')]
                        yy <- data.distr(ys,lp=mc[,'loglike'],plotf=FALSE)[c('med','xminus.1sig','xplus.1sig')]
                        rr <- data.distr(rhos,lp=mc[,'loglike'],plotf=FALSE)[c('med','xminus.1sig','xplus.1sig')]
                        pp <- data.distr(pas,lp=mc[,'loglike'],plotf=FALSE)[c('med','xminus.1sig','xplus.1sig')]
                    }else{
                        xx <- data.distr(xs,lp=mc[,'loglike'],plotf=FALSE)[c('mean','sd')]
                        yy <- data.distr(ys,lp=mc[,'loglike'],plotf=FALSE)[c('mean','sd')]
                        rr <- data.distr(rhos,lp=mc[,'loglike'],plotf=FALSE)[c('mean','sd')]
                        pp <- data.distr(pas,lp=mc[,'loglike'],plotf=FALSE)[c('mean','sd')]
                    }

#####get astrometry
                    ind <- which(info[,'target']==star)
                    type <- "NA"
                    if(length(ind)>0){
                        tmp <- as.numeric(info[ind[1],c('ra','dec','parallax','pmra','pmdec','radial_velocity','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','teff_val','radius_val','lum_val')])
                        type <- info[ind[1],'Type']
                    }else{
                        tmp <- rep(NA,12)
                        fs <- list.files(path=paste0('../data/combined/',star),pattern='hg3$',full.name=TRUE)
                        if(length(fs)>0){
                            cn <- c('ra','dec','parallax','pmra','pmdec','radial_velocity','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','teff_val','radius_val','lum_val')
                            astrometry <- read.table(fs[1],header=TRUE)
                            ii <- match(cn,colnames(astrometry))
                            tmp[!is.na(ii)] <- as.numeric(astrometry[2,ii[!is.na(ii)]])
                        }
                    }
#####
#                   types <- c(types,type)
                                       #	  bd <- rbind(bd,c(star,rho,x,y,Mstar,mm,aa,pp,orb,tmp))
                    bd <- rbind(bd,c(type,star,rr,pp,xx,yy,Mstar,mm,aa,pp,orb,tmp))
                }
            }
        }
        out.bd <- bd
#        rm(e)
#        bd
    }else{
        cat('Not on the AstroCluster!\n')
        out.bd <- c()
    }
}
#bd <- cbind(types,bd)
bd <- out.bd
if(error.type!='meansd'){
    colnames(bd) <- c('StellarType','Star','rho.mas','erho.mas.minor1sig','erho.mas.plus1sig','PA.deg','ePA.deg.minor1sig','ePA.deg.plus1sig','x.mas','ex.mas.minor1sig','ex.mas.plus1sig','y.mas','ey.mas.minor1sig','ey.mas.plus1sig','Mstar.Msun','Mcomp.Mjup','Mcomp.Mjup.minor1sig','Mcomp.Mjup.plus1sig','a.au','a.au.minor1sig','a.au.plus1sig','P.day','P.day.minor1sig','P.day.plus1sig','Period.yr','K.ms','ecc','omega.rad','Tp.BJD','Inc.rad','Omega.rad','ra','dec','parallax','pmra','pmdec','radial_velocity','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','teff_val','radius_val','lum_val')
}else{
    colnames(bd) <- c('StellarType','Star','rho.mas','erho.mas','PA.deg','ePA.deg','x.mas','ex.mas','y.mas','ey.mas','Mstar.Msun','Mcomp.Mjup','Mcomp.sd','a.au','a.au.sd','P.day','P.day.sd','P.yr','K.ms','ecc','omega.rad','Tp.BJD','Inc.rad','Omega.rad','ra','dec','parallax','pmra','pmdec','radial_velocity','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','teff_val','radius_val','lum_val')
}
fout <- paste0('companion_mcomp',mlow,'to',mup,'Mjup_',error.type,'.txt')
cat(fout,'\n')
write.table(bd,file=fout,quote=FALSE,row.names=FALSE)
