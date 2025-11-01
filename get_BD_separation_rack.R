source('mcmc_func.R')
info <- read.csv('../data/code/SigClassify/common_withSig1.csv',sep=',')
tab <- read.table('BD_candidates.txt')[,1]
bd <- c()
types <- c()
targets <- gsub('\\/.+','',tab)
for(i3 in 1:length(tab)){
     star <- targets[i3]
     dir <- paste0('/malibu/ffeng/astro_output/',star)
     f <- tab[i3]
     cat('\n',star,'\n')
     fobj <- gsub('pdf','Robj',f)
     fobj1 <- gsub('.+\\/','',fobj)
     fs <- list.files(path=dir)
     if(any(fs==fobj1)){
	  if(!grepl('malibut',dir)){
          cat('HPC!\n')
          }else{
          cat('RACK!\n')
	  }
#          load(paste0('results/',star,'/',fobj1), envir=e <- new.env())
          load(paste0(dir,'/',fobj1))
	  astro <- astrometry.rel(par.opt,tt=2459488.5%%24e5-tmin%%24e5)#Oct. 1, 2021
	  if(!exists('Mstar')){
	  if(any(names(out)=='Mstar')){
	      Mstar <- out$Mstar
	  }else{
              Mstar <- out$astrometry[1,'Mstar']
	  }
	  }
	  Np <- length(grep('per',names(par.opt)))
	  for(k in 1:Np){
              msini <- K2msini.full(par.opt[paste0('K',k)],exp(par.opt[paste0('per',k)]),par.opt[paste0('e',k)],Ms=Mstar)
              Mp <- msini$ms/sin(par.opt[paste0('Inc',k)])
              Mcomp <- msini$mj/sin(par.opt[paste0('Inc',k)])#companion mass in jupiter mass
	      a <- msini$a#au
####MC -> parameter uncertainty
	      mc <- out$mcmc.opt[[paste0('sig',Np)]]
	      mc <- mc[sample(1:nrow(mc),1e4),]
	      mstars <- rnorm(nrow(mc),out$Mstar,out$eMstar)
      	      mstars[mstars<0] <- median(mstars)
	      ps <- exp(mc[,paste0('per',k)])
              msinis <- K2msini.full(mc[,paste0('K',k)],exp(mc[,paste0('per',k)]),mc[,paste0('e',k)],Ms=mstars)
	      mps <- msinis$mj/sin(mc[,paste0('Inc',k)])
	      as <- msinis$a
	      aa <- data.distr(as,lp=mc[,'loglike'],plotf=FALSE)[c('xopt','xminus.1sig','xplus.1sig')]
      	      pp <- data.distr(ps,lp=mc[,'loglike'],plotf=FALSE)[c('xopt','xminus.1sig','xplus.1sig')]
       	      mm <- data.distr(mps,lp=mc[,'loglike'],plotf=FALSE)[c('xopt','xminus.1sig','xplus.1sig')]
              if(Mcomp>10 & Mcomp<100) break()
	  }
	  orb <- par.opt[paste0(c("per", "K", "e",  "omega", "Mo", "Inc",   "Omega"),k)]
	  Tp <- M02Tp(orb[5],tmin,exp(orb[1]))
#	  cat('Tp=BJD',Tp,'\n')
	  orb[1] <- exp(orb[1])/365.25#yr
	  orb[5] <- Tp
          eta <- Mp/(Mstar+Mp)#
	  rho <- sqrt(astro$ra^2+astro$dec^2)/eta#[mas]
	  x <- astro$ra
  	  y <- astro$dec
	  ind <- which(info[,'target']==star)
	  type <- "NA"
	  if(length(ind)>0){
	      tmp <- as.numeric(info[ind[1],c('ra','dec','parallax','pmra','pmdec','radial_velocity','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','teff_val','radius_val','lum_val')])
	      type <- info[ind[1],'Type']
	  }else{
              tmp <- rep(NA,12)
	  }
	  types <- c(types,type)
#	  bd <- rbind(bd,c(star,rho,x,y,Mstar,Mcomp,a,orb,tmp))
	  bd <- rbind(bd,c(star,rho,x,y,Mstar,mm,aa,pp,orb,tmp))
	  rm('Mstar')
     }else{
	  if(!grepl('malibut',dir)){
          cat('RACK!\n')
          }else{
          cat('HPC!\n')
	  }
     }
}
bd <- cbind(types,bd)
#colnames(bd) <- c('StellarType','Star','rho.mas','x.mas','y.mas','Mstar.Msun','Mcomp.Mjup','a.au','Period.yr','K.ms','ecc','omega.rad','Tp.BJD','Inc.rad','Omega.rad','ra','dec','parallax','pmra','pmdec','radial_velocity','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','teff_val','radius_val','lum_val')
colnames(bd) <- c('StellarType','Star','rho.mas','x.mas','y.mas','Mstar.Msun','Mcomp.Mjup','Mcomp.Mjup.minor1sig','Mcomp.Mjup.plus1sig','a.au','a.au.minor1sig','a.au.plus1sig','P.day','P.day.minor1sig','P.day.plus1sig','Period.yr','K.ms','ecc','omega.rad','Tp.BJD','Inc.rad','Omega.rad','ra','dec','parallax','pmra','pmdec','radial_velocity','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','teff_val','radius_val','lum_val')
write.table(bd,file='bd_hpc.txt',quote=FALSE,row.names=FALSE)