source('mcmc_func.R')
info <- read.csv('../data/code/SigClassify/common_withSig1.csv',sep=',')
tab <- read.table('BD_candidates.txt')[,1]
bd <- c()
types <- c()
targets <- gsub('\\/.+','',tab)
j3 <- which(targets=='HIP111958')
#for(i3 in j3:length(tab)){
for(i3 in 1:length(tab)){
     star <- targets[i3]
     f <- tab[i3]
     cat('\n',star,'\n')
     fobj <- gsub('pdf','Robj',f)
     fobj1 <- gsub('.+\\/','',fobj)
     path <- '/malibu/ffeng/astro_output/'
     fs <- list.files(path=paste0('/malibu/ffeng/astro_output/',star))
     if(!any(fs==fobj1)){
     fs <- list.files(path=paste0('results/',star))
     path <- 'results/'
     }
     if(any(fs==fobj1)){
          cat('HPC!\n')
#          load(paste0('results/',star,'/',fobj1), envir=e <- new.env())
          load(paste0(path,star,'/',fobj1))
          out$Nm <- 0
	  astro <- astrometry.rel(par.opt,tt=2459488.5%%24e5-tmin%%24e5)
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
	  types <- c(types,as.character(type))
	  bd <- rbind(bd,c(star,rho,x,y,Mstar,Mcomp,a,orb,tmp))
	  rm('Mstar')
     }else{
          cat('RACK!\n')
     }
}
bd <- cbind(types,bd)
colnames(bd) <- c('StellarType','Star','rho.mas','x.mas','y.mas','Mstar.Msun','Mcomp.Mjup','a.au','Period.yr','K.ms','ecc','omega.rad','Tp.BJD','Inc.rad','Omega.rad','ra','dec','parallax','pmra','pmdec','radial_velocity','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','teff_val','radius_val','lum_val')
write.table(bd,file='bd_rack.txt',quote=FALSE,row.names=FALSE)