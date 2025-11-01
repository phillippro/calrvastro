fs <- list.files(path='pars',pattern=paste0('^',target,'.+par$'),full.name=TRUE)
replacePar <- FALSE
if(any(grepl(target,fs))){
     ff <- fs[grep(target,fs)][1]
     pars <- read.table(ff)
     paropt <- as.numeric(pars[,2])
     ns <- as.character(pars[,1])
     ns[grep('b\\d',ns)] <- paste0('b_',ins)
                                   names(paropt) <- ns
     ind <- grep('jitter',ns)
     if(length(ind)>0){
         ns[ind] <- gsub('jitter','logJ',ns[ind])
         paropt[ind] <- log(paropt[ind])
         names(paropt) <- ns
     }
     replacePar <- TRUE
     nn <- intersect(names(startvalue),ns)
     cat('replacing parameters by pars/',paste0(target,'.par:'),nn,'\n')
     startvalue[nn] <- paropt[nn]
     startvalue1 <- startvalue
     kep <- RV.kepler(startvalue)
     if(out$Nrv>0){
         for(jj in 1:length(out$ins)){
             i <- out$ins[jj]
             startvalue1[paste0('b',jj)] <- startvalue1[paste0('b',jj)]+mean(out[[i]]$RV[,2]-kep$rv[[i]])
         }
         if(length(out$data.binary)>0){
             for(kk in 1:length(out$ins.binary)){
                 i <- out$ins.binary[kk]
                 inds <- out$ind.binary[[i]]
                 startvalue1[paste0('b_',i)] <- startvalue1[paste0('b_',i)]+mean(out$data.binary[inds,2]-kep$drvC[inds])
             }
         }
     }
     llike <- loglikelihood(startvalue)
     llike1 <- loglikelihood(startvalue1)
     cat('llike=',llike,';llike1=',llike1,'\n')
     if(llike1>llike){
         startvalue <- startvalue1
         llike <- llike1
         cat('update RV offsets!\n')
     }
     if(any('dra'==names(startvalue))){
         dastro <- astrometry.epoch(startvalue,tt=out$astrometry[,1])$epoch
         dastro1 <- dastro[out$iref,]#reflex motion assume the offset is at the reference epoch defined by iref
         startvalue0 <- startvalue1 <- startvalue2 <- startvalue
         if(any('dra'==names(startvalue))) startvalue1['dra'] <- startvalue1['dra']+dastro1$dra
         if(any('ddec'==names(startvalue))) startvalue1['ddec'] <- startvalue1['ddec']+dastro1$ddec
         if(any('dplx'==names(startvalue))) startvalue1['dplx'] <- startvalue1['dplx']+dastro1$dplx
         if(any('dpmra'==names(startvalue))) startvalue1['dpmra'] <- startvalue1['dpmra']+dastro1$dpmra
         if(any('dpmdec'==names(startvalue))) startvalue1['dpmdec'] <- startvalue1['dpmdec']+dastro1$dpmdec
         if(any('drv'==names(startvalue))) startvalue1['drv'] <- startvalue1['drv']+dastro1$drv
         ll1 <- loglikelihood(startvalue1)
         ind.max <- which.max(c(llike,ll1))
         cat('ll1=',ll1,'\n')
         if(length(out$igdr2)>0){
             dastro2 <- dastro[out$igdr2,]#reflex motion assume the offset is at the Gaia DR2 epoch
             dastro2
             if(any('dra'==names(startvalue))) dastro2$dra <- dastro2$dra+startvalue2['dra']
             if(any('ddec'==names(startvalue))) dastro2$ddec <- dastro2$ddec+startvalue2['ddec']
             if(any('dplx'==names(startvalue))) dastro2$dplx <- dastro2$dplx+startvalue2['dplx']
             if(any('dpmra'==names(startvalue))) dastro2$dpmra <- dastro2$dpmra+startvalue2['dpmra']
             if(any('dpmdec'==names(startvalue))) dastro2$dpmdec <- dastro2$dpmdec+startvalue2['dpmdec']
             if(any('drv'==names(startvalue))) dastro2$drv <- dastro2$drv+startvalue2['drv']
             obs <- out$astrometry[out$igdr2,]
             obs['ra'] <- as.numeric(obs['ra']-dastro2$dra/3.6e6/cos(obs['dec']/180*pi))
             obs['dec'] <- as.numeric(obs['dec']-dastro2$ddec/3.6e6)
             obs['parallax'] <- as.numeric(obs['parallax']-dastro2$dplx)
             obs['pmra'] <- as.numeric(obs['pmra']-dastro2$dpmra)
             obs['pmdec'] <- as.numeric(obs['pmdec']-dastro2$dpmdec)
             obs['radial_velocity'] <- as.numeric(obs['radial_velocity']-dastro2$drv)
             barycenter <- obs.lin.prop(obs,out$astrometry[,1] - out$astrometry[out$igdr2,1])
             dastro <- AstroDiff(barycenter[out$iref,],out$astrometry[out$iref,])
             if(any('dra'==names(startvalue))) startvalue2['dra'] <- dastro[1]
             if(any('ddec'==names(startvalue))) startvalue2['ddec'] <- dastro[2]
             if(any('dplx'==names(startvalue))) startvalue2['dplx'] <- dastro[3]
             if(any('dpmra'==names(startvalue))) startvalue2['dpmra'] <- dastro[4]
             if(any('dpmdec'==names(startvalue))) startvalue2['dpmdec'] <- dastro[5]
             if(any('drv'==names(startvalue))) startvalue2['drv'] <- dastro[6]
             ll2 <- loglikelihood(startvalue2)
             ind.max <- which.max(c(llike,ll1,ll2))
             cat('ll2=',ll2,'\n')
         }
         if(ind.max==1){
             startvalue <- startvalue0
             cat('use original par file!')
         }else if(ind.max==2){
             startvalue <- startvalue1
             cat('use updated par file by considering reflex motion correction at reference epoch!\n')
         }else if(ind.max==3){
             startvalue <- startvalue2
             cat('use updated par file by considering reflex motion correction at GDR2 epoch!\n')
         }
     }
     par.hot <- startvalue
}
