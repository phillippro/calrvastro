mstars <- emstars <- par.out <- c()
    ff <- fs[k5]
#    cat('ff:\n',ff,'\n')
    target0 <- targets0[k5]
    target.name <- stars[k5]
     if(!file.exists(ff)){
            cmd <- paste0('scp tdlffb@data.hpc.sjtu.edu.cn:Documents/projects/dwarfs/rvastro/',ff,' results/',target0,'/')
	    cat(cmd,'\n')
	    system(cmd)
     }
#    if(any(target0==exception)){
     if(file.exists(ff)){
    load(ff,env=e0 <- new.env())
#    pp <- any(names(e0$out)=='astrometry')	
    nastro <- e0$out$Nastro
    if(is.null(nastro)) nastro <- 2
    cat('\nwhether astrometry is used in ',ff,'?',nastro>0,'!\n')
    if(nastro>0){
        Nastro <- 2
        astro <- astrometry <- 8
        if(!exists('basis')) basis <- 'natural'
        if(!exists('bases')) bases <- rep('natural',10)
        cat('star=',target.name,'\n')
        cat('ins=',e0$ins,'\n')
        cat('ns=',unlist(lapply(e0$ins,function(i) nrow(e0$out[[i]][['RV']]))),'\n')
        cat('Tp=',unlist(lapply(e0$ins,function(i) max(e0$out[[i]][['RV']][,1])-min(e0$out[[i]][['RV']][,1]))),'\n')
        inss <- c(inss,e0$ins)
        offset <- TRUE
        Nsig <- e0$Nsig
        if(target0=='GJP9066' | target0=='GJ300' | target0=='GJ880' | target0=='GJ251') Nsig <- 1
        trv.all <- e0$trv.all
        basis <- e0$basis
        ins <- e0$ins
        out <- e0$out
        tmin <- e0$tmin
        par.opt <- e0$out$par.stat[[paste0('sig',Nsig)]]['xopt',]
        mcmc <- e0$out$mcmc.opt
        e0$out$mcmc.all <- NULL
#        e0$out$mcmc.opt <- NULL
        out.phase[[ff]] <- list(par.opt=par.opt,ins=ins,Nsig=Nsig,trv.all=trv.all,basis=basis,par.stat=e0$out$par.stat)
        rv.list <- list()
        for(i in ins){
            rv.list[[i]] <- out[[i]]$RV
        }
        out.phase[[ff]]$rv.list <- rv.list

    ins.name <- ins
    if(any(ins=='HARPS')){
        th <- rv.list$HARPS[,1]
        if(min(th)%%2400000>57174.5) ins.name[ins.name=='HARPS'] <- 'HARPSpost'
        if(max(th)%%2400000<=57174.5) ins.name[ins.name=='HARPS'] <- 'HARPSpre'
    }
    if(any(ins=='PFS') & FALSE){
        th <- rv.list$PFS[,1]
        cat('range(th)=',range(th),'\n')
        ins.name[ins=='PFS'] <- 'PFSpre'
        if(min(th)%%2400000>58157) ins.name[ins=='PFS'] <- 'PFSpost'
        if(max(th)%%2400000<=58157) ins.name[ins=='PFS'] <- 'PFSpre'
        cat('ins.name=',ins.name,'\n')
    }
    if(any(target0==exception)) ins.all <- c(ins.all,ins.name)
    source("generate_table.R",local=TRUE)
inss <- unique(inss)
#colnames(par.out) <- c('Star','P','Plow','Pup','K','Klow','Kup','e','elow','eup','omega','omega.low','omega.up','Omega','Omega.low','Omega.up','mp','mp.lower','mp.upper')
   par.out <- cbind(target0,mstars,emstars,par.out)
    par.out[,1] <- gsub(' ','',par.out[,1])
}
}