par <- cn <- c()
par.stat <- e0$out$par.stat[[paste0('sig',Nsig)]]
if(!any(rownames(par.stat)==xlow)){
     par.stat <- sapply(1:(ncol(mc)-2),function(i) data.distr(as.numeric(mc[,i]),lp=mc[,'logpost'],plotf=FALSE))
     colnames(par.stat) <- colnames(mc)[1:(ncol(mc)-2)]
}

#par.kep <- par.stat['xopt',]
if(opt=='map' | !any(row.names(par.stat)==opt)){
   par.kep <- e0$par.opt
}else{
   par.kep <- par.stat[opt,]
   par.opt <- par.stat[opt,]
}
literature <- read.table('../data/code/stellar_mass_manual.txt',header=TRUE)
Ns <- min(nrow(mc),1e5)
ll <- mc[,'loglike']
lp <- mc[,'logpost']
#Mstar <- as.numeric(as.character(tab[k5,'Ms']))
#
if(TRUE){
   f3 <- list.files(paste0('../data/combined/',target0),pattern='htg$',full.name=TRUE)
   f4 <- list.files(paste0('../data/combined/',target0),pattern='astro$',full.name=TRUE)
   if(length(f3)>0){
     e0$out$astrometry <- read.table(f3[1],header=TRUE)
   }else if(length(f4)>0){
     e0$out$astrometry <- read.table(f4[1],header=TRUE)
   }
}
Mstar <- e0$out$astrometry[1,'mass']
if(e0$out$astrometry[1,'mass.upper']<Mstar | e0$out$astrometry[1,'mass.lower']>Mstar){
   e0$out$astrometry[1,'mass.upper'] <- e0$out$astrometry[1,'mass.upper']+Mstar
   e0$out$astrometry[1,'mass.lower'] <- -e0$out$astrometry[1,'mass.lower']+Mstar
}
eMstar <- 0.5*(e0$out$astrometry[1,'mass.upper']-e0$out$astrometry[1,'mass.lower'])
cat('Mstar0=',Mstar,'\n')
cat('eMstar0=',eMstar,'\n')
if(target0=='HD131664'){
    Mstar <- 1.060
    eMstar <- 0.129
}
if(target0=='HD190406'){
    Mstar <- 1.080
    eMstar <- 0.137
}
if(target0=='HD16160'){
    Mstar <- 0.780
    eMstar <- 0.091
}
if(target0=='HD161797'){
##First Results from the Hertzsprung SONG Telescope: Asteroseismology of the G5 Subgiant Star {\ensuremath{\mu}} Herculis by grundahl17
    Mstar <- 1.11
    eMstar <- 0.01
}
if(target0=='HD182488'){
    Mstar <- 0.930
    eMstar <- 0.113
}
if(target0=='HD190360'){
    Mstar <- 0.980
    eMstar <- 0.120
}
if(target0=='HD39587'){
    Mstar <- 1.100
    eMstar <- 0.134
}
if(target0=='HD42581'){
    Mstar <-  0.544
    eMstar <-  0.041
}
if(target0=='HD4747'){
    Mstar <- 0.910
    eMstar <- 0.114
}
if(target0=='HIP2552'){
    Mstar <- 0.530
    eMstar <- 0.083
}
ii <- which(target0==literature[,1])
if(length(ii)>0){
    cat('use stellar mass from stellar_mass_manual.txt\n')
    Mstar <- literature[ii,'m']
    eMstar <- 0.5*(literature[ii,'em.lower']+literature[ii,'em.upper'])
}
cat('Mstar1=',Mstar,'\n')
cat('eMstar1=',eMstar,'\n')
if(any(grepl('Mstar',names(par.opt)))){
   Mstar <- par.stat['mean','Mstar']
   eMstar <- par.stat['sd','Mstar']
}
cat('Mstar2=',Mstar,'\n')
cat('eMstar2=',eMstar,'\n')
ind <- which(!is.na(Mstar))[1]
Mstar <- Mstar[ind]
eMstar <- eMstar[ind]
mstars <- c(mstars,Mstar)
emstars <- c(emstars,eMstar)
Nkep <- 7
Npar <- length(par.opt)
Nval <- Nkep+6
indd <- sort(sample(1:nrow(mc),Ns))
set.seed(99)
ms <- rnorm(2*Ns,Mstar,eMstar)
ms <- ms[ms>0][1:Ns]
###convert period parameters
ind <- grep('per',colnames(par.stat))
inds <- sort(par.stat['mean',ind],index.return=TRUE)$ix
for(i in 1:length(inds)){
    tt <- gsub(' B','B',target.name)
    tt <- gsub(' C','C',tt)
    parshow <- c(paste(tt,pp[i]),paste0('\\colhead{',name.alt[k5],'}'))
    par.kep0 <- par.kep
    par.stat0 <- par.stat
###            ind.out <- c(ind.out,ind2)
    j <- inds[i]
    ind.kep <- (j-1)*Nkep+(1:Nkep)
    if(i==length(inds)){
        ind.val <- c(ind.kep,(Nkep*length(inds)):Npar)
    }else{
        ind.val <- ind.kep
    }
    for(k in ind.val){
        if(k%%Nkep==1 | k%%Nkep==3){
            Ndig <- 3
        }else if(k%%Nkep==4 | k%%Nkep==0){
            Ndig <- 2
        }else{
            Ndig <- Ndig0
        }
        pn <- colnames(mc)
        if(k%%Nkep==1 & grepl('per',pn[k])){
            ps <- exp(mc[indd,k])
            val <- data.distr(ps,plot=FALSE)
            if(length(val)>nrow(par.stat)){
                ind <- match(rownames(par.stat),names(val))
                par.stat[,k] <- val[ind]/yr2d
            }else{
                par.stat[,k] <- val/yr2d
            }
            popt <- par.kep[k] <- exp(par.kep[k])/yr2d
            par.stat[1,k] <- popt
	    cat('popt=',popt,'\n')
        }
        if(k%%Nkep==2 & grepl('^K',pn[k])){
            ks <- mc[indd,k]
            kopt <- par.kep[k]
            cat('kopt=',kopt,'\n')
        }

        if(k%%Nkep==3 & grepl('^e',pn[k])){
            es <- mc[indd,k]
            eopt <- par.kep[k]
        }
        if(k%%Nkep==6 & grepl('^I',pn[k])){
            incs <- mc[indd,k]%%pi
            incopt <- par.kep[k]
            par.stat[,k] <- par.stat[,k]*180/pi
            par.kep[k] <- par.kep[k]*180/pi
        }
        if(any(k%%Nkep==c(4,5,0)) & grepl('^omega|^Omega|^M',pn[k])){
            par.stat[,k] <- par.stat[,k]*180/pi
            par.kep[k] <- par.kep[k]*180/pi
                                        #                out[,k] <- out[,k]*180/pi
        }

####transform from logJ to jitter
                                        #            if(grepl('logJ',names(par.opt)[k])){
                                        #                jitter <- data.distr(exp(mc[,k]),plot=FALSE)
                                        #                par.stat[,k] <- jitter
                                        #            }

                                        #            cat('par:',colnames(par.stat)[k],'\n')
        for(vt in value.type){
            if(vt=='ms'){
                cn <- c(cn,paste0(colnames(par.stat)[k],c('.mean','.sd','.sd')))
                parshow <- c(parshow,paste0('$',show.digit(par.stat['mean',k],Ndig,type),'\\pm',show.digit(par.stat['sd',k],Ndig,type),'$'))
                par <- c(par,par.stat['mean',k],par.stat['sd',k],par.stat['sd',k])
            }else{
	        if(opt=='map' | !any(row.names(par.stat)==opt)){
		    cn <- c(cn,paste0(colnames(par.stat)[k],c('.opt','.lower','.upper')))
		}else{
		    cn <- c(cn,paste0(colnames(par.stat)[k],c('.opt','.lower','.upper')))
#		    cn <- c(cn,paste0(colnames(par.stat)[k],c(paste0('.',opt),'.lower','.upper')))
                }		    
                dpar.lower <- -par.stat[xlow,k]+par.kep[k]
                dpar.upper <- par.stat[xup,k]-par.kep[k]
                Ndig1 <- Ndig
#                if(abs(dpar.lower)< 10^(-Ndig) | abs(dpar.upper) < 10^(-Ndig)){
#                    Ndig1 <- -min(floor(log(c(dpar.lower,dpar.upper))))
#                }
                parshow <- c(parshow,paste0('$',show.digit(par.kep[k],Ndig,type),'_{-',show.digit(dpar.lower,Ndig1,type),'}^{+',show.digit(dpar.upper,Ndig1,type),'}$'))

                par <- c(par,par.kep[k],dpar.lower,dpar.upper)
            }
        }
    }
    Tp <- M02Tp(mc[,'Mo1'],tmin,exp(mc[,'per1']))
    Tp.opt <- M02Tp(par.opt['Mo1'],tmin,exp(par.opt['per1']))
                                        #                Tp.mean <- mean(Tp)
                                        #                Tp.sd <- sd(Tp)
    Ndig2 <- 0
    tps <- data.distr(Tp,plot=FALSE)
    for(vt in value.type){
        if(vt=='ms'){
            cn <- c(cn,paste0('Tp',inds[i],c('.mean','.sd','.sd')))
            parshow <- c(parshow,paste0('$',show.digit(tps['mean'],Ndig2,type),'\\pm',show.digit(tps['sd'],Ndig2,type),'$'))
            par <- c(par,tps['mean'],tps['sd'],tps['sd'])
        }
        if(vt=='mq'){
            if(opt=='map' | !any(row.names(par.stat)==opt)){
	          cn <- c(cn,paste0('Tp',inds[i],c('.opt','.lower','.upper')))
            }else{
	          cn <- c(cn,paste0('Tp',inds[i],c('.opt','.lower','.upper')))
#	          cn <- c(cn,paste0('Tp',inds[i],c(paste0('.',opt),'.lower','.upper')))
            }
            dpar.lower <- -tps[xlow]+Tp.opt
            dpar.upper <- tps[xup]-Tp.opt
            parshow <- c(parshow,paste0('$',show.digit(Tp.opt,Ndig2,type),'_{-',show.digit(dpar.lower,Ndig2,type),'}^{+',show.digit(dpar.upper,Ndig2,type),'}$'))
            par <- c(par,Tp.opt,dpar.lower,dpar.upper)
        }
    }
###calculate semi-major axis
###calculate mass
###https://iopscience.iop.org/article/10.3847/1538-3881/aa6d59/pdf
    cat('k,p,e,m=',kopt,popt,eopt,Mstar,'\n')
    ma.opt <- k2m(kopt,popt*yr2d,eopt,Mstar,Inc=incopt,more=TRUE)
    cat('ma.opt=',ma.opt$mj,'\n')
    ma <- k2m(ks,ps,es,ms,incs,more=TRUE)
    cat('ks[1],ps[1],es[1],ms[1]=',ks[1],ps[1],es[1],'\n')
    mp <- ma$mj
    mp.opt <- ma.opt$mj
    a <- ma$a
    mp.stat <- data.distr(mp,plot=FALSE)
    aopt <- ma.opt$a
    a.stat <- data.distr(a,plot=FALSE)
    for(vt in value.type){
        if(vt=='ms'){
            cn <- c(cn,paste0('mpJ',inds[i],c('.mean','.sd','.sd')))
            str.m <- paste0('$',show.digit(mp.stat['mean'],Nsig),'\\pm',show.digit(mp.stat['sd'],Nsig),'$')
            par <- c(par,mp.stat['mean'],mp.stat['sd'],mp.stat['sd'])
        }else{
            if(opt=='map' | !any(row.names(par.stat)==opt)){
         	    cn <- c(cn,paste0('mpJ',inds[i],c('.opt','.lower','.upper')))
            }else{
         	    cn <- c(cn,paste0('mpJ',inds[i],c('.opt','.lower','.upper')))
#         	    cn <- c(cn,paste0('mpJ',inds[i],c(paste0('.',opt),'.lower','.upper')))
	    }
###                str.m <-  paste0('$',show.digit(ma.opt$me,Ndig,type),'[',show.digit(mp.stat[xlow],Ndig,type),',',show.digit(mp.stat[xup],Ndig,type),']$')
            dm.lower <- -mp.stat[xlow]+mp.stat[opt]
            dm.upper <- mp.stat[xup]-mp.stat[opt]
#           if(dm.lower<0 | dm.upper<0) stop()
            str.m <- paste0('$',show.digit(mp.opt,Ndig0,type),'_{-',show.digit(dm.lower,Ndig0,type),'}^{+',show.digit(dm.upper,Ndig0,type),'}$')
                                        #                str.a <-  paste0('$',show.digit(ma.opt$a,Ndig0+1,type),'[',show.digit(a.stat['x10per'],Ndig0+1,type),',',show.digit(a.stat[xup],Ndig0+1,type),']$')
            par <- c(par,mp.opt,dm.lower,dm.upper)
        }
        cat(str.m,'\n')
        parshow <- c(parshow,str.m)
    }

    for(vt in value.type){
        if(vt=='ms'){
            cn <- c(cn,paste0('a',inds[i],c('.mean','.sd','.sd')))
            str.a <- paste0('$',show.digit(a.stat['mean'],Ndig0),'\\pm',show.digit(a.stat['sd'],Ndig0),'$')
            par <- c(par,a.stat['mean'],a.stat['sd'],a.stat['sd'])
        }else{
            if(opt=='map' | !any(row.names(par.stat)==opt)){
	         cn <- c(cn,paste0('a',inds[i],c('.opt','.lower','.upper')))
            }else{
	         cn <- c(cn,paste0('a',inds[i],c('.opt','.lower','.upper')))
#	         cn <- c(cn,paste0('a',inds[i],c(paste0('.',opt),'.lower','.upper')))
	    }
            da.lower <- -a.stat[xlow]+ma.opt$a
            da.upper <- a.stat[xup]-ma.opt$a
            str.a <- paste0('$',show.digit(ma.opt$a,Ndig0,type),'_{-',show.digit(da.lower,Ndig0,type),'}^{+',show.digit(da.upper,Ndig0,type),'}$')
            par <- c(par,aopt,da.lower,da.upper)
        }
        parshow <- c(parshow,str.a)
    }
}
table.out <- cbind(table.out,parshow)
names(par) <- cn
par.stat <- par.stat0
par.kep <- par.kep0
if(length(par.out)==0){
   par.out <- rbind(par.out,par)
}else{
   par.name <- names(par)
   out.name <- colnames(par.out)
   ind.par <- match(out.name,par.name)
   ind.out <- match(par.name,out.name)
   if(!any(is.na(ind.out))){
     par.out <- rbind(par.out,par[ind.par])
   }else{
     par.name <- names(par)
     out.name <- colnames(par.out)
     ind.out <- match(par.name,out.name)
     if(any(is.na(ind.out))){
         array.add <- array(NA,dim=c(nrow(par.out),length(which(is.na(ind.out)))))
         colnames(array.add) <- par.name[is.na(ind.out)]
         par.out <- cbind(par.out,array.add)
         ind.par <- match(colnames(par.out),names(par))
         par.out <- rbind(par.out,par[ind.par])
     }
  }
}
cat('nrow(par.out)=',nrow(par.out),'\n')
#if(nrow(par.out)>1 & nrow(par.out)!=k5) stop('par.out rows larger than the number of stars, something is wrong!')