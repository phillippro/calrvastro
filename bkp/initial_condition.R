#########################################################################
####this routine is to prepare initial conditions of parameters for DRAM
##########################################################################
###Keplerian parameters
Kmax <- max(100*(rmax-rmin),1e6)
Kmin <- 1e-6
nams <- c()
par.min <- par.max <- par.ini <- c()
if(Np>0){
    kep.ini <- kep.min <- kep.max <- c()
    for(j in 1:(Np+out$Nm)){
        basis <- bases[j]
        per.ini <- runif(1,log(Pmin),log(Pmax))
        per.min <- log(Pmin)
        per.max <- log(Pmax)
        Mo.min <- 0
        Mo.max <- 2*pi
        Mo.ini <- 1
        if(basis=='natural'){
            if(Ntr>0){
                ii <- which(ind.transit==j)
                if(length(ii)>0){
                    per.ini <- log(Ptransit[ii])
                    per.min <- log(Ptransit[ii]-5*ePtransit[ii])
                    per.max <- log(Ptransit[ii]+5*ePtransit[ii])
###if basis is natural, it is not suitable for using Tc because the proposed Mo is always not compatible with Tc
                    if(Ntr>0 & FALSE){
                        Mo.ini <- getM0(e=0,omega=pi/2,P=exp(per.ini),T=Tc[ind.transit[ii]],T0=tmin,type='primary')
                        Mo1 <- 2*pi*((tmin-Tc[ind.transit[ii]]-5*eTc[ind.transit[ii]])%%Ptransit[ind.transit[ii]])/Ptransit[ind.transit[ii]]
                        Mo2 <- 2*pi*((tmin-Tc[ind.transit[ii]]+5*eTc[ind.transit[ii]])%%Ptransit[ind.transit[ii]])/Ptransit[ind.transit[ii]]
                        Mo.min <- min(Mo1,Mo2)
                        Mo.max <- max(Mo1,Mo2)
                    }
                }
            }
            if(!is.null(RV.all)){
                Kini <- max(sd(RV.all),0.1)
            }else{
                Kini <- 10
            }
            kep.ini <- c(kep.ini,c(per.ini,Kini,0.1,pi/2,Mo.ini))#omega=pi/2 -> Tc=Tp (Eastman et al. 2013)
            kep.min <- c(kep.min,c(per.min,Kmin,0,0,Mo.min))
            kep.max <- c(kep.max,c(per.max,Kmax,1,2*pi,Mo.max))
            nams.kep <- c('per','K','e','omega','Mo')
            if(any(j==dP)){
                kep.ini[1] <- 0
                kep.min[1] <- -1000#min
                kep.max[1] <- 1000
                nams.kep[1] <- 'dP'
            }
        }else if(basis=='linear1'){
###lnP, lnK,sqrt(e)*sin(w),sqrt(e)*cos(w),Tc
###lnP, lnK,e*sin(omega),e*cos(omega),Tc; optional
            eini <- 0.1
            omegaini <- 0
            Tc.ini0 <- Tc.ini <- tmin
            Tc.min0 <- Tc.min <- tmin
            Tc.max0 <- Tc.max <- tmax+2*Pmax
            if(Ntr>0){
                ii <- which(ind.transit==j)
                if(length(ii)>0){
                    if(Ntr>0){
                        Tc.ini <- Tc[ii]
                        Tc.min <- Tc[ii]-5*eTc[ii]
                        Tc.max <- Tc[ii]+5*eTc[ii]
                    }
                    per.ini <- log(Ptransit[ind.transit[ii]])
                    per.min <- log(Ptransit[ind.transit[ii]]-5*ePtransit[ind.transit[ii]])
                    per.max <- log(Ptransit[ind.transit[ii]]+5*ePtransit[ind.transit[ii]])
                }
            }
            kep.ini <- c(kep.ini,c(per.ini,log(sd(RV.all)),sqrt(eini)*sin(omegaini),sqrt(eini)*cos(omegaini),Tc.ini))
            kep.min <- c(kep.min,c(per.min,log(Kmin),-1,-1,Tc.min))
            kep.max <- c(kep.max,c(per.max,log(Kmax),1,1,Tc.max))
            nams.kep <- c('per','lnK','sqresinw','sqrecosw','Tc')
        }else if(basis=='linear2'){
###dP, lnK,sqrt(e)*sin(w),sqrt(e)*cos(w),dTc
            eini <- 0.1
            omegaini <- 0
            ii <- which(ind.transit==j)
            dTc.ini <- 0
            dTc.min <- -5*eTc[ii]/alpha.dTc
            dTc.max <- 5*eTc[ii]/alpha.dTc
            dP.ini <- 0
            dP.min <- -5*ePtransit[ind.transit[ii]]/alpha.dP
            dP.max <- 5*ePtransit[ind.transit[ii]]/alpha.dP
            kep.ini <- c(kep.ini,c(dP.ini,log(sd(RV.all)),sqrt(eini)*sin(omegaini),sqrt(eini)*cos(omegaini),dTc.ini))
            kep.min <- c(kep.min,c(dP.min,log(Kmin),-1,-1,dTc.min))
            kep.max <- c(kep.max,c(dP.max,log(Kmax),1,1,dTc.max))
            nams.kep <- c('dP','lnK','sqresinw','sqrecosw','dTc')
        }

###add inclination and Omega
        if(out$astro.par>0){
            kep.ini <- c(kep.ini,c(pi/4,pi/4))
            kep.min <- c(kep.min,c(0,0))
            if(length(out$ins.rv)>0 | length(out$rel)>0 | length(out$timing)>0){
                kep.max <- c(kep.max,c(pi,2*pi))
            }else{
#                kep.max <- c(kep.max,c(pi,pi))#we define 0<Omega<pi for targets without RV data
                kep.max <- c(kep.max,c(pi,2*pi))
            }
            nams.kep <- c(nams.kep,'Inc','Omega')
            if(length(out$timing)>0 & out$Nrv==0){
                kep.ini <- c(kep.ini,1)
                kep.min <- c(kep.min,0)
                kep.max <- c(kep.max,1e6)#au; semi-major axis of photocenter relative to mass center; unit: mas
                nams.kep <- c(nams.kep,'apm')
            }
        }
        if(j>Np){
            nams <- c(nams,paste0('m',nams.kep,j))
        }else{
            nams <- c(nams,paste0(nams.kep,j))
        }
    }
    par.ini <- c(par.ini,as.numeric(kep.ini))
    par.min <- c(par.min,as.numeric(kep.min))
    par.max <- c(par.max,as.numeric(kep.max))
}

####noise models
amax <- (Kmax-Kmin)/((tmax-tmin)/365.25)#m/s/yr
amin <- -amax
bmin <- -max(Kmax,1e6)
bmax <- max(Kmax,1e6)
if(par.global=='trend'){
    par.ini <- c(par.ini,0)
    par.min <- c(par.min,amin)
    par.max <- c(par.max,amax)
    nams <- c(nams,paste0('a',1:Npoly))
}

phi.min <- wmin <- -1
phi.max <- wmax <- 1
smin <- 0
smax <- Kmax
beta.up <- log(tmax-tmin)#time span of the data
if(!is.null(trv.all)){
    beta.low <- log(max(1/24,min(1,min(diff(trv.all)))))#1h or minimum separation
}else{
    beta.low <- 1/24
}
alpha.max <- beta.max <- beta.up#d; limit the range of beta to avoid multimodal or overfitting
alpha.min <- beta.min <- beta.low#24h
ps <- qs <- ns <- c()
if(length(out$tiall) > 0 | length(out$astrometry)>0){
    for(i in out$ins.rv){
###offsets of RV sets
        if(offset){
            par.ini <- c(par.ini,0)
            par.min <- c(par.min,bmin)
            par.max <- c(par.max,bmax)
            nams <- c(nams,paste0('b_',i))
        }
###white jitter
        par.ini <- c(par.ini,median(out[[i]]$RV[,3]))
        par.min <- c(par.min,smin)
        par.max <- c(par.max,smax)
        nams <- c(nams,paste0('J_',i))
###red noise parameters
        nqp <- unlist(out[[i]]$noise$nqp)
        if(nqp[3]>0){
            par.ini <- c(par.ini,rep(0.1,nqp[3]),0)
            par.min <- c(par.min,rep(phi.min,nqp[3]),log(1e-5))
            par.max <- c(par.max,rep(phi.max,nqp[3]),log(1e5))
            nams <- c(nams,paste0('phi',j,'_',i))
            nams <- c(nams,paste0('alpha_',i))
        }
        if(nqp[2]>0){
            par.ini <- c(par.ini,rep(0.1,nqp[2]),0)
            par.min <- c(par.min,rep(wmin,nqp[2]),log(1e-5))
            par.max <- c(par.max,rep(wmax,nqp[2]),log(1e5))
            nams <- c(nams,paste0('w',1:nqp[2],'_',i))
            nams <- c(nams,paste0('beta_',i))
        }
        if(nqp[1]>0){
            par.ini <- c(par.ini,rep(0,nqp[1]))
            par.min <- c(par.min,rep(-Kmax,nqp[1]))
            par.max <- c(par.max,rep(Kmax,nqp[1]))
            nams <- c(nams,paste0('c_',i,1:nqp[1]))
        }
        ps <- c(ps,nqp[3])
        qs <- c(qs,nqp[2])
        ns <- c(ns,nqp[1])
    }
}
###add astrometric parameters
#remove dpmra because it is determined by Mo and other orbital parameters
jitter.max <- 1e5
jitter.min <- 0
logjitter.max <- log(1e5)
logjitter.min <- log(1e-5)
Jmin <- 1#error inflation
Jmax <- 100
#jitter.sd <- 10
jitter.sd <- 0.2
dpos.max <- 1e6
dpm.max <- 1e3
dplx.max <- 100
Nins <- length(out$astro.ins)
if(Nins>0 & out$astro.par>0){
###jitter as a fraction of error
#    if(length(out$ins.rv)>0){
#    if(TRUE){
    if(jitter!='fixed'){
        par.ini <- c(par.ini,rep(1,Nins))
        par.min <- c(par.min,rep(Jmin,Nins))
        par.max <- c(par.max,rep(Jmax,Nins))
        nams <- c(nams,paste0('J_',out$astro.ins))
    }
    if(out$astro.par==5){
        pp <- c(rep(dpos.max,2),dplx.max,rep(dpm.max,2))
        nams <- c(nams,c('dra','ddec','dplx','dpmra','dpmdec'))
    }
    if(out$astro.par==4){
        pp <- c(rep(dpos.max,2),rep(dpm.max,2))
        nams <- c(nams,c('dra','ddec','dpmra','dpmdec'))
    }
    if(out$astro.par==2){
        pp <- rep(dpos.max,2)
        nams <- c(nams,c('dra','ddec'))
    }
    par.ini <- c(par.ini,rep(0,out$astro.par))
    par.min <- c(par.min,-pp)
    par.max <- c(par.max,pp)

#####epoch absolute astrometry data
    if(length(out$data.epoch)>0){
        ind <- which(!grepl('hip|gaia',out$ins.epoch))
        if(length(ind)>0){
            ##dra
            par.ini <- c(par.ini,rep(0,length(ind)))
            par.min <- c(par.min,rep(-1e6,length(ind)))
            par.max <- c(par.max,rep(1e6,length(ind)))
            ##ddec
            par.ini <- c(par.ini,rep(0,length(ind)))
            par.min <- c(par.min,rep(-1e6,length(ind)))
            par.max <- c(par.max,rep(1e6,length(ind)))
            nams <- c(nams,paste0('dra_',out$ins.epoch[ind]))
            nams <- c(nams,paste0('ddec_',out$ins.epoch[ind]))
        }
    }
}
#!any(target==c('HD30219'))
if(any(names(out)=='rel') & jitter!='fixed'){
#        if(length(out$ins.rv)>0){
#        if(TRUE){
    par.ini <- c(par.ini,rep(1,length(rel.name)))
    par.min <- c(par.min,rep(Jmin,length(rel.name)))
    par.max <- c(par.max,rep(Jmax,length(rel.name)))
    nams <- c(nams,paste0('J_',rel.name))
}
###add Mstar
#        if(!any(par.fix=='Mstar')){
if(any(names(out)=='rel') | any(out$comp.epoch==2)){
    par.ini <- c(par.ini,out$Mstar)
    par.min <- c(par.min,out$Mstar/4)
    par.max <- c(par.max,4*out$Mstar)
    nams <- c(nams,'Mstar')
}

if(any(target==c('CPD-632495','TYC3588-11669-1')) & !any(names(out)=='rel')){
    par.ini <- c(par.ini,20)
    par.min <- c(par.min,0.1)
    par.max <- c(par.max,100)
    nams <- c(nams,'Mstar')
}

###companion RV
    if(out$Nrvc>0){
        if(rvc.type=='relative'){
            par.ini <- c(par.ini,rep(0,out$Nrvc))
            par.min <- c(par.min,rep(-1e5,out$Nrvc))
            par.max <- c(par.max,rep(1e5,out$Nrvc))
            nams <- c(nams,paste0('bc',out$Irvc))
        }else if(rvc.type=='barycentric'){
            par.ini <- c(par.ini,-10)#HR8799
            par.min <- c(par.min,-1e3)
            par.max <- c(par.max,1e3)
            nams <- c(nams,'rvb')
        }
    }

if(out$relativity){
    if(gammaf){
        par.ini <- c(par.ini,1)
        par.min <- c(par.min,0)
        par.max <- c(par.max,10)
        nams <- c(nams,'gamma')
    }

    nadd <- 6-out$astro.par
    par.ini <- c(par.ini,rep(0,nadd))
    par.min <- c(par.min,rep(-1e3,nadd))
    par.max <- c(par.max,rep(1e3,nadd))
    nn <- c('dra','ddec','dplx','dpmra','dpmdec','drv')
    ii <- which(is.na(match(nn,nams)))
    if(length(ii)>0) nams <- c(nams,nn[ii])
}

#####reflex motion; jitter
#if(length(out$data.ref)>0){
#    par.ini <- c(par.ini,rep(0,length(out$ins.ref)))
#    par.min <- c(par.min,rep(0,length(out$ins.ref)))
#    par.max <- c(par.max,rep(10,length(out$ins.ref)))
#}

if(photovary){
    par.ini <- c(par.ini,0)
    par.min <- c(par.min,0)
    par.max <- c(par.max,1)
    nams <- c(nams,'eta')
}

####timing offset
if(length(out$timing)>0){
    if(!any(grepl('pulsation|eb',out$timetypes))){
        par.ini <- c(par.ini,0)#min
        par.min <- c(par.min,-1e3)
        par.max <- c(par.max,1e3)
        nams <- c(nams,'t0')
    }else if(grepl('pulsation',out$timetypes)){
        nn <- gsub('pulsation','',out$timetypes)
        for(n in nn){
            par.ini <- c(par.ini,0)#min
            par.min <- c(par.min,-1e3)
            par.max <- c(par.max,1e3)
            nams <- c(nams,paste0('t',n))

            par.ini <- c(par.ini,0)#min/yr
            par.min <- c(par.min,-1e3)
            par.max <- c(par.max,1e3)
            nams <- c(nams,paste0('g',n))

            par.ini <- c(par.ini,0)#min
            par.min <- c(par.min,-1e3)
            par.max <- c(par.max,1e3)
            nams <- c(nams,paste0('h',n))
        }
    }else if(grepl('eb',out$timetypes)){
        nn <- gsub('pulsation','',out$timetypes)
        for(n in nn){
            par.ini <- c(par.ini,0)#min
            par.min <- c(par.min,-1e3)
            par.max <- c(par.max,1e3)
            nams <- c(nams,'t0')

            par.ini <- c(par.ini,0)#min/yr
            par.min <- c(par.min,-1e3)
            par.max <- c(par.max,1e3)
            nams <- c(nams,'g0')

            par.ini <- c(par.ini,0)#min/yr
            par.min <- c(par.min,-1e3)
            par.max <- c(par.max,1e3)
            nams <- c(nams,'h0')
        }
    }
}

#####binary component
if(length(out$data.binary)>0){
    ins.binary <- unique(out$data.binary[,'ins'])
    db <- max(out$data.binary[,2])-min(out$data.binary[,2])
    if(offset){
###offsets of RV sets
        par.ini <- c(par.ini,rep(0,length(ind.binary)))
        par.min <- c(par.min,min(db,rep(bmin,length(ind.binary))))
        par.max <- c(par.max,max(db,rep(bmax,length(ind.binary))))
        nams <- c(nams,paste0('b_',ins.binary))
    }
###white jitter
    par.ini <- c(par.ini,rep(0,length(ins.binary)))
    par.min <- c(par.min,rep(0,length(ins.binary)))
    par.max <- c(par.max,rep(db,length(ins.binary)))
    nams <- c(nams,paste0('s_',out$ins.binary))
}

###if timing but RV data
if(length(out$timing)>0 & out$Nrv==0){
    indk <- grep('^K\\d',nams)
    if(length(indk)>0){
        nams[indk] <- gsub('K','arc',nams[indk])
        par.min[indk] <- 1e-1#mas
        par.max[indk] <- 1e6#mas
    }
}

##assign names
names(par.max) <- names(par.min) <- names(par.ini) <- nams
#ftxt <- paste0('par_',target,'.txt')
#cat(ftxt,'\n')
#write.table(rbind(par.ini,par.min,par.max),ftxt,quote=FALSE,row.names=FALSE)
#startvalue <- assign.names(,Np=Np,p=ps,q=qs,n=ns,bases=bases)
#names(startvalue)
startvalue <- par.ini
Npar <- length(startvalue)
Sd <- 2.4^2/Npar#hyper par of s prior
epsilon <- 1e-6
#if(!is.null(out$timing)) epsilon <- 1e-9
#if(target=='CPD-632495') epsilon <- 1e-10
cov.start <- epsilon*diag(Npar)

