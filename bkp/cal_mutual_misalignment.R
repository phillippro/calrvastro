source('mcmc_func.R')
library(plotrix)
ms <- c()
ems <- c()
pars <- t(rep(NA,203))#40*5+3
#fs <- read.table('companion_files.txt')[,1]
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    fs <- args[1]
}else{
    fs <- 'HR8799/HR8799_noPAgost_photovaryFALSE_relativityFALSE_Niter3100000_Ncores16_ofac2_Nset0_hg123+_transit0_P157582d148183d40880d22305_acc0.015_sinI_lnlmax-1585.pdf'
#    fs <- 'UCAC4569-026385/UCAC4569-026385_noPAgost_photovaryFALSE_relativityFALSE_Niter2200000_Ncores16_ofac2_Nset1_hg123_transit0_P886_acc0.31_sinI_lnlmax-196.pdf'
}
#stars <- c('HD39091')
#stars <- c('HR8799')
stars <- gsub('\\/.+','',fs)
ss <- gsub('\\/.+','',fs)
ind <- match(stars,ss)
fs <- fs[ind]
fs <- paste0('results/',gsub('.pdf','.Robj',fs))
for(k in 1:length(stars)){
star <- stars[k]
                                        #    fs <- list.files(paste0('results/',stars[k]),pattern='hg123.+Robj',full.name=TRUE)
                                        #    if(length(fs)>0){
                                        #    lls <- gsub('.+lnlmax|\\.Robj','',fs)
                                        #    ind <- which.max(as.numeric(lls))
                                        #    per <- gsub('.+_P|_acc.+','',fs[ind])
                                        #    if(grepl('d',per)){
#    tmp <- gsub('\\/.+','',fs[ind])
#    cat(tmp,';',k,'/',length(stars),'\n')
#    load(fs[ind],env=e0 <- new.env())
#    load(fs[k],env=e0 <- new.env())
    load(fs[k])
    par.stat <- out$par.stat[[paste0('sig',Nsig)]]
    cn <- colnames(par.stat)
    ii <- grep('Inc',cn)
    jj <- grep('per',cn)
    kk <- grep('Omega',cn)
    ll <- grep('^e',cn)
    nn <- grep('K',cn)
    oo <- grep('omega',cn)
    rr <- grep('Mo',cn)
mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
#    ss <- grep('Omega',cn)
    if(!any(names(par.opt)=='Mstar')){
        Mstar <- out$Mstar
        eMstar <- out$eMstar
        if(star=='UCAC4569-026385'){
            Mstar <- 1.7
            eMstar <- 0.15
        }
        Ms <- rnorm(nrow(mc)*2,Mstar,eMstar)
        Ms <- Ms[Ms>0][1:nrow(mc)]
    }else{
        Ms <- mc[,'Mstar']
        Mstar <- mean(Ms)
        eMstar <- sd(Ms)
    }
    tmp <- c(stars[k],Mstar,eMstar)
Nsig <- length(ii)
    for(i in 1:Nsig){
        inc.opt <- par.stat['med',ii[i]]
        inc.map <- par.stat['xopt',ii[i]]
        inc.minus <- par.stat['xminus.1sig',ii[i]]
        inc.plus <- par.stat['xplus.1sig',ii[i]]
        e.opt <- par.stat['med',ll[i]]
        e.map <- par.stat['xopt',ll[i]]
        e.minus <- par.stat['xminus.1sig',ll[i]]
        e.plus <- par.stat['xplus.1sig',ll[i]]
        K.opt <- par.stat['med',nn[i]]
        K.map <- par.stat['xopt',nn[i]]
        K.minus <- par.stat['xminus.1sig',nn[i]]
        K.plus <- par.stat['xplus.1sig',nn[i]]
        Ks <- out$mcmc.opt[[paste0('sig',Nsig)]][,nn[i]]
        es <- out$mcmc.opt[[paste0('sig',Nsig)]][,ll[i]]
        incs <- out$mcmc.opt[[paste0('sig',Nsig)]][,ii[i]]
        ps <- exp(out$mcmc.opt[[paste0('sig',Nsig)]][,jj[i]])#day
        p.map <- exp(par.stat['xopt',jj[i]])
        omegas <- convert.angle(out$mcmc.opt[[paste0('sig',Nsig)]][,oo[i]])
        omega.map <- par.stat['xopt',oo[i]]
        Mos <- convert.angle(out$mcmc.opt[[paste0('sig',Nsig)]][,rr[i]])
        Mo.map <- par.stat['xopt',rr[i]]
        Omegas <- convert.angle(out$mcmc.opt[[paste0('sig',Nsig)]][,kk[i]])
        Omega.map <- par.stat['xopt',kk[i]]
        tmp3 <- as.numeric(quantile(ps,c(0.16,0.5,0.84)))
        p.opt <- tmp3[1]
        p.minus <- tmp3[2]
        p.plus <- tmp3[3]
        tmp1 <- as.numeric(quantile(omegas,c(0.16,0.5,0.84)))
        omega.opt <- tmp1[1]
        omega.minus <- tmp1[2]
        omega.plus <- tmp1[3]
        tmp2 <- as.numeric(quantile(Mos,c(0.16,0.5,0.84)))
        Mo.opt <- tmp2[1]
        Mo.minus <- tmp2[2]
        Mo.plus <- tmp2[3]
        tmp4 <- as.numeric(quantile(Omegas,c(0.16,0.5,0.84)))
        Omega.opt <- tmp4[1]
        Omega.minus <- tmp4[2]
        Omega.plus <- tmp4[3]

        msini <- K2msini.full(Ks,ps,es,Ms)
        mps <- msini$mj/sin(incs)#Mjup
        as <- ((mps/1047+Ms)*(ps/365.25)^2)^(1/3)#au
        aa <- as.numeric(quantile(as,c(0.16,0.5,0.84)))
        a.opt <- aa[2]
        a.minus <- aa[1]
        a.plus <- aa[3]
        mm <- as.numeric(quantile(mps,c(0.16,0.5,0.84)))
        mp.opt <- mm[2]
        mp.minus <- mm[1]
        mp.plus <- mm[3]
        msini.map <- K2msini.full(K.opt,p.opt,e.opt,Mstar)$mj
        mp.map <- msini.map/sin(inc.map)
        a.map <- ((mp.map/1047+Mstar)*(p.map/365.25)^2)^(1/3)#au
        inc1 <- incs
        Omega1 <- Omegas
        if(Nsig>1){
            psis <- c()
            for(j in 1:length(ii)){
                if(j!=i){
                    inc2 <- out$mcmc.opt[[paste0('sig',Nsig)]][,ii[j]]
                    inc2.map <- out$mcmc.opt[[paste0('sig',Nsig)]][which.max(out$mcmc.opt[[paste0('sig',Nsig)]][,'logpost']),ii[j]]
                    Omega2 <- out$mcmc.opt[[paste0('sig',Nsig)]][,kk[j]]
                    Omega2.map <- out$mcmc.opt[[paste0('sig',Nsig)]][which.max(out$mcmc.opt[[paste0('sig',Nsig)]][,'logpost']),kk[j]]
                    psi <- acos(cos(inc1)*cos(inc2)+cos(Omega2-Omega1)*sin(inc1)*sin(inc2))*180/pi#mutual inclination
                    psi.map <- acos(cos(inc.map)*cos(inc2.map)+cos(Omega2.map-Omega.map)*sin(inc.map)*sin(inc2.map))*180/pi#mutual inclination
                    psis <- rbind(psis,c(quantile(psi,c(0.16,0.5,0.84)),psi.map))
                }
            }
            psi.opt <- paste(psis[,2],collapse='d')
            psi.minus <- paste(psis[,1],collapse='d')
            psi.plus <- paste(psis[,3],collapse='d')
            psi.map <- paste(psis[,4],collapse='d')
        }else{
            psi.opt <- psi.minus <- psi.plus <- psi.map <- NA
        }
        tmp0 <- c(p.opt,p.minus,p.plus,p.map,inc.opt,inc.minus,inc.plus,inc.map,Omega.opt,Omega.minus,Omega.plus,Omega.map,psi.opt,psi.minus,psi.plus,psi.map,e.opt,e.minus,e.plus,e.map,K.opt,K.minus,K.plus,K.map,mp.opt,mp.minus,mp.plus,mp.map,omega.opt,omega.minus,omega.plus,omega.map,Mo.opt,Mo.minus,Mo.plus,Mo.map,a.opt,a.minus,a.plus,a.map)
                                        #	cat('length(tmp0)=',length(tmp0),'\n')
        tmp <- c(tmp,tmp0)
                                        #	cat('length(tmp)=',length(tmp),'\n')
    }
    if(length(tmp)<ncol(pars)) tmp  <- c(tmp,rep(NA,ncol(pars)-length(tmp)))
                                        #	pars <- rbind(pars,tmp)
    pars <- data.frame(rbind(pars,tmp))

}

colnames(pars) <- c('name','Mstar','eMstar',as.vector(outer(c('Per.opt','Per.lower','Per.upper','Per.map','Inc.opt','Inc.lower','Inc.upper','Inc.map','Omega.opt','Omega.lower','Omega.upper','Omega.map','psi.opt','psi.lower','psi.upper','psi.map','e.opt','e.lower','e.upper','e.map','K.opt','K.lower','K.upper','K.map','mp.opt','mp.lower','mp.upper','mp.map','omega.opt','omega.lower','omega.upper','omega.map','Mo.opt','Mo.lower','Mo.upper','Mo.map','a.opt','a.lower','a.upper','a.map'),1:5,paste0)))
#####units
#Mstar: Msun
#per: day
#Inc: rad
#Omega: rad
#psi: deg
#omega: rad
#mo: rad
#a: au
#fout <- paste0('results/misaligned_starpar_N',length(stars),'.txt')
fout <- gsub('Robj','txt',fs[1])
cat(fout,'\n')
pars <- pars[-1,]
write.table(pars,file=fout,quote=FALSE,row.names=FALSE)

for(j in 1:nrow(pars)){
####plot mutual misalignment
#    fpdf <- paste0(stars[j],'_polar.pdf')
    fpdf <- gsub('.Robj','_polar.pdf',fs[1])
    cat(fpdf,'\n')
    pdf(fpdf,6,6)
    indI <- grep('Inc.opt',colnames(pars))
    indO <- grep('Omega.opt',colnames(pars))
    Iopt <- as.numeric(pars[j,indI])*180/pi
    Oopt <- as.numeric(pars[j,indO])*180/pi
    ii <- which(!is.na(Iopt)&!is.na(Oopt))
    polar.plot(Iopt[ii],Oopt[ii],main=stars[k],lwd=3,line.col=4,rp.type='s',point.symbols=20,point.col='steelblue',radial.lim=c(0,max(Iopt[ii])))
    dev.off()
}
