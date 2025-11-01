source('mcmc_func.R')
library(plotrix)
ds  <- gsub('results/','',list.dirs('results'))
stars <- ds[ds!='results']
#stars  <- gsub('\\./','',ds)
ms <- c()
ems <- c()
nc <- 204
pars <- t(rep(NA,nc))#40*5+3+1
ff  <- c()
ss <- c()
fe <- read.table('except_files.txt')[,1]
se <- gsub('\\/.+','',fe)
for(k in 1:length(stars)){
    star <- stars[k]
    fs <- list.files(paste0('results/',star),pattern='hg123.+Robj',full.name=TRUE)
    if(length(fs)>0){
        lls <- gsub('.+lnlmax|\\.Robj','',fs)
        ind <- which.max(as.numeric(lls))
        kk <- which(star==se)
        if(length(kk)>0){
            ind <- which(gsub('.+\\/|pdf|Robj','',fs)==gsub('.+\\/|pdf|Robj','',fe[kk]))
            cat('new file used:',fs[ind],';ind=',ind,'\n')
        }
        per <- gsub('.+_P|_acc.+','',fs[ind])
        if(grepl('d',per) & !grepl('pimen',fs[ind])){
#            f <- gsub('Robj','pdf',fs[ind])
            f <- fs[ind]
            ff <- c(ff,f)
            ss <- c(ss,star)
        }
    }
}
Ncores <- 16
#Ncores <- 1
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
ind <- which(ss=='HD25311')
#for(k in 1:length(ss)){
#for(k in ind+(-1:1)){
pars <- foreach(k = 1:length(ss), .combine='rbind',.errorhandling="pass") %dopar% {
#pars <- foreach(k = 1:10, .combine='rbind',.errorhandling="pass") %dopar% {
    star <- ss[k]
    f5 <- ff[k]
    cat('k=',k,';',f5,'\n')
    load(f5,env= e0 <- new.env())
    out <- e0$out
    Nsig <- e0$Nsig
    par.opt <- e0$par.opt
    par.stat <- out$par.stat[[paste0('sig',Nsig)]]
    cn <- colnames(par.stat)
    indi <- grep('Inc',cn)
    indp <- grep('per',cn)
    indO <- grep('Omega',cn)
    inde <- grep('^e',cn)
    indk <- grep('K',cn)
    indo <- grep('omega',cn)
    indm <- grep('Mo',cn)
#    ss <- grep('Omega',cn)
    mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
    npar <- ncol(mc)-2
    par.stat <- c()
    par.name <- colnames(mc)[1:npar]
    for(j in 1:npar){
        pp <- mc[,j]
        if(grepl('per',par.name[j])){
            pp <- exp(pp)
        }
        if(grepl('omega|Omega|Mo',par.name[j])){
            pp <- convert.angle(pp)
        }
        par.stat <- cbind(par.stat,data.distr(as.numeric(pp),lp=mc[,'logpost'],plotf=FALSE))
    }
    colnames(par.stat) <- par.name
    if(!any(names(par.opt)=='Mstar')){
        Ms <- rnorm(nrow(mc)*2,out$Mstar,out$eMstar)
        Ms <- Ms[Ms>0][1:nrow(mc)]
        Mstar <- out$Mstar
        eMstar <- out$eMstar
    }else{
        Ms <- mc[,'Mstar']
        Mstar <- mean(Ms)
        eMstar <- sd(Ms)
    }
    tmp <- c(star,Mstar,eMstar)
    for(i in 1:length(indp)){
        inc.med <- par.stat['med',indi[i]]
        inc.map <- par.stat['xopt',indi[i]]
        inc.minus <- par.stat['xlow',indi[i]]
        inc.plus <- par.stat['xup',indi[i]]

        p.med <- par.stat['med',indp[i]]
        p.map <- par.stat['xopt',indp[i]]
        p.minus <- par.stat['xlow',indp[i]]
        p.plus <- par.stat['xup',indp[i]]

        e.med <- par.stat['med',inde[i]]
        e.map <- par.stat['xopt',inde[i]]
        e.minus <- par.stat['xlow',inde[i]]
        e.plus <- par.stat['xup',inde[i]]

        omega.med <- par.stat['med',indo[i]]
        omega.map <- par.stat['xopt',indo[i]]
        omega.minus <- par.stat['xlow',indo[i]]
        omega.plus <- par.stat['xup',indo[i]]

        Omega.med <- par.stat['med',indO[i]]
        Omega.map <- par.stat['xopt',indO[i]]
        Omega.minus <- par.stat['xlow',indO[i]]
        Omega.plus <- par.stat['xup',indO[i]]

        K.med <- par.stat['med',indk[i]]
        K.map <- par.stat['xopt',indk[i]]
        K.minus <- par.stat['xlow',indk[i]]
        K.plus <- par.stat['xup',indk[i]]

        Mo.med <- par.stat['med',indm[i]]
        Mo.map <- par.stat['xopt',indm[i]]
        Mo.minus <- par.stat['xlow',indm[i]]
        Mo.plus <- par.stat['xup',indm[i]]

        ps <- exp(mc[,indp[i]])
        Ks <- mc[,indk[i]]
        es <- mc[,inde[i]]
        incs <- mc[,indi[i]]
        Omegas <- mc[,indO[i]]

        mps <- k2m(Ks,ps,es,Ms,Inc=incs)$mj
        mm <- data.distr(mps,lp=mc[,'logpost'],plotf=FALSE)
        mp.med <- mm['med']
        mp.minus <- mm['xlow']
        mp.plus <- mm['xup']
        mp.map <- mm['xopt']

        as <- ((mps/1047+Ms)*(ps/365.25)^2)^(1/3)#au
        aa <- data.distr(as,lp=mc[,'logpost'],plotf=FALSE)
#as.numeric(quantile(as,c(0.16,0.5,0.84)))
        a.med <- aa['med']
        a.minus <- aa['xlow']
        a.plus <- aa['xup']
        a.map <- aa['xopt']

        inc1 <- incs
        Omega1 <- Omegas
        psis <- c()
        for(j in 1:length(indi)){
            if(j!=i){
                inc2 <- mc[,indi[j]]
                Omega2 <- mc[,indO[j]]
                psi <- acos(cos(inc1)*cos(inc2)+cos(Omega2-Omega1)*sin(inc1)*sin(inc2))*180/pi#mutual inclination in unit of degrees
                psi.out <- data.distr(as.numeric(psi),lp=mc[,'logpost'],plotf=FALSE)
                psis <- rbind(psis,psi.out[c('xlow','med','xup','xopt')])
            }
        }
        psi.med <- paste(psis[,2],collapse='d')
        psi.minus <- paste(psis[,1],collapse='d')
        psi.plus <- paste(psis[,3],collapse='d')
        psi.map <- paste(psis[,4],collapse='d')
        tmp0 <- c(p.med,p.minus,p.plus,p.map,inc.med,inc.minus,inc.plus,inc.map,Omega.med,Omega.minus,Omega.plus,Omega.map,psi.med,psi.minus,psi.plus,psi.map,e.med,e.minus,e.plus,e.map,K.med,K.minus,K.plus,K.map,mp.med,mp.minus,mp.plus,mp.map,omega.med,omega.minus,omega.plus,omega.map,Mo.med,Mo.minus,Mo.plus,Mo.map,a.med,a.minus,a.plus,a.map)
        tmp <- c(tmp,tmp0)
    }
    if(length(tmp)<(nc-1)) tmp  <- c(tmp,rep(NA,nc-1-length(tmp)))
    cat('star:',star,';file:',f5,'\n')
#    pars <- data.frame(rbind(pars,c(tmp,f5)))
    t(c(tmp,f5))
}

colnames(pars) <- c('name','Mstar','eMstar',as.vector(outer(c('Per.med','Per.lower','Per.upper','Per.map','Inc.med','Inc.lower','Inc.upper','Inc.map','Omega.med','Omega.lower','Omega.upper','Omega.map','psi.med','psi.lower','psi.upper','psi.map','e.med','e.lower','e.upper','e.map','K.med','K.lower','K.upper','K.map','mp.med','mp.lower','mp.upper','mp.map','omega.med','omega.lower','omega.upper','omega.map','Mo.med','Mo.lower','Mo.upper','Mo.map','a.med','a.lower','a.upper','a.map'),1:5,paste0)),'file')

#####units
#Mstar: Msun
#per: day
#Inc: rad
#Omega: rad
#omega: rad
#mo: rad
#a: au
#psi: deg

fout <- paste0('results/misaligned_starpar_N',length(ss),'.txt')
cat(fout,'\n')
if(all(is.na(pars[1,]))) pars <- pars[-1,]
write.csv(pars,file=fout,quote=FALSE,row.names=FALSE)

if(FALSE){
for(j in 1:nrow(pars)){
####plot mutual misalignment
    fpdf <- paste0(ss[j],'_polar.pdf')
    cat(fpdf,'\n')
    pdf(fpdf,6,6)
    indI <- grep('Inc.map',colnames(pars))
    indO <- grep('Omega.map',colnames(pars))
    Iopt <- as.numeric(pars[j,indI])*180/pi
    Oopt <- as.numeric(pars[j,indO])*180/pi
    ii <- which(!is.na(Iopt)&!is.na(Oopt))
    polar.plot(Iopt[ii],Oopt[ii],main=stars[k],lwd=3,line.col=4,rp.type='s',point.symbols=20,point.col='steelblue',radial.lim=c(0,max(Iopt[ii])))
    dev.off()
}
}
