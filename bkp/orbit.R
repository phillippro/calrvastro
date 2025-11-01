library(magicaxis)
library(graphics)
library(RColorBrewer)
source('mcmc_func.R')
options(scipen = 0)
library("plotrix")
source('OrbitFunction.R')
source('timing_function.R')
source('general_function.R')
source('sofa_function.R')
set.seed(9999)
Ndig <- 2
modifypar <- function(par){
    Nsig <- length(grep('^per',names(par)))
    ps <- exp(par[paste0('per',1:Nsig)])
    ind <- which(ps<1000)
    if(length(ind)>0){
        par[paste0('K',ind)] <- 0
    }
    par
}
pmfun <- pemfun <- NULL
if(file.exists('plx_mass.txt')){
    tab <- read.table('plx_mass.txt',header=TRUE)
    pmfun.spec <- approxfun(tab[,'plx'],tab[,'spec_mass'])
    pmfun.evo <- approxfun(tab[,'plx'],tab[,'evo_mass'])
    pemfun.spec <- approxfun(tab[,'plx'],0.5*(tab[,'esmass1']+tab[,'esmass2']))
    pemfun.evo <- approxfun(tab[,'plx'],0.5*(tab[,'eemass1']+tab[,'eemass2']))
}

m22 <- read.table('../data/combined/mamajek22.txt',header=TRUE)
ms <- as.numeric(m22[,'Msun'])
mg <- as.numeric(m22[,'M_G'])
ind <- which(!is.na(ms) & !is.na(mg))
mrl.m22 <- approxfun(ms[ind],mg[ind])
mlow.m22 <- min(ms[ind])
mup.m22 <- max(ms[ind])
eta0 <- rep(NA,10)

Nmc <- 0
#Nmc <- 1000
#Nmc <- 100
star <- 'UCAC4569-026385'
#star <- 'HD209100'
comptype <- 'companion'
version <- 1
set.seed(100)
#ff <- read.table('companion_files.txt')[,1]
#ff <- read.table('2companion_files.txt')[,1]
ff <- read.table('1companion_files.txt')[,1]
ff <- ff[grepl(star,ff)]

#nn <- c('GL676A','HD145675','HD149806','HD62364','HD115404A','HD42581')
#nn <- c('HD42581')
#nn <- c('HD209100','HD22049')
#nn <- c('UCAC4569-026385')
nn <- stars <- gsub('\\/.+','',ff)
ii <- match(nn,stars)
#ii <- as.integer(sapply(nn,function(i) ))
ff <- ff[ii]
stars <- stars[ii]
host <- read.table('host_par.txt',header=TRUE)
rds <- mpacs <- c()

#stars <- gsub('\\/.+','',ff)
#stars <- c('HD42581','HD39060','HD182488','HD215257')
#stars <- c('HD39060','HD182488','HD215257','HD39091')
#stars <- c('HD39060','HD182488')
#stars <- c('HD7449')
Nstar <- length(stars)
fpdf <- paste0('prediction_orbit_Nmc',Nmc,'.pdf')
#cat('output:\n',fpdf,'\n')
etaH <- etaG <- 1
###modifiy input file names
#Nsig <- length(fs)
planets <- read.csv('../data/code/SigClassify/ranking_complex2.csv')
n0 <- gsub(' ','',as.character(planets[,'Name']))
n1 <- gsub(' ','',as.character(planets[,'ID']))
n2 <- gsub(' ','',as.character(planets[,'StarKnown']))
nfig <- 0
npdf <- 0

####pdf figure
#f1 <- gsub('mas',paste0('mas_p',npdf),fpdf)
f1 <- fpdf
cat(f1,'\n')
pdf(f1,16,16)
size <- 1.2
if(Nmc>0){
    size <- 1.2
    par(mfrow=c(4,4),mar=c(4,4,2,1),cex=size,cex.lab=size,cex.axis=1,oma=c(0,0,1,0),mgp=c(2,1,0))
    ##par(mfrow=c(4,4),cex=size,cex.lab=size,cex.axis=size,mar=c(4,4,2,1),mgp=c(2,0.5,0),cex.main=1.2*size,oma=c(0,0,1,0),cex.axis=0.8)
}else{
    par(mfrow=c(4,4),cex=size,cex.lab=size,cex.axis=size,mar=c(4,4,2,1),mgp=c(2,0.5,0),cex.main=1.2*size)
}

###collect all data needed for the plot
#out.orbit <- c()
for(j in 1:Nstar){
    star <- stars[j]
    star.name <- gsub('HD','HD ',star)
    star.name <- gsub('HIP','HIP ',star.name)
    star.name <- gsub('GJ|GL','GJ ',star.name)
###check whether to show
    ind <- which(host[,1]==gsub('GL','GJ',star))
    plx <- host[ind,'plxG']
#    out.orbit <- rbind(out.orbit,c(pp[pp[,1]==star,],host[ind,]))
    cat('plotting!\n')
    ii <- grep(paste0('^',star,'\\/'),ff)
    fobj <- gsub('pdf','Robj',ff[ii])
    fs <- paste0('results/',fobj)
    cat('load ',fs,'\n')
    ind <- which(stars==star)
    load(fs)
    out0 <- out
#    if(star=='UCAC4569-026385'){
#        out$Mstar <- 2.07
#        out$eMstar <- (0.19+0.26)/2
#    }
    source('mcmc_func.R')
    popts <- Popt
    Nkep <- length(Popt)
    mc <- out$mcmc.opt[[paste0('sig',Nkep)]]
    plxs <- out$astrometry[out$iref,'parallax']-mc[,'dplx']
    Mstars <- pmfun.spec(plxs)
    eMstars <- pemfun.spec(plxs)
    out$Mstar <- mean(Mstars)
    out$eMstar <- sd(Mstars)
    if(Nmc==0){
        nepoch <- 2
                                        #        astro.type <- 'astro'
        Nastro <- 2
                                        #ell.col <- rep('darkgrey',2)
                                        #ell.col <- c('darkblue','darkgreen')
        ell.col <- c('darkblue','darkgreen')
                                        #ell.pch <- c(16,17)
                                        #ell.pch <- c(20,20)
        ell.pch <- 18
        rel.pch <- 15
        cols <- tcol(c('black','blue','green','purple','cyan','brown','orange','steelblue','yellow',brewer.pal(8,'Accent')),50)
                                        #cc <- c('black','blue','orange','steelblue','green','darkgreen',rainbow(10))
                                        #        cc <- brewer.pal(8,'Accent')
        size <- 0.9
        show.center <- FALSE
        margin <- 1.3
        pch.size <- 1
        if(Nmc>0){
            line.size <- 2
        }else{
                                        #    line.size <- 3
            line.size <- 1
        }
        fit.col <- tcol('red',90)
                                        #fit.opt <- tcol('red',50)
                                        #fit.opt <- c('red','blue','green')
        fit.opt <- tcol(c('red','red'),20)
                                        #layout(matrix(c(1, 1,2,2), nrow=2, byrow=TRUE))
        Nstar <- length(stars)
                                        #Nstar <- Nsig
####RV
        imaging <- FALSE
        if(grepl('hg3\\+|astro\\+',fs) & !grepl('HD39091',fs)) imaging <- TRUE
                                        #        if(grepl('hg3\\+|astro\\+',fs)) imaging <- TRUE
        ins <- out$ins.rv
        Nset <- length(out$ins.rv)
        ind <- sort(trv.all,index.return=TRUE)$ix
        tmin <- min(trv.all)
###split the RV plot into two sub-panel
        instr <- ins
###simulation for ylim
        jitter <- TRUE
        parallax <- FALSE
        if(!exists('bases')) bases <- rep('natural',5)
        astrometry <- 8
        kep.type <- 'pure'
        yr2d <- 365.25
        period.par <- 'logP'
        offset <- TRUE
        time.unit <- 1
        Np <- length(ind)
        prior.type <- 'mt'
        pars <- par.opt
        ind <- grep('per',names(pars))
        Popt <- exp(par.opt[ind])
        indI <- grep('Inc',names(pars))
        par.opt[indI] <- par.opt[indI]%%pi
        par.opt1 <- par.opt
###
        ind.show <- which.max(par.opt[grep('per',names(par.opt))])
        yr2d <- 365.25
        ##astrometry fit
        data.astrometry <- out$astrometry
        if(!any(names(out)=='shortP')){
            shortP <- FALSE
        }else{
            shortP <- out$shortP
        }
        if(!any(names(out)=='plx')) out$plx <- out$astrometry$parallax[2]
        astro0 <- astrometry.kepler(par.opt1,bases=bases,Pmin=1000)$barycenter
        Dt <- max(popts)/2
        tsim0 <- seq(data.astrometry[1,1]-Dt,data.astrometry[2,1]+Dt,length.out=1e4)
        tsim <- tsim0-tmin
        Nt <- length(tsim)

        if(!any(names(out)=='rel')){
            out$rel <- c()
            out$Nrel <- rep(10,10)
        }
        planet <- astrometry.rel(par.opt,tt=tsim0,bases=bases)
        Nsamp <- 100

        Popt <- exp(par.opt[paste0('per',ind.show)])
        Kopt <- par.opt[paste0('K',ind.show)]
        eopt <- par.opt[paste0('e',ind.show)]
        Iopt <- par.opt[paste0('Inc',ind.show)]%%(pi)
        Mopt <- par.opt[paste0('Mo',ind.show)]
        Oopt <- par.opt[paste0('Omega',ind.show)]
        oopt <- par.opt[paste0('omega',ind.show)]
##        msini <- K2msini.full(Kopt,Popt,eopt,Ms=Mstars)
###revised function
#        tmp <- k2m(Kopt,Popt,eopt,Ms=out$Mstar,Inc=Iopt,more=TRUE)
#        mp.opt <- Mp <- tmp$mj
        Mo <- (Mopt%%(2*pi))*180/pi#deg
        inc <- Iopt*180/pi#deg
        ##        cc <- c('black','darkgrey','grey')
        cc <- c('black','red','blue')

####plotting

#        j2 <- which(popts>1000)
        j2 <- which(popts>1)
        jj <- j2[sort(popts[j2],index.return=TRUE,decreasing=TRUE)$ix]
        for(j3 in 1:length(jj)){
            j1 <- jj[j3]
            n <- paste0('p',j1)
            popt <- exp(par.opt[paste0('per',j1)])
            eopt <- par.opt[paste0('e',j1)]
            show.epoch <- FALSE
            if(popt==max(popts)){
                show.epoch <- TRUE
                yrs <- time_Jd2yr(cbind(tsim0,0))
                dt <- max(popts)/365.25
                ##                tts <- c(2020,2020+min(dt/2,max(5,round(dt/10))),2020+min(dt/2,max(10,round(dt/10))))
                tts <- 2020+round(seq(0,dt/2,length.out=4))
                if(dt<5){
                    tts <- seq(2022,2025,by=1)
                }else if(dt<10){
                    tts <- seq(2022,2028,by=2)
                }else if(dt<15){
                    tts <- seq(2022,2031,by=3)
                }else if(star=='HD39060'){
                    tts <- c(2020,2024)
                }
                ii <- unlist(lapply(tts, function(tt) which.min(abs(yrs-tt))))
            }
            n <- paste0('p',j1)
            if(FALSE){
                draS <- planet[[n]][,'ra']
                ddecS <-planet[[n]][,'dec']
                xlim <- range(unlist(sapply(names(planet), function(i) planet[[i]][,'ra'])))
                ylim <- range(unlist(sapply(names(planet), function(i) planet[[i]][,'dec'])))
            }else{
                draS <- planet$rel$dra[[n]]
                ddecS <-planet$rel$ddec[[n]]
                xlim <- range(unlist(planet$rel$dra))
                ylim <- range(unlist(planet$rel$ddec))
            }
            nc <- length(popts)
            if(j3==1){
                plot(draS,ddecS,xlab=expression(Delta*alpha[c]*"*[mas]"),ylab=expression(Delta*delta[c]*"[mas]"),xlim=rev(xlim),ylim=ylim,type='l',col=cc[j3],main=star.name)
                points(0,0,pch='+',col='black')
            }else{
                lines(draS,ddecS,col=cc[j3])
            }
            if(show.epoch){
                ##                points(draS[ii],ddecS[ii],pch='+',col=tcol('red',50))
                pos <- rep(4,length(ii))
                pos[draS[ii]<0] <- 2
                ##                jj <- which(sqrt(diff(draS[ii]/(max(draS)-min(draS)))^2+diff(ddecS[ii]/(max(ddecS)-min(ddecS)))^2)<1/20)
                j3 <- c()
                if(any(star==c('HD71881','GJ494','GJ864','HD104263','HD154682','HD156279','HD16905','HD181234','HD188641','HD221420','HD24633','HD25015','HD56957','HD66428','HD71881','HD80883','HD80913','HD85472','HIP116838','HD14729')) | eopt>0.3) j3 <- 4
                if(any(star==c('HIP26196','HIP6344','GJ234','HD167665','HD65430','HD87883')) | eopt>0.5) j3 <- 3:4
                if(dt<5 & eopt>0.3) j3 <- 2:4
                if(length(j3)==0){
                    text(draS[ii],ddecS[ii],labels=tts,pos=pos,xpd=NA,col='darkgrey',cex=0.8)
                    points(draS[ii],ddecS[ii],pch='+',col='darkgrey')
                }else{
                    text(draS[ii[-j3]],ddecS[ii[-j3]],labels=tts[-j3],pos=pos,xpd=NA,col='darkgrey',cex=0.8)
                    points(draS[ii[-j3]],ddecS[ii[-j3]],pch='+',col='darkgrey')
                }
            }
        }
    }else if(Nmc>0){
        source('orbit_predict.R')
    }
}
if(Nmc>0) title(main=paste('Reference epoch: ',paste(as.numeric(cal),collapse='-')),line=0,cex=0.8,xpd=NA,outer=TRUE)
dev.off()
fout <- paste0('targets_par_with_predicted_orbits.txt')
cat(fout,'\n')
##write.table(out.orbit,file=fout,quote=FALSE,row.names=FALSE)

if(Nmc>0){
#    fout <- paste0('companion_mpa.txt')
    fout <- paste0(star,'_mpa.txt')
    cat(fout,'\n')
    mpacs <- cbind(mpacs,rds)
    colnames(mpacs) <- c('star','planet','p','e_p','E_p','K','e_K','E_K','e','e_e','E_e','I','e_I','E_I','omega','e_omega','E_omega','Omega','e_Omega','E_Omega','M0','e_M0','E_M0','dra','e_dra','E_dra','ddec','e_ddec','E_ddec','dplx','e_dplx','E_dplx','dpmra','e_dpmra','E_dpmra','dpmdec','e_dpmdec','E_dpmdec','pyr','e_pyr','E_pyr','a','e_a','E_a','m','e_m','E_m','Tp','e_Tp','E_Tp','mstar','e_mstar','E_mstar','plx','e_plx','E_plx','rho','e_rho','E_rho','theta','e_theta','E_theta')
    write.table(mpacs,file=fout,quote=FALSE,row.names=FALSE)
    dat <- read.table(fout,header=TRUE)
    source('produce_latex_table.R')
}
