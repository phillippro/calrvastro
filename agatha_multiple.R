options(warn=2)
library(magicaxis)
source('periodograms.R')
source('periodoframe.R')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    ind.min <- as.integer(args[1])
    ind.max <- as.integer(args[2])
}else{
    ind.min <- 1
    ind.max <- 1
}

###global parameters
fmin <- 1/1000
fmax <- 1/1.01
Nmax <- 2
ofac <- 1
#noise.types0 <- c('white','MA','AR','GP011')
noise.types0 <- c('white','MA','AR')
#noise.types0 <- c('white','MA')
#noise.types0 <- c('white')
#noise.types0 <- c('GP011')
pq <- 1
Nsamp <- 1
#refine <- TRUE
refine <- FALSE
sampling <- 'combined'
ind.proxy <- 0#
noise.only <- FALSE

###data files
dir <- '/Users/phillippro/Documents/projects/dwarfs/data/pfs/pfs_vels_unbinned'
flist <- list.files(path=dir)[ind.min:ind.max]
targets <- gsub('_.+','',flist)
ins <- gsub('.+_|\\.vels|\\.dat','',flist)
solutions <- list()

for(k3 in 1:length(targets)){
    target <- targets[k3]
##load data
    tab <- read.table(paste0(dir,'/',flist[k3]))#jd,rv,erv,sindex,Halpha,photon count,observation time
##remove outliers
    ind.ref <- 2:3
    if(ncol(tab)<4){
        Sindex <- Halpha <- FALSE
    }else{
        if(all(tab[,4]==-1)){
            Sindex <- FALSE
        }else{
            Sindex <- TRUE
        }
        if(all(tab[,5]==-1)){
            Halpha <- FALSE
        }else{
            Halpha <- TRUE
        }
    }

    sigmas <- sapply(ind.ref,function(j) sd(tab[,j]))
    means <- sapply(ind.ref,function(j) mean(tab[,j]))
    medians <- sapply(ind.ref,function(j) median(tab[,j]))
    ind <- which(sapply(1:nrow(tab),function(i) (any(abs(tab[i,ind.ref]-medians) > 3*sigmas))))
    if(length(ind)>0){
        tab <- tab[-ind,]
    }
    trv <- tab[,1]
    ind <- sort(trv,index.return=TRUE)$ix
    trv <- trv[ind]
    tt <- trv-min(trv)
    tsim <- seq(min(tt),max(tt),length.out=10001)
    tab <- tab0 <- tab[ind,]
    tab[,1] <- trv-min(trv)
    if(Sindex){
        iss <- which(abs(tab[,4]-median(tab[,4]))<3*sd(tab[,4]))
        ss <- as.numeric(scale(tab[iss,4]))#scale Halpha value 
        dss <- sd(ss)/sd(tab[iss,2])*tab[iss,3]
    }
    if(Halpha){
        ihh <- which(abs(tab[,5]-median(tab[,5]))<3*sd(tab[,5]))
        hh <- as.numeric(scale(tab[ihh,5]))#scale Halpha value 
        dhh <- sd(hh)/sd(tab[ihh,2])*tab[ihh,3]
    }
    if(!Sindex & !Halpha & any(grepl('GP',noise.types0))){
        noise.types <- noise.types0[-grep('GP',noise.types0)]
    }else{
        noise.types <- noise.types0
    }
    Ntype <- length(noise.types)
    DT <- max(trv)-min(trv)

###plot setting up
    if(!file.exists(paste0('results/',target))){
        system(paste0('mkdir results/',target))
    }
    fname <- paste0('results/',target,'/',target,'_BFP_',ins[k3],'_ofac',ofac)
    fpdf <- paste0(fname,'.pdf')
    fobj <- paste0(fname,'.Robj')
    cat('output pdf:\n',fpdf,'\n')
    pdf(fpdf,16,16)
    par(mfrow=c(4,4))
    plot(tab0[,1],tab0[,2],xlab='JD',ylab='RV [m/s]')
    arrows(tab0[,1],tab0[,2]+tab0[,3],tab0[,1],tab0[,2]-tab0[,3],length=0.05,angle=90,code=3,col='grey')
    if(Sindex){
        plot(tab0[iss,1],ss,xlab='JD',ylab='Sindex')
    }
    if(Halpha){
        plot(tab0[ihh,1],hh,xlab='JD',ylab='Halpha')
    }
    if(ncol(tab)>5){
        plot(tab0[,1],tab0[,6],xlab='JD',ylab='Photon Count')
    }

    Nmas <- Nars <- durs <- c()
    cnames <- c()
    out <- c()
    opt.par <- ws <- sigs <- list()
    Popt <- logBFs <- c()
    for(j1 in 1:length(noise.types)){
        noise.type <- noise.types[j1]
        cat('noise.type=',noise.type,'\n')
        start <- c(sj=0,sigmaGP=0,logProt=1,logtauGP=1)
        ll <- pp <- c()
        if(grepl('GP',noise.type)){
###fit the GP kernel to Halpha
###Determine error of Halpha by scaling the RV error
            if(Halpha){
                ind.act <- ihh
                activity <- hh
                Dactivity <- dhh
            }else if(Sindex){
                ind.act <- iss
                activity <- ss
                Dactivity <- dss
            }
            data <- cbind(tab[ind.act,1],activity,Dactivity)

            par.up <- c(sj=sd(activity),sigmaGP=sd(activity),logProt=log(2*DT),logtauGP=log(2*DT))
            par.low <- c(sj=0,sigmaGP=0,logProt=-log(2*DT),logtauGP=-log(2*DT))
            Ntry <- 100
            for(j3 in 1:Ntry){
                if(Ntry>1){
                    start <- runif(length(par.low),par.low,par.up)
                    names(start) <- names(par.low)
                }
                tmp <- try(sopt(omega=NA,phi=NA,Nma=0,Nar=0,Indices=NULL,data=data,type='noise',par.low=par.low,par.up=par.up,start=as.list(start),noise.only=TRUE,GP=TRUE,gp.par=rep(NA,3),Nrep=Nsamp),TRUE)
                if(class(tmp)!='try-error'){
                    ll <- c(ll,tmp$logLmax)
                    pp <- rbind(pp,unlist(tmp$par0))
                }
            }
            ind <- which.max(ll)
            logProt <- pp[ind,'logProt']
            logtauGP <- pp[ind,'logtauGP']
            sigmaGP <- pp[ind,'sigmaGP']
        }
#        per <- BFP(tab[,1],tab[,5],dHalpha,Nma=0,Nar=0,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=TRUE,gp.par=rep(NA,3),noise.only=TRUE,Nsamp=Nsamp,sampling=sampling)

        if(grepl('MA',noise.type)){
            Nma <- pq
        }else{
            Nma <- 0
                                        #        Nma <- 1
        }
        if(grepl('AR',noise.type)){
            Nar <- pq
        }else{
            Nar <- 0
        }
        Nars <- c(Nars,Nar)
        Nmas <- c(Nmas,Nma)
        gp.par <- rep(NA,3)
        if(grepl('GP',noise.type)){
            GP <- TRUE
            tmp <- as.integer(gsub('GP|\\+.+','',noise.type))
            if(floor(tmp/100)){
                gp.par[1] <- sigmaGP
            }
            if(floor((tmp%%100)/10)){
                gp.par[2] <- logProt
            }
            if(floor(tmp%%10)){
                gp.par[3] <- logtauGP
            }
            if(grepl('HD10700_Sindex',target)){
                tab[,2] <- tab[,2]*1e3
                tab[,3] <- tab[,3]*1e3
            }
        }else{
            GP <- FALSE
        }
        par.all <- ws.all <- sig.all <- c()
        for(k in 1:Nmax){
            t1 <- proc.time()
            if(k==1){
                res <- tab[,2]
            }else{
                res <- per$res.s
            }
            Indices <- NULL
            if(!all(ind.proxy==0)) Indices <- tab[,ind.proxy+3]
###first find the baseline optimal parameters
            if(refine){
                per <- BFP(tt,res,tab[,3],Nma=Nma,Nar=Nar,Indices=Indices,ofac=ofac/10,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=GP,gp.par=gp.par,noise.only=noise.only,Nsamp=Nsamp,sampling=sampling)
####Then find the 1-planet solution
                par.opt <- per$par[1,!grepl('^A$|^B$|^gamma$|^beta$',colnames(per$par))]
            }else{
                par.opt <- NULL
            }
            per <- BFP(tt,res,tab[,3],Nma=Nma,Nar=Nar,Indices=Indices,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=GP,gp.par=gp.par,noise.only=noise.only,Nsamp=Nsamp,sampling=sampling,par.opt=par.opt)
            ws.all <- cbind(ws.all,per$res.nt)#for phase curve
            sig.all <- cbind(sig.all,per$ysig)#for phase curve

            dur <- (proc.time()-t1)[3]
            cat('computation time for ',noise.type,' for ',k,'signal is ',dur,'s\n')
            durs <- c(durs,dur)
            if(k==1 & j1==1){
                out <- cbind(out,per$P,per$power)
                cnames <- c(cnames,'period',paste0('noise',noise.type,'Nsig',k))
            }else{
                out <- cbind(out,per$power)
                cnames <- c(cnames,paste0('noise',noise.type,'Nsig',k))
            }
            par.all <- rbind(par.all,per$par[1,])
            Popt <- c(Popt,per$P[which.max(per$power)])
            logBFs <- c(logBFs,max(per$logBF))
            source('agatha_plot.R')
        }
        opt.par[[j1]] <- par.all
        ws[[j1]] <- ws.all
        sigs[[j1]] <- sig.all
    }
####plot activity white noise periodogram
    Pact <- c()
    if(Sindex){
       acts <- cbind(tab[iss,1],ss,dss)
       act.name <- 'Sindex'
       source('activity.R')
       Pact <- rbind(Pact,per$Popt)
    }
    if(Halpha){
       acts <- cbind(tab[iss,1],hh,dhh)
       act.name <- 'Halpha'
       source('activity.R')
       Pact <- rbind(Pact,per$Popt)
    }
    if(ncol(tab)>5){
        acts <- cbind(tab[,1],tab[,6],sqrt(tab[,6]))
        act.name <- 'Photon Count'
        source('activity.R')
        Pact <- rbind(Pact,per$Popt)
    }

    colnames(out) <- cnames
    cat('output Robj:\n',fobj,'\n')
    LBF <- Nact <- Psig <- array(NA,dim=c(Nmax,length(noise.types)))
    Ndet <- c()
    for(j in 1:Nmax){
        Ps <- Popt[Nmax*((1:length(noise.types))-1)+j]
        LBF[j,] <- logBFs[Nmax*((1:length(noise.types))-1)+j]
        Psig[j,] <- Ps
        ind <- which(abs(Ps-median(Ps))/median(Ps)<0.1)
        Ndet <- c(Ndet,length(ind))
        for(m in 1:length(Ps)){
            Nact[j,m] <- length(which(abs(Ps[m]-Pact)/Ps[m]<0.1))
        }
    }
    save(out,Popt,opt.par,Psig,Ndet,Nact,durs,file=fobj)
    source('phase_plot.R')
    dev.off()
    solutions[[target]]$par <- par.shows
    solutions[[target]]$P <- Pshows
    solutions[[target]]$Psig <- Psig
    solutions[[target]]$Ksig <- Ksig
    solutions[[target]]$logBF <- LBF
    solutions[[target]]$Ndet <- Ndet
    solutions[[target]]$Nact <- Nact
    solutions[[target]]$fobj <- fobj
    solutions[[target]]$fpdf <- fpdf
}
##final output
ff <- paste0('results/solutions',ind.min,'to',ind.max,'.Robj')
cat('output solutions:\n',ff,'\n')
save(solutions,file=ff)
