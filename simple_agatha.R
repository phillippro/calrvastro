library(magicaxis)
source('periodograms.R')
source('periodoframe.R')
options(warn=2)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    fin <- args[1]
    Nmax <- as.integer(args[2])
    ofac <- as.numeric(args[3])
    noise.types <- args[4]
    pq <- as.integer(args[5])
    Nsamp <- as.integer(args[6])
}else{
    fin <- '/Users/ffeng/Documents/projects/dwarfs/bary/data/HD128621_HARPSbin.rv'
    Nmax <- 3
    ofac <- 1
#    noise.types <- c('white','MA','AR','GP011')
#    noise.types <- c('white','MA','AR')
#    noise.types <- c('MA')
#    noise.types <- 'white'
#    noise.types <- 'GP111'
    pq <- 1
    Nsamp <- 1
noise.types <- 'AR'
#noise.types <- c('white')
}
#refine <- TRUE
refine <- FALSE
sampling <- 'combined'
#folder <- '../data/aperture/'
#folder <- '../auto/data/'
#folder <- paste0('../data/',ins,'/')
target <- star <- gsub('.+\\/|_.+','',fin)
tab <- read.table(fin)
ind.proxy <- 0#
cat('input data file:\n',fin,'\n')
if(class(tab[,1])=='factor'){
    tab <- read.table(fin,header=TRUE)
}
if(any(duplicated(tab[,1]))){
    tab <- tab[-which(duplicated(tab[,1])),]
}
tab[,1] <- tab[,1]%%2400000
trv <- tab[,1]
ind <- sort(trv,index.return=TRUE)$ix
trv <- trv[ind]
tab <- tab[ind,]
#noise.only <- TRUE
noise.only <- FALSE
sigmaGP <- sd(tab[,2])#m/s
logProt <- log(39)
logtauGP <- log(100)

Nmas <- Nars <- durs <- c()
out <- c()
opt.par <- list()
#fmin <- 1/(4*(max(trv)-min(trv)))
fmin <- 1/1000
fmax <- 1/1.1
Popt <- c()
for(j in 1:length(noise.types)){
    noise.type <- noise.types[j]
    cat('noise.type=',noise.type,'\n')
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
    for(k in 1:Nmax){
        t1 <- proc.time()
        if(k==1){
            res <- tab[,2]
        }else{
            res <- per$res.s
#            Nma <- 0
#            Nar <- 0
        }
        Indices <- NULL
        if(!all(ind.proxy==0)) Indices <- tab[,ind.proxy+3]
###first find the baseline optimal parameters
        if(refine){
            per <- BFP(tab[,1],res,tab[,3],Nma=Nma,Nar=Nar,Indices=Indices,ofac=ofac/10,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=GP,gp.par=gp.par,noise.only=noise.only,Nsamp=Nsamp,sampling=sampling,renew=TRUE)
####Then find the 1-planet solution
            par.opt <- per$par[1,!grepl('^A$|^B$|^gamma$|^beta$',colnames(per$par))]
        }else{
            par.opt <- NULL
        }
        per <- BFP(tab[,1]-min(tab[,1]),res,tab[,3],Nma=Nma,Nar=Nar,Indices=Indices,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=GP,gp.par=gp.par,noise.only=noise.only,Nsamp=Nsamp,sampling=sampling,par.opt=par.opt,renew=TRUE)
        dur <- (proc.time()-t1)[3]
        cat('computation time for ',noise.type,' for ',k,'signal is ',dur,'s\n')
        durs <- c(durs,dur)
        if(j==1 & k==1){
            out <- cbind(out,per$P,per$power)
        }else{
            out <- cbind(out,per$power)
        }
        Popt <- c(Popt,per$P[which.max(per$power)])
        opt.par[[j]] <- per$par
    }
}
####plot
fname <- paste0('results/',target,'_Nsamp',Nsamp,'_',paste(noise.types,collapse='_'),'_p',Nar,'q',Nma,'_ind',paste(ind.proxy,collapse='.'),'Nmax',Nmax,'_noise',noise.only,'_ofac',ofac,'logtauGP',round(logtauGP),'_P',paste(round(Popt,1),collapse='d'),'.pdf')
cat('output pdf:\n',fname,'\n')
source('plot_agatha.R')
save(list=ls(all=TRUE),file=gsub('.pdf','.Robj',fname))
