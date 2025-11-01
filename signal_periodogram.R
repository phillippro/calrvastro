plotf <- TRUE
Nmas <- Nars <- durs <- c()
cnames <- c()
out <- c()
opt.par <- ws <- sigs <- list()
Aopt <- Bopt <- Popt <- logBFs <- array(NA,dim=c(Nmax,length(noise.types)))
for(j1 in 1:length(noise.types)){
    noise.type <- noise.types[j1]
    cat('noise.type=',noise.type,'\n')
    start <- c(sj=0,sigmaGP=0,logProt=log(Protation),logtauGP=1)
    if(GPtype=='MA+GP'){
        start <- c(start,m1=0.1,logtau=0)
    }
    ll <- prot <- c()
    if(grepl('GP',noise.type)){
###fit the GP kernel to Halpha
###Determine error of Halpha by scaling the RV error
        Prot <- mean(Protations)
        Prot.min <- max(Prot-3*sd(Protations),0)
        Prot.max <- Prot+3*sd(Protations)
        n1 <- sapply(1:nrow(Pact),function(k) any(abs(Pact[k,]-Prot)/Prot<0.1))
        activity.name <- rownames(Pact)[n1]
        if(any(activity.name=='Halpha') & Halpha){
            ind.act <- ihh
            activity <- hh
            Dactivity <- dhh
        }else if(any(activity.name=='Sindex') & Sindex){
            ind.act <- iss
            activity <- ss
            Dactivity <- dss
        }else if(Halpha){
            ind.act <- ihh
            activity <- hh
            Dactivity <- dhh
        }else if(Sindex){
            ind.act <- iss
            activity <- ss
            Dactivity <- dss
        }else{
            ind.act <- ipp
            activity <- pp
            Dactivity <- dpp
        }

        data <- cbind(tab[ind.act,1]-min(tab[ind.act,1]),activity,Dactivity)

        par.up <- c(sj=10,sigmaGP=10,logProt=log(Prot.max),logtauGP=log(2*DT))
        par.low <- c(sj=0,sigmaGP=0,logProt=log(Prot.min),logtauGP=0)
        if(GPtype=='MA+GP'){
            par.up <- c(par.up,m1=1,logtau=log(2*DT))
            par.low <- c(par.low,m1=0,logtau=0)
        }

        Ntry <- 100
        for(j3 in 1:Ntry){
            if(Ntry>1){
                start <- runif(length(par.low),par.low,par.up)
                names(start) <- names(par.low)
            }
            tmp <- try(sopt(omega=NA,phi=NA,Nma=1,Nar=0,Indices=NULL,data=data,type='noise',par.low=par.low,par.up=par.up,start=as.list(start),noise.only=TRUE,GP=TRUE,gp.par=rep(NA,3),Nrep=Nsamp),TRUE)
            if(class(tmp)!='try-error'){
                ll <- c(ll,tmp$logLmax)
                prot <- rbind(prot,unlist(tmp$par0))
                break()
            }
        }
        cat('GP parameter resampling times=',j3,'\n')
        ind <- which.max(ll)
        logProt <- prot[ind,'logProt']
        logtauGP <- prot[ind,'logtauGP']
        sigmaGP <- prot[ind,'sigmaGP']
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
        Popt[k,j1] <- per$P[which.max(per$power)]
        logBFs[k,j1] <- max(per$logBF)
        source('agatha_plot.R')
    }
    opt.par[[j1]] <- par.all
    ws[[j1]] <- ws.all
    sigs[[j1]] <- sig.all
}
