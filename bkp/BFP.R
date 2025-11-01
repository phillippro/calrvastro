BFP <- function(t, y, dy, Nma=0, Nar=0,Inds=Inds,Indices=Indices,sj=0,logtau=NULL,ofac=1, fmax=NULL,fmin=NA,tspan=NULL,sampling='combined',model.type='man',progress=TRUE,quantify=FALSE,dP=0.1,section=1,GP=TRUE,gp.par=rep(NA,3),noise.only=FALSE){
###gp.par is the free parameters of SHO-GP
    if(Nma==0 & Nar==0 & !GP) noise.only <- FALSE
    if(noise.only) quantify <- FALSE
    unit <- 1
    if(length(Inds)>0){
        if(all(Inds==0)){
            NI <- 0
        }else{
            Inds <- Inds[Inds>0]
            NI <- length(Inds)
            Indices <- as.matrix(Indices[,Inds,drop=FALSE])
            for(j in 1:ncol(Indices)){
                Indices[,j] <- scale(Indices[,j])
            }
        }
    }else{
        Inds <- 0
    }
    data <- cbind(t,y,dy)
    if(NI==0) d <- 0
    if(Nar==0) l <- 0
    if(Nma==0) m <- 0
    if(is.null(tspan)){
        tspan <- max(t)-min(t)
    }
    step <- 1/(tspan*ofac)
    if(is.na(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    Nrep <- 1
    if(quantify) Nrep <- 2
    f <- fsample(fmin,fmax,sampling,section,ofac,unit)
    nout <- length(f)
    t <- (t-min(t))/unit#rescale time
    Ndata <- length(t)
    omegas <- 2*pi*f
    phi <- 0
    if(is.null(logtau)){
        logtau <- log(median(diff(t)))
    }
    #######define notations and variables
    vars1 <- global.notation(t,y,dy,Indices,Nma,Nar,NI,GP,gp.par)
    if(noise.only){
        vars0 <- global.notation(t,y,dy,Indices,Nma=0,Nar=0,NI=0,GP=FALSE,gp.par=gp.par)
    }else{
        vars0 <- vars1
    }
#    var <- names(vars)
####shor
#    for(k in 1:length(var)){
#        assign(var[k],vars[[var[k]]])
#    }
##    
##########################################
#####optimizing the noise parameters
#########################################
    if(model.type=='auto'){
        t1 <- proc.time()
        out <- bfp.inf(vars1,Indices,logtau=NULL,m=NULL)
        t2 <- proc.time()
        dur <- format((t2-t1)[3],digit=3)
        cat('model comparison computation time:',dur,'s\n\n')
        NI <- out$NI
        Nma <- out$Nma
        Nar <- out$Nar
        vars1 <- out$vars1
        Indices <- out$Indices
    }
#    m <- rep(vars$mini,Nma)
#    l <- rep(vars$lini,Nar)
#    logtau <- vars$logtau.ini
#    d <- rep(vars$dini,NI)
########################################################
##############baseline model
########################################################    
    tmp <- sopt(omega=NA,phi=NA,Nma=vars0$Nma,Nar=vars0$Nar,NI=vars0$NI,data=data,type='noise',par.low=vars0$par.low,par.up=vars0$par.up,start=vars0$start,noise.only=noise.only,GP=vars0$GP,gp.par=vars0$gp.par)
    logLmax0 <- tmp$logLmax
##
    logLmax <- rep(NA,length(omegas))
#    cat('names(tmp$par)=',names(tmp$par),'\n')
    Npar <- length(tmp$par)
    if(noise.only){
        Nextra <- 0
        if(GP){
            ind <- which(is.na(gp.par[-2]))
            Nextra <- Nextra+length(ind)
        }
        if(Nma>0){
            Nextra <- Nextra+Nma
        }
        if(Nar>0){
            Nextra <- Nextra+Nar
        }
    }else{
        Nextra <- 2
    }
    omegas.all <- c()
    opt.pars <- array(NA,dim=c(length(f),Npar+Nextra))
    P <- c()
    for(nn in 1:Nrep){
        P <- c(unit/f,P)
        if(progress){
            withProgress(message = 'Calculating BFP', value = 0, {
              ###later add on
            })
        }else{
            for(kk in 1:length(f)){
#####for periodogram, set sj=0
#                cat(kk,'/',length(f),'\n')
                omega <- omegas[kk]
################################################
####I optimization
####################################################
                t0 <- proc.time()
                if(noise.only){
                    if(GP){
                        gp.par[2] <- log(2*pi/omega)
                        vars1$gp.par <- gp.par
                    }
                    opt <- sopt(omega=omega,phi=phi,Nma=Nma,Nar=Nar,NI=NI,data=data,type='noise',par.low=vars1$par.low,par.up=vars1$par.up,start=vars1$start,noise.only=noise.only,GP=GP,gp.par=gp.par)
                }else{
                    tmp <- local.notation(t,y,dy,Indices,NI,omega,phi)
                    vars1 <- c(vars1,tmp)
                    opt <- sopt(omega=omega,phi=phi,Nma=Nma,Nar=Nar,NI=NI,data=data,type='period',par.low=vars1$par.low,par.up=vars1$par.up,start=vars1$start,noise.only=noise.only,GP=GP,gp.par=gp.par)
                }
                vars1$start <- opt$par0
                logLmax[kk] <- opt$logLmax
                chi2 <- opt$chi2
                opt.pars[kk,] <- unlist(opt$par)
            }
        }
####signals
        Pmax <- P[which.max(logLmax)]
        if(quantify & nn==1){
            P1 <- P
            logLmax1 <- logLmax
            opt.pars1 <- opt.pars
            ind.max <- which.max(logLmax)
            fmin <- (1-dP)*f[ind.max]
            fmax <- (1+dP)*f[ind.max]
###oversampling
            f <- seq(fmin,fmax,by=step/20)*unit
                                        #        f <- seq(fmin,fmax,length.out=1000)*unit
            Nf <- length(f)
            omegas <- f*2*pi
            logLmax <- c(rep(NA,length(f)),logLmax)
            opt.pars <- rbind(array(data=NA,dim=c(length(f),ncol(opt.pars))),opt.pars)

        }else if(quantify & nn==2){
            ind.max <- which.max(logLmax)
            ind1 <- which.max(logLmax1)
            if(ind.max<Nf){
                P1[ind1] <- P[ind.max]
                logLmax1[ind1] <- logLmax[ind.max]
                opt.pars1[ind1,] <- opt.pars[ind.max,]
            }
            P <- P1
            logLmax <- logLmax1
            opt.pars <- opt.pars1
        }
        omegas.all <- c(omegas.all,omegas)
    }
    colnames(opt.pars) <- names(opt$par)
    if(noise.only){
        name0 <- names(opt$par)
        cat('names(opt$par)=',names(opt$par),'\n')
        if(Nma>0){ 
            opt.pars <- cbind(opt.pars,log(2*pi/omegas.all))
            name0 <- colnames(opt.pars) <- c(name0,'logtau')            
        }
        if(Nar>0){
            opt.pars <- cbind(opt.pars,log(2*pi/omegas.all))
            colnames(opt.pars) <- c(name0,'logtauAR')
        }
    }
    logBF <- logLmax-logLmax0-Nextra/2*log(length(y))#BIC-estimated BF; the extra free parameter n=2, which are {A, B}
####signals
    ind.max <- which.max(logBF)
    inds <- sort(logBF,decreasing=TRUE,index.return=TRUE)$ix[1:10]
    Popt <- P[inds]
    opt.par <- opt.pars[inds,]
#    rv.red()
    logBF.opt <- logBF[inds]
####calculate the residual
    par.fix <- list(omega=2*pi/Popt[1],phi=0)
#    par.fix <- list(omega=2*pi/131,phi=0)
    df <- list(data=cbind(t,y,dy),Indices=Indices,par.fix=par.fix,Nma=Nma,Nar=Nar,NI=NI,GP=GP)
    if(is.matrix(opt.par) | is.data.frame(opt.par)){
        pp <- as.list(opt.par[1,])
    }else{
        pp <- as.list(opt.par)
    }
    names(pp) <- colnames(opt.pars)
    yall <- RV.model(pp,data=df)
    if(!noise.only){
        ysig <- pp$A*cos(2*pi/Popt[1]*t)+pp$B*sin(2*pi/Popt[1]*t)
        ysigt <- ysig+pp$gamma+pp$beta*t
    }else{
        ysig <- ysigt <- rep(0,length(t))
    }
    yred <- yall-ysigt
    yredt <- yall-ysig
    res.nst <- y-yall
    res.n <- y-yred
    res.nt <- y-yredt
    res.s <- y-ysig
    res.st <- y-ysigt
    inds <- sort(logBF,decreasing=TRUE,index.return=TRUE)$ix
    ps <- P[inds[1:5]]
    power.opt <- logBF[inds[1:5]]
    return(list(data=data,logBF=logBF,P=P,Popt=Popt,logBF.opt=logBF.opt,par=opt.par,res.nst=res.nst,res=res.nst,res.n=res.n,res.s=res.s,res.st=res.st,res.nt=res.nt,sig.level=log(150), power=logBF,ps=ps,power.opt=power.opt,ysig=ysig,pars=opt.pars))
}
