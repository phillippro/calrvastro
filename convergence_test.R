#####This file is to judge whether the chain is converged or not. I divide the long chain into equal-length subchains and apply the Gelman-Rubin criteria to them. also see Ford 2005

###deal with angular variables
N2 <- length(out$mcmc.opt)
sigs <- names(out$mcmc.opt)
#Nsigs <- sort(as.integer(gsub('sig','',names(out$mcmc.opt))))
for(n in sigs){
    coln <- colnames(out[['mcmc.opt']][[n]])
#    angs <- c('omega','Mo','Omega','Inc')
    angs <- c('omega','Mo','Omega')
    for(ang in angs){
        ind <- grep(ang,coln)
	rad <- 2*pi
#	if(ang=='Inc') rad <- pi
        if(length(ind)>0){
            for(j in ind){
                x <- out[['mcmc.opt']][[n]][,j]%%rad
                x[(x-median(x))> rad/2] <- x[(x-median(x))> rad/2]-rad
                x[(x-median(x))< -rad/2] <- x[(x-median(x)) < -rad/2]+rad
                out[['mcmc.opt']][[n]][,j] <- x
            }
        }
    }
    ind.Tc <- grep('Tc',coln)
    if(length(ind.Tc)>0){
        for(j in ind.Tc){
            x <- out[['mcmc.opt']][[n]][,j]
            P <- exp(out[['mcmc.opt']][[n]][,paste0('per',j)])
            x[(x-median(x))>P] <- x[(x-median(x))>P]-P*floor((x[(x-median(x))>P]-median(x))/P)
            x[(x-median(x))< -P] <- x[(x-median(x))< -P]+P*floor((median(x)-x[(x-median(x))< -P])/P)
        }
    }
}

####calculate Gelman-Rubin criterion Rhat for the final solution
n <- paste0('sig',Nsig)
Npar <- ncol(out[['mcmc.opt']][[n]])-2
if(!save.memory){
    Nchain <- length(out[['mcmc.all']][[n]])
    ind.chain <- which.max(sapply(1:Nchain,function(j) max(out[['mcmc.all']][[n]][[j]][,Npar+2])))
}else{
    Nchain <- 1
    ind.chain <- 1
}
#ind.chain <- 1:Nchain
if(Nsig>0 & chain.type=='multiple' & !save.memory){
    dlnPs <- lnPs <- c()
    lls <- c()
    for(j in 1:Nchain){
        ind <- which.max(out[['mcmc.all']][[n]][[j]][,Npar+2])
        lnPs <- cbind(lnPs,out[['mcmc.all']][[n]][[j]][ind,5*(1:Nsig)-4])
        dlnP <- sapply(1:Nsig,function(i) sd(out[['mcmc.all']][[n]][[j]][,5*i-4]))
        dlnPs <- cbind(dlnPs,dlnP)
        lls <- c(lls,)
    }
    ind.max <- which.max(lls)
    delta.lnP <- t(sapply(1:nrow(lnPs),function(k) abs(lnPs[k,]-lnPs[k,ind.max])))
    sigma.lnP <- t(sapply(1:nrow(lnPs),function(k) sqrt(dlnPs[k,]^2+dlnPs[k,ind.max]^2)))
    ind.chain <- which(sapply(1:Nchain,function(i2) all(delta.lnP[,i2]<sigma.lnP[,i2])))
}

###calculate Rhat
meanvar <- rep(NA,Npar)
mean.all <- rep(NA,Npar)
var.single.est <- rep(NA,Npar)
par.name <- colnames(out[['mcmc.opt']][[n]][,1:Npar])
if(Nchain>1 & length(ind.chain)>1 & !save.memory){
    subchain.len <- nrow(out[['mcmc.all']][[n]][[1]])
    var.sub <- array(data=NA,dim=c(length(ind.chain),Npar))
    mean.sub <- array(data=NA,dim=c(length(ind.chain),Npar))
    for(j in 1:Npar){
        for(i1 in 1:length(ind.chain)){
            i <- ind.chain[i1]
            var.sub[i1,j] <- var(out[['mcmc.all']][[n]][[i]][,j])
            mean.sub[i1,j] <- mean(out[['mcmc.all']][[n]][[i]][,j])
        }
        meanvar[j] <- mean(var.sub[,j])
        mean.all[j] <- mean(mean.sub[,j])
        var.single.est[j] <- subchain.len/Nchain*sum((mean.sub[,j]-mean.all[j])^2)
    }
}else{
    Nsub <- 4
    var.sub <- array(data=NA,dim=c(Nsub,Npar))
    mean.sub <- array(data=NA,dim=c(Nsub,Npar))
    len <- nrow(out[['mcmc.opt']][[n]])
    dN <- floor(len/Nsub)
    subchain.len <- dN
    for(j in 1:Npar){
        for(i in 1:Nsub){
            inds <- (i-1)*dN+(1:dN)
            var.sub[i,j] <- var(out[['mcmc.opt']][[n]][inds,j])
            mean.sub[i,j] <- mean(out[['mcmc.opt']][[n]][inds,j])
        }
        meanvar[j] <- mean(var.sub[,j])
        mean.all[j] <- mean(mean.sub[,j])
        var.single.est[j] <- subchain.len/Nsub*sum((mean.sub[,j]-mean.all[j])^2)
    }
}
var.est <- (subchain.len-1)/subchain.len*meanvar+1/subchain.len*var.single.est
Rhat <- sqrt(var.est/meanvar)
names(Rhat) <- par.name
ind <- grep('Tc|Mo|omega',par.name)
if(length(ind)>0) Rhat <- Rhat[-ind]
if(any(is.na(Rhat))) Rhat[is.na(Rhat)] <- 1
if(any(Rhat>1.1)){
    conv <- FALSE
}else{
    conv <- TRUE
}
cat('convergence parameter Rhat=',Rhat,'\n')
cat('conv=',conv,'\n')
