#####This file is to judge whether the chain is converged or not. I divide the long chain into equal-length subchains and apply the Gelman-Rubin criteria to them. also see Ford 2005
###deal with angular variables
Nchain <- length(out[['mcmc.all']][[paste0('sig',Nsig)]])
if(Nsig>0){
    for(k in 1:Nchain){
        coln <- colnames(out[['mcmc.all']][[paste0('sig',Nsig)]][[k]])
        omega = omega1 = omega2 = omega3 = out[['mcmc.all']][[paste0('sig',Nsig)]][[k]][,grep('omega([[:digit:]]{1})',coln)]
        if(any(grepl('Mo',coln))){
            Mo = Mo1 = Mo2 = Mo3 = out[['mcmc.all']][[paste0('sig',Nsig)]][[k]][,grep('Mo([[:digit:]]{1})',coln)]
        }
        if(Nsig==1){
            omega = omega1 = omega2 = omega3 = matrix(omega1)
            if(any(grepl('Mo',coln))){
                Mo = Mo1 = Mo2 = Mo3 = matrix(Mo1)
            }
        }
        for(j in 1:Nsig){
            ind.o2 <- which(omega1[,j]>pi)
            ind.o3 <- which(omega1[,j]<pi)
            omega2[ind.o2,j] <- omega2[ind.o2,j]-2*pi
            omega3[ind.o3,j] <- omega3[ind.o3,j]+2*pi
            k <- which.min(c(sd(omega1[,j]),sd(omega2[,j]),sd(omega3[,j])))
            if(k==1){
                omega[,j] <- omega1[,j]
            }else if(k==2){
                omega[,j] <- omega2[,j]
            }else{
                omega[,j] <- omega3[,j]
            }
            if(any(grepl('Mo',coln))){
                ind.o2 <- which(Mo1[,j]>pi)
                ind.o3 <- which(Mo1[,j]<pi)
                Mo2[ind.o2,j] <- Mo2[ind.o2,j]-2*pi
                Mo3[ind.o3,j] <- Mo3[ind.o3,j]+2*pi
                k <- which.min(c(sd(Mo1[,j]),sd(Mo2[,j]),sd(Mo3[,j])))
                if(k==1){
                    Mo[,j] <- Mo1[,j]
                }else if(k==2){
                    Mo[,j] <- Mo2[,j]
                }else{
                    Mo[,j] <- Mo3[,j]
                }
                out[['mcmc.all']][[paste0('sig',Nsig)]][[k]][,grep('Mo([[:digit:]]{1})',coln)] <- Mo
            }
        }
        out[['mcmc.all']][[paste0('sig',Nsig)]][[k]][,grep('omega([[:digit:]]{1})',coln)] <- omega
    }
}

####calculate Gelman-Rubin criterion
Npar <- ncol(out[['mcmc.all']][[paste0('sig',Nsig)]][[1]])
subchain.len <- nrow(out[['mcmc.all']][[paste0('sig',Nsig)]][[1]])
var.sub <- array(data=NA,dim=c(Nchain,Npar))
mean.sub <- array(data=NA,dim=c(Nchain,Npar))
meanvar <- rep(NA,Npar)
mean.all <- rep(NA,Npar)
var.single.est <- rep(NA,Npar)
for(j in 1:Npar){
    for(i in 1:Nchain){
        var.sub[i,j] <- var(out[['mcmc.all']][[paste0('sig',Nsig)]][[i]][,j])
        mean.sub[i,j] <- mean(out[['mcmc.all']][[paste0('sig',Nsig)]][[i]][,j])
    }
    meanvar[j] <- mean(var.sub[,j])
    mean.all[j] <- mean(mean.sub[,j])
    var.single.est[j] <- subchain.len/Nchain*sum((mean.sub[,j]-mean.all[j])^2)
}
var.est <- (subchain.len-1)/subchain.len*meanvar+1/subchain.len*var.single.est
Rhat <- sqrt(var.est/meanvar)
if(any(is.na(Rhat))) Rhat[is.na(Rhat)] <- 1
if(any(Rhat>1.1)){
    conv <- FALSE
}else{
    conv <- TRUE
}
cat('convergence parameter Rhat=',Rhat,'\n')
cat('conv=',conv,'\n')
