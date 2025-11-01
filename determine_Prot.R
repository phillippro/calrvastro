source('periodograms.R')
source('periodoframe.R')
kappa <- 1
periodogram.type <- 'GLST'
Nsamp <- 1
power.act <- Psig <- Pact <- c()
act.data <- sig.name <- act.names <- c()
Nact <- 3
for(i in ins){
    p <- names(out$BFP[[i]])
    ind1 <- which(!grepl('sig|window',p))
    for(j in p[ind1]){
        q <- names(out$BFP[[i]][[j]])
        for(k in q){
            per <- out$BFP[[i]][[j]][[k]]
            if(!is.null(per)){
                inds <- sort(per$power,decreasing=TRUE,index.return=TRUE)$ix[1:Nact]
                Pact <- rbind(Pact,as.numeric(per$P[inds]))
                power.act <- rbind(power.act,as.numeric(per$power[inds]))
                act.name <- paste0(i,'_',j,'_',k)
                act.names <- c(act.names,act.name)
                act.data[[act.name]] <- out[[i]][['activity']][[j]]
            }
        }
    }
}
colnames(Pact) <- paste0('peak',1:Nact)
rownames(Pact) <- act.names

conservative <- FALSE
###find rotation period
Ntry <- 10
Nboot <- 1e2
Nrot <- array(NA,dim=c(nrow(Pact),ncol(Pact)))
NN <- rep(NA,nrow(Pact))
power.rot <- array(NA,dim=c(nrow(Pact),ncol(Pact)))
for(j in 1:nrow(Pact)){
    for(k in 1:ncol(Pact)){
        index <- c(1:nrow(Pact))[-j]
#        inds <- which(abs(Pact[index,]-Pact[j,k])/Pact[j,k] < 0.1,arr.ind=TRUE)
        inds <- which(abs(Pact[index,]-Pact[j,k])/Pact[j,k] < dnorm(log(Pact[j,k]),log(100),log(100)),arr.ind=TRUE)
        Nrot[j,k] <- nrow(inds)
        NN[j] <- nrow(act.data[[act.names[j]]])
        power.rot[j,k] <- sum(power.act[index,][inds])
    }
}
#ind.rot <- which(Nrot==max(Nrot),arr.ind=TRUE)#number of activity indices (including window function and photon count) showing rotation signals
ind.rot <- which(Nrot==max(Nrot),arr.ind=TRUE)
if(nrow(ind.rot)>1){
    Ns <- c()
    for(k in 1:nrow(ind.rot)){
        j <- ind.rot[k,]
        inds <- (1:nrow(ind.rot))[-k]
        Ns <- c(Ns,length(which(abs(Pact[inds,]-Pact[j[1],j[2]])/Pact[j[1],j[2]]<0.1)))
    }
    ind.rot <- ind.rot[Ns==max(Ns),]
}
#ind.rot <- which(Nrot[,1]==max(Nrot[,1]))
#ind.rot <- ind.rot[which.max(power.rot[ind.rot]),]
#ind.rot <- ind.rot[which.max(power.rot[ind.rot,1]),]
Protation <- Pact[ind.rot[1],ind.rot[2]]
###find the activity indicies which has the most data points
inds <- which(abs(Pact-Protation)/Protation<0.1,arr.ind=TRUE)
nn <- unique(inds[,1])
ii <- nn[which.max(NN[nn])]
rot.name <- act.names[ii]
rot.names <- act.names[nn]
cat('activity indices:',rot.name,'\n')
acts0 <- act.data[[act.names[ii]]]
plotf <- FALSE
Protations <- power.rots <- c()
for(j in 1:Nboot){
    fmin1 <- 0.9/Protation
    fmax1 <- 1.1/Protation
    if(periodogram.type=='BFP'){
        ofac <- 0.01
    }else{
        ofac <- 1
    }
    acts <- acts0[sort(sample(1:nrow(acts0),nrow(acts0),replace=TRUE)),]
    if(periodogram.type=='BFP'){
        per <- BFP(acts[,1],acts[,2],acts[,3],Nma=1,Nar=0,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin1,fmax=fmax1,quantify=TRUE,progress=FALSE,GP=FALSE,noise.only=FALSE,Nsamp=Nsamp,sampling='random')
    }else if(periodogram.type=='GLST'){
        per <- glst(t=acts[,1],y=acts[,2],err=acts[,2],ofac=ofac,fmin=fmin1,fmax=fmax1,sampling='random')
    }else if(periodogram.type=='LS'){
        per <- lsp(times=acts[,1],x=acts[,2],ofac=ofac,from=fmin1,to=fmax1,alpha=c(0.1,0.01,0.001))
    }
    if(periodogram.type=='GLST'){
        power.rots <- c(power.rots,max(per$lnBFs))
    }else{
        power.rots <- c(power.rots,max(per$power))
    }
    Protations <- c(Protations,per$P[which.max(per$power)])
}
plotf <- TRUE
rownames(Pact) <- act.names
colnames(Pact) <- paste(1:Nact,'sig')
