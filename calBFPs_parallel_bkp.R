###list the names of all periodograms
ns0 <- c()
for(i in ins){
    acts <- names(out[[i]][['activity']])
    for(act in acts){
        if(grepl('window',act)){
            ns0 <- rbind(ns0,c(i,act,'white'))
        }else if(!grepl('AP',act)){
            ns0 <- rbind(ns0,c(i,act,'MA'))
        }
    }
}

if(any(grepl('photometry',names(out)))){
    ns0 <- rbind(ns0,c('photometry','photometry','MA'))
    ns0 <- rbind(ns0,c('photometry','photometry','GP'))
}

ns1 <- c()
for(j in 1:Nset){
    for(k in 1:length(noise.types)){
            ns1 <- rbind(ns1,c(ins[j],noise.types[k]))
    }
}

###activity periodograms
bfp.act <- foreach(j = 1:nrow(ns0)) %dopar% {
    if(ns0[j,1]=='photometry'){
        data <- out[['photometry']]$data
    }else{
        data <- out[[ns0[j,1]]][['activity']][[ns0[j,2]]]
    }
    Nma <- Nar <- 0
    GP <- FALSE
    gp.par <- rep(NA,3)
    noise.only <- FALSE
    ofac1 <- ofac
    if(ns0[j,3]=='MA'){
        Nma <- 1
    }
    if(ns0[j,3]=='AR'){
        Nar <- 1
    }
    if(ns0[j,3]=='GP'){
        GP <- TRUE
        gp.par <- rep(NA,3)
        noise.only <- TRUE
        ofac1 <- min(1,ofac)
    }
    cat('calulate BFP for',ns0[j,2],'for',ns0[j,3],'for',ns0[j,1],'RV set of',target,'\n')
    if(ns0[j,2]=='window'){
        per <- lsp(times=data[,1],x=data[,2],ofac=ofac,from=fmin,to=fmax,alpha=c(0.1,0.01,0.001))
    }else{
        per <- BFP(data[,1],data[,2],data[,3],Nma=0,Nar=0,Indices=NULL,ofac=ofac1,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=GP,noise.only=noise.only,Nsamp=Nsamp,gp.par=gp.par,renew=TRUE)
    }
    c(list(per=per),list(name=ns0[j,]))
}
out$BFP <- list()
for(j in 1:nrow(ns0)){
    tmp <- list(bfp.act[[j]]$per)
    nn <- bfp.act[[j]]$name
    names(tmp) <-nn[3]
    out$BFP[[nn[1]]][[nn[2]]] <- tmp
}

###signal periodograms
bfp.sig <-foreach(j = 1:nrow(ns1)) %dopar% {
    data <- out[[ns1[j,1]]]$RV
    bfp.all <- list()
    for(i in 1:min(Nmax,4)){
        Nma <- Nar <- 0
        GP <- FALSE
        gp.par <- rep(NA,3)
        if(ns1[j,2]=='MA'){
            Nma <- 1
        }
        if(ns1[j,2]=='AR'){
            Nar <- 1
        }
        if(ns1[j,2]=='GP'){
            GP <- TRUE
            if(any('photometry'==names(out$BFP))){
                per <- out$BFP$photometry$photometry$GP$
                gp.par <- c(NA,per$P[1],per$par.opt[1,'logtauGP'])
            }
        }
        cat('calulate BFP for',i,'signal for',ns1[j,2],'for',ns1[j,1],'RV set of',target,'\n')
        per <- BFP(data[,1],data[,2],data[,3],Nma=Nma,Nar=Nar,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=GP,noise.only=FALSE,Nsamp=Nsamp, gp.par=gp.par,renew=TRUE)
        bfp.all[[paste0('sig',i)]] <- per
        data[,2] <- per$res.s
    }
    c(list(bfp=bfp.all),list(name=ns1[j,]))
}

for(j in 1:nrow(ns1)){
                                        #    tmp <- list()
    nn <- bfp.sig[[j]]$name
                                        #    names(tmp) <-nn[2]
    Nsig <- length(bfp.sig[[j]])
    for(k in 1:Nsig){
                                        #        tmp <- list(bfp.sig[[j]]$bfp[[k]])
                                        #        names(tmp) <- nn[2]
        name.new <- paste0('sig',k)
        tmp <- list(bfp.sig[[j]]$bfp[[k]])
        names(tmp) <- nn[2]
        if(!any(name.new==names(out$BFP[[nn[1]]]))){
            out$BFP[[nn[1]]][[name.new]] <- tmp
        }else{
            out$BFP[[nn[1]]][[name.new]][[nn[2]]] <- bfp.sig[[j]]$bfp[[k]]
        }
    }
}
