###list the names of all periodograms
#Nmin0 <- 1
Nmin0 <- Nmin
ns0 <- c()
for(i in out$ins.rv){
    acts <- names(out[[i]][['activity']])
    for(act in acts){
        if(grepl('window',act)){
            ns0 <- rbind(ns0,c(i,act,'white'))
        }else if(!grepl('AP',act)){
            ns0 <- rbind(ns0,c(i,act,'MA'))
        }
    }
}
ns0 <- rbind(ns0,c('all','window','white'))

if(any(grepl('photometry',names(out)))){
    ns0 <- rbind(ns0,cbind('photometry',names(out$photometry),'MA'))
#    ns0 <- rbind(ns0,c('photometry','photometry','GP'))
}
colnames(ns0) <- NULL

ns1 <- c()
for(j in 1:(Nset+1)){
    for(k in 1:length(noise.types)){
        if(j<=Nset){
            ns1 <- rbind(ns1,c(out$ins.rv[j],noise.types[k]))
        }else{
            ns1 <- rbind(ns1,c('all',noise.types[k]))
        }
    }
}

###activity periodograms
bfp.act <- foreach(j = 1:nrow(ns0)) %dopar% {
#for(j in 1:nrow(ns0)){
    if(ns0[j,1]=='photometry'){
        if(grepl('GP',ns0[j,3])){
            ind <- grep('bin',names(out[['photometry']]))
            if(length(ind)==0){
                data <- out[['photometry']][[ns0[j,2]]]
            }else{
                data <- out[['photometry']][[ind[1]]]
            }
        }else{
            data <- out[['photometry']][[ns0[j,2]]]
        }
    }else if(ns0[j,1]=='all'){
        data <- out$all
        data[,2] <- 1
    }else{
        data <- out[[ns0[j,1]]][['activity']][[ns0[j,2]]]
    }
    if(nrow(data)>2){
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
                                        #        noise.only <- TRUE
            ofac1 <- min(1,ofac)
        }
        cat('calulate BFP for',ns0[j,2],'for',ns0[j,3],'for',ns0[j,1],'set of',target,'\n')
        if(ns0[j,2]=='window'){
            per <- try(lsp(times=data[,1],x=data[,2],ofac=ofac,from=fmin,to=fmax,alpha=c(0.1,0.01,0.001)),TRUE)
        }else{
            per <- try(BFP(data[,1],data[,2],data[,3],Nma=0,Nar=0,Indices=NULL,ofac=ofac1,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=GP,noise.only=noise.only,Nsamp=Nsamp,gp.par=gp.par,renew=TRUE),TRUE)
                                        #            per <- BFP(data[,1],data[,2],data[,3],Nma=0,Nar=0,Indices=NULL,ofac=ofac1,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=GP,noise.only=noise.only,Nsamp=Nsamp,gp.par=gp.par,renew=TRUE)
            if(class(per)=='try-error'){
                per <- try(glst(t=data[,1],y=data[,2],err=data[,3],ofac=ofac1,fmin=fmin,fmax=fmax),TRUE)
            }
        }
        if(class(per)=='try-error'){
            c(list(per=NULL),list(name=ns0[j,]))
        }else{
            c(list(per=per),list(name=ns0[j,]))
        }
    }else{
#        per <- glst(t=data[,1],y=data[,2],err=data[,3],ofac=ofac,fmin=fmin,fmax=fmax)
        c(list(per=NULL),list(name=ns0[j,]))
#        c(list(per=per),list(name=ns0[j,]))
    }
}

out$BFP <- list()
for(j in 1:length(bfp.act)){
#    if(nn[2]==nn[1]){
#        out$BFP[[nn[1]]][[nn[3]]] <- bfp.act[[j]]$per
#    }else{
#    if(!any(grepl(nn[1],names(out$BFP)))){
#        out$BFP[[nn[1]]] <- list()
#    }
    tmp <- list(bfp.act[[j]]$per)
    nn <- bfp.act[[j]]$name
    names(tmp) <-nn[3]
    if(nn[1]!=nn[2]){
    out$BFP[[nn[1]]][[nn[2]]] <- tmp
}else{
    out$BFP[[nn[1]]][[nn[3]]] <- bfp.act[[j]]$per
}
#    }
}

###signal periodograms
Nbfp <- max(min(Nsig+1,Nmax),Nmin0+1)
bfp.sig <-foreach(j = 1:nrow(ns1)) %dopar% {
#for(j in 1:nrow(ns1)){
    if(ns1[j,1]=='all'){
        data <- out$all
    }else{
        data <- out[[ns1[j,1]]]$RV
    }
    if(nrow(data)>1){
        bfp.all <- list()
        for(i in (Nmin0+1):Nbfp){
#        for(i in (1+1):Nbfp){
#            cat('i=',i,'\n')
            if(ns1[j,1]=='all'){
#                if(i==1){
                    data[,2] <- out$res.comb.sigtrend[[paste0('sig',i-1)]]
#                    data[,2] <- out$res.comb.all[[paste0('sig',i-1)]]
#                }else{
#                    data[,2] <- out$res.comb.all[[paste0('sig',i-1)]]
#                }
            }else{
                data[,2] <- out$res.sig[[paste0('sig',i-1)]][[ns1[j,1]]]
            }

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
                    per <- out$BFP$photometry$photometry$MA
                                        #                gp.par <- c(NA,per$Popt[1],per$par.opt[1,'logtauGP'])
                    gp.par <- c(NA,per$Popt[1],log(100))
                }
            }
            cat('calulate BFP for',i,'signal for',ns1[j,2],'for',ns1[j,1],'RV set of',target,'\n')
            per <- try(BFP(data[,1],data[,2],data[,3],Nma=Nma,Nar=Nar,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=GP,noise.only=FALSE,Nsamp=Nsamp, gp.par=gp.par,renew=TRUE),TRUE)
            if(class(per)=='try-error'){
                bfp.all[[paste0('sig',i)]] <- NA
            }else{
                bfp.all[[paste0('sig',i)]] <- per
            }
        }
        c(list(bfp=bfp.all),list(name=ns1[j,]))
    }
}

for(j in 1:length(bfp.sig)){
    nn <- bfp.sig[[j]]$name
    sigs <- names(bfp.sig[[j]]$bfp)
    for(sig in sigs){
        if(length(bfp.sig[[j]]$bfp[[sig]])>0 & grepl('sig',sig)){
            tmp <- try(list(bfp.sig[[j]]$bfp[[sig]]),TRUE)
            if(class(tmp)!='try-error'){
                names(tmp) <- nn[2]
                if(!any(sig==names(out$BFP[[nn[1]]]))){
                    out$BFP[[nn[1]]][[sig]] <- tmp
                }else{
                    out$BFP[[nn[1]]][[sig]][[nn[2]]] <- bfp.sig[[j]]$bfp[[sig]]
                }
            }else{
                out$BFP[[nn[1]]][[sig]] <- NA
            }
        }else{
            out$BFP[[nn[1]]][[sig]] <- NA
        }
    }
}
