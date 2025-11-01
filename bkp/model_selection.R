##model.out <- list()
#target.except <- read.table('target_noMA.txt')[,1]
target.except <- ''
model.out <- foreach(i = out$ins.rv) %dopar% {
####without including differential RVs
    cat(i,'RV set!\n')
#    if(nrow(out[[i]]$RV)<10 | target=='HD56380' | target=='GJ288'){
    if(nrow(out[[i]]$RV)<10 | any(target==target.except)){
#    if(TRUE){
        list(mc=list(nqp=rep(0,3)),instr=i)
    }else{
        ns <- names(out[[i]][['activity']])
        ind <- grep('AP',ns)
###sort dRVs in order of decreasing correlation with RVs
        if(length(ind)>0){
            rs <- c()
            for(k3 in ind){
                x <- out[[i]]$RV
                y <- out[[i]][['activity']][[ns[k3]]]
                if(nrow(x)!=nrow(y)){
                    pp <- x
                    ind1 <- match(x[,1]%%2400000,y[,1]%%2400000)
                    indx.na <- which(is.na(ind1))
                    indx <- which(!is.na(ind1))
                    indy <- ind1[!is.na(ind1)]
                    pp[indx,] <- y[indy,]
                    pp[indx.na,2] <- mean(y[indy,2])
                    pp[indx.na,3] <- mean(y[indy,3])
                    y <- pp
                    out[[i]][['activity']][[ns[k3]]] <- y
                }
                rs <- c(rs,cor(x[,2],y[,2]))
            }
            ii <- sort(rs,decreasing=TRUE,index.return=TRUE)$ix
            ind <- ind[ii]
            ##        ns[ind[ii][1:min(5,length(rs))]]
            proxy <- sapply(ind,function(j) out[[i]]$activity[[ns[j]]][,2])
            name.proxy <- sapply(ind,function(j) names(out[[i]]$activity)[j])
            colnames(proxy) <- name.proxy
            proxy <- proxy[,1:min(ncol(proxy),5),drop=FALSE]
            NI.inds <- list(0)
            for(j in 1:ncol(proxy)){
                NI.inds[[j+1]] <- 1:j
            }
            ##        mc <- bfp.inf.norm(data=cbind(out[[i]]$RV,proxy),Nmas=0:5,Nars=0,NI.inds=NI.inds)
            mc <- bfp.inf.combined(data=cbind(out[[i]]$RV,proxy),Nmas=0:5,Nars=0,NI.inds=NI.inds,Nrep=10)
        }else{
            proxy <- NULL
            mc <- try(bfp.inf.combined(data=out[[i]]$RV,Nmas=0:5,Nars=0,NI.inds=list(0),Nrep=10),TRUE)
        }
##        model.out[[i]] <- list(mc=mc,instr=i)
        list(mc=mc,instr=i)
    }
}
for(j in 1:length(model.out)){
    out[[model.out[[j]]$instr]][['noise']] <- model.out[[j]]$mc
}

nqp <- sapply(out$ins.rv, function(i) out[[i]]$noise$nqp)
noise.models <- sapply(1:ncol(nqp),function(i) paste0('AR',nqp[3,i],'MA',nqp[2,i]))
###
#mc <- model.infer(data=out[[i]]$RV,proxy=proxy,Nma.max=5,Nar.max=0,Nrep=5,GP=FALSE)
##    noise.models <- c(noise.models,paste0('ARMA',mc$Nar,mc$Nma,'ind',mc$Nind))
#out[[i]][['noise']] <- mc
