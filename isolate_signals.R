#plows0 <- plows
#pups0 <- pups
#pups<- pups0
cat('old plows=',plows,'\n')
cat('old pups=',pups,'\n')
#plow.new <- min(max(P.bin[,ncol(P.bin)],1.03*Popt),max(pups))
#pup.new <- max(min(P.bin[,ncol(P.bin)],0.97*Popt),min(plows))
for(i3 in 1:ncol(P.bin)){
    if(any(plows==pups)){
        ind <- which(plows==pups)
        plows <- plows[-ind]
        pups <- pups[-ind]
    }
    if(all(min(P.bin[,i3])>pups) | all(max(P.bin[,i3])<plows)){
        plow.new <- pup.new <- c()
    }else{
        plow.new <- min(max(P.bin[,i3]),max(pups))
        pup.new <- max(min(P.bin[,i3]),min(plows))
    }
    cat('pup.new=',pup.new,'\n')
    cat('plow.new=',plow.new,'\n\n')
    if(identical(plows,plow.new) & identical(pups,pup.new)){
    if(length(plows)>1 & length(pup.new)>0){
    pf.low <- pups[-length(pups)]
    pf.up <- plows[-1]
    ind.interval <- which(pup.new>plows & pup.new<pups & plow.new>plows & plow.new<pups)
    if(length(ind.interval)==0){
        ind.up <- which(pup.new>=plows & pup.new<=pups)
        ind.low <- which(plow.new>=plows & plow.new<=pups)
        ind.fup <- which(pup.new>=pf.low & pup.new<=pf.up)
        ind.flow <- which(plow.new>=pf.low & plow.new<=pf.up)
        if(length(ind.low)>0 & length(ind.flow)>0){
            ind.low <- c()
        }
        if(length(ind.up)>0 & length(ind.fup)>0){
            ind.up <- c()
        }
        if(length(ind.up)>0){
            pups.close <- pups[sort(abs(pups-pup.new),index.return=TRUE)$ix[1:2]]
            pup.close <- pups.close[pups.close>pup.new][1]
            ind <- which(pups==pup.close)
            pups[ind] <- pup.new
            if(length(ind.low)>0){
                plows.close <- plows[sort(abs(plows-plow.new),index.return=TRUE)$ix[1:2]]
                plow.close <- plows.close[plows.close<plow.new][1]
                ind <- which(plows==plow.close)
                plows[ind] <- plow.new
                if(ind.low!=ind.up){
                    plow1 <- plow.new
                    pup1 <- pup.new
                    indc <- which(pups>pup1 & pups<plow1)
                    if(length(indc)>0){
                        pups <- pups[-indc]
                        plows <- plows[-which(plows>pup1 & plows<plow1)]
                    }
                }
            }
        }else if(length(ind.low)==0){
            if(ind.fup!=ind.flow){
                pup1 <- pf.low[ind.fup]
                plow1 <- pf.up[ind.flow]
                indc <- which(pups>pup1 & pups<plow1)
                if(length(indc)>0){
                    pups <- pups[-indc]
                    plows <- plows[-which(plows>pup1 & plows<plow1)]
                }
            }
        }else if(length(ind.low)==1){
            plows.close <- plows[sort(abs(plows-plow.new),index.return=TRUE)$ix[1:2]]
            plow.close <- plows.close[plows.close<plow.new][1]
            ind <- which(plows==plow.close)
            plows[ind] <- plow.new
            if(ind.low!=(ind.fup+1)){
                pup1 <- pf.low[ind.fup]
                plow1 <- plow.new
                indc <- which(pups>pup1 & pups<plow1)
                if(length(indc)>0){
                    pups <- pups[-indc]
                    plows <- plows[-which(plows>pup1 & plows<plow1)]
                }
            }
        }
    }else{
        if(!any(plow.new==plows)){
            plows <- sort(c(plows,plow.new))
        }
        if(!any(pup.new==pups)){
            pups <- sort(c(pups,pup.new))
        }
    }
}else{
    if(!any(plow.new==plows)){
        plows <- sort(c(plows,plow.new))
    }
    if(!any(pup.new==pups)){
        pups <- sort(c(pups,pup.new))
    }
}
ind1 <- c()
for(j1 in 1:length(plows)){
    if(plows[j1]==pups[j1]) ind1 <- c(ind1,j1)
}
if(length(ind1)>0){
    plows <- plows[-ind1]
    pups <- pups[-ind1]
}
}
}
if(any(plows==pups)){
    ind <- which(plows==pups)
    plows <- plows[-ind]
    pups <- pups[-ind]
}
cat('plows=',plows,'\n')
cat('pups=',pups,'\n')
#if(length(plows)==0 | length(pups)==0){
#   break()
#}
