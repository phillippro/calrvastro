show.table <- function(table,alpha=0.3,b=0){
    fonts <- rep(1,length(table))
    plot.new()
    cn <- c('',sapply(1:length(table),function(i) table[[i]]['name']),'')
    rn1 <- names(table[[1]][-1])
    rn2 <- names(table[[length(table)]][-1])
#    legend('topleft',inset=c(0,b-0.1),legend=rn,bty='n',xpd=NA,horiz=TRUE)
#    legend('topleft',inset=c(0,b),legend=cn,bty='n',xpd=NA)
    legend('topleft',inset=c(0,b),legend=c('',rn1),bty='n',xpd=NA,text.font=2)
    for(j in 1:length(table)){
        if(length(table[[j]])>0){
            legend('topleft',inset=c(alpha*j,b),legend=table[[j]],bty='n',xpd=NA,text.font=c(2,rep(fonts[j],length(table[[j]])-1)))
        }
    }
    legend('topleft',inset=c((length(table)+1)*alpha,b),legend=c('',rn2),bty='n',xpd=NA,text.font=2)
}

###activity signals
Psig <- Pact <- c()
sig.name <- act.name <- c()
if(out$Nrv>0){
    for(i in out$ins){
        p <- names(out$BFP[[i]])
        ind0 <- grep('sig',p)
        ind1 <- which(!grepl('sig',p))
        for(j in p[ind0]){
            q <- names(out$BFP[[i]][[j]])
            for(k in q){
                per <- out$BFP[[i]][[j]][[k]]
                if(!is.null(per)){
                    if(!all(is.na(per))){
                        inds <- sort(per$power,decreasing=TRUE,index.return=TRUE)$ix[1:2]
                        pp <- as.numeric(per$P[inds])
                        if(is.null(pp)){
                            pp <- rep(NA,2)
                        }
                        Psig <- rbind(Psig,pp)
                        sig.name <- c(sig.name,paste0(i,'_',j,'_',k))
                    }
                }
            }
        }
        for(j in p[ind1]){
            q <- names(out$BFP[[i]][[j]])
            for(k in q){
                per <- out$BFP[[i]][[j]][[k]]
                if(!is.null(per)){
                    inds <- sort(per$power,decreasing=TRUE,index.return=TRUE)$ix[1:2]
                    Pact <- rbind(Pact,as.numeric(per$P[inds]))
                    act.name <- c(act.name,paste0(i,'_',j,'_',k))
                }
            }
        }
    }
    colnames(Psig) <- colnames(Pact) <- c('peak1','peak2')
    rownames(Pact) <- act.name
    rownames(Psig) <- sig.name
}
lnbf5 <- lnbf3 <- Nsets <- Nmodel <- Noverlap <- quality <- 0
###diagnostics;
if(any(grepl('sig0',names(out$mcmc.opt)))){
    mcmc0 <- out$mcmc.opt[['sig0']]
    l0 <- max(mcmc0[,'loglike'])
}else{
    l0 <- -1e4
}
if(Nmin==Nmax){
    inds <- Nmin
}else{
    inds <- 1:Nsig
}
for(i in inds){
    score <- 0
    score.sub <- 0
###criterion 1: statistical significant
    nn <- paste0('sig',i)
    if(any(names(out$mcmc.opt)==nn)){
        mcmc <- out$mcmc.opt[[nn]]
        l1 <- max(mcmc[,'loglike'])
    }else{
        l1 <- -1e4
    }
    bf3 <- l1-l0-1.5*log(length(trv.all))
    bf5 <- l1-l0-2.5*log(length(trv.all))
    l0 <- l1
    lnbf3 <- c(lnbf3,bf3)
    lnbf5 <- c(lnbf5,bf5)
    if(bf5>=5)  score <- score+1

###mcmc signals overlap with BFP signals
    if(out$Nrv>0){
        ind0 <- grep(paste0('sig',i),rownames(Psig))
        pp <- Psig[ind0,]
        ind00 <- which(abs(pp-Popt[i])<0.2*Popt[i],arr.ind=TRUE)
        n0 <- nrow(ind00)
        set.noise <- rownames(ind00)
        xi <- sapply(out$ins,function(i) sum(1/out[[i]]$RV[,3]^2))#fitting ability
        ind <- grep(out$ins[which.max(xi)],set.noise)

####criterion 2: check noise model dependency
        Nmodel <- c(Nmodel,length(ind))
        if(length(ind)==length(noise.types)) score <- score+1

####criterion 3: check data set dependency
        ind <- sapply(out$ins,function(i) length(grep(i,set.noise)))
        ns <- length(which(ind>0))
        Nsets <- c(Nsets,ns)
        if(ns>=max(2,ceiling(length(out$ins)/2))) score <- score+1

###criterion 4: activity signal dependency
        ind1 <- which(abs(Pact-Popt[i])<0.1*Popt[i])
        n1 <- length(ind1)
        Noverlap <- c(Noverlap,n1)
        if(n1==0)   score <- score+1

###criterion 5: time dependency
        if(!is.na(gammas[i])){
            if(gammas[i]<0.2){
                score.sub <- score+1
            }else{
                score.sub <- score-1
            }
        }

###criterion 6: wavelength dependency
        if(!all(is.na(etas[,i]))){
            if(min(etas[,i]<0.2)){
                score.sub <- score+1
            }else{
                score.sub <- score-1
            }
        }

###overall quality assessment
        if(score==4) quality[i] <- 'A'
        if(score==3) quality[i] <- 'B'
        if(score==2) quality[i] <- 'C'
        if(score==1) quality[i] <- 'D'
        if(score==0) quality[i] <- 'E'
        if(score.sub==2) quality[i] <- paste0(quality[i],'++')
        if(score.sub==1) quality[i] <- paste0(quality[i],'+')
        if(score.sub==-1) quality[i] <- paste0(quality[i],'-')
        if(score.sub==-2) quality[i] <- paste0(quality[i],'--')
    }
}
if(!exists('gammas')) gammas <- rep(NA,Nsig)
if(!exists('etas')) etas <- array(NA,dim=c(10,Nsig))

#####show derived parameters of each fit
par.out <- list()
par.conv <- extract.par(par.opt,bases=bases)

for(j in 1:(Nsig+1)){
    if(j<=Nsig){
        msini <- K2msini.full(par.conv$K[j],par.conv$P[j],par.conv$e[j],Ms=out$Mstar)
        ni <- paste0('Inc',j)
        if(any(names(par.opt)==ni)){
            Mp <- msini$mj/sin(par.opt[ni])
            par.show <- c(Mp,par.conv$P[j],par.conv$K[j],par.conv$e[j],par.conv$omega[j],par.conv$Mo[j],par.conv$Inc[j],par.conv$Omega[j])
            ns <- c('name','mass[Mjup]','Period[day]','K[m/s]','e','omega[rad]','Mo[rad]','I[rad]','Omega[rad]','lnBF3','lnBF5','Nset','Nmodel','Noverlap','gamma','eta','quality')
        }else{
            Mp <- msini$mj
            par.show <- c(Mp,par.conv$P[j],par.conv$K[j],par.conv$e[j],par.conv$omega[j],par.conv$Mo[j])
            ns <- c('name','mass[Mjup]','Period[day]','K[m/s]','e','omega[rad]','Mo[rad]','lnBF3','lnBF5','Nset','Nmodel','Noverlap','gamma','eta','quality')
        }
        tmp <- c(paste(j,'signal'),round(par.show,3),round(lnbf3[j],1),round(lnbf5[j],3),Nset[j],Nmodel[j],Noverlap[j],round(gammas[j],3),round(min(etas[,j]),3),quality[j])
        names(tmp) <- ns
        par.out[[j]] <- tmp
    }else{
        par.out[[j]] <- c('noise',round(par.opt[(Nkeppar*Nsig+1):length(par.opt)],3))
    }
}
#fonts <- rep(1,ncol(par.out))
#inds <- ncol(par.shows)+Ntype*(0:5)+ind.show+1
#fonts[inds] <- 2
#par(mfrow=c(4,4),mar=c(5,5,5,1),cex.axis=1.5,cex.lab=1.5)
par(mar=c(1,1,1,1),mfrow=c(4,4))
show.table(par.out,alpha=0.2)
mtext("Model Parameters and Quality Flags for Signals", outer=TRUE,  cex=1.5, line=-0.5)
#fonts=fonts
