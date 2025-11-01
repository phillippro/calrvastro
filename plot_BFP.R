#plotbfp <- TRUE
plotbfp <- FALSE
if(plotbfp){
#if(TRUE){
pdf('test.pdf',16,16)
par(mfrow=c(4,4))
}
####activity
Nmin0 <- 1
Popt0 <- Popt
if(target=='HD86226') Popt <- c(Popt,6.4)
par(mfrow=c(4,4))
panel <- 0
inss <- c(out$ins.rv,'all')
for(i in inss){
    ns0 <- names(out$BFP[[i]])
    ind <- which(!grepl('sig',ns0))
    for(j in ind){
        ns <- names(out$BFP[[i]][[j]])
        for(n in ns){
            bfpf <- FALSE
            if(!is.null(unlist(out$BFP[[i]][[j]]))){
                if(any(names(out$BFP[[i]][[j]])==n)){
                    if(!is.null(out$BFP[[i]][[j]][[n]])) bfpf <- TRUE
                }
            }
            if(bfpf){
                bfp <-out$BFP[[i]][[j]][[n]]
                lnbf <- bfp$power
                plot(bfp$P,lnbf,xlab='Period [day]',ylab='ln(BF)',type='l',ylim=range(min(median(lnbf),0),max(max(lnbf),5)),log='x',main=paste0(i,';',ns0[j],';',n),xaxt='n')
                magaxis(side=1,tcl=-0.5)
                if(target=='HIP47103') abline(v=c(75,78),col='blue',lty=3,lwd=2)
                panel <- panel+1
                if(panel==1){
                    mtext("BFPs for Noise Proxies", outer=TRUE,  cex=1.5, line=-0.5)
                }
                abline(h=c(0,5),lty=3:2,col='darkgrey')
                abline(v=Popt,col=tcol('red',70),lty=1,lwd=2)
                lines(bfp$P,lnbf)
            }else{
	        plot(0,pch='',xlab='',ylab='',xaxt='n',yaxt='n')
                legend('center',legend=paste0('No plot for ',i,';',ns0[j],';',n,'!\n'))
	    }
        }
    }
}
if(any(grepl('photometry',names(out$BFP)))){
    photo <- out$BFP$photometry
    sets <- names(photo)
    for(set in sets){
        noises <- names(photo[[set]])
        for(noise in noises){
            bfpf <- FALSE
            if(any(names(photo[[set]])==noise)){
                if(!is.null(unlist(photo[[set]][[noise]]))) bfpf <- TRUE
            }
            if(bfpf){
                lnbf <- photo[[set]][[noise]]$power
                plot(photo[[set]][[noise]]$P,lnbf,xlab='Period [day]',ylab='ln(BF)',type='l',ylim=range(min(median(lnbf),0),max(max(lnbf),5)),log='x',main=paste0('photometry;',set,';',noise),xaxt='n')
                magaxis(side=1,tcl=-0.5)
                if(target=='HIP47103') abline(v=c(75,78),col='blue',lty=3,lwd=2)
                abline(h=c(0,5),lty=3:2,col='darkgrey')
                abline(v=Popt,col=tcol('red',70),lty=1,lwd=2)
                lines(photo[[set]][[noise]]$P,photo[[set]][[noise]]$power)
            }else{
                plot(0,pch='',xlab='',ylab='',xaxt='n',yaxt='n')
                legend('center',legend=paste0('No plot for photometry;',set,';',noise))
            }
        }
    }
}

par(mfrow=c(length(noise.types),length(noise.types)))
panel <- 0
Nbfp <- length(grep('sig',names(out$BFP[[1]])))
inss <- c(out$ins,'all')
for(i in inss){
    sigs <- names(out$BFP[[i]])
    sigs <- sigs[grepl('sig',sigs)]
#    for(j in 1:(Nbfp-Nmin0)){
    for(sig in sigs){
        for(noise in noise.types){
	if(any(names(out$BFP[[i]][[sig]])==noise)){
            bfp <- out$BFP[[i]][[sig]][[noise]]
            bfpf <- FALSE
            if(!is.null(unlist(bfp))){
                if(!is.null(names(bfp))){
                    if(any(names(bfp)=='power')) bfpf <- TRUE
                }
            }
            if(bfpf){
                lnbf <- bfp$power
                plot(bfp$P,lnbf,xlab='Period [day]',ylab='ln(BF)',type='l',ylim=range(min(median(lnbf),0),max(max(lnbf),5)),log='x',main=paste0(i,';',noise,';',sig),xaxt='n')
            magaxis(side=1,tcl=-0.5)
            if(target=='HIP47103') abline(v=c(75,78),col='blue',lty=3,lwd=2)
                panel <- panel+1
                if(panel==1){
                    mtext("BFPs for RV Data Sets", outer=TRUE,  cex=1.5, line=-0.5)
                }
                abline(h=c(0,5),lty=3:2,col='darkgrey')
                abline(v=Popt,col=tcol('red',70),lty=1,lwd=2)
                lines(bfp$P,lnbf)
            }else{
                plot(0,pch='',xlab='',ylab='',xaxt='n',yaxt='n')
                legend('center',legend=paste0('No plot for ',i,';',noise,';sig'))
            }
        }
        }
    }
}
Popt <- Popt0
if(plotbfp){
    dev.off()
}
