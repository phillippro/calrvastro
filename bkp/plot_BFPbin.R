if(period.par=='P'){
    Ps <- out.all[,1]
}else if(period.par=='logP'){
    Ps <- exp(out.all[,1])
}
ll <- out.all[,ncol(out.all)]
###binning
Nbin <-1000
logPbin <- seq(log(min(Ps)),log(max(Ps)),length.out=Nbin+1)
Pmids <-exp((logPbin[-1]+logPbin[-length(Pbin)])/2)
lls <- c()
for(j in 1:(length(Pbin)-1)){
    inds <-which(Ps>=Pbin[j] & Ps<Pbin[j+1])
    if(length(inds)==0){
        ll.max <-NA
    }else{
        ll.max <- max(out.all[inds,ncol(out.all)])
    }
    lls <- c(lls,ll.max)
}
logbf <- lls-ll0-Nkeppar/2*log(nrow(data))
plot(Pmids,logbf,xlab='Period [day]',ylab='ln(BF)',type='l')#,ylim=range(median(logbf),max(logbf)),log='x')
abline(h=c(0,5),lty=3:2,col='grey')
inds <-sort(logbf,index.return=TRUE,decreasing=TRUE)$ix
Popts <- Pmids[inds]
abline(v=Popts[1:3],col=tcol('red',70),lty=1,lwd=3)
legend('topright',legend=paste0('P=',round(Popt[1:3],2),'d'),bty='n',col='red')
