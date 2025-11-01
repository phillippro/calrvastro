#RV3 <- RV2[ind.ok,k1]
xlim <- range(trv)
####RV data
plot(trv,RV,xlab=expression((t-t[0])*'[days]'),ylab='RV [m/s]',ylim=range(RV+eRV,RV-eRV),main=ids[k1],col=cols[k1])
arrows(trv,RV-eRV,trv,RV+eRV,length=0.05,angle=90,code=3,col=cols[k1])
if(length(epochs)>2){
    abline(v=epochs,col='blue')
}
points(trv,RV2,pch=20,col='grey')
lines(ts,RV.kepler(pars.kep=par.opt.post,tt=ts,kep.only=TRUE)[[1]],col="red")
a <- par.opt.post[grep(paste0('a\\d_',k1),names(par.opt.post))]
b <- par.opt.post[grep(paste0('\\db',k1),names(par.opt.post))]
lines(ts,cal.trend(a=a,b=b,t=(ts-min(ts))/time.unit),col=cols[k1])
#lines(ts,cal.trend(a=a,b=b,t=poly((ts-min(ts))/time.unit,degree=Npoly)),col='green')
legend('topright',legend=paste0('RMS=',format(sd(RV),digit=3),'m/s'),bty='n')
######phase-folded figure
if(Np>0){
data.combine <- FALSE
source('phase_fold.R')
}
###RV periodogram
if(Niter<1e5){
    ofac <- 0.5
}else if((tmax-tmin)<100){
    ofac <- 10
}else{
    ofac <- 1
}
ofac <- 1
fmax <- max(1/pmin,1)
#if(length(trv)>20){
rv.glst <- try(glst(t=trv-min(trv),y=RV,err=eRV,ofac=ofac,fmax=fmax),TRUE)
if(class(rv.glst)!='error'){
plot(rv.glst$P,rv.glst$power,xlab='Period[day]',ylab='Power of RV',log='x',type='l',main='')
#paste0('GLST; p value:',format(rv.glst$pvalue,digit=3))
abline(h=rv.glst$sig.level[1],lty=2)
abline(h=rv.glst$sig.level[2],lty=3)
abline(h=rv.glst$sig.level[3],lty=4)
abline(v=Popt,lty=2)
}
if(k1==1){
#####add universal title
title(main=paste(length(trv),' data points; acceptance rate: ',format(acceptance*100,digit=3),'%; Nsamp=',Niter,'; Pmin=',format(Pmin,digit=3),'; Pmax=',format(Pmax,digit=3),'; tempering: ',format(tem,digit=3),';\n duration=',time.consumed,'; Nsep=',Nsep,'; fchop=',fchop,'; ',Ncores,'cores; max(log(post))=',format(max(post.out),digits=3),'; max(log(like))=',format(max(loglike.out),digits=3),sep=""),outer=TRUE,line=-3,cex=1.5)
}
tit <- ''
###residuals
plot(trv,res,xlab=expression((t-t[0])*'[days]'),ylab='Residuals [m/s]',main='')
legend('topright',legend=paste0('RMS=',format(sd(res),digit=3),'m/s'),bty='n')
###residual periodogram
#if(length(trv)>20){
res.glst <- try(glst(t=trv-min(trv),y=res,err=eRV,ofac=ofac,fmax=fmax),TRUE)
if(class(res.glst)!='error'){
plot(res.glst$P,res.glst$power,xlab='Period[day]',ylab='Power of Residuals',log='x',type='l',main='')
ind <- which.max(res.glst$power)
if(length(ind)>0){
    pmax <- res.glst$P[ind]
    power.max <- res.glst$power[ind]
    par(xpd=TRUE)
    text(pmax,1.1*power.max,pos=3,labels=format(pmax,digit=3),col='red',cex=1.0)
    arrows(pmax,1.1*power.max,pmax,1.05*power.max,col='red',length=0.05)
    par(xpd=FALSE)
}
abline(h=res.glst$sig.level[1],lty=2)
abline(h=res.glst$sig.level[2],lty=3)
abline(h=res.glst$sig.level[3],lty=4)
}
######differential radial velocity and periodograms
if(!all(Inds==0)){
    for(j1 in 1:ncol(proxies)){
        if(sd(eproxies[,j1])>0){
            ylim <- range(proxies[,j1]+eproxies[,j1],proxies[,j1]-eproxies[,j1])
        }else{
            ylim <- range(proxies[,j1])
        }
        plot(trv,proxies[,j1],xlab=expression((t-t[0])*'[days]'),ylab=paste0(proxy.names[j1],'[m/s]'),ylim=ylim,main=proxy.names[j1])
        if(sd(eproxies[,j1])>0){
            arrows(trv,proxies[,j1]-eproxies[,j1],trv,proxies[,j1]+eproxies[,j1],length=0.05,angle=90,code=3)
        }
        legend('topright',legend=paste0('RMS=',format(sd(proxies[,j1]),digit=3),'m/s'),bty='n')
###residual periodogram
        if(length(trv)>20){
            if(sd(eproxies[,j1])>0){
                proxy.glst <- glst(t=trv-min(trv),y=proxies[,j1],err=eproxies[,j1],ofac=ofac,fmax=fmax)
            }else{
                proxy.glst <- lsp(times=trv-min(trv),x=proxies[,j1],ofac=ofac,to=fmax)
            }
            plot(proxy.glst$P,proxy.glst$power,xlab='Period[day]',ylab=paste0('Power of ',proxy.names[j1]),log='x',type='l',main=proxy.names[j1])
            abline(v=Popt,lty=2)
        }
    }
}
#win.bgls <- bgls(t=trv-min(trv),y=rep(2,length(trv)),err=rnorm(length(trv),1,0.01), ofac=ofac,fmax=fmax)
win.bgls <- bgls(t=trv-min(trv),y=rep(2,length(trv)),err=rep(1,length(trv)), ofac=ofac,fmax=fmax)
#win.per <- lsp(times=trv-min(trv),x=rep(2,length(trv)), ofac=ofac,to=fmax)
plot(win.bgls$P,win.bgls$power,xlab='Period[day]',ylab='Power of window function',log='x',type='l',main=tit)
abline(v=Popt,lty=2)
