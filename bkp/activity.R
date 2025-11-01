if(periodogram.type=='BFP'){
    per <- BFP(acts[,1],acts[,2],acts[,3],Nma=1,Nar=0,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=FALSE,noise.only=FALSE,Nsamp=Nsamp)
}else if(periodogram.type=='LS'){
    per <- lsp(times=acts[,1],x=acts[,2],ofac=ofac,from=fmin,to=fmax,alpha=c(0.1,0.01,0.001))
}
if(plotf){
    plot(per$P,per$power,xlab='P [d]',ylab='ln(BF)',main=act.name,log='x',type='l')
#    magaxis(side=1,tcl=-0.5)
    abline(h=c(0,5),lty=3:2,col='grey')
    Popt <- per$P[which.max(per$power)]
    abline(v=Popt,col='red',lty=3,lwd=3)
    legend('topright',legend=paste0('P=',round(Popt,2),'d'),bty='n',col='red')
}
