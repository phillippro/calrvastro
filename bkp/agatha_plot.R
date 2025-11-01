ylim <- range(out[,ncol(out)])
if(grepl('GP',noise.type)){
    mm <- paste0('RV-',k-1,'signal; ',noise.type,'; Prot=',round(exp(logProt)),'d;tau=',round(exp(logtauGP)),'d;dt=',round(durs[length(durs)]),'s')
}else if(grepl('AR|MA',noise.type)){
    mm <- paste0('RV-',k-1,'signal; ',noise.type,'; p',Nar,'q',Nma,';dt=',round(durs[length(durs)]),'s')
}else{
    mm <- paste0('RV-',k-1,'signal ',noise.type,'; dt=',round(durs[length(durs)]),'s')
}
plot(out[,1],out[,ncol(out)],xlab='P [d]',ylab='ln(BF)',main=mm,log='x',type='l',ylim=ylim)
#magaxis(side=1,tcl=-0.5)
abline(h=c(0,5),lty=3:2,col='grey')
if(noise.only & noise.type!='white'){
    abline(h=-(Nma+Nar+GP*length(which(is.na(gp.par[-2]))))/2*log(nrow(tab)),lty=2,col='red')
}
Po <- out[which.max(out[,ncol(out)]),1]
abline(v=Po,col='red',lty=3,lwd=3)
legend('topright',legend=paste0('P=',round(Po,2),'d'),bty='n',col='red')
