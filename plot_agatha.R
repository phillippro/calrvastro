pdf(fname,12,12)
#pdf(fname,8,8)
par(mfrow=c(3,3))
for(j in 1:(length(noise.types)*Nmax)){
    noise.type <- noise.types[ceiling(j/Nmax)]
#    if(!grepl('GP',noise.type)){
#        ylim <- range(median(out[,j+1]),max(5,out[,j+1]))
#        ylim <- range(min(out[,j+1]),max(5,out[,j+1]))
#    }else{
#    }
    ylim <- range(out[,j+1])
    plot(out[,1],out[,j+1],xlab='P [d]',ylab='ln(BF)',main=paste0('RV-',(j-1)%%Nmax,'signal; ',noise.type,'; p=',Nars[ceiling(j/Nmax)],';q=',Nmas[ceiling(j/Nmax)],';dt=',round(durs[j]),'s'),log='x',type='l',ylim=ylim,xaxt='n')
    magaxis(side=1,tcl=-0.5)
    abline(h=c(0,5),lty=3:2,col='grey')
    if(noise.only & noise.type!='white'){
        abline(h=-(Nma+Nar+GP*length(which(is.na(gp.par[-2]))))/2*log(nrow(tab)),lty=2,col='red')
    }
    Pmax <- out[which.max(out[,j+1]),1]
    abline(v=Pmax,col='red',lty=3,lwd=3)
    legend('topright',legend=paste0('P=',round(Pmax,2),'d'),bty='n',col='red')
}
dev.off()