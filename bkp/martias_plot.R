load('results/HD95338/HD95338_natural_Nmax6_Niter1100000_Ncores8_ofac2_Nset2_Esd0.2_transit0_P55_acc1.5_lnlmax-199.Robj')
colors <- c('black','blue','orange','cyan','brown','green','steelblue','purple','pink',rainbow(10))
if(!exists('bases')) bases <- rep('natural',10)

pdf('matias_MA.pdf',16,16)
par(mfrow=c(4,4))
##without MA
par.opt0 <- par.opt
par.opt['w11'] <- 0
source('phase_plot.R',local=TRUE)
par.opt <- par.opt0

load('results/HD95338/HD95338_natural_Nmax6_Niter1100000_Ncores8_ofac2_Nset2_Esd0.2_transit0_P55_acc1.5_lnlmax-199.Robj')
##with MA(1)
source('phase_plot.R',local=TRUE)

##MA(1) component
for(i in ins[1]){
    tt <- out[[i]]$RV[,1]
    rv <- rv.arma[[i]]
    plot(tt,rv,xlab='JD',ylab='RV',main=paste('MA for',i,';RMS=',round(sd(rv),2),'m/s'))
    plot(tt,out$res.all$sig1[[i]]+rv.arma[[i]],xlab='JD',ylab='RV',main=paste('residual for',i),ylim=1.1*range(out$res.all$sig1[[i]],out$res.all$sig1[[i]]+rv.arma[[i]]))
    points(tt,out$res.all$sig1[[i]],col='red')
    legend('topleft',bty='n',legend=c(paste0('RMS(MA(0))=',round(sd(out$res.all$sig1[[i]]+rv.arma[[i]]),2),'m/s'),paste0('RMS(MA(1))=',round(sd(out$res.all$sig1[[i]]),2),'m/s')),text.col=c('black','red'))
##periodogram difference
    plot(out$BFP$PFS$sig1$white$P,out$BFP$PFS$sig1$white$logBF,xlab='Period [day]',ylab='ln(BF)',main='MA(0)',log='x',type='l')
    abline(v=Popt,col=tcol('red',50))
    plot(out$BFP$PFS$sig1$MA$P,out$BFP$PFS$sig1$MA$logBF,xlab='Period [day]',ylab='ln(BF)',main='MA(1)',log='x',type='l')
    abline(v=Popt,col=tcol('red',50))
    plot(out$BFP$PFS$sig1$MA$P,out$BFP$PFS$sig1$MA$logBF-out$BFP$PFS$sig1$white$logBF,xlab='Period [day]',ylab='relative ln(BF)',main='MA(1)-MA(0)',log='x',type='l')
    abline(v=Popt,col=tcol('red',50))
#    abline(v=Popt,col=tcol('red',50))
}
dev.off()
