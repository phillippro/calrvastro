res <- out$res.comb.all$sig1
t <- out$all[,1]
rv <- out$all[,2]
erv <- out$all[,3]
pdf('einstein_fit.pdf',16,16)
par(mfrow=c(4,4))
plot(t,rv,xlab='JD',ylab='RV [m/s]',main=paste0('RMS=',round(sd(rv)),'\n'))
plot(t,res,xlab='JD',ylab='O-C [m/s]',main=paste0('RMS=',round(sd(res)),'\n'))
par.opt <- out$par.stat$sig1[1,]
Popt <- exp(par.opt['per1'])
fit <- lm(res~t)
res1 <- fit$residuals
phase <- (t-min(t))%%Popt
plot(phase,res1,xlab='Phase [day]',ylab='O-C [m/s]', main='phase-folded residual')
dev.off()
