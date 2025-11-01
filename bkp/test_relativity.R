load('results/HD56380/HD56380_fix1_relativityFALSE_Niter4000000_Ncores8_ofac2_Nset2_Esd1_transit0_P3254_acc3.2_lnlmax-171.Robj')
source('mcmc_func.R')
kep0 <- RV.kepler(par.opt)
par.opt[c('w11','w21')] <- 0
ll0 <- loglikelihood(par.opt)
rv0 <- c(kep0$rv$HARPSpre,kep0$rv$HARPSpost)
res0 <- out$all[,2]-rv0
llike0 <- sum(dnorm(res0,0,out$all[,3],log=TRUE))

load('results/HD56380/HD56380_fix1_relativityTRUE_Niter4000000_Ncores8_ofac2_Nset2_Esd1_transit0_P3255_acc3_lnlmax-171.Robj')
source('mcmc_func.R')
kep1 <- RV.kepler(par.opt)
par.opt[c('w11','w21')] <- 0
ll1 <- loglikelihood(par.opt)
rv1 <- c(kep1$rv$HARPSpre,kep1$rv$HARPSpost)
res1 <- out$all[,2]-rv1
llike1 <- sum(dnorm(res1,0,out$all[,3],log=TRUE))
gs <- kep1$rvs+kep1$rvg

out$relativity <- FALSE
kep2 <- RV.kepler(par.opt)
par.opt[c('w11','w21')] <- 0
ll2 <- loglikelihood(par.opt)
rv2 <- c(kep2$rv$HARPSpre,kep2$rv$HARPSpost)
res2 <- out$all[,2]-rv2
llike2 <- sum(dnorm(res2,0,out$all[,3],log=TRUE))

cat('ll0=',ll0,'\n')
cat('llike0=',llike0,'\n')

cat('ll1=',ll1,'\n')
cat('llike1=',llike1,'\n')

pdf('relativity.pdf',12,12)
par(mfrow=c(3,3))
plot(trv.all,rv0,ylim=range(rv0,rv1))
plot(trv.all,out$all[,2],ylim=range(rv0,rv1))
points(trv.all,rv0,col='red')
points(trv.all,rv1,col='blue')
plot(trv.all,rv1,ylim=range(rv0,rv1))
plot(trv.all,rv1-rv0)
plot(trv.all,rv1-gs-rv0)
plot(trv.all,rv2-rv0)
plot(trv.all,gs)
dev.off()
