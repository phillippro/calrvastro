Prot <- Protation
tau <- NA
###white noise model
cat('\nlikelihood for the white noise model for the rotation proxy:\n')
tem <- 1
GP <- FALSE
Nma.opt <- Nar.opt <- 0
data <-tab
Np <- 0
source('setting_up.R')
source('initial_condition.R')
out <- AMH(round(Niter/2/Ncores),Niter/Ncores)$out
logLmax.white <-max(out[,ncol(out)])

cat('\nlikelihood for the GP model for the rotation proxy:\n')
GP <- TRUE
Nma.opt <- 0
Nar.opt <- 0
Np <- 0
source('setting_up.R')
source('initial_condition.R')
out <- AMH(round(Niter/2/Ncores),Niter/Ncores)$out
logLmax.GP <-max(out[,ncol(out)])

logBF.GP <- logLmax.GP-logLmax.white-1.5*log(nrow(data))
if(logBF.GP>5){
    Prots <- exp(mcmc.out[,'1logProt1'])
    taus <- exp(mcmc.out[,'1logl1'])*Prots/pi
    logtau.stat <- data.distr(log(taus),plotf=FALSE)
    tau <- exp(logtau.stat['mean'])
    tau.sd <- exp(logtau.stat['sd'])
    tau.opt <- taus[ind.opt]
    logProt.stat <- data.distr(log(Prots),plotf=FALSE)
    Prot.opt <- Prots[ind.opt]
    Prot <- exp(logProt.stat['mean'])
    Prot.sd <-exp(logProt.stat['sd'])
}
