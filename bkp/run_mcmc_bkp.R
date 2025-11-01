##########################################################################
####this code is for parameter inference
##########################################################################
Np <- j3-1
if(Np==0){
    iniref <- c()
    par.min0 <- par.min <- c()
    par.max0 <- par.max <- c()
}else{
    if(period.par=='P'){
        iniref <- c(Popt,Kopt,0,0,0)
        par.min0 <- par.min <- c(Pmin,Kmin,0,0,0)
        par.max0 <- par.max <- c(Pmax,Kmax,1,2*pi,2*pi)
    }else{
        iniref <- c(log(Popt),Kopt,0,0,0)
        par.min0 <- par.min <- c(log(Pmin),Kmin,0,0,0)
        par.max0 <- par.max <- c(log(Pmax),Kmax,1,2*pi,2*pi)
    }
}
logPmaxs <- rep(log(Pmax),Np)
logPmins <- rep(log(Pmin),Np)
Sd <- 0.2
if(Np>1){
    iniref <- c(par.opt[1:(Nkeppar*(Np-1))],iniref)
    par.min <- c(par.min,par.min0)
    par.max <- c(par.max,par.max0)
}

amax <- (Kmax-Kmin)/DT
amin <- -amax
bmin <- Kmin
bmax <- Kmax
phi.min <- wmin <- -1
phi.max <- wmax <- 1
smin <- 0
smax <- Kmax
beta.up <- log(tmax-tmin)#time span of the data
beta.low <- log(max(1/24,min(1,min(diff(trv.all)))))#1h or minimum separation
alpha.max <- beta.max <- beta.up#d; limit the range of beta to avoid multimodal or overfitting
alpha.min <- beta.min <- beta.low#24h
if(par.global=='trend'){
#linear trend
#    p <- lm(tab[,2] ~ tab[,1])
#    iniref <- c(iniref,c(p$coeff[2],p$coeff[1]))
    iniref <- c(iniref,c(0,0))
    par.min <- c(par.min,c(amin,bmin))
    par.max <- c(par.max,c(amax,bmax))
}else{
    iniref <- c(iniref,-mean(tab[,2]))#intersect
    par.min <- c(par.min,Kmin)
    par.max <- c(par.max,Kmax)
}

###jitter and red noise parameters
iniref <- c(iniref,0.1*sd(tab[,2]))
par.min <- c(par.min,0)
par.max <- c(par.max,Kmax)
if(Nma.opt>0){
    iniref <- c(iniref,rep(0.1,Nma.opt),0)
    par.min <- c(par.min,rep(0,Nma.opt),log(1e-5))
    par.max <- c(par.max,rep(1,Nma.opt),log(1e5))
}

startvalue <- assign.names(iniref,Np=Np,p=Nar.opt,q=Nma.opt)
Npar <- length(startvalue)
Sd <- 2.4^2/Npar#hyper par of s prior
cov.start <- 1e-6*diag(Npar)

if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
chain.type <- 'adapt'
source('gen_mcmc.R',local=TRUE)

cat('\n convergency test!\n\n')
source('convergence_test.R',local=TRUE)
#source('figure_gen.R',local=TRUE)
###find optimal parameters and their uncertainty
ind.opt <- which.max(post.out)
logL.max <- max(loglike.out)
logP.max <- max(post.out)
par.opt <- mcmc.out[ind.opt,1:Npar]
par.stat <- c()
for(j in 1:Npar){
    tmp <- data.distr(mcmc.out[,j],xlab=colnames(mcmc.out)[j],ylab='Freq.',plotf=FALSE)
    par.stat <- cbind(par.stat,tmp)
}
colnames(par.stat) <- names(startvalue)
