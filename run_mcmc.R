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
    x <- mcmc.out[,j]
    if(grepl('per',colnames(mcmc.out)[j]) & period.par=='logP'){
        x <- exp(x)
    }
    tmp <- data.distr(x,loglike=loglike.out,xlab=colnames(mcmc.out)[j],ylab='Freq.',plotf=FALSE)
    par.stat <- cbind(par.stat,tmp)
}
colnames(par.stat) <- names(startvalue)
