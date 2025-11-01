#temps <- tem.min*10^(log10(seq(10^(0),10^(-log10(tem.min)),length.out=10)))
temps <- exp(seq(log(tem.min),log(1),length.out=10))
acceptance <- 100
Niter <- 1000
for(tem in temps){
    Nsamp <- Niter/Ncores
    amh <- AMH(round(Nsamp/2),Nsamp)
    out <- amh$out
    startvalue <- out[which.max(out[,Npar+2]),1:Npar]
    cov <- amh$cov
    mcmc.out <- out[,1:Npar]
    post.out <- out[,Npar+1]
    like.out <- out[,Npar+2]
    acceptance  <-  (1-mean(duplicated(mcmc.out)))*100
    cat('\ntem=',tem,'\n')
    cat('Nsamp=',Nsamp,'\n')
    cat('acceptance:',acceptance,'%\n')
    if(acceptance<25) break()
    if(acceptance>=50) Niter <- 1000
    if(acceptance<50 & acceptance>=25) Niter <- 10000
}
