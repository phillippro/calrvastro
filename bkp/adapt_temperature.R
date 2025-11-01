Ntem <- round(log(1/tem.min)/log(2))#tem.min*2^n=tem; tem<=1
Ntmp <- Nexp/10
i1 <- 0
for(i0 in 0:Ntem){
    tem <- min(tem.min*2^i0,1)
    cat('Adaptive tempering burning: tem=',format(tem,digit=3),'\n')
    nburn <- round(Ntmp/2)
    amh <- AMH(nburn,Ntmp)
    out <- amh$out
    cov <- amh$cov
    mcmc.out <- out[,1:Npar]
    acceptance = 1-mean(duplicated(mcmc.out))
    cat('acceptance:',acceptance*100,'\n')
    if((acceptance*100)<(tem.up+5)){
        cat('Finding the optimal temperture with longer chain!\n')
        Ntmp <- Nexp
        if(i1>0 & (acceptance*100)<(tem.up-5)) break()
#        if((acceptance*100)<(tem.up-5)) break()
        i1 <- i1+1
    }
}
cat('run chains using tem=',tem,' and initial period of ',Pini,'\n')
