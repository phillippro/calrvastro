load('results/HD222237/HD222237_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1100000_sppriorFALSE_dr1FALSE_230617_Nset4_hg123_Nsig1_P21646_Esd1_astro5TRUE_acc0.81_rvcbarycentric_lnpmax-705_lnlmax-677.Robj')
mc <- out$mcmc.opt$sig1
npar <- ncol(mc)-2
out <- array(NA,dim=c(npar,npar))
mc[,1] <- exp(mc[,1])#logP to P
for(i in 1:npar){
    for(j in i:npar){
        out[i,j] <- out[j,i] <- cor(mc[,i],mc[,j])
    }
}
colnames(out) <- rownames(out) <- colnames(mc)[1:npar]
