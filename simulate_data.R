source('mcmc_func.R')
#load('../output/HD128621/keppure_priormt_poly2_Ndata16868_quantifyTRUE_1per1_HD128621_TERRA_1AP1_erv_ind0.2.5.3.1.14.12.18.9_0planet_ARMA09_Nsamp4000000_tem1_acc23_pretem1Pd_negLmax27042.Robj',envir= e <- new.env())
load('../output/LHS1140/keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_2planet_GPqp_Nsamp4000000_tem1_acc0.67_pretem1P93.1d24.5d_negLmax400.Robj',envir= e <- new.env())
par.data <- e$par.data
ids <- e$ids
tmin <- e$tmin
kep.type <- e$kep.type
period.par <- e$period.par
res <- e$RV-e$RV.kep[[1]]
tab <- e$data
tab[,2] <- res
#out <- tab <- read.table('../data/aperture/HD128621/HD128621_HARPS.dat',header=TRUE,check.names=FALSE)
#out <- tab <- read.table('../data/aperture/LHS1140/LHS1140.dat',header=TRUE,check.names=FALSE)
out <- tab
###planetary signals to be injected
Ks <- c(0.4,0.8,1.6,3.2,6.4,12.8,5.6)
es <- c(0,0.2)
Ps <- c(3,5, 12,17, 59, 109, 130, 310)
#####
for(i in length(Ks)){
    for(j in 1:length(es)){
        for(k in 1:length(Ps)){
            rv.kep <- RV.kepler(pars.kep=c('per1'=log(Ps[k]),'K1'=Ks[i],'e1'=es[j],'omega1'=0,'Mo1'=0),tt=tab[,1]%%2400000,Np.kep=1,prior.kep='mt',kep.only=TRUE)[[1]]
            out[,2] <- tab[,2]+rv.kep
            cat('rv.kep[1:10]=',rv.kep[1:10],'\n')
            cat('tab[1:10,2]=',tab[1:10,2],'\n')
            cat('out[1:10,2]=',out[1:10,2],'\n')
            fout <- paste0('../data/aperture/LHS1140/LHS1140_resGP2_sim',i,j,k,'.dat')
            cat('output file:',fout,'\n')
            write.table(out,file=fout,row.names=FALSE,quote=FALSE)
        }
    }
}
