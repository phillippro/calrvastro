source('mcmc_func.R')
#target <- 'HD26965_CHI'
target <- 'HD26965_KECK_noout'
star <- gsub('_.+','',target)
file <- paste0('../data/aperture/',star,'/',target,'.vels')
cat('read file:\n')
cat(file,'\n')
tab <- read.table(file)
#Ks <- c(0.4,0.8,1.6,3.2,6.4,12.8,5.6)
#Ps <- c(42,25)
Ps <- c(42,25,60)
es <- 0
#c(0.0,0.07)
Ks <- c(1.5,1.83)
#c(1.83,1.5)
#es <- c(0,0.2)
#Ps <- c(3,5, 12,17, 59, 109, 130, 310)
#####
out <- tab
for(i in length(Ks)){
    for(j in 1:length(es)){
        for(k in length(Ps)){
            rv.kep <- RV.kepler(pars.kep=c('per1'=log(Ps[k]),'K1'=Ks[i],'e1'=es[j],'omega1'=0,'Mo1'=0),tt=tab[,1]%%2400000,Np.kep=1,prior.kep='mt',kep.only=TRUE)[[1]]
            out[,2] <- tab[,2]+rv.kep
            fout <- paste0('../data/aperture/',star,'/',target,'_sim',i,j,k,'.dat')
            cat('output file:\n',fout,'\n')
            write.table(out,file=fout,row.names=FALSE,quote=FALSE,col.names=FALSE)
        }
    }
}
