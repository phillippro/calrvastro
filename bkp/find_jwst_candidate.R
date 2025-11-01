tab <- read.table('companion_par0_med.txt',header=TRUE)
host <- read.table('host_par_m0_1000.txt',header=TRUE)
ind <- match(tab[,1],host[,1])
host <- host[ind,]

mu <- 5
ind <- which(tab[,'mpJ1.opt']<mu | tab[,'mpJ2.opt']<mu | tab[,'mpJ3.opt']<mu | tab[,'mpJ4.opt']<mu)
jup <- tab[ind,]
host <- host[ind,]

nn <- sapply(1:nrow(jup),function(i) length(which(!is.na(jup[i,grep('mpJ\\d.opt',colnames(jup))]))))
#ind1 <- which(nn>1)
ind1 <- which(nn>0)
out <- c()
for(j in ind1){
    ms <- jup[j,grep('mpJ\\d.mean',colnames(jup))]
    dms <- jup[j,grep('mpJ\\d.sd$',colnames(jup))]
    as <- jup[j,grep('^a\\d.mean',colnames(jup))]
    das <- jup[j,grep('^a\\d.sd$',colnames(jup))]
    ems <- dms/ms
    rho <- as*host[j,'plxG']/1e3#arcsec
    drho <- das*host[j,'plxG']/1e3#arcsec
    N <- length(which(ms<mu & ms>0.2 & ems<0.2 & rho>0.4))
    if(N>0){
        out <- rbind(out,c(jup[j,1],rho,drho,ms,dms,as,das))
        cat(jup[j,1],'\n')
    }
}
nn <- 1:4
colnames(out) <- c('star',paste0('rho.as',nn),paste0('erho.as',nn),paste0('mjup',nn),paste0('emjup',nn),paste0('a',nn),paste0('ea',nn))
#jup[ind1,1:3]
