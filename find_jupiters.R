tab <- read.table('companion_par0_med.txt',header=TRUE)
ind <- which(tab[,'mpJ1.opt']<5 | tab[,'mpJ2.opt']<5 | tab[,'mpJ3.opt']<5 | tab[,'mpJ4.opt']<5)
jup <- tab[ind,]
nn <- sapply(1:nrow(jup),function(i) length(which(!is.na(jup[i,grep('mpJ\\d.opt',colnames(jup))]))))
ind1 <- which(nn>1)
for(j in ind1){
    ms <- jup[j,grep('mpJ\\d.mean',colnames(jup))]
    dms <- jup[j,grep('mpJ\\d.sd$',colnames(jup))]
    ems <- dms/ms
    N <- length(which(ms<5 & ms>0.2 & ems<0.2))
    if(N>1){
        cat(jup[j,1],'\n')
    }
}
#jup[ind1,1:3]
