#fs <- read.table('companion_files.txt')[,1]
tab <- read.table('companion_m5_120.txt',header=TRUE)
index <- sapply(1:nrow(tab),function(i) which(!is.na(tab[i,grepl('^mpJ\\d.opt',colnames(tab))])))
Ns <- sapply(1:nrow(tab),function(i) length(which(!is.na(tab[i,grepl('^mpJ\\d.opt',colnames(tab))]))))
inds <- which(Ns>1)
cat('Number of multiple-companion systems:',length(inds),'\n')
fout <- 'companion_multi.txt'
write.table(tab[inds,],file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
