if(!exists('fs')){
   fs <- list.files('results',pattern='Robj',recursive=TRUE,full.name=TRUE)
   fs <- fs[grep('_cal|lpspm|photovary',fs)]
}
targets  <- gsub('_.+','',gsub('.+\\/','',fs))
lls <- as.numeric(gsub('.Robj','',gsub('.+_lnlmax','',fs)))
stars <- unique(targets)

out <- c()
for(s in stars){
      ii <- which(targets==s)
      out <- c(out,fs[ii[which.max(lls[ii])]])
}
write.table(out,file='all_companions.txt',quote=FALSE,col.names=FALSE,row.names=FALSE)
