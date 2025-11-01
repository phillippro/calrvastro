ds <- list.dirs('.')[-1]
ff <- c()
N <- length(ds)
for(i in 1:N){
    cat(i,'/',N,'\n')
    d <- ds[i]
   fs <- list.files(d,pattern='pdf',full.name=TRUE)
   dlnl <- as.numeric(gsub('.+_lnlmax|.pdf','',fs))
   f <- fs[which.max(dlnl)]
    cat(f,'\n')
   ff <- c(ff,f)
}
ff1 <- gsub('\\.\\/','',ff)
stars <- gsub('\\/.+','',ff1)
#write.csv(cbind(stars,ff1),file='solution_list.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(ff1,file='solution_list.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(stars,file='solution_ids.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
