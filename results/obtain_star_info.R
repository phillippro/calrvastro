ds  <- list.dirs('.')
ds <- ds[ds!='.']
ds  <- gsub('\\./','',ds)
stars <- c()
for(d in ds){
    fs <- list.files(d,pattern='hg123')
    if(length(fs)>0) stars <- c(stars,d)
}
fout <- 'stars_with_hg123.txt'
cat(fout,'\n')
write.table(stars,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
