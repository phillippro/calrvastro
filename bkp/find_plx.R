tab  <- read.table('all_par.txt',header=TRUE)
plxs <- c()
for(j in 1:nrow(tab)){
    star <- tab[j,1]
    f <- paste0('../data/combined/',star,'/',star,'_hipgaia.hg123')
    if(!file.exists(f)){
        f <- paste0('../data/combined/',gsub('J','L',star),'/',gsub('J','L',star),'_hipgaia.hg123')
    }
    if(!file.exists(f)){
        plxs <- c(plxs,NA)
    }else{
        dat <- read.table(f,header=TRUE)
        plxs <- c(plxs,dat[nrow(dat),'parallax'])
    }
}
write.table(plxs,file='plx.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
