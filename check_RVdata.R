tab <- read.table('companion_m5_120.txt',header=TRUE)
stars <- tab[,1]
ins <- c()
Nrv <- c()
for(star in stars){
   fs <- list.files(paste0('../data/combined/',star),pattern='rv|dat|vels',full.name=TRUE)
   for(f in fs){
        i1 <- gsub('.+\\/','',f)
        i2 <- gsub('.+_','',i1)
        i3 <- gsub('\\..+','',i2)
        ins <- c(ins,i3)
        tab <- read.table(f,header=TRUE)
        Nrv <- c(Nrv,nrow(tab))
   }
}
ii <- unique(ins)
nn <- c()
for(i in ii){
     index <- which(ins==i)
     n <- sum(Nrv[index])
     nn <- c(nn,n)
}
jj <- sort(nn,index.return=TRUE,decreasing=TRUE)$ix
print(cbind(ii[jj],nn[jj]))
#     cat('Number of RVs contributed by',i,'is',n,'\n')
