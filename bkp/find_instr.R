targets <- read.table('BD_hosts.txt')[,1]
stars <- ii <- c()
info <- list()
for(star in targets){
    dir.in <- paste0('../data/combined/',star)
    fs <- list.files(dir.in,pattern='rv$|dat$|vels$')
    fs <- fs[!grepl('photo|ASAS',fs)]
    ins <- gsub('.+_|\\..+','',fs)
    ii <- c(ii,ins)
    stars <- c(stars,rep(star,length(ins)))
    Ns <- c()
    for(f in fs){
        Ns <- c(Ns,nrow(read.table(paste0(dir.in,'/',f),header=TRUE)))
    }
    info[[star]] <- list()
    info[[star]][['ins']] <- ins
    info[[star]][['Nrv']] <- Ns
}
out <- cbind(ii,stars)
inss <- c('CORALIE','LICK','ELODIE','SOPHIE')
for(i in inss){
ind <- grep(i,out[,1])
cat(i,'data:',unique(out[ind,2]),'\n')
}
