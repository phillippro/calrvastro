red3.dir <- list.dirs('../red3/results/',recursive=FALSE)
hg2.dir <- list.dirs('../rvastro/results/',recursive=FALSE)
hg123.dir <- list.dirs('results/',recursive=FALSE)
red3.dir <- gsub('.+\\/','',red3.dir)
hg2.dir <- gsub('.+\\/','',hg2.dir)
hg123.dir <- gsub('.+\\/','',red3.dir)
red3.dir <- red3.dir[grepl('^[A-Z]|^[a-z]',red3.dir)]
hg2.dir <- hg2.dir[grepl('^[A-Z]|^[a-z]',hg2.dir)]
hg123.dir <- hg123.dir[grepl('^[A-Z]|^[a-z]',hg123.dir)]
red3.dir <- red3.dir[grepl('^HD|^GJ|^GL|^XO|^HIP|^BD|^WASP|^TYC|^KOI|^K2|^Kepler|^HTR|^EPIC|^CoRoT|^CLT|^CD|^G\\d|^PDS|^LHS|^LSPM|teegarden|^RXJ|^HR|^KELT|^TIC|^TOI|^TRES',red3.dir)]
hg2.dir <- hg2.dir[grepl('^HD|^GJ|^GL|^XO|^HIP|^BD|^WASP|^TYC|^KOI|^K2|^Kepler|^HTR|^EPIC|^CoRoT|^CLT|^CD|^G\\d|^PDS|^LHS|^LSPM|teegarden|^RXJ|^HR|^KELT|^TIC|^TOI|^TRES',hg2.dir)]
hg123.dir <- hg123.dir[grepl('^HD|^GJ|^GL|^XO|^HIP|^BD|^WASP|^TYC|^KOI|^K2|^Kepler|^HTR|^EPIC|^CoRoT|^CLT|^CD|^G\\d|^PDS|^LHS|^LSPM|teegarden|^RXJ|^HR|^KELT|^TIC|^TOI|^TRES',hg123.dir)]
out <- c()
for(d in red3.dir){
    fs <- list.files(paste0('../red3/results/',d),pattern='pdf')
    lls <- as.numeric(gsub('.+_lnlmax|\\.pdf','',fs))
    if(length(lls)>0){
        if(!all(is.na(lls))){
            f <- fs[which.max(lls)]
            out <- rbind(out,c(d,'red3',f))
        }
    }
}
for(d in hg2.dir){
    fs <- list.files(paste0('../rvastro/results/',d),pattern='pdf')
    lls <- as.numeric(gsub('.+_lnlmax|\\.pdf','',fs))
    if(length(lls)>0){
        if(!all(is.na(lls))){
            f <- fs[which.max(lls)]
            out <- rbind(out,c(d,'rvastro',f))
        }
    }
}
for(d in hg123.dir){
    fs <- list.files(paste0('results/',d),pattern='pdf')
    lls <- as.numeric(gsub('.+_lnlmax|\\.pdf','',fs))
    if(length(lls)>0){
        if(!all(is.na(lls))){
            f <- fs[which.max(lls)]
            out <- rbind(out,c(d,'calrvastro',f))
        }
    }
}
write.table(out,file='all_targets.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
