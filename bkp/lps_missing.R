tab <- read.table('companion_files.txt')[,1]
#tab <- read.table('all_targets.txt.txt')[,1]
stars <- gsub('\\/.+','',tab)
missing <- c()
for(star in stars){
    fs <- list.files(paste0('results/',star),pattern='_hg123_.+pdf')
    cat('star:',star,'\n')
    cat(length(fs),'\n')
    if(length(fs)==0){
          missing <- c(missing,star)
    }
}
write.table(missing,file='candidates1.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
