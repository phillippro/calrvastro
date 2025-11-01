arg <- commandArgs(trailingOnly=TRUE)
#arg <- 'HD85440'
f1 <- 'companion_files.txt'
files <- read.table(f1)[,1]
f2 <- 'new_companion.txt'
new <- readLines(f2)
f3 <- 'companion_m5_120.txt'
trim <- read.table(f3,header=TRUE)
cat(f1,'\n')
cat(f2,'\n')
cat(f3,'\n')
#cat('^',arg,'\/\n')
if(any(grepl(paste0('^',arg,'\\/'),files))){
    files <- files[-grep(paste0('^',arg,'\\/'),files)]
    cat('rm',arg,'entry from companion_files.txt!\n')
    write.table(files,file=f1,quote=FALSE,row.names=FALSE,col.names=FALSE)
}
if(any(grepl(paste0('^',arg,'\\/|',arg,' '),new))){
    new <- new[-grep(paste0('^',arg,'\\/|',arg,' '),new)]
    cat('rm',arg,'entry from new_companion.txt!\n')
    write.table(new,file=f2,quote=FALSE,row.names=FALSE,col.names=FALSE)
}
if(any(trim[,1]==arg)){
    trim <- trim[-which(trim[,1]==arg),]
    cat('rm',arg,'entry from companion_m5_120.txt!\n')
    write.table(trim,file=f3,quote=FALSE,row.names=FALSE)
}
