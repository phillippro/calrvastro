#targets <- as.character(unlist(read.table('miss_list.txt')))
targets <- as.character(unlist(read.table('pfs2_targets.txt')))
for(target in targets){
#    fs <- list.files(path=paste0('../data/combined/',target,'/'),pattern='')
#    cat('fs:',fs,'\n')
    fs <- list.files(path=paste0('../data/combined/',target,'/'),pattern='HARPS|harps')
    if(length(fs)==0){
        cat(target,'\n')
    }
}
