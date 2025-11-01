fs <- read.table('BD_files.txt')[,1]
for(f in fs){
system(paste0('cp results/',f,' BD_targets/'))
}
