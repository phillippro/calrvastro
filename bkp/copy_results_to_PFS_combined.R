fs <- commandArgs(trailingOnly=TRUE)
targets <- gsub('.+/|_.+','',fs)
fs <- gsub('.+/','',fs)
for(k in 1:length(fs)){
      system(paste0('cp results/',targets[k],'/',fs[k],'*pdf ../PFS_combined/results/'))
}