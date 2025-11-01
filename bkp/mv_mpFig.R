args <- commandArgs(trailingOnly=TRUE)
target <- gsub('b|c|d|_.+','',args)
system(paste0('cp results/',target,'/',args,' ../paper/earth_auto2/v1/'))