args <- commandArgs(trailingOnly=TRUE)
tab <- read.table('targets_m0_1e+05_wc.txt',header=TRUE)
print(tab[tab[,1]==args[1],c('Star','sep_AU', 'binary_type')])