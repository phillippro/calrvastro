hpc <- read.table('bd_hpc.txt',header=TRUE)
rack <- read.table('bd_rack.txt',header=TRUE)
out <- rbind(hpc,rack)
write.table(out,file='../../../BDproject/bd_list.txt',quote=FALSE,row.names=FALSE)
