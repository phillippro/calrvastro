red <- readLines('new_red.txt')
jupiter <- readLines('new_jupiter.txt')
bd <- readLines('new_bd.txt')
fout <- 'new_companion.txt'
cat(fout,'\n')
out <- c(red,jupiter,bd)
out <- out[!grepl('list|List',out)]
write.table(out,file=fout,row.names=FALSE,quote=FALSE,col.names=FALSE)

