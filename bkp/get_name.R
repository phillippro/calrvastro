#type <- 'caltrue'
#type <- 'calfalse'
type <- commandArgs(trailingOnly=TRUE)
tab <- read.table(paste0(type,'.txt'))[,1]
stars <- c()
for(f in tab){
     stars <- c(stars,gsub('.+\\/|_.+','',f)) 
}
fout <- paste0('candidates_',type,'.txt')
cat(fout,'\n')
write.table(stars,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)