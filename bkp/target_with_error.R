fs <- list.files(path='logs',pattern='log$')
targets <- c()
for(f in fs){
      s <- readLines(paste0('logs/',f))
      if(grepl('Execution halted',s[length(s)])){
	  targets <- c(targets,gsub('_PFS.+','',f))
      }
}
write.table(targets,file='PFS_kesey.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
