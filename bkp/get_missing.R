stars <- read.table('candidates_more.txt')[,1]
missing <- c()
for(star in stars){
     fs <- list.files(paste0('../red3/results/',star),pattern='*Robj')
     fs <- c(fs,list.files('optpar',pattern=paste0(star,'_')))
     if(length(fs)==0){
	 missing <- c(missing,star)
     }else{
         fout <- paste0('pars/',star,'.par')
         if(!file.exists(fout)){

             write.table(optpar,file=fout)
             stop()
         }
     }
}
write.table(missing,file='star_with_missing_red3_solutions.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
