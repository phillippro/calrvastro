targets <- read.table('PFS_targets.txt')
servers <- c('intrepid','paradise','marx','benway','ginsberg','babbs','moriarty','further')
cores <- c(16,8,24,8,24,24,24,24)
N <- 4
for(j in 1:length(servers)){
      if(j==1){
           ind.min <- 1
      }else{
           ind.min <- N*sum(cores[1:(j-1)])+1
      }
      ind.max <- min(N*sum(cores[1:j]),nrow(targets))
      cat('ind.min=',ind.min,';ind.max=',ind.max,'\n')
      fout <- paste0('PFS_',servers[j],'.txt')
      cat('fout:',fout,'\n')
      write.table(targets[ind.min:ind.max,1],file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
}

