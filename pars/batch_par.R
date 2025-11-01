library(foreach)
library(doMC)
library(doParallel)
library(parallel)
Ncores <- 16
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}

dir.in <- '../../data/combined'
dirs <- list.dirs(dir.in)
stars <- unique(gsub('.+\\/','',dirs))
stars <- stars[grepl('GL|GJ|HD|HIP|WASP|XO|TYC|K2|BD|CoRoT|CD',stars)]
#ds <- list.dirs('../../data/combined')
p <- foreach(star = stars) %dopar% {
#for(star in stars){
    fpar <- paste0(star,'.par')
    	 cat(fpar,'\n')
    if(!file.exists(fpar)){
#      fs <- list.files(paste0('../../rvastro/results/',star),pattern='Robj')
      fs <- list.files(paste0('../../red3/results/',star),pattern='Robj')
      if(length(fs)>0){
         lls <- as.numeric(gsub('.+_lnlmax|\\.Robj','',gsub('.Robj','',fs)))
         fs <- fs[which.max(lls)]
         cmd <- paste('Rscript get_par.R',fs)
         cat(cmd,'\n')
         system(cmd)
      }
    }
}
