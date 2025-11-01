stars <- read.table('candidates0.txt')[,1]
for(star in stars){
f <- paste0('logs/',star,'_Nmax1_fix1_Niter1.4e6_Ncore16_natural_relTRUE_hg3.log')
cat('\n',f,'\n')
system(paste('tail',f))
}
